// parallel_bench - the same Kerr image rendered with different parallel
// programming techniques, plus a false-sharing experiment.
//
//   parallel_bench [--res 256] [--threads N] [--type double|float] [--reps 3]
//
// Every technique renders the identical image (a=0.9375, i=60 deg, thin disk,
// RK4 STEPSIZE 0.03) and the result is checked against the serial render, so
// only the scheduling differs.  The cost of a ray varies a lot (rays that skim
// the photon ring take ~10x more steps), which is what separates static from
// dynamic scheduling.
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <execution>
#include <functional>
#include <numeric>
#include <string>
#include <thread>
#include <vector>

#include <omp.h>

#include "kerr/render.hpp"

using namespace kerr;
namespace stdx = std::experimental;
using clk = std::chrono::steady_clock;

template <class T> struct Scene {
    Spacetime<T> st{T(0.9375)};
    Camera<T> cam;
    TraceConfig<T> cfg;
    explicit Scene(int res) : cam{T(1e4), T(60) * m::pi<T>() / T(180), T(40), res, res} {
        cfg.method = Method::RK4;
        cfg.step_frac = T(0.03);
        cfg.disk.enabled = true;
        cfg.disk.r_in = st.r_isco();
        cfg.disk.r_out = T(20);
    }
};

// Per-thread busy time -> load imbalance = max / mean
struct Busy {
    std::vector<double> sec;
    explicit Busy(int n) : sec(n, 0.0) {}
    double imbalance() const {
        const double mx = *std::max_element(sec.begin(), sec.end());
        const double mean = std::accumulate(sec.begin(), sec.end(), 0.0) / sec.size();
        return mean > 0 ? mx / mean : 1.0;
    }
};

// ---------------------------------------------------------------------------
// SIMD packet tracer: W neighbouring pixels advance in lock-step.  Finished
// lanes are frozen with dl = 0 (they still cost arithmetic -> "lane
// utilisation" measures the SIMD divergence overhead).
// ---------------------------------------------------------------------------
template <class T> struct PacketStats {
    long vector_steps = 0, lane_steps = 0;
};

template <class T>
void trace_packet(const Scene<T>& sc, int first, int count, Image& img, PacketStats<T>& ps) {
    using V = stdx::native_simd<T>;
    using M = typename V::mask_type;
    constexpr int W = int(V::size());
    const Spacetime<V> sv{V(sc.st.a)};

    alignas(64) T buf[8][W];
    for (int l = 0; l < W; ++l) {
        const int idx = first + std::min(l, count - 1);  // pad with a duplicate ray
        T al, be;
        sc.cam.impact(idx % sc.cam.width, idx / sc.cam.width, al, be);
        const auto y0 = sc.cam.initial_state(sc.st, al, be);
        for (int k = 0; k < 8; ++k) buf[k][l] = y0[k];
    }
    State<V> y, k1, k2, k3, k4, tmp;
    for (int k = 0; k < 8; ++k) y[k].copy_from(buf[k], stdx::vector_aligned);

    M active(true);
    for (int l = count; l < W; ++l) active[l] = false;
    V status(T(0)), inten(T(0)), nsteps(T(0));
    const T r_stop = sc.st.r_plus * (T(1) + sc.cfg.horizon_margin);
    const T frac = sc.cfg.step_frac;
    const V pi(m::pi<T>());

    for (long s = 0; s < sc.cfg.max_steps && stdx::any_of(active); ++s) {
        sv.rhs(y, k1);
        // RAPTOR step heuristic, vectorised; inactive lanes get dl = 0
        const V dist_pole = stdx::min(y[2], pi - y[2]);
        const V inv = stdx::abs(k1[1]) / (y[1] * frac) + stdx::abs(k1[2]) / (dist_pole * frac) + stdx::abs(k1[3]) / frac;
        V dl = stdx::max(V(T(1)) / inv, V(T(1e-12)));
        where(!active, dl) = V(T(0));

        const V r_old = y[1], c_old = stdx::cos(y[2]);
        const V half = dl * T(0.5);
        for (int k = 0; k < 8; ++k) tmp[k] = y[k] + half * k1[k];
        sv.rhs(tmp, k2);
        for (int k = 0; k < 8; ++k) tmp[k] = y[k] + half * k2[k];
        sv.rhs(tmp, k3);
        for (int k = 0; k < 8; ++k) tmp[k] = y[k] + dl * k3[k];
        sv.rhs(tmp, k4);
        const V h6 = dl / T(6);
        for (int k = 0; k < 8; ++k) y[k] += h6 * (k1[k] + T(2) * k2[k] + T(2) * k3[k] + k4[k]);

        ps.vector_steps++;
        ps.lane_steps += stdx::popcount(active);
        where(active, nsteps) += V(T(1));

        // disk crossing
        const V c_new = stdx::cos(y[2]);
        const V f = c_old / (c_old - c_new);
        const V rc = r_old + f * (y[1] - r_old);
        const M hit = active && (c_old * c_new < T(0)) && (rc > sc.cfg.disk.r_in) && (rc < sc.cfg.disk.r_out);
        if (stdx::any_of(hit)) {
            for (int l = 0; l < W; ++l)
                if (hit[l]) inten[l] = sc.cfg.disk.intensity(sc.st, rc[l], y[7][l]);
            where(hit, status) = V(T(int(Status::Disk)));
            active = active && !hit;
        }
        const M cap = active && (y[1] < r_stop);
        where(cap, status) = V(T(int(Status::Captured)));
        active = active && !cap;
        sv.rhs(y, k1);  // dr/dl sign for the escape test
        const M esc = active && (y[1] > sc.cfg.r_out) && (k1[1] > T(0));
        where(esc, status) = V(T(int(Status::Escaped)));
        active = active && !esc;
    }
    where(active, status) = V(T(int(Status::MaxSteps)));
    for (int l = 0; l < count; ++l) {
        const int idx = first + l;
        img.intensity[idx] = double(inten[l]);
        img.status[idx] = std::uint8_t(status[l]);
        img.steps[idx] = long(nsteps[l]);
    }
}

// ---------------------------------------------------------------------------

struct Result {
    std::string name;
    double sec, imbalance;
    long mismatches;
    std::string note;
};

template <class T> long compare(const Image& a, const Image& b) {
    long bad = 0;
    for (size_t i = 0; i < a.status.size(); ++i)
        if (a.status[i] != b.status[i] || std::fabs(a.intensity[i] - b.intensity[i]) > 1e-6 * (1e-12 + std::fabs(a.intensity[i]))) ++bad;
    return bad;
}

template <class T> void run(int res, int nthreads, int reps) {
    Scene<T> sc(res);
    const int n = res * res;
    const int W = int(stdx::native_simd<T>::size());
    std::vector<Result> results;
    Image ref(res, res);

    auto time_it = [&](const std::string& name, auto&& body, Image& img, Busy& busy) {
        double best = 1e30, imb = 1;
        for (int r = 0; r < reps; ++r) {
            std::fill(busy.sec.begin(), busy.sec.end(), 0.0);
            const auto t0 = clk::now();
            body(img, busy);
            const double s = std::chrono::duration<double>(clk::now() - t0).count();
            if (s < best) best = s, imb = busy.imbalance();
        }
        return std::pair{best, imb};
    };
    auto add = [&](const std::string& name, auto&& body, const std::string& note = "") {
        Image img(res, res);
        Busy busy(nthreads);
        auto [s, imb] = time_it(name, body, img, busy);
        results.push_back({name, s, imb, compare<T>(ref, img), note});
        std::fprintf(stderr, "  %-34s %.3f s\n", name.c_str(), s);
    };
    auto pixel = [&](Image& img, int i) { trace_pixel(sc.st, sc.cam, sc.cfg, i % res, i / res, img); };

    // 1. serial
    {
        Busy busy(1);
        auto [s, imb] = time_it("serial", [&](Image& img, Busy&) { for (int i = 0; i < n; ++i) pixel(img, i); }, ref, busy);
        results.push_back({"serial", s, 1.0, 0, ""});
    }

    // 2. std::thread, contiguous blocks of rows
    add("std::thread static block", [&](Image& img, Busy& busy) {
        std::vector<std::jthread> pool;
        for (int t = 0; t < nthreads; ++t)
            pool.emplace_back([&, t] {
                const auto t0 = clk::now();
                const int lo = int(long(n) * t / nthreads), hi = int(long(n) * (t + 1) / nthreads);
                for (int i = lo; i < hi; ++i) pixel(img, i);
                busy.sec[t] = std::chrono::duration<double>(clk::now() - t0).count();
            });
    });

    // 3. std::thread, cyclic (row-interleaved) distribution
    add("std::thread static cyclic rows", [&](Image& img, Busy& busy) {
        std::vector<std::jthread> pool;
        for (int t = 0; t < nthreads; ++t)
            pool.emplace_back([&, t] {
                const auto t0 = clk::now();
                for (int row = t; row < res; row += nthreads)
                    for (int x = 0; x < res; ++x) pixel(img, row * res + x);
                busy.sec[t] = std::chrono::duration<double>(clk::now() - t0).count();
            });
    });

    // 4. std::thread + atomic work counter (dynamic self-scheduling)
    add("std::thread dynamic (atomic, 16)", [&](Image& img, Busy& busy) {
        std::atomic<int> next{0};
        std::vector<std::jthread> pool;
        for (int t = 0; t < nthreads; ++t)
            pool.emplace_back([&, t] {
                const auto t0 = clk::now();
                for (int lo; (lo = next.fetch_add(16, std::memory_order_relaxed)) < n;)
                    for (int i = lo; i < std::min(lo + 16, n); ++i) pixel(img, i);
                busy.sec[t] = std::chrono::duration<double>(clk::now() - t0).count();
            });
    });

    // 5-7. OpenMP schedules
    auto omp_variant = [&](const char* name, auto sched) {
        add(name, [&](Image& img, Busy& busy) {
#pragma omp parallel num_threads(nthreads)
            {
                const auto t0 = clk::now();
                sched(img);
                busy.sec[omp_get_thread_num()] = std::chrono::duration<double>(clk::now() - t0).count();
            }
        });
    };
    omp_variant("OpenMP schedule(static)", [&](Image& img) {
#pragma omp for schedule(static) nowait
        for (int i = 0; i < n; ++i) pixel(img, i);
    });
    omp_variant("OpenMP schedule(static,1)", [&](Image& img) {
#pragma omp for schedule(static, 1) nowait
        for (int i = 0; i < n; ++i) pixel(img, i);
    });
    omp_variant("OpenMP schedule(dynamic,16)", [&](Image& img) {
#pragma omp for schedule(dynamic, 16) nowait
        for (int i = 0; i < n; ++i) pixel(img, i);
    });
    omp_variant("OpenMP schedule(guided)", [&](Image& img) {
#pragma omp for schedule(guided) nowait
        for (int i = 0; i < n; ++i) pixel(img, i);
    });

    // 8. C++17 parallel algorithm (libstdc++ -> oneTBB work stealing)
    add("std::for_each(execution::par)", [&](Image& img, Busy&) {
        std::vector<int> idx(n);
        std::iota(idx.begin(), idx.end(), 0);
        std::for_each(std::execution::par, idx.begin(), idx.end(), [&](int i) { pixel(img, i); });
    }, "TBB, busy time n/a");

    // 9-10. SIMD packets, serial and with OpenMP
    PacketStats<T> ps;
    add("SIMD packets (1 thread)", [&](Image& img, Busy& busy) {
        ps = {};
        const auto t0 = clk::now();
        for (int i = 0; i < n; i += W) trace_packet(sc, i, std::min(W, n - i), img, ps);
        busy.sec[0] = std::chrono::duration<double>(clk::now() - t0).count();
        for (int t = 1; t < nthreads; ++t) busy.sec[t] = busy.sec[0];
    });
    const double util = double(ps.lane_steps) / double(ps.vector_steps * W);
    results.back().note = "lanes=" + std::to_string(W) + ", lane utilisation " + std::to_string(int(util * 100 + 0.5)) + "%";
    add("SIMD packets + OpenMP dynamic", [&](Image& img, Busy& busy) {
#pragma omp parallel num_threads(nthreads)
        {
            PacketStats<T> local;
            const auto t0 = clk::now();
#pragma omp for schedule(dynamic, 2) nowait
            for (int i = 0; i < n; i += W) trace_packet(sc, i, std::min(W, n - i), img, local);
            busy.sec[omp_get_thread_num()] = std::chrono::duration<double>(clk::now() - t0).count();
        }
    });

    const double serial = results[0].sec;
    std::printf("\n%s, %dx%d rays, %d threads (hardware_concurrency=%u), best of %d\n", type_name<T>(), res, res,
                nthreads, std::thread::hardware_concurrency(), reps);
    std::printf("%-34s %9s %8s %10s %10s %8s  %s\n", "technique", "time [s]", "speedup", "efficiency", "imbalance",
                "mismatch", "note");
    for (const auto& r : results)
        std::printf("%-34s %9.3f %8.2f %9.0f%% %10.2f %8ld  %s\n", r.name.c_str(), r.sec, serial / r.sec,
                    100 * serial / r.sec / (r.name.rfind("SIMD packets (1", 0) == 0 || r.name == "serial" ? 1 : nthreads),
                    r.imbalance, r.mismatches, r.note.c_str());

    // cost distribution: how uneven is the work?
    std::vector<long> st(ref.steps.begin(), ref.steps.end());
    std::sort(st.begin(), st.end());
    std::printf("steps per ray: min %ld  median %ld  p99 %ld  max %ld\n", st.front(), st[st.size() / 2],
                st[size_t(0.99 * (st.size() - 1))], st.back());
}

// ---------------------------------------------------------------------------
// False sharing: every thread increments its own counter; the counters are
// either packed into one cache line or padded to 64 bytes each.
// ---------------------------------------------------------------------------
struct alignas(64) Padded {
    std::atomic<long> v{0};
};
void false_sharing(int nthreads) {
    const long iters = 20'000'000;
    auto bench = [&](auto& counters) {
        const auto t0 = clk::now();
        {
            std::vector<std::jthread> pool;
            for (int t = 0; t < nthreads; ++t)
                pool.emplace_back([&, t] {
                    for (long i = 0; i < iters; ++i) counters[t].fetch_add(1, std::memory_order_relaxed);
                });
        }
        return std::chrono::duration<double>(clk::now() - t0).count();
    };
    std::vector<std::atomic<long>> packed(nthreads);
    std::vector<Padded> padded(nthreads);
    const double tp = bench(packed);
    struct Wrap {
        std::vector<Padded>& p;
        std::atomic<long>& operator[](int i) { return p[i].v; }
    } w{padded};
    const double tq = bench(w);
    std::printf("\nFalse sharing, %d threads x %ld relaxed atomic increments on per-thread counters\n", nthreads, iters);
    std::printf("  packed  (adjacent, share one 64-byte line): %.3f s  (%.2f ns/inc)\n", tp, tp * 1e9 / iters);
    std::printf("  padded  (alignas(64), one line per counter): %.3f s  (%.2f ns/inc)   -> %.1fx faster\n", tq,
                tq * 1e9 / iters, tp / tq);
}

int main(int argc, char** argv) {
    int res = 256, reps = 3;
    int nthreads = int(std::thread::hardware_concurrency());
    std::string type = "double";
    bool fs = true;
    for (int i = 1; i + 1 < argc; i += 2) {
        if (!std::strcmp(argv[i], "--res")) res = std::atoi(argv[i + 1]);
        else if (!std::strcmp(argv[i], "--threads")) nthreads = std::atoi(argv[i + 1]);
        else if (!std::strcmp(argv[i], "--reps")) reps = std::atoi(argv[i + 1]);
        else if (!std::strcmp(argv[i], "--type")) type = argv[i + 1];
        else if (!std::strcmp(argv[i], "--false-sharing")) fs = std::atoi(argv[i + 1]);
    }
    if (type == "float") run<float>(res, nthreads, reps);
    else run<double>(res, nthreads, reps);
    if (fs) false_sharing(nthreads);
}
