// cache_sweep - how the memory hierarchy changes the float-vs-double picture.
//
// (1) streaming sum over an array of W bytes (W = 4 KiB ... 1 GiB), float vs
//     double: inside L1/L2 the loop is compute/SIMD bound, beyond the LLC it is
//     bandwidth bound - there float wins simply because it moves half the bytes.
// (2) dependent random walk (pointer chasing) over the same working sets: the
//     load-to-use latency of L1 / L2 / L3 / DRAM.
//
//   cache_sweep [max_MiB=1024]
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <numeric>
#include <random>
#include <vector>

using clk = std::chrono::steady_clock;
inline void clobber() { asm volatile("" ::: "memory"); }

template <class T> double stream_gbs(size_t bytes, double& gelem) {
    const size_t n = bytes / sizeof(T);
    std::vector<T> v(n, T(1));
    const size_t total = std::max<size_t>(size_t(1) << 31, bytes * 4);  // bytes to stream per measurement
    const size_t reps = total / bytes;
    double best = 1e30;
    T s = 0;
    for (int k = 0; k < 3; ++k) {
        const auto t0 = clk::now();
        for (size_t r = 0; r < reps; ++r) {
            // L independent partial sums (= 4 x 512-bit registers of T) so that
            // neither the add latency nor missing vectorisation is the bottleneck
            constexpr size_t L = 4 * 64 / sizeof(T);
            T acc[L] = {};
            const T* p = v.data();
            for (size_t i = 0; i + L <= n; i += L)
                for (size_t j = 0; j < L; ++j) acc[j] += p[i + j];
            for (size_t j = 0; j < L; ++j) s += acc[j];
            clobber();
        }
        best = std::min(best, std::chrono::duration<double>(clk::now() - t0).count());
    }
    if (s == T(-1)) std::puts("");
    gelem = double(n) * reps / best * 1e-9;
    return double(bytes) * reps / best * 1e-9;
}

double chase_ns(size_t bytes) {
    // one pointer per 64-byte line, random cyclic permutation
    const size_t lines = bytes / 64;
    struct alignas(64) Line { Line* next; char pad[56]; };
    std::vector<Line> mem(lines);
    std::vector<size_t> perm(lines);
    std::iota(perm.begin(), perm.end(), 0);
    std::shuffle(perm.begin() + 1, perm.end(), std::mt19937_64(42));
    for (size_t i = 0; i < lines; ++i) mem[perm[i]].next = &mem[perm[(i + 1) % lines]];
    const size_t hops = 20'000'000;
    Line* p = &mem[perm[0]];
    const auto t0 = clk::now();
    for (size_t i = 0; i < hops; ++i) p = p->next;
    const double s = std::chrono::duration<double>(clk::now() - t0).count();
    if (p == nullptr) std::puts("");
    return s * 1e9 / hops;
}

int main(int argc, char** argv) {
    const size_t max_mib = argc > 1 ? std::strtoul(argv[1], nullptr, 10) : 1024;
    std::printf("%10s | %10s %10s %10s %10s | %12s\n", "working set", "float GB/s", "double GB/s", "float Gel/s",
                "double Gel/s", "chase ns/load");
    for (size_t kib = 4; kib <= max_mib * 1024; kib *= 2) {
        const size_t bytes = kib * 1024;
        double gf, gd;
        const double bf = stream_gbs<float>(bytes, gf);
        const double bd = stream_gbs<double>(bytes, gd);
        const double lat = chase_ns(bytes);
        if (kib < 1024) std::printf("%7zu KiB | %10.1f %10.1f %10.2f %10.2f | %12.2f\n", kib, bf, bd, gf, gd, lat);
        else std::printf("%7zu MiB | %10.1f %10.1f %10.2f %10.2f | %12.2f\n", kib / 1024, bf, bd, gf, gd, lat);
        std::fflush(stdout);
    }
}
