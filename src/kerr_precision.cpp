// kerr_precision - Q1 benchmark: the same Kerr ray-tracing workload in every
// floating-point type.  Speed is measured single-threaded; accuracy is the
// round-off error relative to the identical algorithm run in __float128
// (same scheme, same step rule), so truncation error cancels out.
//
//   kerr_precision [--res 48] [--method rk4|dp45] [--reps 3]
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "kerr/tracer.hpp"

using namespace kerr;
using Q = __float128;

struct PixelOut {
    int status;
    Q intensity;
    State<Q> y;
    long steps;
};

struct Setup {
    int res = 48;
    Method method = Method::RK4;
    int reps = 3;
    bool exact_ktheta = false;
};

template <class T>
std::vector<PixelOut> run(const Setup& s, double& sec_best, long& steps_total, long& evals_total, double& maxH,
                          double& maxQ) {
    Spacetime<T> st{T(0.9375)};
    Camera<T> cam{T(1e4), T(60) * m::pi<T>() / T(180), T(40), s.res, s.res};
    cam.exact_ktheta = s.exact_ktheta;
    TraceConfig<T> cfg;
    cfg.method = s.method;
    cfg.step_frac = T(0.03);
    // tolerance that every type can honour, so that all types take ~ the same steps
    cfg.rtol = T(1e-6);
    cfg.atol = T(1e-8);
    cfg.disk.enabled = true;
    cfg.disk.r_in = st.r_isco();
    cfg.disk.r_out = T(20);
    cfg.track_invariants = true;
    const int n = s.res * s.res;
    std::vector<PixelOut> out(n);
    sec_best = 1e30;
    for (int rep = 0; rep < s.reps; ++rep) {
        steps_total = evals_total = 0;
        maxH = maxQ = 0;
        const auto t0 = std::chrono::steady_clock::now();
        for (int i = 0; i < n; ++i) {
            T al, be;
            cam.impact(i % s.res, i / s.res, al, be);
            const auto r = trace(st, cam.initial_state(st, al, be), cfg);
            steps_total += r.steps;
            evals_total += r.rhs_evals;
            maxH = std::max(maxH, double(r.max_abs_H));
            maxQ = std::max(maxQ, double(r.max_rel_dQ));
            State<Q> yq;
            for (int k = 0; k < 8; ++k) yq[k] = Q(r.y[k]);
            out[i] = {int(r.status), Q(r.intensity), yq, r.steps};
        }
        sec_best = std::min(sec_best, std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
    }
    return out;
}

static double angle(const State<Q>& u, const State<Q>& v) {
    const Q c = sinq(u[2]) * sinq(v[2]) * cosq(u[3] - v[3]) + cosq(u[2]) * cosq(v[2]);
    const Q s = hypotq(sinq(u[2] - v[2]), sinq(u[2]) * sinq(v[2]) * sinq(u[3] - v[3]));
    return double(atan2q(s, c));
}

template <class T> void report(const Setup& s, const std::vector<PixelOut>& ref) {
    double sec, maxH, maxQ;
    long steps, evals;
    const auto out = run<T>(s, sec, steps, evals, maxH, maxQ);
    const int n = int(out.size());
    long mismatch = 0;
    std::vector<double> dI, dAng;
    double flux = 0, flux_ref = 0;
    for (int i = 0; i < n; ++i) {
        flux += double(out[i].intensity), flux_ref += double(ref[i].intensity);
        if (out[i].status != ref[i].status) { ++mismatch; continue; }
        if (out[i].status == int(Status::Disk))
            dI.push_back(double(fabsq(out[i].intensity - ref[i].intensity) / ref[i].intensity));
        if (out[i].status == int(Status::Escaped) && out[i].steps == ref[i].steps) dAng.push_back(angle(out[i].y, ref[i].y));
    }
    auto q = [](std::vector<double> v, double p) {
        if (v.empty()) return 0.0;
        std::sort(v.begin(), v.end());
        return v[size_t(p * (v.size() - 1))];
    };
    std::printf("%-12s %8.1f %9.1f %8.2f %6ld | %9.2e %9.2e %9.2e | %9.2e %9.2e | %9.2e %9.2e\n", type_name<T>(),
                sec * 1e9 / double(steps), sec * 1e9 / double(evals), n / sec, mismatch, q(dI, 0.5), q(dI, 1.0),
                std::fabs(flux - flux_ref) / flux_ref, q(dAng, 0.5), q(dAng, 1.0), maxH, maxQ);
    std::fflush(stdout);
}

int main(int argc, char** argv) {
    Setup s;
    for (int i = 1; i < argc; ++i) {
        if (!std::strcmp(argv[i], "--res") && i + 1 < argc) s.res = std::atoi(argv[++i]);
        else if (!std::strcmp(argv[i], "--reps") && i + 1 < argc) s.reps = std::atoi(argv[++i]);
        else if (!std::strcmp(argv[i], "--method") && i + 1 < argc)
            s.method = std::string(argv[++i]) == "dp45" ? Method::DP45 : Method::RK4;
        else if (!std::strcmp(argv[i], "--exact-ktheta")) s.exact_ktheta = true;
    }
    std::printf("Kerr ray tracing, a=0.9375 i=60, %dx%d rays with thin disk, %s, single thread, best of %d\n", s.res,
                s.res, s.method == Method::DP45 ? "DP45 rtol=1e-6" : "RK4 STEPSIZE=0.03", s.reps);
    std::printf("reference = identical algorithm in __float128; errors are therefore pure round-off%s\n",
                s.exact_ktheta ? " (k_theta = beta)" : " (RAPTOR k_theta formula)");
    std::printf("%-12s %8s %9s %8s %6s | %9s %9s %9s | %9s %9s | %9s %9s\n", "type", "ns/step", "ns/rhs", "rays/s",
                "status", "dI med", "dI max", "d flux", "dir med", "dir max", "max|2H|", "max dQ");
    std::printf("%-12s %8s %9s %8s %6s | %9s %9s %9s | %9s %9s | %9s %9s\n", "", "", "", "", "diff", "(rel)", "(rel)",
                "(rel)", "[rad]", "[rad]", "", "(rel)");
    double sec, maxH, maxQ;
    long steps, evals;
    Setup sq = s;
    sq.reps = 1;
    const auto ref = run<Q>(sq, sec, steps, evals, maxH, maxQ);
    report<float>(s, ref);
    report<double>(s, ref);
    report<long double>(s, ref);
    report<Q>(sq, ref);
#ifdef __STDCPP_FLOAT16_T__
    report<std::float16_t>(s, ref);
#endif
}
