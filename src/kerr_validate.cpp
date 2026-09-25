// kerr_validate - validate the engine against the analytic Kerr shadow.
//
// The edge of the black-hole shadow seen by a distant observer is Bardeen's
// critical curve (spherical photon orbits):
//   xi(r)  = [r^2 (3 - r) - a^2 (r + 1)] / [a (r - 1)]
//   eta(r) = r^3 [4 a^2 - r (r - 3)^2] / [a^2 (r - 1)^2]
//   alpha  = -xi / sin i,   beta = +-sqrt(eta + a^2 cos^2 i - xi^2 cot^2 i)
// For N polar angles psi we locate the capture boundary rho(psi) by bisection
// with the geodesic integrator and compare with the analytic curve.
//
//   kerr_validate [--spin 0.9375] [--incl 60] [--angles 24]
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "kerr/critical_curve.hpp"
#include "kerr/tracer.hpp"

using namespace kerr;

template <class T>
bool captured(const Spacetime<T>& st, const Camera<T>& cam, const TraceConfig<T>& cfg, double rho, double psi,
              long& steps) {
    const T alpha = T(rho * std::cos(psi)), beta = T(rho * std::sin(psi));
    auto res = trace(st, cam.initial_state(st, alpha, beta), cfg);
    steps += res.steps;
    return res.status == Status::Captured;
}

template <class T>
void run(double spin, double incl, int nang, double rtol, const CriticalCurve& cc, Method method, double step) {
    Spacetime<T> st{T(spin)};
    Camera<T> cam{T(1e4), T(incl) * m::pi<T>() / T(180), T(40), 1, 1};
    TraceConfig<T> cfg;
    cfg.method = method;
    cfg.step_frac = T(step);
    cfg.rtol = T(rtol);
    cfg.atol = T(rtol * 1e-2);
    cfg.horizon_margin = T(0.01);  // must stay below the prograde photon orbit (r ~ 1.43 for a = 0.9375)
    cfg.max_steps = 2000000;

    double max_err = 0, sum_err = 0;
    long steps = 0;
    const auto t0 = std::chrono::steady_clock::now();
    for (int k = 0; k < nang; ++k) {
        const double psi = 2 * M_PI * (k + 0.5) / nang;
        const double ref = cc.at(psi);
        double lo = ref * 0.9, hi = ref * 1.1;  // captured at lo, escapes at hi
        if (!captured(st, cam, cfg, lo, psi, steps) || captured(st, cam, cfg, hi, psi, steps)) {
            std::printf("  bracket failed at psi=%g\n", psi);
            continue;
        }
        for (int it = 0; it < 64 && hi - lo > 1e-15 * ref; ++it) {
            const double mid = 0.5 * (lo + hi);
            (captured(st, cam, cfg, mid, psi, steps) ? lo : hi) = mid;
        }
        const double err = std::fabs(0.5 * (lo + hi) - ref) / ref;
        max_err = std::max(max_err, err);
        sum_err += err;
    }
    const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    const char* mname = method == Method::DP45 ? "dp45" : method == Method::RK4 ? "rk4" : "rk2";
    std::printf("%-12s %-5s %8.0e %6d | rel. edge error: mean %.3e  max %.3e | %8.2f s  %ld steps\n",
                type_name<T>(), mname, method == Method::DP45 ? rtol : step, nang, sum_err / nang, max_err,
                sec, steps);
    std::fflush(stdout);
}

int main(int argc, char** argv) {
    double spin = 0.9375, incl = 60;
    int nang = 24;
    for (int i = 1; i + 1 < argc; i += 2) {
        if (!std::strcmp(argv[i], "--spin")) spin = std::atof(argv[i + 1]);
        else if (!std::strcmp(argv[i], "--incl")) incl = std::atof(argv[i + 1]);
        else if (!std::strcmp(argv[i], "--angles")) nang = std::atoi(argv[i + 1]);
    }
    CriticalCurve cc(spin, incl);
    std::printf("Shadow edge vs Bardeen critical curve  (a=%g, i=%g deg, %d angles)\n", spin, incl, nang);
    std::printf("%-12s %-5s %8s %8s\n", "type", "meth", "tol/h", "angles");

    // Fixed-step RK4 with RAPTOR's heuristic: truncation error dominates
    for (double h : {0.03, 0.01, 0.003})
        run<double>(spin, incl, nang, 0, cc, Method::RK4, h);
    // Adaptive DP45 at the tightest tolerance each type supports
    run<float>(spin, incl, nang, 1e-6, cc, Method::DP45, 0.03);
    run<float>(spin, incl, nang, 1e-7, cc, Method::DP45, 0.03);
    run<double>(spin, incl, nang, 1e-9, cc, Method::DP45, 0.03);
    run<double>(spin, incl, nang, 1e-12, cc, Method::DP45, 0.03);
    run<long double>(spin, incl, nang, 1e-14, cc, Method::DP45, 0.03);
    run<__float128>(spin, incl, nang, 1e-15, cc, Method::DP45, 0.03);
}
