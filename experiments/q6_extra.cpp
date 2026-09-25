// q6_extra - three floating-point pitfalls not covered by Q1-Q5, all taken
// from the Kerr ray tracer.
//
//  A. Conditioning: near the shadow edge the map (alpha, beta) -> final ray
//     direction is chaotic.  A 1e-15 relative change of the input is amplified
//     by ~1/delta (delta = distance to the critical curve).  No floating type
//     can fix an ill-conditioned problem - it only moves the wall.
//  B. Catastrophic cancellation: RAPTOR's k_theta formula is algebraically
//     equal to beta but loses digits when |alpha| >> |beta|.
//  C. Subnormal numbers: correct, but they can be very slow.
#include <chrono>
#include <cmath>
#include <cstdio>

#include "kerr/critical_curve.hpp"
#include "kerr/tracer.hpp"

using namespace kerr;

static double angle(const State<double>& u, const State<double>& v) {
    const double c = std::sin(u[2]) * std::sin(v[2]) * std::cos(u[3] - v[3]) + std::cos(u[2]) * std::cos(v[2]);
    const double s = std::hypot(std::sin(u[2] - v[2]), std::sin(u[2]) * std::sin(v[2]) * std::sin(u[3] - v[3]));
    return std::atan2(s, c);
}

void part_a() {
    const double spin = 0.9375, incl = 60, psi = 0.3;  // off-axis (alpha = 0 would send the ray through the coordinate pole)
    CriticalCurve cc(spin, incl);
    const double rc = cc.at(psi);
    Spacetime<double> st{spin};
    Camera<double> cam{1e4, incl * M_PI / 180, 40, 1, 1};
    TraceConfig<double> cfg;
    cfg.method = Method::DP45;
    cfg.rtol = 1e-12, cfg.atol = 1e-14;
    cfg.horizon_margin = 0.01;
    cfg.max_steps = 1000000;
    std::printf("A. Conditioning near the photon ring (a=%g, i=%g, psi=0.3 rad, rho_c=%.12f)\n", spin, incl, rc);
    std::printf("   rays at rho = rho_c (1 + delta) and rho (1 + 1e-13); angle between final directions\n");
    std::printf("   %-10s %-10s %-14s %-14s\n", "delta", "steps", "d(direction)", "amplification");
    const double eps = 1e-13;
    for (double delta : {1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10}) {
        const double rho = rc * (1 + delta);  // just outside the shadow: the ray escapes
        auto r1 = trace(st, cam.initial_state(st, rho * std::cos(psi), rho * std::sin(psi)), cfg);
        auto r2 = trace(st, cam.initial_state(st, rho * (1 + eps) * std::cos(psi), rho * (1 + eps) * std::sin(psi)), cfg);
        const double d = angle(r1.y, r2.y);
        std::printf("   %-10.0e %-10ld %-14.3e %-14.3e %s\n", delta, r1.steps, d, d / eps,
                    r1.status == Status::Escaped && r2.status == Status::Escaped ? "" : "(status differs)");
    }
}

void part_b() {
    std::printf("\nB. Catastrophic cancellation in the camera set-up (float, i = 60 deg)\n");
    std::printf("   k_theta = sqrt(Q - L^2 cot^2 i + cos^2 i) with Q = beta^2 + cos^2 i (alpha^2 - 1), L = -alpha sin i\n");
    std::printf("   %-8s %-8s %-16s %-16s %s\n", "alpha", "beta", "RAPTOR formula", "exact (=beta)", "rel. error");
    for (auto [al, be] : {std::pair{0.5f, 3.0f}, {5.0f, 3.0f}, {20.0f, 3.0f}, {20.0f, 0.3f}, {20.0f, 0.03f}, {20.0f, 0.003f}}) {
        Spacetime<float> st{0.9375f};
        Camera<float> cam{1e4f, 60.0f * m::pi<float>() / 180, 40, 1, 1};
        const float k = -cam.initial_state(st, al, be)[6];  // reversed momentum
        std::printf("   %-8g %-8g %-16.9g %-16.9g %.2e\n", al, be, k, be, std::fabs(k - be) / be);
    }
}

void part_c() {
    std::printf("\nC. Subnormal operands (x *= 0.999999 in a dependent chain, 1e7 iterations)\n");
    for (double start : {1.0, 1e-300, 1e-310}) {
        volatile double v = start;
        double x = v;
        const double f = 0.999999;
        const auto t0 = std::chrono::steady_clock::now();
        for (int i = 0; i < 10'000'000; ++i) x = x * f + (x == 0 ? start : 0.0);
        const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        v = x;
        std::printf("   start %-8g (%s): %.2f ns per multiply\n", start,
                    std::fpclassify(start) == FP_SUBNORMAL ? "subnormal" : "normal   ", sec * 1e9 / 1e7);
    }
}

int main() {
    part_a();
    part_b();
    part_c();
}
