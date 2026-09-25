// raptor_compare - ray-by-ray comparison of this engine against RAPTOR.
//
// For every pixel of an N x N camera (RAPTOR's model.in: 40 x 40 M field of
// view, a = 0.9375, i = 60 deg) the same ray is traced by
//   (R) RAPTOR's initialize_photon + stepsize + rk4_step/rk2_step (double),
//   (E) this engine with the same step heuristic and scheme (double),
//   (X) a reference: this engine, long double, adaptive DP45 at rtol 1e-13.
// Escaping rays are then propagated to exactly r = 1.01e4 so that the final
// sky coordinates (theta, phi) can be compared, and the constants of motion
// (E, L, Carter Q, H) of the final states are evaluated with the same
// long-double routine for all three.
//
//   raptor_compare [--res 100] [--method rk4|rk2] [--step 0.03] [--margin 0.1] [--dump file]
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <algorithm>

#include "kerr/tracer.hpp"
#include "raptor_bridge.h"

using namespace kerr;
using LD = long double;

// Propagate (x, p) to exactly r = R in nearly flat space (RK4, |dr| <= 10 per step).
static State<LD> to_radius(const Spacetime<LD>& st, State<LD> y, LD R) {
    State<LD> k;
    for (int it = 0; it < 10000; ++it) {
        st.rhs(y, k);
        const LD dr = R - y[1];
        if (std::fabs(dr) < 1e-13L * R) break;
        LD dl = dr / k[1];
        if (std::fabs(dl * k[1]) > 10) dl = std::copysign(10.0L, dr) / k[1];
        detail::rk4(st, y, k, dl);
    }
    return y;
}

struct Invariants {
    LD E, L, Q, H;
};
static Invariants invariants(const Spacetime<LD>& st, const State<LD>& y) {
    return {y[4], y[7], st.carter(y), st.two_H(y)};
}

// RAPTOR final state (contravariant k, physical direction) -> engine state (covariant, reversed)
static State<LD> from_raptor(const Spacetime<LD>& st, const double o[8]) {
    const LD ku[4] = {o[4], o[5], o[6], o[7]};
    LD kd[4];
    st.lower(o[1], o[2], ku, kd);
    return {LD(o[0]), LD(o[1]), LD(o[2]), LD(o[3]), -kd[0], -kd[1], -kd[2], -kd[3]};
}

template <class T> static State<LD> widen(const State<T>& y) {
    State<LD> o;
    for (int i = 0; i < 8; ++i) o[i] = LD(y[i]);
    return o;
}

int main(int argc, char** argv) {
    int res = 100;
    std::string method = "rk4", dump;
    double step = 0.03, margin = 0.1, spin = 0.9375, incl = 60, fov = 40;
    for (int i = 1; i + 1 < argc; i += 2) {
        if (!std::strcmp(argv[i], "--res")) res = std::atoi(argv[i + 1]);
        else if (!std::strcmp(argv[i], "--method")) method = argv[i + 1];
        else if (!std::strcmp(argv[i], "--step")) step = std::atof(argv[i + 1]);
        else if (!std::strcmp(argv[i], "--margin")) margin = std::atof(argv[i + 1]);
        else if (!std::strcmp(argv[i], "--dump")) dump = argv[i + 1];
    }
    const bool rk4 = method == "rk4";
    raptor_setup(spin, incl, step, margin);

    Spacetime<double> st{spin};
    Camera<double> cam{1e4, incl * M_PI / 180, fov, res, res};
    TraceConfig<double> cfg;
    cfg.method = rk4 ? Method::RK4 : Method::RK2;
    cfg.step_frac = step;
    cfg.horizon_margin = margin;

    Spacetime<LD> stx{LD(spin)};
    Camera<LD> camx{1e4L, LD(incl) * m::pi<LD>() / 180, LD(fov), res, res};
    TraceConfig<LD> cfgx;
    cfgx.method = Method::DP45;
    cfgx.rtol = 1e-13L;
    cfgx.atol = 1e-15L;
    cfgx.horizon_margin = LD(margin);
    cfgx.max_steps = 1000000;

    const int n = res * res;
    std::vector<int> sR(n), sE(n), sX(n);
    std::vector<State<LD>> yR(n), yE(n), yX(n), y0(n);
    std::vector<long> nR(n), nE(n);
    double max_init_diff = 0;

    // --- RAPTOR (single thread, timed) ---
    auto t0 = std::chrono::steady_clock::now();
    for (int i = 0; i < n; ++i) {
        double al, be, o[8];
        cam.impact(i % res, i / res, al, be);
        sR[i] = raptor_trace(al, be, rk4 ? RAPTOR_RK4 : RAPTOR_RK2, o, &nR[i]);
        yR[i] = from_raptor(stx, o);
    }
    const double tR = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

    // --- this engine, same scheme (single thread, timed) ---
    t0 = std::chrono::steady_clock::now();
    for (int i = 0; i < n; ++i) {
        double al, be;
        cam.impact(i % res, i / res, al, be);
        auto r = trace(st, cam.initial_state(st, al, be), cfg);
        sE[i] = int(r.status);
        nE[i] = r.steps;
        yE[i] = widen(r.y);
    }
    const double tE = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

    // --- reference (parallel) + initial-condition check ---
#pragma omp parallel for schedule(dynamic, 4) reduction(max : max_init_diff)
    for (int i = 0; i < n; ++i) {
        LD al, be;
        camx.impact(i % res, i / res, al, be);
        y0[i] = camx.initial_state(stx, al, be);
        auto r = trace(stx, y0[i], cfgx);
        sX[i] = int(r.status);
        yX[i] = r.y;
        double o[8];
        raptor_initial(double(al), double(be), o);
        const State<LD> yr = from_raptor(stx, o);
        const State<double> ye = cam.initial_state(st, double(al), double(be));
        for (int k : {1, 2, 5, 6, 7})  // r, theta, p_r, p_theta, p_phi
            max_init_diff = std::max(max_init_diff, double(std::fabs(yr[k] - LD(ye[k])) / (1 + std::fabs(yr[k]))));
    }

    // --- statistics ---
    long agreeRE = 0, agreeRX = 0, agreeEX = 0, esc = 0, capR = 0, capE = 0, capX = 0;
    std::vector<double> dRX, dEX, dRE;  // angular distance on the sky between final directions
    double dE_R = 0, dL_R = 0, dQ_R = 0, H_R = 0, dE_E = 0, dL_E = 0, dQ_E = 0, H_E = 0;
    std::vector<float> map_dth(n, NAN);
    for (int i = 0; i < n; ++i) {
        agreeRE += sR[i] == sE[i];
        agreeRX += sR[i] == sX[i];
        agreeEX += sE[i] == sX[i];
        capR += sR[i] == RAPTOR_CAPTURED, capE += sE[i] == int(Status::Captured), capX += sX[i] == int(Status::Captured);
        const auto I0 = invariants(stx, y0[i]);
        const auto IR = invariants(stx, yR[i]), IE = invariants(stx, yE[i]);
        const LD q0 = std::max(1.0L, std::fabs(I0.Q));
        dE_R = std::max(dE_R, double(std::fabs(IR.E - I0.E))), dE_E = std::max(dE_E, double(std::fabs(IE.E - I0.E)));
        dL_R = std::max(dL_R, double(std::fabs(IR.L - I0.L))), dL_E = std::max(dL_E, double(std::fabs(IE.L - I0.L)));
        dQ_R = std::max(dQ_R, double(std::fabs(IR.Q - I0.Q) / q0)), dQ_E = std::max(dQ_E, double(std::fabs(IE.Q - I0.Q) / q0));
        H_R = std::max(H_R, double(std::fabs(IR.H))), H_E = std::max(H_E, double(std::fabs(IE.H)));
        if (sR[i] == 2 && sE[i] == int(Status::Escaped) && sX[i] == int(Status::Escaped)) {
            ++esc;
            const LD R = 1.01e4L;
            const auto a = to_radius(stx, yR[i], R), b = to_radius(stx, yE[i], R), c = to_radius(stx, yX[i], R);
            auto ang = [](const State<LD>& u, const State<LD>& v) {  // great-circle distance
                const LD c = std::sin(u[2]) * std::sin(v[2]) * std::cos(u[3] - v[3]) + std::cos(u[2]) * std::cos(v[2]);
                const LD s = std::hypot(std::sin(u[2] - v[2]), std::sin(u[2]) * std::sin(v[2]) * std::sin(u[3] - v[3]));
                return double(std::atan2(s, c));
            };
            dRX.push_back(ang(a, c)), dEX.push_back(ang(b, c)), dRE.push_back(ang(a, b));
            map_dth[i] = float(ang(a, c));
        }
    }
    long stepsR = 0, stepsE = 0;
    for (int i = 0; i < n; ++i) stepsR += nR[i], stepsE += nE[i];

    std::printf("RAPTOR vs this engine: a=%g i=%g fov=%g res=%d method=%s STEPSIZE=%g margin=%g\n", spin, incl, fov,
                res, method.c_str(), step, margin);
    std::printf("  max rel. difference of initial (r,th,p_r,p_th,p_phi): %.3e\n", max_init_diff);
    std::printf("  %-26s %12s %12s %12s\n", "", "RAPTOR", "engine", "reference");
    std::printf("  %-26s %12ld %12ld %12ld\n", "captured rays", capR, capE, capX);
    std::printf("  %-26s %12.3f %12.3f\n", "time [s] (1 thread)", tR, tE);
    std::printf("  %-26s %12.1f %12.1f\n", "mean steps / ray", double(stepsR) / n, double(stepsE) / n);
    std::printf("  %-26s %12.1f %12.1f\n", "ns / step", tR * 1e9 / stepsR, tE * 1e9 / stepsE);
    std::printf("  %-26s %12.3e %12.3e\n", "max |dE| (E=p_t)", dE_R, dE_E);
    std::printf("  %-26s %12.3e %12.3e\n", "max |dL| (L=p_phi)", dL_R, dL_E);
    std::printf("  %-26s %12.3e %12.3e\n", "max |dQ|/max(1,Q)", dQ_R, dQ_E);
    std::printf("  %-26s %12.3e %12.3e\n", "max |g^ab p_a p_b|", H_R, H_E);
    std::printf("  status agreement: RAPTOR-engine %.4f%%  RAPTOR-ref %.4f%%  engine-ref %.4f%%\n",
                100.0 * agreeRE / n, 100.0 * agreeRX / n, 100.0 * agreeEX / n);
    auto stats = [](const char* name, std::vector<double> v) {
        std::sort(v.begin(), v.end());
        auto q = [&](double p) { return v.empty() ? 0.0 : v[size_t(p * (v.size() - 1))]; };
        std::printf("    %-20s median %.3e   p99 %.3e   max %.3e\n", name, q(0.5), q(0.99), q(1.0));
    };
    std::printf("  escaping rays (%ld): angle between final directions at r = 1.01e4 [rad]\n", esc);
    stats("RAPTOR - reference", dRX);
    stats("engine - reference", dEX);
    stats("RAPTOR - engine", dRE);

    if (!dump.empty()) {
        FILE* f = std::fopen(dump.c_str(), "wb");
        std::fwrite(&res, sizeof res, 1, f);
        std::fwrite(sR.data(), sizeof(int), n, f);
        std::fwrite(sE.data(), sizeof(int), n, f);
        std::fwrite(sX.data(), sizeof(int), n, f);
        std::fwrite(map_dth.data(), sizeof(float), n, f);
        std::fclose(f);
    }
}
