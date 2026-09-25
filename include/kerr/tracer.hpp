// tracer.hpp - integrators and the per-ray driver.
#pragma once

#include "kerr.hpp"

namespace kerr {

enum class Method : std::uint8_t {
    RK4,   // classic 4th-order Runge-Kutta with RAPTOR's step-size heuristic
    RK2,   // RAPTOR's default midpoint scheme with the same heuristic
    DP45,  // Dormand-Prince 5(4) with embedded error control
};

template <class T>
struct TraceConfig {
    Method method = Method::RK4;
    T step_frac = T(0.03);        // RAPTOR "STEPSIZE"
    T rtol = T(1e-10), atol = T(1e-12);  // DP45 only
    T r_out = T(1.01e4);          // escape radius (RAPTOR: 1.01 * rcam)
    T horizon_margin = T(0.1);    // stop at r_+ (1 + margin), RAPTOR uses 0.1
    long max_steps = 100000;      // RAPTOR: 1e5
    ThinDisk<T> disk{};
    bool track_invariants = false;
};

template <class T>
struct RayResult {
    Status status = Status::Running;
    State<T> y{};
    long steps = 0;
    long rhs_evals = 0;
    T intensity = T(0);
    T max_abs_H = T(0);       // max |2H| along the ray
    T max_rel_dQ = T(0);      // max |Q - Q0| / max(1, |Q0|)
};

namespace detail {

template <class T> inline void axpy(State<T>& out, const State<T>& y, T h, const State<T>& k) {
    for (int i = 0; i < 8; ++i) out[i] = y[i] + h * k[i];
}

// RAPTOR stepsize(): dl = f / (|r'|/r + |th'|/min(th, pi-th) + |phi'|), bounded below by 1e-12.
template <class T> inline T raptor_step(const State<T>& y, const State<T>& dy, T frac) {
    const T tiny = T(1e-40f);
    const T dist_pole = m::min(y[2], m::pi<T>() - y[2]);
    const T inv = m::fabs(dy[1]) / (y[1] * frac) + m::fabs(dy[2]) / (dist_pole * frac + tiny) +
                  m::fabs(dy[3]) / frac;
    return m::max(T(1) / (inv + tiny), T(1e-12));
}

template <class T> inline void rk4(const Spacetime<T>& st, State<T>& y, const State<T>& k1, T h) {
    State<T> k2, k3, k4, tmp;
    axpy(tmp, y, h * T(0.5), k1);
    st.rhs(tmp, k2);
    axpy(tmp, y, h * T(0.5), k2);
    st.rhs(tmp, k3);
    axpy(tmp, y, h, k3);
    st.rhs(tmp, k4);
    const T h6 = h / T(6);
    for (int i = 0; i < 8; ++i) y[i] += h6 * (k1[i] + T(2) * k2[i] + T(2) * k3[i] + k4[i]);
}

template <class T> inline void rk2(const Spacetime<T>& st, State<T>& y, const State<T>& k1, T h) {
    State<T> k2, tmp;
    axpy(tmp, y, h * T(0.5), k1);
    st.rhs(tmp, k2);
    for (int i = 0; i < 8; ++i) y[i] += h * k2[i];
}

// Dormand-Prince 5(4) coefficients
template <class T> struct DP {
    static constexpr double c2 = 1. / 5, c3 = 3. / 10, c4 = 4. / 5, c5 = 8. / 9;
    static constexpr double a21 = 1. / 5;
    static constexpr double a31 = 3. / 40, a32 = 9. / 40;
    static constexpr double a41 = 44. / 45, a42 = -56. / 15, a43 = 32. / 9;
    static constexpr double a51 = 19372. / 6561, a52 = -25360. / 2187, a53 = 64448. / 6561, a54 = -212. / 729;
    static constexpr double a61 = 9017. / 3168, a62 = -355. / 33, a63 = 46732. / 5247, a64 = 49. / 176,
                            a65 = -5103. / 18656;
    static constexpr double b1 = 35. / 384, b3 = 500. / 1113, b4 = 125. / 192, b5 = -2187. / 6784, b6 = 11. / 84;
    static constexpr double e1 = 71. / 57600, e3 = -71. / 16695, e4 = 71. / 1920, e5 = -17253. / 339200,
                            e6 = 22. / 525, e7 = -1. / 40;
};

}  // namespace detail

// Trace a single ray until it is captured, escapes, hits the disk or runs out of steps.
template <class T>
RayResult<T> trace(const Spacetime<T>& st, State<T> y, const TraceConfig<T>& cfg) {
    RayResult<T> res;
    const T r_stop = st.r_plus * (T(1) + cfg.horizon_margin);
    const T Q0 = cfg.track_invariants ? st.carter(y) : T(0);
    State<T> k1;
    st.rhs(y, k1);
    res.rhs_evals = 1;
    T h = T(0);  // DP45 current step
    if (cfg.method == Method::DP45) h = detail::raptor_step(y, k1, cfg.step_frac);

    while (res.steps < cfg.max_steps) {
        const T r_old = y[1];
        const T c_old = m::cos(y[2]);

        if (cfg.method == Method::DP45) {
            using D = detail::DP<T>;
            State<T> k2, k3, k4, k5, k6, k7, tmp, y5;
            for (;;) {
                for (int i = 0; i < 8; ++i) tmp[i] = y[i] + h * (T(D::a21) * k1[i]);
                st.rhs(tmp, k2);
                for (int i = 0; i < 8; ++i) tmp[i] = y[i] + h * (T(D::a31) * k1[i] + T(D::a32) * k2[i]);
                st.rhs(tmp, k3);
                for (int i = 0; i < 8; ++i)
                    tmp[i] = y[i] + h * (T(D::a41) * k1[i] + T(D::a42) * k2[i] + T(D::a43) * k3[i]);
                st.rhs(tmp, k4);
                for (int i = 0; i < 8; ++i)
                    tmp[i] = y[i] + h * (T(D::a51) * k1[i] + T(D::a52) * k2[i] + T(D::a53) * k3[i] +
                                         T(D::a54) * k4[i]);
                st.rhs(tmp, k5);
                for (int i = 0; i < 8; ++i)
                    tmp[i] = y[i] + h * (T(D::a61) * k1[i] + T(D::a62) * k2[i] + T(D::a63) * k3[i] +
                                         T(D::a64) * k4[i] + T(D::a65) * k5[i]);
                st.rhs(tmp, k6);
                for (int i = 0; i < 8; ++i)
                    y5[i] = y[i] + h * (T(D::b1) * k1[i] + T(D::b3) * k3[i] + T(D::b4) * k4[i] +
                                        T(D::b5) * k5[i] + T(D::b6) * k6[i]);
                st.rhs(y5, k7);
                res.rhs_evals += 6;
                // error norm over r, theta, phi, p_r, p_theta (t, p_t, p_phi excluded)
                T err = T(0);
                for (int i : {1, 2, 3, 5, 6}) {
                    const T e = h * (T(D::e1) * k1[i] + T(D::e3) * k3[i] + T(D::e4) * k4[i] +
                                     T(D::e5) * k5[i] + T(D::e6) * k6[i] + T(D::e7) * k7[i]);
                    const T sc = cfg.atol + cfg.rtol * m::max(m::fabs(y[i]), m::fabs(y5[i]));
                    err = m::max(err, m::fabs(e) / sc);
                }
                if (!m::isfinite(err)) err = T(1e10);
                if (err <= T(1)) {
                    y = y5;
                    k1 = k7;  // FSAL
                    T fac = err > T(0) ? T(0.9) * T(std::pow(double(err), -0.2)) : T(5);
                    h *= m::min(T(5), m::max(T(0.2), fac));
                    // never step further than the heuristic allows (keeps the disk crossing resolvable)
                    h = m::min(h, T(20) * detail::raptor_step(y, k1, cfg.step_frac));
                    break;
                }
                h *= m::max(T(0.1), T(0.9) * T(std::pow(double(err), -0.25)));
                if (h < T(1e-14)) break;
            }
        } else {
            const T dl = detail::raptor_step(y, k1, cfg.step_frac);
            if (cfg.method == Method::RK4) {
                detail::rk4(st, y, k1, dl);
                res.rhs_evals += 3;
            } else {
                detail::rk2(st, y, k1, dl);
                res.rhs_evals += 1;
            }
            st.rhs(y, k1);
            res.rhs_evals += 1;
        }
        ++res.steps;

        if (cfg.track_invariants) {
            res.max_abs_H = m::max(res.max_abs_H, m::fabs(st.two_H(y)));
            const T dQ = m::fabs(st.carter(y) - Q0) / m::max(T(1), m::fabs(Q0));
            res.max_rel_dQ = m::max(res.max_rel_dQ, dQ);
        }

        const T r = y[1];
        if (!m::isfinite(r) || !m::isfinite(y[2])) { res.status = Status::NotFinite; break; }

        if (cfg.disk.enabled) {
            const T c_new = m::cos(y[2]);
            if (c_old * c_new < T(0)) {
                const T f = c_old / (c_old - c_new);
                const T rc = r_old + f * (r - r_old);
                if (rc > cfg.disk.r_in && rc < cfg.disk.r_out) {
                    res.status = Status::Disk;
                    res.intensity = cfg.disk.intensity(st, rc, y[7]);
                    break;
                }
            }
        }
        if (r < r_stop) { res.status = Status::Captured; break; }
        if (r > cfg.r_out && k1[1] > T(0)) { res.status = Status::Escaped; break; }
    }
    if (res.status == Status::Running) res.status = Status::MaxSteps;
    res.y = y;
    return res;
}

}  // namespace kerr
