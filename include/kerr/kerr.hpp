// kerr.hpp - null geodesics in the Kerr spacetime (G = c = M = 1).
//
// Formulation
// -----------
// Coordinates : ingoing Kerr-Schild (t, r, theta, phi); regular on the future
//               horizon.  r and theta coincide with Boyer-Lindquist.
// Equations   : Hamiltonian form, H = 1/2 g^{ab}(x) p_a p_b = 0,
//                   dx^a/dl = g^{ab} p_b,   dp_a/dl = -dH/dx^a.
//               The inverse KS metric only depends on (r, theta), so p_t (=-E)
//               and p_phi (=L) are *exactly* constant - the integrator can only
//               drift in H and in the Carter constant Q, which we use as
//               accuracy diagnostics.  (RAPTOR instead integrates the
//               second-order geodesic equation with 64 Christoffel symbols in
//               modified KS coordinates.)
// Camera      : distant observer, impact parameters (alpha, beta) mapped to
//               (E, L, Q) exactly as RAPTOR / Cunningham & Bardeen (1973).
//               Rays are traced backwards by flipping the momentum p = -k.
#pragma once

#include <array>
#include <cstdint>

#include "real_math.hpp"

namespace kerr {

template <class T> using State = std::array<T, 8>;  // x^0..3 then p_0..3

enum class Status : std::uint8_t { Running, Captured, Escaped, Disk, MaxSteps, NotFinite };

inline const char* status_name(Status s) {
    switch (s) {
        case Status::Running: return "running";
        case Status::Captured: return "captured";
        case Status::Escaped: return "escaped";
        case Status::Disk: return "disk";
        case Status::MaxSteps: return "maxsteps";
        case Status::NotFinite: return "notfinite";
    }
    return "?";
}

template <class T>
struct Spacetime {
    T a;   // dimensionless spin, |a| < 1
    T a2;
    T r_plus, r_minus;  // outer / inner horizon

    explicit Spacetime(T spin) : a(spin), a2(spin * spin) {
        const T s = m::sqrt(T(1) - a2);
        r_plus = T(1) + s;
        r_minus = T(1) - s;
    }

    // Prograde innermost stable circular orbit (Bardeen, Press & Teukolsky 1972)
    T r_isco() const {
        using std::cbrt;
        const double ad = double(a);
        const double z1 = 1 + std::cbrt(1 - ad * ad) * (std::cbrt(1 + ad) + std::cbrt(1 - ad));
        const double z2 = std::sqrt(3 * ad * ad + z1 * z1);
        return T(3 + z2 - std::sqrt((3 - z1) * (3 + z1 + 2 * z2)));
    }

    // dy/dlambda for y = (x^a, p_a).  2H = -p_t^2 + M / Sigma with
    //   M = -2r p_t^2 + 4r p_t p_r + Delta p_r^2 + 2a p_r p_phi + p_th^2 + p_phi^2 / sin^2
    void rhs(const State<T>& y, State<T>& dy) const {
        const T r = y[1], th = y[2];
        const T pt = y[4], pr = y[5], pth = y[6], pph = y[7];
        const T s = m::sin(th), c = m::cos(th);
        const T s2 = s * s;
        const T r2 = r * r;
        const T sigma = r2 + a2 * c * c;
        const T isig = T(1) / sigma;
        const T delta = r2 - T(2) * r + a2;
        const T is2 = T(1) / s2;

        // x-dot = g^{ab} p_b
        dy[0] = -(T(1) + T(2) * r * isig) * pt + T(2) * r * isig * pr;
        dy[1] = isig * (T(2) * r * pt + delta * pr + a * pph);
        dy[2] = isig * pth;
        dy[3] = isig * (a * pr + pph * is2);

        const T M = -T(2) * r * pt * pt + T(4) * r * pt * pr + delta * pr * pr + T(2) * a * pr * pph +
                    pth * pth + pph * pph * is2;
        const T dM_dr = -T(2) * pt * pt + T(4) * pt * pr + (T(2) * r - T(2)) * pr * pr;
        const T dM_dth = -T(2) * pph * pph * c * is2 / s;
        const T dsig_dr = T(2) * r;
        const T dsig_dth = -T(2) * a2 * s * c;

        dy[4] = T(0);
        dy[5] = -T(0.5f) * (dM_dr - M * dsig_dr * isig) * isig;
        dy[6] = -T(0.5f) * (dM_dth - M * dsig_dth * isig) * isig;
        dy[7] = T(0);
    }

    // 2H = g^{ab} p_a p_b   (exactly 0 for a null geodesic)
    T two_H(const State<T>& y) const {
        const T r = y[1], th = y[2];
        const T pt = y[4], pr = y[5], pth = y[6], pph = y[7];
        const T s = m::sin(th), c = m::cos(th);
        const T sigma = r * r + a2 * c * c;
        const T delta = r * r - T(2) * r + a2;
        const T M = -T(2) * r * pt * pt + T(4) * r * pt * pr + delta * pr * pr + T(2) * a * pr * pph +
                    pth * pth + pph * pph / (s * s);
        return -pt * pt + M / sigma;
    }

    // Lower an index with the covariant KS metric: k_a = g_ab k^b
    void lower(T r, T th, const T ku[4], T kd[4]) const {
        const T s = m::sin(th), c = m::cos(th), s2 = s * s;
        const T sigma = r * r + a2 * c * c;
        const T z = T(2) * r / sigma;
        const T gtt = -(T(1) - z), gtr = z, gtp = -a * z * s2, grr = T(1) + z, grp = -a * s2 * (T(1) + z),
                gthth = sigma, gpp = s2 * (sigma + a2 * s2 * (T(1) + z));
        kd[0] = gtt * ku[0] + gtr * ku[1] + gtp * ku[3];
        kd[1] = gtr * ku[0] + grr * ku[1] + grp * ku[3];
        kd[2] = gthth * ku[2];
        kd[3] = gtp * ku[0] + grp * ku[1] + gpp * ku[3];
    }

    // Carter constant Q = p_th^2 + cos^2(th) (-a^2 p_t^2 + p_phi^2 / sin^2(th))
    T carter(const State<T>& y) const {
        const T s = m::sin(y[2]), c = m::cos(y[2]);
        return y[6] * y[6] + c * c * (-a2 * y[4] * y[4] + y[7] * y[7] / (s * s));
    }
};

// Distant camera at (r_cam, inclination) with RAPTOR's impact-parameter mapping.
template <class T>
struct Camera {
    T r_cam;         // camera radius [M]
    T incl;          // inclination [rad]
    T fov;           // image covers alpha, beta in [-fov/2, fov/2]  [M]
    int width, height;
    // RAPTOR computes k_theta = sgn(beta) sqrt(|Q - L^2 cot^2 i + cos^2 i|), which is
    // algebraically just beta but cancels catastrophically when |alpha| >> |beta|.
    bool exact_ktheta = false;

    // pixel centre -> impact parameters (identical to RAPTOR LINEAR_IMPACT_CAM)
    void impact(int x, int y, T& alpha, T& beta) const {
        const T stepx = fov / T(width), stepy = fov / T(height);
        alpha = -fov * T(0.5) + (T(x) + T(0.5)) * stepx;
        beta = -fov * T(0.5) + (T(y) + T(0.5)) * stepy;
    }

    // Initial state (x^a, p_a) in KS coordinates for impact parameters (alpha, beta).
    // p is the *time-reversed* photon momentum, so integrating forward in lambda
    // traces the ray from the camera back towards the black hole.
    State<T> initial_state(const Spacetime<T>& st, T alpha, T beta) const {
        const T a = st.a, a2 = st.a2;
        const T mu0 = m::cos(incl);
        const T th = m::acos(mu0);
        const T s = m::sin(th), c = m::cos(th);
        const T s2 = s * s, c2 = c * c;
        const T r = r_cam;

        // Constants of motion (E = 1)
        const T kt = T(-1);
        const T L = -alpha * m::sqrt(T(1) - mu0 * mu0);
        const T Q = beta * beta + mu0 * mu0 * (alpha * alpha - T(1));
        const T sgn = beta < T(0) ? T(-1) : (beta > T(0) ? T(1) : T(0));
        const T kth = exact_ktheta ? beta : sgn * m::sqrt(m::fabs(Q - L * L * c2 / s2 + c2));

        // Boyer-Lindquist inverse metric at the camera, solve H = 0 for k_r
        const T sigma = r * r + a2 * c2;
        const T delta = r * r - T(2) * r + a2;
        const T A = (r * r + a2) * (r * r + a2) - delta * a2 * s2;
        const T gtt = -A / (sigma * delta);
        const T gtp = -T(2) * a * r / (sigma * delta);
        const T gpp = (delta - a2 * s2) / (sigma * delta * s2);
        const T gthth = T(1) / sigma;
        const T grr = delta / sigma;
        const T rest = gtt * kt * kt + T(2) * gtp * kt * L + gpp * L * L + gthth * kth * kth;
        const T kr_bl = m::sqrt(-rest / grr);  // outgoing photon at the camera

        // BL -> KS: covariant k_r picks up the dt/dr and dphi/dr terms
        const T kr_ks = kr_bl - T(2) * r / delta * kt - a / delta * L;

        // Same KS phi offset as RAPTOR's BL_to_KS_u
        const T phi0 = m::pi<T>() / T(2) +
                       a / (st.r_plus - st.r_minus) * m::log((r - st.r_plus) / (r - st.r_minus));

        return State<T>{T(0), r, th, phi0, -kt, -kr_ks, -kth, -L};
    }
};

// Geometrically thin, optically thick Keplerian disk in the equatorial plane.
template <class T>
struct ThinDisk {
    bool enabled = false;
    T r_in = T(6), r_out = T(20);

    // Observed bolometric intensity for a ray with p_phi (reversed momentum) that
    // hits the disk at radius r: I = g^4 * r^-3 (1 - sqrt(r_in/r)).
    T intensity(const Spacetime<T>& st, T r, T p_phi) const {
        const T a = st.a;
        const T L = -p_phi;                              // physical k_phi, E = 1
        const T omega = T(1) / (r * m::sqrt(r) + a);     // prograde Keplerian
        const T gtt = -(T(1) - T(2) / r);
        const T gtp = -T(2) * a / r;
        const T gpp = r * r + st.a2 + T(2) * st.a2 / r;
        const T ut = T(1) / m::sqrt(-(gtt + T(2) * gtp * omega + gpp * omega * omega));
        const T g = T(1) / (ut * (T(1) - omega * L));    // redshift factor nu_obs/nu_em
        const T emiss = (T(1) - m::sqrt(r_in / r)) / (r * r * r);
        const T g2 = g * g;
        return g2 * g2 * emiss;
    }
};

}  // namespace kerr
