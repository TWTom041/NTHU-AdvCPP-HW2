// critical_curve.hpp - Bardeen's analytic shadow boundary (binary128).
//   xi(r)  = [r^2 (3 - r) - a^2 (r + 1)] / [a (r - 1)]
//   eta(r) = r^3 [4 a^2 - r (r - 3)^2] / [a^2 (r - 1)^2]
//   alpha  = -xi / sin i,   beta = +-sqrt(eta + a^2 cos^2 i - xi^2 cot^2 i)
#pragma once

#include <cmath>
#include <quadmath.h>

namespace kerr {
using Q = __float128;

// Analytic critical curve, evaluated in binary128.  rho(psi) is found by
// bisection in the photon-orbit radius r on the branch with sign(beta) = sign(sin psi).
struct CriticalCurve {
    Q a, i, rp, rm;
    CriticalCurve(double spin, double incl_deg) : a(spin), i(Q(incl_deg) * M_PIq / 180) {
        rp = 2 * (1 + cosq(2 * acosq(-a) / 3));  // prograde photon orbit
        rm = 2 * (1 + cosq(2 * acosq(a) / 3));   // retrograde photon orbit
    }
    // point on the curve for orbit radius r; returns false where beta^2 < 0
    bool point(Q r, int sgn, Q& al, Q& be) const {
        const Q xi = (r * r * (3 - r) - a * a * (r + 1)) / (a * (r - 1));
        const Q eta = r * r * r * (4 * a * a - r * (r - 3) * (r - 3)) / (a * a * (r - 1) * (r - 1));
        const Q b2 = eta + a * a * cosq(i) * cosq(i) - xi * xi * cosq(i) * cosq(i) / (sinq(i) * sinq(i));
        if (b2 < 0) return false;
        al = -xi / sinq(i);
        be = sgn * sqrtq(b2);
        return true;
    }
    double at(double psi) const {
        const int sgn = std::sin(psi) >= 0 ? 1 : -1;
        const Q target = psi > M_PI ? Q(psi) - 2 * M_PIq : Q(psi);  // atan2 range
        auto f = [&](Q r, bool& ok) {
            Q al, be;
            ok = point(r, sgn, al, be);
            return ok ? remainderq(atan2q(be, al) - target, 2 * M_PIq) : Q(0);
        };
        // bracket on a coarse grid, then bisect
        const int n = 20000;
        Q r0 = 0, f0 = 0;
        bool have = false;
        for (int k = 0; k <= n; ++k) {
            const Q r = rp + (rm - rp) * Q(k) / n;
            bool ok;
            const Q fr = f(r, ok);
            if (!ok) { have = false; continue; }
            if (have && f0 * fr <= 0 && fabsq(f0 - fr) < 1) {
                Q lo = r0, hi = r, flo = f0;
                for (int it = 0; it < 200; ++it) {
                    const Q mid = (lo + hi) / 2;
                    const Q fm = f(mid, ok);
                    if ((fm <= 0) == (flo <= 0)) lo = mid, flo = fm; else hi = mid;
                }
                Q al, be;
                point((lo + hi) / 2, sgn, al, be);
                return double(hypotq(al, be));
            }
            r0 = r, f0 = fr, have = true;
        }
        return NAN;
    }
};

}  // namespace kerr
