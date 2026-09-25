// q2_sin - how accurate and how fast is sin(x), and why.
//
// Compares against MPFR (256-bit, correctly rounded reference):
//   glibc sin (double) / sinf (float) / sinl (x87 long double) / sinq (binary128)
//   the x87 hardware instruction FSIN
//   mini_sin: a 60-line fdlibm-style implementation (Cody-Waite reduction +
//             minimax polynomial) to show what the library has to get right.
// Inputs are split into the same ranges that glibc's dbl-64/s_sin.c branches on.
#include <mpfr.h>

#include <chrono>
#include <cfenv>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <quadmath.h>
#include <random>
#include <vector>

// ---------------------------------------------------------------------------
// mini_sin: fdlibm (Sun, 1993) algorithm, simplified
// ---------------------------------------------------------------------------
namespace mini {
constexpr double S1 = -1.66666666666666324348e-01, S2 = 8.33333333332248946124e-03,
                 S3 = -1.98412698298579493134e-04, S4 = 2.75573137070700676789e-06,
                 S5 = -2.50507602534068634195e-08, S6 = 1.58969099521155010221e-10;
constexpr double C1 = 4.16666666666666019037e-02, C2 = -1.38888888888741095749e-03,
                 C3 = 2.48015872894767294178e-05, C4 = -2.75573143513906633035e-07,
                 C5 = 2.08757232129817482790e-09, C6 = -1.13596475577881948265e-11;
// pi/2 split into 33-bit pieces so that n * piece is exact for |n| < 2^20
constexpr double invpio2 = 6.36619772367581382433e-01, pio2_1 = 1.57079632673412561417e+00,
                 pio2_2 = 6.07710050630396597660e-11, pio2_2t = 2.02226624879595063154e-21,
                 pio2_3 = 2.02226624871116645580e-21, pio2_3t = 8.47842766036889956997e-32;

inline uint32_t hi(double x) { uint64_t u; std::memcpy(&u, &x, 8); return uint32_t(u >> 32); }

// sin(x + y) on [-pi/4, pi/4], y = tail of the reduced argument
inline double ksin(double x, double y) {
    const double z = x * x, v = z * x;
    const double r = S2 + z * (S3 + z * (S4 + z * (S5 + z * S6)));
    return x - ((z * (0.5 * y - v * r) - y) - v * S1);
}
inline double kcos(double x, double y) {
    const double z = x * x;
    const double r = z * (C1 + z * (C2 + z * (C3 + z * (C4 + z * (C5 + z * C6)))));
    const uint32_t ix = hi(x) & 0x7fffffff;
    if (ix < 0x3FD33333) return 1.0 - (0.5 * z - (z * r - x * y));  // |x| < 0.3
    double qx;
    if (ix > 0x3fe90000) qx = 0.28125;
    else { const uint64_t u = uint64_t(ix - 0x00200000) << 32; std::memcpy(&qx, &u, 8); }  // x/4 truncated
    const double hz = 0.5 * z - qx, a = 1.0 - qx;
    return a - (hz - (z * r - x * y));
}
// Cody-Waite reduction x = n pi/2 + (y0 + y1), three pieces (~118 bits of pi/2)
inline int reduce(double x, double& y0, double& y1) {
    const double fn = std::nearbyint(x * invpio2);
    double t = x - fn * pio2_1;
    double w = fn * pio2_2;
    double r = t - w;
    w = fn * pio2_2t - ((t - r) - w);
    t = r;
    w = fn * pio2_3;
    r = t - w;
    w = fn * pio2_3t - ((t - r) - w);
    y0 = r - w;
    y1 = (r - y0) - w;
    return int(int64_t(fn) & 3);
}
double sin(double x) {
    if (std::fabs(x) < 0.7853981633974483) return ksin(x, 0.0);
    double y0, y1;
    switch (reduce(x, y0, y1)) {
        case 0: return ksin(y0, y1);
        case 1: return kcos(y0, y1);
        case 2: return -ksin(y0, y1);
        default: return -kcos(y0, y1);
    }
}
}  // namespace mini

// x87 FSIN instruction (valid for |x| < 2^63), 64-bit significand
inline long double x87_fsin(long double x) {
    long double r;
    asm("fsin" : "=t"(r) : "0"(x));
    return r;
}

// ---------------------------------------------------------------------------
// ULP error of y (a value of precision `prec` bits) w.r.t. the exact sin(x)
// ---------------------------------------------------------------------------
struct Ulp {
    mpfr_t ex, diff, xx;
    Ulp() { mpfr_inits2(320, ex, diff, xx, (mpfr_ptr)0); }
    ~Ulp() { mpfr_clears(ex, diff, xx, (mpfr_ptr)0); }
    void exact(double x) { mpfr_set_d(xx, x, MPFR_RNDN), mpfr_sin(ex, xx, MPFR_RNDN); }
    void exactf(float x) { exact(double(x)); }
    double err(int prec) {  // diff already holds y
        if (mpfr_zero_p(ex)) return mpfr_zero_p(diff) ? 0 : 1e30;
        const long e = mpfr_get_exp(ex);  // ex = m * 2^e, 0.5 <= |m| < 1
        mpfr_sub(diff, diff, ex, MPFR_RNDN);
        mpfr_mul_2si(diff, diff, -(e - prec), MPFR_RNDN);  // one ulp = 2^(e - prec)
        return std::fabs(mpfr_get_d(diff, MPFR_RNDN));
    }
    double of_double(double y, double x) { exact(x); mpfr_set_d(diff, y, MPFR_RNDN); return err(53); }
    double of_float(float y, float x) { exactf(x); mpfr_set_d(diff, double(y), MPFR_RNDN); return err(24); }
    double of_ld(long double y, double x) { exact(x); mpfr_set_ld(diff, y, MPFR_RNDN); return err(64); }
    double of_q(__float128 y, double x) {
        exact(x);
        char buf[64];
        quadmath_snprintf(buf, sizeof buf, "%.40Qa", y);
        mpfr_set_str(diff, buf, 0, MPFR_RNDN);
        return err(113);
    }
};

struct Range { const char* name; double lo, hi; bool log; };

struct Acc {
    double max = 0, sum = 0;
    long n = 0, correct = 0;
    void add(double u) { max = std::max(max, u); sum += u; ++n; correct += u <= 0.5; }
};

int main(int argc, char** argv) {
    const long samples = argc > 1 ? std::atol(argv[1]) : 200000;
    const Range ranges[] = {
        {"[2^-26, 0.855469)", 0x1p-26, 0.855469, true},
        {"[0.855469, 2.426265)", 0.855469, 2.426265, false},
        {"[2.426265, 105414350)", 2.426265, 105414350.0, true},
        {"[105414350, 1e300]", 105414350.0, 1e300, true},
    };
    std::mt19937_64 rng(2026);
    Ulp u;

    std::printf("ULP error vs MPFR, %ld random inputs per range  (max ulp / %% correctly rounded)\n", samples);
    std::printf("%-22s | %-17s | %-17s | %-17s | %-17s | %-17s | %-17s\n", "range (glibc branch)", "glibc sin",
                "mini_sin (fdlibm)", "x87 fsin", "glibc sinl", "glibc sinf", "libquadmath sinq");
    for (const auto& r : ranges) {
        Acc a_sin, a_mini, a_fsin, a_sinl, a_sinf, a_sinq;
        std::uniform_real_distribution<double> U(r.log ? std::log(r.lo) : r.lo, r.log ? std::log(r.hi) : r.hi);
        for (long i = 0; i < samples; ++i) {
            double x = U(rng);
            if (r.log) x = std::exp(x);
            a_sin.add(u.of_double(std::sin(x), x));
            a_mini.add(u.of_double(mini::sin(x), x));
            if (x < 0x1p63) a_fsin.add(u.of_ld(x87_fsin(x), x));
            a_sinl.add(u.of_ld(sinl((long double)x), x));
            if (i % 8 == 0) a_sinq.add(u.of_q(sinq(__float128(x)), x));
            const float xf = float(x);
            if (std::isfinite(xf)) a_sinf.add(u.of_float(sinf(xf), xf));
        }
        auto cell = [](const Acc& a) {
            static char buf[8][32];
            static int k = 0;
            char* b = buf[k++ & 7];
            if (a.n == 0) std::snprintf(b, 32, "%17s", "n/a");
            else if (a.max > 1e6) std::snprintf(b, 32, "%8.1e / %5.1f%%", a.max, 100.0 * a.correct / a.n);
            else std::snprintf(b, 32, "%8.3f / %5.1f%%", a.max, 100.0 * a.correct / a.n);
            return b;
        };
        std::printf("%-22s | %s | %s | %s | %s | %s | %s\n", r.name, cell(a_sin), cell(a_mini), cell(a_fsin),
                    cell(a_sinl), cell(a_sinf), cell(a_sinq));
    }

    // Hard cases for argument reduction
    std::printf("\nHard inputs (value printed with 17 significant digits, error in ulp)\n");
    const double hard[] = {M_PI, 2 * M_PI, 355.0, 103993.0, 1e22,
                           6381956970095103.0 * 0x1p797 /* Kahan-McDonald: closest double to a multiple of pi/2 */};
    std::printf("%-26s %-26s %10s %10s %10s\n", "x", "sin(x) glibc", "glibc", "mini_sin", "x87 fsin");
    for (double x : hard) {
        const double g = std::sin(x);
        const double eg = u.of_double(g, x), em = u.of_double(mini::sin(x), x);
        const double ef = x < 0x1p63 ? u.of_double(double(x87_fsin(x)), x) : NAN;
        std::printf("%-26.17g %-26.17g %10.3g %10.3g %10.3g\n", x, g, eg, em, ef);
    }

    // Timing per range
    std::printf("\nTime per call [ns] (glibc sin, 1M calls per range, independent inputs)\n");
    for (const auto& r : ranges) {
        std::vector<double> xs(1 << 20);
        std::uniform_real_distribution<double> U(r.log ? std::log(r.lo) : r.lo, r.log ? std::log(r.hi) : r.hi);
        for (auto& x : xs) x = r.log ? std::exp(U(rng)) : U(rng);
        auto time = [&](auto f) {
            double best = 1e30, s = 0;
            for (int k = 0; k < 3; ++k) {
                const auto t0 = std::chrono::steady_clock::now();
                for (double x : xs) s += f(x);
                best = std::min(best, std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
            }
            if (s == 12345.678) std::puts("");
            return best * 1e9 / double(xs.size());
        };
        std::printf("  %-22s glibc %6.1f   mini_sin %6.1f   x87 fsin %6.1f   sinl %6.1f\n", r.name,
                    time([](double x) { return std::sin(x); }), time([](double x) { return mini::sin(x); }),
                    time([](double x) { return double(x87_fsin(x)); }), time([](double x) { return double(sinl(x)); }));
    }

    // Special values (C Annex F / IEEE 754)
    std::printf("\nSpecial values\n");
    auto special = [](const char* name, double x) {
        errno = 0;
        std::feclearexcept(FE_ALL_EXCEPT);
        const double y = std::sin(x);
        std::printf("  sin(%-8s) = %-12g  errno=%-6s FE_INVALID=%d FE_UNDERFLOW=%d FE_INEXACT=%d\n", name, y,
                    errno == EDOM ? "EDOM" : "0", !!std::fetestexcept(FE_INVALID), !!std::fetestexcept(FE_UNDERFLOW),
                    !!std::fetestexcept(FE_INEXACT));
    };
    special("+0", 0.0);
    special("-0", -0.0);
    special("1e-310", 1e-310);
    special("+inf", INFINITY);
    special("NaN", NAN);
    special("1", 1.0);
}
