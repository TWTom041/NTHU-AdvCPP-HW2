// q3_fastmath - the same source compiled with and without -ffast-math.
//
//   g++ -O3            q3_fastmath.cpp -o q3_strict
//   g++ -O3 -ffast-math q3_fastmath.cpp -o q3_fast
//
// All inputs come from argv / volatile so that nothing is constant-folded; the
// differences are caused only by what -ffast-math allows the optimiser to do.
#include <cerrno>
#include <cfenv>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <chrono>
#include <vector>
#include <xmmintrin.h>

#ifdef __FAST_MATH__
static const char* mode = "-O3 -ffast-math";
#else
static const char* mode = "-O3";
#endif

__attribute__((noinline)) bool contains_nan(const double* v, int n) {
    for (int i = 0; i < n; ++i)
        if (std::isnan(v[i])) return true;
    return false;
}

__attribute__((noinline)) bool is_inf(double x) { return std::isinf(x); }

__attribute__((noinline)) float naive_sum(const float* v, int n) {
    float s = 0;
    for (int i = 0; i < n; ++i) s += v[i];
    return s;
}

// Kahan (compensated) summation: c captures the rounding error of each add.
// Algebraically c == 0, so -fassociative-math may delete it.
__attribute__((noinline)) float kahan_sum(const float* v, int n) {
    float s = 0, c = 0;
    for (int i = 0; i < n; ++i) {
        const float y = v[i] - c;
        const float t = s + y;
        c = (t - s) - y;
        s = t;
    }
    return s;
}

__attribute__((noinline)) double dot(const double* a, const double* b, int n) {
    double s = 0;
    for (int i = 0; i < n; ++i) s += a[i] * b[i];
    return s;
}

__attribute__((noinline)) double add_zero(double x) { return x + 0.0; }

__attribute__((noinline)) void divide(double* out, const double* in, double d, int n) {
    for (int i = 0; i < n; ++i) out[i] = in[i] / d;  // -freciprocal-math: in[i] * (1/d)
}

__attribute__((noinline)) double my_sqrt(double x) { return std::sqrt(x); }

int main(int argc, char** argv) {
    volatile double vzero = 0.0, vone = 1.0, vbig = 1e308, vtiny = 1e-308, vten = 10.0, vthree = 3.0;
    std::printf("=== build: %s ===\n", mode);
    const unsigned csr = _mm_getcsr();
    std::printf("MXCSR = 0x%04x  (FTZ=%d, DAZ=%d)  <- crtfastmath.o sets these at program start\n", csr,
                !!(csr & 0x8000), !!(csr & 0x40));

    // 1. -ffinite-math-only: isnan / isinf may be assumed false
    std::vector<double> v(1000, 1.0);
    v[500] = vzero / vzero;  // NaN created at run time
    std::printf("[1] contains_nan(v) with v[500] = 0.0/0.0      -> %s\n", contains_nan(v.data(), 1000) ? "true" : "false");
    std::printf("    is_inf(1e308 * 10)                          -> %s\n", is_inf(vbig * vten) ? "true" : "false");

    // 2. -fassociative-math: Kahan summation silently becomes naive summation
    const int n = 10'000'000;
    std::vector<float> f(n, 0.1f);
    std::printf("[2] sum of 1e7 x 0.1f: naive = %.6f   kahan = %.6f   (exact 1000000.0149...)\n", naive_sum(f.data(), n),
                kahan_sum(f.data(), n));

    // 3. reassociation -> vectorised reduction in a different order
    std::vector<double> a(1 << 20), b(1 << 20);
    for (int i = 0; i < (1 << 20); ++i) a[i] = std::sin(double(i)) * std::pow(10.0, i % 17 - 8), b[i] = std::cos(double(i));
    std::printf("[3] dot(a, b) over 2^20 terms of mixed magnitude = %.17g\n", dot(a.data(), b.data(), 1 << 20));

    // 4. -fno-signed-zeros: x + 0.0 -> x
    const double mz = -vzero;
    std::printf("[4] add_zero(-0.0) = %g,  1/add_zero(-0.0) = %g\n", add_zero(mz), vone / add_zero(mz));

    // 5. -freciprocal-math: x / d -> x * (1/d)
    const int m = 1'000'000;
    std::vector<double> in(m), out(m);
    for (int i = 0; i < m; ++i) in[i] = 1.0 + i * 1e-6;
    divide(out.data(), in.data(), vthree, m);
    long wrong = 0;
    for (int i = 0; i < m; ++i) {  // reference: the hardware DIVSD instruction (correctly rounded)
        double q = in[i];
        const double d = vthree;
        asm("divsd %1, %0" : "+x"(q) : "x"(d));
        wrong += out[i] != q;
    }
    uint64_t h = 1469598103934665603ull;
    for (double y : out) { uint64_t u; std::memcpy(&u, &y, 8); h = (h ^ u) * 1099511628211ull; }
    std::printf("[5] x/3 for 1e6 values: %ld results differ from correctly-rounded x/3 (hash %016llx)\n", wrong,
                (unsigned long long)h);
    for (int i = 0; i < m; ++i) {
        double q = in[i];
        const double d = vthree;
        asm("divsd %1, %0" : "+x"(q) : "x"(d));
        if (out[i] != q) {
            std::printf("    e.g. %.17g / 3 -> %.17g  (correctly rounded: %.17g)\n", in[i], out[i], q);
            break;
        }
    }

    // 6. flush-to-zero / denormals-are-zero (MXCSR set by crtfastmath.o)
    const double sub = vtiny * 1e-5;
    std::printf("[6] 1e-308 * 1e-5 = %g   (subnormal %s)\n", sub, sub == 0 ? "flushed to zero" : "kept");
    volatile double vsub = 4.9e-324;
    std::printf("    4.9e-324 (smallest subnormal) * 1e300 = %g\n", vsub * 1e300);

    // 7. -fno-math-errno: sqrt of a negative number no longer sets errno
    errno = 0;
    const double r = my_sqrt(-vone);
    std::printf("[7] sqrt(-1) = %g, errno = %s\n", r, errno == EDOM ? "EDOM" : "0");

    // 8. speed of the reduction (vectorisable only with reassociation)
    auto t0 = std::chrono::steady_clock::now();
    double s = 0;
    for (int k = 0; k < 200; ++k) {
        asm volatile("" : : "r"(a.data()) : "memory");  // forbid hoisting the call out of the loop
        s += dot(a.data(), b.data(), 1 << 20);
    }
    const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::printf("[8] dot() speed: %.3f ns per element  (checksum %.6g)\n", sec * 1e9 / (200.0 * (1 << 20)), s);
    return 0;
}
