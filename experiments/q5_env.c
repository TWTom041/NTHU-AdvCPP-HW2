/* q5_env.c - identical source, identical compiler version (GCC 13.3), identical
 * flags (-O2), no undefined behaviour - but different results depending on
 * the target ISA / ABI / C library / CPU the program runs on.
 *
 *   x86-64 glibc :  gcc -O2 q5_env.c -lm
 *   x86-64 musl  :  musl-gcc -O2 -static q5_env.c -lm        (musl-gcc wraps the same gcc 13.3)
 *   aarch64 glibc:  aarch64-linux-gnu-gcc -O2 q5_env.c -lm   (run with qemu-aarch64)
 *
 *   q5_env [out.bin]   -> writes libm results for later cross-platform diffing
 */
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#if defined(__x86_64__)
#define ARCH "x86-64"
#elif defined(__aarch64__)
#define ARCH "aarch64"
#else
#define ARCH "other"
#endif
#if defined(__GLIBC__)
#include <gnu/libc-version.h>
#define LIBC "glibc"
#else
#define LIBC "musl"
#endif

/* a*b + c: GCC contracts this into one fused multiply-add (-ffp-contract=fast is
 * the GNU default) *if the target has an FMA instruction in its baseline ISA*.
 * aarch64: always.  x86-64 baseline (SSE2): no FMA -> separate mul + add. */
__attribute__((noinline)) double mul_sub(double a, double b, double c) { return a * b - c; }
__attribute__((noinline)) double discriminant(double a, double b, double c) { return b * b - 4 * a * c; }

/* u*k + c computed WITHOUT contraction (the volatile store forces the rounding
 * of u*k), so that every platform feeds bit-identical inputs to libm */
static double lin(double u, double k, double c) {
    volatile double t = u * k;
    return t + c;
}

static uint64_t fnv(uint64_t h, double y) {
    uint64_t u;
    memcpy(&u, &y, 8);
    return (h ^ u) * 1099511628211ull;
}

/* deterministic input generator, identical on every platform */
static uint64_t s = 88172645463325252ull;
static double next_uniform(void) {
    s ^= s << 13, s ^= s >> 7, s ^= s << 17;
    return (double)(s >> 11) * 0x1p-53;
}

int main(int argc, char** argv) {
    volatile double one = 1.0, tenth = 0.1, ten = 10.0;
#if defined(__GLIBC__)
    printf("platform: %s / %s %s, sizeof(long double)=%zu, LDBL_MANT_DIG=%d, FLT_EVAL_METHOD=%d\n", ARCH, LIBC,
           gnu_get_libc_version(), sizeof(long double), LDBL_MANT_DIG, (int)FLT_EVAL_METHOD);
#else
    printf("platform: %s / %s, sizeof(long double)=%zu, LDBL_MANT_DIG=%d, FLT_EVAL_METHOD=%d\n", ARCH, LIBC,
           sizeof(long double), LDBL_MANT_DIG, (int)FLT_EVAL_METHOD);
#endif

    /* 1. FMA contraction */
    printf("[1] 0.1*10 - 1          = %.17g\n", mul_sub(tenth, ten, one));
    const double x = 0.1 * 3;  /* double root of t^2 - 2x t + x^2 */
    const double d = discriminant(one, -2 * x, x * x);
    printf("    b*b-4ac (double root) = %.17g  -> sqrt = %g\n", d, sqrt(d));

    /* 2. long double is a different type on different ABIs */
    volatile long double tiny = 1e-17L;
    printf("[2] (1.0L + 1e-17L) - 1.0L = %.6Lg\n", (1.0L + tiny) - 1.0L);
    long double acc = 0;
    for (int i = 0; i < 10; ++i) acc += 0.1L;
    printf("    sum of ten 0.1L == 1.0L ? %s  (diff %.3Lg)\n", acc == 1.0L ? "yes" : "no", acc - 1.0L);

    /* 3. libm: the C standard does not require correct rounding */
    enum { N = 1000000 };
    const char* names[] = {"sin", "cos", "tan", "exp", "log", "pow", "atan2", "cbrt", "sinh", "lgamma"};
    const int nf = 10;
    FILE* out = argc > 1 ? fopen(argv[1], "wb") : NULL;
    printf("[3] libm hashes over %d deterministic inputs:\n    ", N);
    for (int f = 0; f < nf; ++f) {
        uint64_t h = 1469598103934665603ull;
        s = 88172645463325252ull + f;
        for (int i = 0; i < N; ++i) {
            const double u = next_uniform(), w = next_uniform();
            double y;
            switch (f) {
                case 0: y = sin(lin(u, 200, -100)); break;
                case 1: y = cos(lin(u, 200, -100)); break;
                case 2: y = tan(lin(u, 20, -10)); break;
                case 3: y = exp(lin(u, 1400, -700)); break;
                case 4: y = log(lin(u, 1e6, 0)); break;
                case 5: y = pow(lin(u, 10, 0), lin(w, 60, -30)); break;
                case 6: y = atan2(lin(u, 1, -0.5), lin(w, 1, -0.5)); break;
                case 7: y = cbrt(lin(u, 2e6, -1e6)); break;
                case 8: y = sinh(lin(u, 40, -20)); break;
                default: y = lgamma(lin(u, 100, 0.01)); break;
            }
            h = fnv(h, y);
            if (out) fwrite(&y, 8, 1, out);
        }
        printf("%s=%04x ", names[f], (unsigned)(h & 0xffff));
    }
    printf("\n");
    if (out) fclose(out);

    return 0;
}
