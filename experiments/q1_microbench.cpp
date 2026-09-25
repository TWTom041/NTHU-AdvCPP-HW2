// q1_microbench - latency and throughput of basic operations per floating type.
//
// latency    : one dependent chain  x = x OP c   -> time of one operation
//              (what scalar code with data dependencies, e.g. an ODE step, sees)
// throughput : c[i] = a[i] OP b[i] over 2048-element arrays (L1 resident),
//              which the compiler may auto-vectorise -> best case per element
//
// Build twice to see the effect of the ISA:
//   -O2                         (x86-64 baseline: SSE2, 128-bit vectors)
//   -O3 -march=native           (here: AVX-512 incl. AVX512-FP16)
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <quadmath.h>
#include <vector>
#if __has_include(<stdfloat>)
#include <stdfloat>
#endif

using clk = std::chrono::steady_clock;

template <class T> inline void escape(T& x) { asm volatile("" : "+m"(x)); }
inline void clobber() { asm volatile("" ::: "memory"); }

// ---- per-type math -------------------------------------------------------
inline float msqrt(float x) { return std::sqrt(x); }
inline double msqrt(double x) { return std::sqrt(x); }
inline long double msqrt(long double x) { return std::sqrt(x); }
inline __float128 msqrt(__float128 x) { return sqrtq(x); }
inline float msin(float x) { return std::sin(x); }
inline double msin(double x) { return std::sin(x); }
inline long double msin(long double x) { return std::sin(x); }
inline __float128 msin(__float128 x) { return sinq(x); }
#ifdef __STDCPP_FLOAT16_T__
using f16 = std::float16_t;
inline f16 msqrt(f16 x) { return f16(std::sqrt(float(x))); }
inline f16 msin(f16 x) { return f16(std::sin(float(x))); }
#endif

enum Op { ADD, MUL, DIV, SQRT, SIN, NOPS };
const char* op_name[] = {"add", "mul", "div", "sqrt", "sin"};

template <class T> double latency(Op op) {
    const long n = op == SIN ? 2'000'000 : 20'000'000;
    T x = T(1.5f), c = T(1e-3f), c1 = T(1.0009765625f), c2 = T(0.9990234375f), cd = T(2.25f);
    escape(c), escape(c1), escape(c2), escape(cd);
    double best = 1e30;
    for (int rep = 0; rep < 3; ++rep) {
        const auto t0 = clk::now();
        switch (op) {
            case ADD: for (long i = 0; i < n; i += 2) { x = x + c; x = x - c; } break;
            case MUL: for (long i = 0; i < n; i += 2) { x = x * c1; x = x * c2; } break;
            case DIV: for (long i = 0; i < n; ++i) x = cd / x; break;
            case SQRT: for (long i = 0; i < n; ++i) x = msqrt(x); break;
            case SIN: for (long i = 0; i < n; ++i) x = msin(x); break;
            default: break;
        }
        escape(x);
        best = std::min(best, std::chrono::duration<double>(clk::now() - t0).count() * 1e9 / double(n));
    }
    return best;
}

template <class T> double throughput(Op op) {
    constexpr int N = 2048;
    std::vector<T> a(N), b(N), c(N);
    for (int i = 0; i < N; ++i) a[i] = T(1.0f + float(i % 97) * 0.01f), b[i] = T(2.0f - float(i % 89) * 0.01f);
    const long reps = op == SIN ? 200 : 4000;
    T* __restrict pa = a.data();
    T* __restrict pb = b.data();
    T* __restrict pc = c.data();
    double best = 1e30;
    for (int rep = 0; rep < 3; ++rep) {
        const auto t0 = clk::now();
        for (long r = 0; r < reps; ++r) {
            switch (op) {
                case ADD: for (int i = 0; i < N; ++i) pc[i] = pa[i] + pb[i]; break;
                case MUL: for (int i = 0; i < N; ++i) pc[i] = pa[i] * pb[i]; break;
                case DIV: for (int i = 0; i < N; ++i) pc[i] = pa[i] / pb[i]; break;
                case SQRT: for (int i = 0; i < N; ++i) pc[i] = msqrt(pa[i]); break;
                case SIN: for (int i = 0; i < N; ++i) pc[i] = msin(pa[i]); break;
                default: break;
            }
            clobber();
        }
        best = std::min(best, std::chrono::duration<double>(clk::now() - t0).count() * 1e9 / double(reps * N));
    }
    return best;
}

template <class T> void row(const char* name) {
    std::printf("%-12s %3zu B |", name, sizeof(T));
    for (int op = 0; op < NOPS; ++op) std::printf(" %7.2f", latency<T>(Op(op)));
    std::printf(" |");
    for (int op = 0; op < NOPS; ++op) std::printf(" %7.3f", throughput<T>(Op(op)));
    std::printf("\n");
    std::fflush(stdout);
}

int main() {
    std::printf("ns per operation (lower is better).  Build: %s\n", BUILD_FLAGS);
    std::printf("%-12s %5s |%-40s |%-40s\n", "type", "size", "  latency (dependent chain)", "  throughput (independent, arrays)");
    std::printf("%-12s %5s |", "", "");
    for (int k = 0; k < 2; ++k) {
        for (int op = 0; op < NOPS; ++op) std::printf(" %7s", op_name[op]);
        std::printf(" |");
    }
    std::printf("\n");
#ifdef __STDCPP_FLOAT16_T__
    row<std::float16_t>("float16_t");
#endif
    row<float>("float");
    row<double>("double");
    row<long double>("long double");
    row<__float128>("__float128");
}
