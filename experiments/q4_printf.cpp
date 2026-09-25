/* q4_printf - is printf("%.9f") / std::fixed << setprecision(9) really
 * "round half up at the 9th decimal"?
 *
 * Compiles as C++ (g++) and as C (gcc / musl-gcc -x c), the iostream part is
 * C++ only.  Every value is also printed with %.40f, which glibc prints
 * *exactly* (the binary value has a finite decimal expansion).
 */
#include <fenv.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#ifdef __cplusplus
#include <iomanip>
#include <iostream>
#include <sstream>
#endif

/* "school" rounding of the *decimal literal* the programmer had in mind */
static void show(const char* literal, double x, const char* naive) {
    char buf[64];
    snprintf(buf, sizeof buf, "%.9f", x);
#ifdef __cplusplus
    std::ostringstream os;
    os << std::fixed << std::setprecision(9) << x;
    const char* same = os.str() == buf ? "same" : os.str().c_str();
#else
    const char* same = "(C build)";
#endif
    printf("%-16s %%.9f=%-16s naive=%-14s %-4s iostream:%-5s exact=%.40f\n", literal, buf, naive,
           strcmp(buf, naive) ? "DIFF" : "", same, x);
}

int main(void) {
    volatile double v;
    printf("--- 1. the binary value, not the decimal literal, is rounded ---\n");
    show("1.0000000005", 1.0000000005, "1.000000001");
    show("2.0000000005", 2.0000000005, "2.000000001");
    show("0.1234567895", 0.1234567895, "0.123456790");
    show("1.0000000015", 1.0000000015, "1.000000002");

    printf("--- 2. exact ties are rounded half-to-even (IEEE default), not half-up ---\n");
    show("1/1024", 1.0 / 1024, "0.000976563");       /* 0.0009765625 exactly */
    show("5/1024", 5.0 / 1024, "0.004882813");       /* 0.0048828125 exactly */
    show("3/1024", 3.0 / 1024, "0.002929688");       /* 0.0029296875 exactly */
    printf("    %%.0f: 0.5 -> %.0f   1.5 -> %.0f   2.5 -> %.0f   3.5 -> %.0f\n", 0.5, 1.5, 2.5, 3.5);

    printf("--- 3. negative values that round to zero keep their sign ---\n");
    show("-1e-10", -1e-10, "0.000000000");

    printf("--- 4. float arguments are promoted to double and printed exactly ---\n");
    show("0.1f", (double)0.1f, "0.100000000");
    show("16777217.0f", (double)16777217.0f, "16777217.000000000");

    printf("--- 5. printf honours the dynamic rounding mode (fesetround) ---\n");
    v = 0.1;
    const int modes[] = {FE_TONEAREST, FE_UPWARD, FE_DOWNWARD, FE_TOWARDZERO};
    const char* names[] = {"FE_TONEAREST", "FE_UPWARD", "FE_DOWNWARD", "FE_TOWARDZERO"};
    for (int i = 0; i < 4; ++i) {
        fesetround(modes[i]);
        char a[32], b[32], c[32];
        snprintf(a, sizeof a, "%.9f", v);
        snprintf(b, sizeof b, "%.9f", 0.3);
        snprintf(c, sizeof c, "%.9f", -0.3);
#ifdef __cplusplus
        std::ostringstream os;
        os << std::fixed << std::setprecision(9) << double(v);
        printf("  %-14s 0.1 -> %s   0.3 -> %s   -0.3 -> %s   (iostream 0.1 -> %s)\n", names[i], a, b, c,
               os.str().c_str());
#else
        printf("  %-14s 0.1 -> %s   0.3 -> %s   -0.3 -> %s\n", names[i], a, b, c);
#endif
    }
    fesetround(FE_TONEAREST);
    return 0;
}
