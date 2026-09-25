/* q2_ifunc.c - which implementation of sin() did glibc's IFUNC resolver pick?
 * Prints the offset of the resolved function inside libm.so.6; compare it with
 * the section table (readelf -S): glibc builds the FMA variant (__sin_fma)
 * into a separate section ".text.fma".
 *   gcc -O2 q2_ifunc.c -o q2_ifunc -lm -ldl */
#define _GNU_SOURCE
#include <dlfcn.h>
#include <math.h>
#include <stdio.h>

int main(void) {
    double (*volatile f)(double) = sin; /* resolved through the GOT */
    Dl_info info;
    dladdr((void*)f, &info);
    printf("sin resolved to %p = %s + 0x%lx\n", (void*)f, info.dli_fname,
           (unsigned long)((char*)f - (char*)info.dli_fbase));
    return 0;
}
