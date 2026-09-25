// q5_threads - a well-defined OpenMP reduction whose result depends on how many
// cores the machine has (OpenMP's default thread count = number of CPUs).
//
//   g++ -O2 -fopenmp q5_threads.cpp && for t in 1 2 3 4; do OMP_NUM_THREADS=$t ./a.out; done
#include <cmath>
#include <cstdio>
#include <vector>

#include <omp.h>

int main() {
    const int n = 10'000'000;
    std::vector<double> v(n);
    for (int i = 0; i < n; ++i) v[i] = std::sin(i * 0.001) * std::pow(10.0, i % 9 - 4);  // mixed magnitudes

    double s = 0;
#pragma omp parallel for reduction(+ : s) schedule(static)
    for (int i = 0; i < n; ++i) s += v[i];

    double serial = 0;
    for (int i = 0; i < n; ++i) serial += v[i];
    std::printf("threads=%d  parallel sum = %.17g   serial sum = %.17g   diff = %.3g\n", omp_get_max_threads(), s,
                serial, s - serial);
}
