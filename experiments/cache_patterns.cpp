// cache_patterns - access patterns whose cache hit rates we measure with
// valgrind --tool=cachegrind (hardware PMU counters are not exposed in the VM).
//
//   cache_patterns row|col float|double [N=2048]
// Sums an N x N matrix either along rows (unit stride, 16 floats / 8 doubles
// per 64-byte line) or along columns (stride N*sizeof(T): every access touches
// a new line).  Also prints wall time so the effect can be seen without the
// simulator.
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

template <class T> double run(bool row, size_t n) {
    std::vector<T> a(n * n);
    for (size_t i = 0; i < n * n; ++i) a[i] = T(i % 7);
    T s = 0;
    const auto t0 = std::chrono::steady_clock::now();
    if (row)
        for (size_t i = 0; i < n; ++i)
            for (size_t j = 0; j < n; ++j) s += a[i * n + j];
    else
        for (size_t j = 0; j < n; ++j)
            for (size_t i = 0; i < n; ++i) s += a[i * n + j];
    const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::printf("%s %s n=%zu sum=%.1f time=%.3f s (%.2f ns/element)\n", row ? "row" : "col",
                sizeof(T) == 4 ? "float" : "double", n, double(s), sec, sec * 1e9 / double(n * n));
    return sec;
}

int main(int argc, char** argv) {
    if (argc < 3) { std::fprintf(stderr, "usage: %s row|col float|double [N]\n", argv[0]); return 2; }
    const bool row = !std::strcmp(argv[1], "row");
    const size_t n = argc > 3 ? std::strtoul(argv[3], nullptr, 10) : 2048;
    if (!std::strcmp(argv[2], "float")) run<float>(row, n);
    else run<double>(row, n);
}
