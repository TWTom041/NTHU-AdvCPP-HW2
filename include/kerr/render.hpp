// render.hpp - image-level helpers shared by the command-line tools.
#pragma once

#include <cstdio>
#include <string>
#include <vector>

#include "tracer.hpp"

namespace kerr {

struct Image {
    int width = 0, height = 0;
    std::vector<double> intensity;     // row-major, y = 0 is the bottom row (beta < 0)
    std::vector<std::uint8_t> status;  // Status per pixel
    std::vector<double> final_theta, final_phi;
    std::vector<long> steps;

    Image(int w, int h)
        : width(w), height(h), intensity(size_t(w) * h), status(size_t(w) * h), final_theta(size_t(w) * h),
          final_phi(size_t(w) * h), steps(size_t(w) * h) {}

    void store(int idx, const auto& res) {
        intensity[idx] = double(res.intensity);
        status[idx] = std::uint8_t(res.status);
        final_theta[idx] = double(res.y[2]);
        final_phi[idx] = double(res.y[3]);
        steps[idx] = res.steps;
    }

    // Binary dump: header "KIMG" w h, then intensity[], status[], theta[], phi[], steps[]
    bool save(const std::string& path) const {
        FILE* f = std::fopen(path.c_str(), "wb");
        if (!f) return false;
        std::fwrite("KIMG", 1, 4, f);
        std::fwrite(&width, sizeof width, 1, f);
        std::fwrite(&height, sizeof height, 1, f);
        const size_t n = intensity.size();
        std::fwrite(intensity.data(), sizeof(double), n, f);
        std::fwrite(status.data(), 1, n, f);
        std::fwrite(final_theta.data(), sizeof(double), n, f);
        std::fwrite(final_phi.data(), sizeof(double), n, f);
        std::fwrite(steps.data(), sizeof(long), n, f);
        std::fclose(f);
        return true;
    }
};

template <class T>
void trace_pixel(const Spacetime<T>& st, const Camera<T>& cam, const TraceConfig<T>& cfg, int x, int y,
                 Image& img) {
    T alpha, beta;
    cam.impact(x, y, alpha, beta);
    const auto res = trace(st, cam.initial_state(st, alpha, beta), cfg);
    img.store(y * cam.width + x, res);
}

// Default parallel renderer: OpenMP, dynamic schedule (see parallel_bench for alternatives).
template <class T>
void render(const Spacetime<T>& st, const Camera<T>& cam, const TraceConfig<T>& cfg, Image& img) {
    const int n = cam.width * cam.height;
#pragma omp parallel for schedule(dynamic, 16)
    for (int i = 0; i < n; ++i) trace_pixel(st, cam, cfg, i % cam.width, i / cam.width, img);
}

}  // namespace kerr
