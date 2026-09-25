// kerr_render - render a Kerr black hole (shadow + optional thin disk).
//
//   kerr_render [--type float|double|ldouble|quad|f16] [--method rk4|rk2|dp45]
//               [--spin 0.9375] [--incl 60] [--fov 40] [--res 256] [--step 0.03]
//               [--rtol 1e-10] [--disk] [--rcam 1e4] [--out img.bin]
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

#include "kerr/render.hpp"

using namespace kerr;

struct Args {
    std::string type = "double", method = "rk4", out = "";
    double spin = 0.9375, incl = 60, fov = 40, step = 0.03, rtol = 1e-10, rcam = 1e4;
    int res = 256;
    bool disk = false;
};

template <class T> int run(const Args& a) {
    Spacetime<T> st{T(a.spin)};
    Camera<T> cam{T(a.rcam), T(a.incl) * m::pi<T>() / T(180), T(a.fov), a.res, a.res};
    TraceConfig<T> cfg;
    cfg.method = a.method == "rk2" ? Method::RK2 : a.method == "dp45" ? Method::DP45 : Method::RK4;
    cfg.step_frac = T(a.step);
    cfg.rtol = T(a.rtol);
    cfg.atol = T(a.rtol * 1e-2);
    cfg.r_out = T(a.rcam * 1.01);
    cfg.disk.enabled = a.disk;
    cfg.disk.r_in = st.r_isco();
    cfg.disk.r_out = T(20);

    Image img(a.res, a.res);
    const auto t0 = std::chrono::steady_clock::now();
    render(st, cam, cfg, img);
    const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

    long counts[6] = {}, steps = 0;
    for (size_t i = 0; i < img.status.size(); ++i) counts[img.status[i]]++, steps += img.steps[i];
    double flux = 0;
    for (double v : img.intensity) flux += v;
    std::printf("type=%s method=%s spin=%g incl=%g res=%d  time=%.3f s  (%.0f rays/s, %.1f ns/step)\n",
                type_name<T>(), a.method.c_str(), a.spin, a.incl, a.res, sec, img.status.size() / sec,
                sec * 1e9 / double(steps));
    std::printf("captured=%ld escaped=%ld disk=%ld maxsteps=%ld notfinite=%ld  mean steps=%.1f  flux=%.10e\n",
                counts[1], counts[2], counts[3], counts[4], counts[5], double(steps) / img.status.size(), flux);
    if (!a.out.empty() && !img.save(a.out)) {
        std::fprintf(stderr, "cannot write %s\n", a.out.c_str());
        return 1;
    }
    return 0;
}

int main(int argc, char** argv) {
    Args a;
    for (int i = 1; i < argc; ++i) {
        auto next = [&] { return i + 1 < argc ? argv[++i] : (std::fprintf(stderr, "missing value\n"), std::exit(2), ""); };
        if (!std::strcmp(argv[i], "--type")) a.type = next();
        else if (!std::strcmp(argv[i], "--method")) a.method = next();
        else if (!std::strcmp(argv[i], "--spin")) a.spin = std::atof(next());
        else if (!std::strcmp(argv[i], "--incl")) a.incl = std::atof(next());
        else if (!std::strcmp(argv[i], "--fov")) a.fov = std::atof(next());
        else if (!std::strcmp(argv[i], "--res")) a.res = std::atoi(next());
        else if (!std::strcmp(argv[i], "--step")) a.step = std::atof(next());
        else if (!std::strcmp(argv[i], "--rtol")) a.rtol = std::atof(next());
        else if (!std::strcmp(argv[i], "--rcam")) a.rcam = std::atof(next());
        else if (!std::strcmp(argv[i], "--disk")) a.disk = true;
        else if (!std::strcmp(argv[i], "--out")) a.out = next();
        else { std::fprintf(stderr, "unknown option %s\n", argv[i]); return 2; }
    }
    if (a.type == "float") return run<float>(a);
    if (a.type == "double") return run<double>(a);
    if (a.type == "ldouble") return run<long double>(a);
#ifdef KERR_HAS_FLOAT128
    if (a.type == "quad") return run<__float128>(a);
#endif
#ifdef __STDCPP_FLOAT16_T__
    if (a.type == "f16") return run<std::float16_t>(a);
#endif
    std::fprintf(stderr, "unknown type %s\n", a.type.c_str());
    return 2;
}
