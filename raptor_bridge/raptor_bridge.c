/* raptor_bridge.c - drive RAPTOR's own geodesic code ray by ray.
 *
 * RAPTOR's sources define their globals inside headers, which modern GCC
 * (-fno-common) refuses to link across translation units.  Including every
 * RAPTOR .c file (except main.c) into this single translation unit side-steps
 * that without modifying RAPTOR.
 *
 * Coordinates: RAPTOR's default "MKS" metric with hslope = 1 and R0 = 0, i.e.
 *   x1 = log r,  x2 = theta / pi,  x3 = phi   (plain Kerr-Schild, rescaled).
 * The ray loop below is RAPTOR's integrate_geodesic() without the radiative
 * transfer part.
 */
#include "core.c"
#include "GRmath.c"
#include "integrator.c"
#include "metric.c"
#include "radiative_transfer.c"
#include "raptor_harm_model.c"
/* utilities.c defines write_VTK_image with a signature that contradicts functions.h */
#define write_VTK_image write_VTK_image_unused
#include "utilities.c"
#undef write_VTK_image
#include "j_nu.c"
#include "rcarry.c"
#include "newtonraphson.c"
#include "grmhd.c"

#include "raptor_bridge.h"

void raptor_setup(double spin, double incl_deg, double step, double horizon_margin_) {
    a = spin;
    hslope = 1.0;
    R0 = 0.0;
    INCLINATION = incl_deg;
    STEPSIZE = step;
    MBH = 1.0;
    CAM_SIZE_X = CAM_SIZE_Y = 40;
    IMG_WIDTH = IMG_HEIGHT = 1;
    set_constants(0.0);  /* sets cutoff_inner = r_+ (1 + horizon_marg) */
    cutoff_inner = (1. + sqrt(1. - a * a)) * (1. + horizon_margin_);
    initrcarry(982451653);
}

void raptor_initial(double alpha, double beta, double out_ks[8]) {
    real photon_u[8];
    initialize_photon(alpha, beta, photon_u, 0.0);
    const double r = exp(photon_u[1]);
    const double o[8] = {photon_u[0], r, M_PI * photon_u[2], photon_u[3],
                         photon_u[4], r * photon_u[5], M_PI * photon_u[6], photon_u[7]};
    for (int i = 0; i < 8; i++) out_ks[i] = o[i];
}

int raptor_trace(double alpha, double beta, int method, double out_ks[8], long* nsteps) {
    real photon_u[8];
    initialize_photon(alpha, beta, photon_u, 0.0);

    real r_current = exp(photon_u[1]);
    long steps = 0;
    while (r_current < cutoff_outer && r_current > cutoff_inner && steps < max_steps) {
        real X_u[4] = {photon_u[0], photon_u[1], photon_u[2], photon_u[3]};
        real k_u[4] = {photon_u[4], photon_u[5], photon_u[6], photon_u[7]};
        real dl = stepsize(X_u, k_u);
        if (method == RAPTOR_RK4)
            rk4_step(photon_u, dl);
        else
            rk2_step(photon_u, dl);
        r_current = exp(photon_u[1]);
        steps++;
    }
    /* convert (log r, theta/pi, phi) and contravariant U to plain KS */
    const double r = exp(photon_u[1]);
    out_ks[0] = photon_u[0];
    out_ks[1] = r;
    out_ks[2] = M_PI * photon_u[2];
    out_ks[3] = photon_u[3];
    out_ks[4] = photon_u[4];
    out_ks[5] = r * photon_u[5];
    out_ks[6] = M_PI * photon_u[6];
    out_ks[7] = photon_u[7];
    *nsteps = steps;
    if (!(r_current == r_current)) return RAPTOR_NOTFINITE;
    if (r_current <= cutoff_inner) return RAPTOR_CAPTURED;
    if (r_current >= cutoff_outer) return RAPTOR_ESCAPED;
    return RAPTOR_MAXSTEPS;
}
