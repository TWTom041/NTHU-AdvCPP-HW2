/* raptor_bridge.h - C interface to RAPTOR's geodesic integrator. */
#ifndef RAPTOR_BRIDGE_H
#define RAPTOR_BRIDGE_H

#ifdef __cplusplus
extern "C" {
#endif

enum { RAPTOR_RK2 = 0, RAPTOR_RK4 = 1 };
enum { RAPTOR_CAPTURED = 1, RAPTOR_ESCAPED = 2, RAPTOR_MAXSTEPS = 4, RAPTOR_NOTFINITE = 5 };

/* spin a, inclination [deg], RAPTOR STEPSIZE, inner cutoff r_+ (1 + margin) */
void raptor_setup(double spin, double incl_deg, double step, double horizon_margin);

/* Initial ray (same layout as raptor_trace's output) from RAPTOR's initialize_photon. */
void raptor_initial(double alpha, double beta, double out_ks[8]);

/* Trace one ray with RAPTOR's initialize_photon / stepsize / rk{2,4}_step.
 * out_ks = final (t, r, theta, phi, k^t, k^r, k^theta, k^phi) in plain KS,
 * contravariant k (RAPTOR integrates backwards with dlambda < 0). */
int raptor_trace(double alpha, double beta, int method, double out_ks[8], long* nsteps);

#ifdef __cplusplus
}
#endif
#endif
