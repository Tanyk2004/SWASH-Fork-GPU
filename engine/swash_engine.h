/* C API of libswash_engine: control layer around the SWASH solver.
 *
 * One SWASH instance per process (all state is global).  The working
 * directory passed to swash_init must contain the INPUT command file, exactly
 * as for swash.exe.  Single MPI rank only.
 *
 * Dense fields are float arrays of nx*ny values in C row-major order
 * ([ny][nx], x fastest) where nx, ny come from swash_grid_info and include
 * SWASH's virtual boundary cells.  Inactive points carry -99.
 *
 * Return codes: 0 ok, negative = error (see the Fortran source for details).
 */
#ifndef SWASH_ENGINE_H
#define SWASH_ENGINE_H

#ifdef __cplusplus
extern "C" {
#endif

/* lifecycle */
int    swash_init(const char *workdir);      /* 0 ready, 1 nothing to compute, 2 error, <0 init failure */
int    swash_step(int nsteps);               /* steps taken (0 when the run has ended), -1 on error */
int    swash_status(void);                   /* -1 not initialised, 0 stepping, 1 finished, 2 error */
int    swash_finish(void);                   /* output collection, timings, cleanup, MPI shutdown */
typedef void (*swash_step_cb)(double t);
void   swash_set_step_callback(swash_step_cb fn);   /* called after every completed step */

/* queries */
double swash_time(void);                     /* current simulation time [s] */
double swash_dt(void);                       /* current time step [s] (adaptive in explicit mode) */
double swash_end_time(void);                 /* end time of the current COMPUTE block [s] */
void   swash_grid_info(int *nx, int *ny, int *nk, int *npts,
                       double *x0, double *y0, double *dx, double *dy);
/* name: WATL BOTL DEP HS U V W Q HK ZK WETS WETU WETV BRKS PRESP PATM SRCM XP YP
 * k: layer index for U V Q HK (1..nk), W ZK (0..nk); ignored otherwise.
 * U and V are given in their staggered (u-/v-point) positions. */
int    swash_get_field(const char *name, float *out, int k);

/* forcing (dense layout, wl-points) */
int    swash_set_pressure_field(const float *pa);   /* surface pressure [Pa]; hull = rho*g*draft */
int    swash_clear_pressure_field(void);
int    swash_set_mass_source(const float *srcm);    /* [m/s], added to the water level each step */
int    swash_set_bathymetry(const float *depth);    /* depth below still water level [m], + down */

/* in-memory checkpoints, slots 1..8 */
int    swash_save_state(int slot);
int    swash_restore_state(int slot);

int    swash_acc_enabled(void);              /* 1 when kernels run on the GPU */

#ifdef __cplusplus
}
#endif
#endif
