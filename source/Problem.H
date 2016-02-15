// Preprocessor directive for allowing us to do a post-initialization update.

#if (BL_SPACEDIM==3)
#ifndef do_problem_post_init
#define do_problem_post_init
#endif
#endif

// Preprocessor directive for allowing us to do a post-timestep update.

#ifndef do_problem_post_timestep
#define do_problem_post_timestep
#endif

// Routine for determining white dwarf masses, centers of mass and velocities.

void wdCOM (Real time, Real& mass_p, Real& mass_s, Real* com_p, Real* com_s, Real* vel_p, Real* vel_s);

// Volume within white dwarfs at a given density cutoff.

void volInBoundary (Real               time,
		    Real&              vol_p,
		    Real&              vol_s,
                    Real               rho_cutoff);

// Calculate gravitational wave signal.

void gwstrain (Real time,
	       Real& h_plus_1, Real& h_cross_1,
	       Real& h_plus_2, Real& h_cross_2,
	       Real& h_plus_3, Real& h_cross_3);

// Computes standard dot product of two three-vectors.

Real dot_product(const Real a[], const Real b[]);

// Computes standard cross-product of two three-vectors.

void cross_product(const Real a[], const Real b[], Real c[]);

// Computes norm of a three vector.

Real norm(const Real a[]);

// Problem post-initialization routine.

void problem_post_init();

// End of timestep analytics on the white dwarf properties.

void problem_post_timestep();




// Data we want to store, for interfacing with the Fortran

static bool relaxation_is_done;
static int problem;
static int accurate_IC_frame;
