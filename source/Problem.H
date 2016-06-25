// Preprocessor directive for allowing us to do a post-initialization update.

#ifndef do_problem_post_init
#define do_problem_post_init
#endif

// Same for post-restart.

#ifndef do_problem_post_restart
#define do_problem_post_restart
#endif

// Preprocessor directive for allowing us to do a post-timestep update.

#ifndef do_problem_post_timestep
#define do_problem_post_timestep
#endif

// Routine for determining white dwarf masses, centers of mass and velocities.

void wdCOM (Real time, Real& mass_p, Real& mass_s, Real* com_p, Real* com_s, Real* vel_p, Real* vel_s, bool local = false);

// Volume within white dwarfs at a given density cutoff.

void volInBoundary (Real time, Real& vol_p, Real& vol_s, Real rho_cutoff, bool local = false);

// Calculate gravitational wave signal.

void gwstrain (Real time,
	       Real& h_plus_1, Real& h_cross_1,
	       Real& h_plus_2, Real& h_cross_2,
	       Real& h_plus_3, Real& h_cross_3,
	       bool local = false);

// Computes standard dot product of two three-vectors.

Real dot_product(const Real a[], const Real b[]);

// Computes standard cross-product of two three-vectors.

void cross_product(const Real a[], const Real b[], Real c[]);

// Computes norm of a three vector.

Real norm(const Real a[]);

// Problem post-initialization routine.

void problem_post_init();

// Problem post-restart routine.

void problem_post_restart();

// End of timestep analytics on the white dwarf properties.

void problem_post_timestep();

// Write out the git hashes for the various parts of the code.

void writeGitHashes(std::ostream& log);



// Data we want to store, for interfacing with the Fortran

static bool relaxation_is_done;
static int problem;

// Radii of the WDs at various density thresholds.

Real rad_p[7];
Real rad_s[7];

// Rate of change of stellar masses

Real mdot_p;
Real mdot_s;

// Global extrema of various quantities over the whole simulation

Real T_global_max;
Real rho_global_max;
Real ts_te_global_max;

// Current values of the above extrema

Real T_curr_max;
Real rho_curr_max;
Real ts_te_curr_max;

// Value of the total energy on the domain over the last several timesteps

static const int num_previous_ener_timesteps = 5;
Real total_ener_array[num_previous_ener_timesteps];
