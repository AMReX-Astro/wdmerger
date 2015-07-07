// Define a preprocessor directive that determines if we're doing a full merger simulation.
// This is useful for sharing C++ files like the diagnostics between multiple test problems.

#define merger

// Preprocessor directive for allowing us to do a post-initialization step.

#define do_problem_post_init

// Routine for determining white dwarf masses, centers of mass and velocities.

void wdCOM (Real time, Real& mass_p, Real& mass_s, Real* com_p, Real* com_s, Real* vel_p, Real* vel_s);

// Volume within white dwarfs at a given density cutoff.

void volInBoundary (Real               time,
		    Real&              vol_p,
		    Real&              vol_s,
                    Real               rho_cutoff);

// Calculate gravitational wave signal.

void gwstrain (Real time, Real& h_plus, Real& h_cross);

// Problem post-initialization routine.

void problem_post_init();

// Parameter that controls whether we do a post-initialization relaxation.

static int do_relax;

