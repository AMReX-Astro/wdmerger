// Define a preprocessor directive that determines if we're doing a full merger simulation.
// This is useful for sharing C++ files like the diagnostics between multiple test problems.

#define merger
#define do_problem_post_init

void wdCOM (Real time, Real& mass_p, Real& mass_s, Real* com_p, Real* com_s, Real* vel_p, Real* vel_s);

void volInBoundary (Real               time,
		    Real&              vol_p,
		    Real&              vol_s,
                    Real               rho_cutoff);

void problem_post_init();
