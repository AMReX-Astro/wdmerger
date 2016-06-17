#include "Castro.H"
#include "Castro_F.H"

#include "Gravity.H"
#include <Gravity_F.H>

#include "ParmParse.H"

#include "buildInfo.H"

bool Castro::relaxation_is_done = false;
int Castro::problem = -1;

#ifdef do_problem_post_timestep
void
Castro::problem_post_timestep()
{

    if (level != 0) return;

    int jobDoneStatus;

    int finest_level = parent->finestLevel();
    Real time = state[State_Type].curTime();
    Real dt = parent->dtLevel(0);

    if (time == 0.0) dt = 0.0; // dtLevel returns the next timestep for t = 0, so overwrite    

    Real mass_p       = 0.0;
    Real mass_s       = 0.0;

    Real com_p[3]     = { 0.0 };
    Real com_s[3]     = { 0.0 };

    Real vel_p[3]     = { 0.0 };
    Real vel_s[3]     = { 0.0 };

    Real old_mass_p   = 0.0;
    Real old_mass_s   = 0.0;

    Real old_com_p[3] = { 0.0 };
    Real old_com_s[3] = { 0.0 };

    Real old_vel_p[3] = { 0.0 };
    Real old_vel_s[3] = { 0.0 };

    // Effective volume of the stars at various density cutoffs.

    Real vol_p[7] = { 0.0 };
    Real vol_s[7] = { 0.0 };

    // Average density of the stars.

    Real rho_avg_p = 0.0;
    Real rho_avg_s = 0.0;

    // Gravitational free-fall timescale of the stars.

    Real t_ff_p       = 0.0;
    Real t_ff_s       = 0.0;

    bool local_flag = true;


    // Get the current job done status.
    get_job_status(&jobDoneStatus);

    // Get the current stellar data
    get_star_data(com_p, com_s, vel_p, vel_s, &mass_p, &mass_s, &t_ff_p, &t_ff_s);

    // Update the problem center using the system bulk velocity
    update_center(&time);

    for ( int i = 0; i < 3; i++ ) {
      com_p[i] += vel_p[i] * dt;
      com_s[i] += vel_s[i] * dt;
    }

    // Now send this first estimate of the COM to Fortran, and then re-calculate
    // a more accurate result using it as a starting point.

    set_star_data(com_p, com_s, vel_p, vel_s, &mass_p, &mass_s, &t_ff_p, &t_ff_s);

    old_mass_p = mass_p;
    old_mass_s = mass_s;

    mass_p = 0.0;
    mass_s = 0.0;

    for ( int i = 0; i < 3; i++ ) {

      old_com_p[i] = com_p[i];
      old_com_s[i] = com_s[i];

      old_vel_p[i] = vel_p[i];
      old_vel_s[i] = vel_s[i];

      com_p[i] = 0.0;
      com_s[i] = 0.0;

      vel_p[i] = 0.0;
      vel_s[i] = 0.0;

    }

    // Compute white dwarf centers of mass

    Castro::wdCOM(time, mass_p, mass_s, com_p, com_s, vel_p, vel_s, local_flag);

    // Compute effective radii of stars at various density cutoffs

    for (int i = 0; i <= 6; ++i)
        Castro::volInBoundary(time, vol_p[i], vol_s[i], pow(10.0,i), local_flag);

    // Compute extrema

    T_curr_max     = 0.0;
    rho_curr_max   = 0.0;
    ts_te_curr_max = 0.0;

    for (int lev = 0; lev <= finest_level; lev++) {

      MultiFab& S_new = parent->getLevel(lev).get_new_data(State_Type);

      T_curr_max = std::max(T_curr_max, S_new.max(Temp, 0, local_flag));
      rho_curr_max = std::max(rho_curr_max, S_new.max(Density, 0, local_flag));

      if (lev == finest_level) {

        MultiFab* ts_te_MF = parent->getLevel(lev).derive("t_sound_t_enuc", time, 0);
	ts_te_curr_max = std::max(ts_te_curr_max, ts_te_MF->max(0,0,local_flag));
	delete ts_te_MF;

      }

    }

    // Do all of the reductions in a single step.

    const int nfoo_sum = 28;
    Real foo_sum[nfoo_sum] = { 0.0 };

    for (int i = 0; i <= 6; ++i) {
      foo_sum[i  ] = vol_p[i];
      foo_sum[i+7] = vol_s[i];
    }

    foo_sum[14] = mass_p;
    foo_sum[15] = mass_s;

    for (int i = 0; i <= 2; ++i) {
      foo_sum[i+16] = com_p[i];
      foo_sum[i+19] = com_s[i];
      foo_sum[i+22] = vel_p[i];
      foo_sum[i+25] = vel_s[i];
    }

    ParallelDescriptor::ReduceRealSum(foo_sum, nfoo_sum);

    for (int i = 0; i <= 6; ++i) {
      vol_p[i] = foo_sum[i  ];
      vol_s[i] = foo_sum[i+7];
    }

    mass_p = foo_sum[14];
    mass_s = foo_sum[15];

    for (int i = 0; i <= 2; ++i) {
      com_p[i] = foo_sum[i+16];
      com_s[i] = foo_sum[i+19];
      vel_p[i] = foo_sum[i+22];
      vel_s[i] = foo_sum[i+25];
    }

    if (mass_p > 0.0 && dt > 0.0)
      mdot_p = (mass_p - old_mass_p) / dt;
    else
      mdot_p = 0.0;

    if (mass_s > 0.0 && dt > 0.0)
      mdot_s = (mass_s - old_mass_s) / dt;
    else
      mdot_s = 0.0;

    // Max reductions

    const int nfoo_max = 3;

    Real foo_max[3];

    foo_max[0] = T_curr_max;
    foo_max[1] = rho_curr_max;
    foo_max[2] = ts_te_curr_max;

    ParallelDescriptor::ReduceRealSum(foo_max, nfoo_max);

    T_curr_max     = foo_max[0];
    rho_curr_max   = foo_max[1];
    ts_te_curr_max = foo_max[2];

    T_global_max     = std::max(T_global_max, T_curr_max);
    rho_global_max   = std::max(rho_global_max, rho_curr_max);
    ts_te_global_max = std::max(ts_te_global_max, ts_te_curr_max);

    // Send extrema data to Fortran

    set_extrema(&T_global_max, &rho_global_max, &ts_te_global_max);

    // Compute effective WD radii

    for (int i = 0; i <= 6; ++i) {

        rad_p[i] = std::pow(vol_p[i] * 3.0 / 4.0 / M_PI, 1.0/3.0);
        rad_s[i] = std::pow(vol_s[i] * 3.0 / 4.0 / M_PI, 1.0/3.0);

    }

     // Complete calculations for center of mass quantities

    for ( int i = 0; i < 3; i++ ) {

      if ( mass_p > 0.0 ) {
        com_p[i] = com_p[i] / mass_p;
        vel_p[i] = vel_p[i] / mass_p;
      }

      if ( mass_s > 0.0 ) {
        com_s[i] = com_s[i] / mass_s;
        vel_s[i] = vel_s[i] / mass_s;
      }

    }

    // Free-fall timescale ~ 1 / sqrt(G * rho_avg}

    Real Gconst;

    get_grav_const(&Gconst);

    if (mass_p > 0.0 && vol_p[2] > 0.0) {
      rho_avg_p = mass_p / vol_p[2];
      t_ff_p = sqrt(3.0 * M_PI / (32.0 * Gconst * rho_avg_p));
    }

    if (mass_s > 0.0 && vol_s[2] > 0.0) {
      rho_avg_s = mass_s / vol_s[2];
      t_ff_s = sqrt(3.0 * M_PI / (32.0 * Gconst * rho_avg_s));
    }

    // Send this updated information back to the Fortran probdata module

    set_star_data(com_p, com_s, vel_p, vel_s, &mass_p, &mass_s, &t_ff_p, &t_ff_s);

    // If we are doing problem 3, which has an initial relaxation step,
    // determine whether the criterion for terminating the relaxation
    // has been satisfied.
    // Note that at present the following code is only done on the
    // coarse grid but if we wanted more accuracy we could do a loop
    // over levels as above.

    Real L1[3] = { -1.0e200 };
    Real L2[3] = { -1.0e200 };
    Real L3[3] = { -1.0e200 };

    if (problem == 3 && !relaxation_is_done && mass_p > 0.0 && mass_s > 0.0)
    {

      // First, calculate the location of the L1 Lagrange point.

      get_lagrange_points(mass_p, mass_s, com_p, com_s, L1, L2, L3);

      // Then, figure out the effective potential corresponding to that
      // Lagrange point.

      Real potential = 0.0;

      MultiFab* mfphieff = derive("phiEff", time, 0);

#ifdef _OPENMP
#pragma omp parallel reduction(+:potential)
#endif
      for (MFIter mfi(*mfphieff,true); mfi.isValid(); ++mfi) {

	const Box& box = mfi.tilebox();

	const int* lo  = box.loVect();
	const int* hi  = box.hiVect();

	get_critical_roche_potential(BL_TO_FORTRAN_3D((*mfphieff)[mfi]),
				     lo, hi, L1, &potential);

      }

      ParallelDescriptor::ReduceRealSum(potential);

      // Now cycle through the grids and determine if any zones
      // have crossed the density threshold outside the critical surface.

      MultiFab& S_new = get_new_data(State_Type);

      int is_done = 0;

#ifdef _OPENMP
#pragma omp parallel reduction(+:is_done)
#endif
      for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi) {

	  const Box& box  = mfi.tilebox();

	  const int* lo   = box.loVect();
	  const int* hi   = box.hiVect();

	  check_relaxation(BL_TO_FORTRAN_3D(S_new[mfi]),
			   BL_TO_FORTRAN_3D((*mfphieff)[mfi]),
			   ARLIM_3D(lo),ARLIM_3D(hi),
			   &potential,&is_done);

      }

      ParallelDescriptor::ReduceIntSum(is_done);

      if (is_done > 0)
	relaxation_is_done = true;

      delete mfphieff;

      if (relaxation_is_done) {

	turn_off_relaxation(&time);

      }

      // Given the new estimate of the velocities relative to the old,
      // estimate the acceleration of each star. This is then used to
      // construct a new estimate for the rotation frequency.

      set_stellar_acceleration(&time, &dt, old_com_p, old_com_s, old_vel_p, old_vel_s);

    }

    // Some of the problems might have stopping conditions that depend on
    // the state of the simulation; those are checked here.

    // For the collision problem, we know we are done when the total energy
    // is positive (indicating that we have become unbound due to nuclear
    // energy release) and when it is decreasing in magnitude (indicating
    // all of the excitement is done and fluid is now just streaming off
    // the grid). We don't need to be super accurate for this, so let's check
    // on the coarse grid only. It is possible that a collision could not
    // generate enough energy to become unbound, so possibly this criterion
    // should be expanded in the future to cover that case.

    if (problem == 0) {

      Real rho_E_old = 0.0;
      Real rho_E_new = 0.0;

      Real rho_phi_old = 0.0;
      Real rho_phi_new = 0.0;

      // Note that we'll define the total energy using only
      // gas energy + gravitational. Rotation is never on
      // for the collision problem so we can ignore it.

      Real E_tot_old = 0.0;
      Real E_tot_new = 0.0;

      Real prevTime  = state[State_Type].prevTime();
      Real curTime   = state[State_Type].curTime();

      bool local_flag = true;

      rho_E_old += volWgtSum("rho_E", prevTime, local_flag);
      rho_E_new += volWgtSum("rho_E", curTime,  local_flag);

#ifdef GRAVITY
      if (do_grav) {
        rho_phi_old += volProductSum("density", "phiGrav", prevTime, local_flag);
        rho_phi_new += volProductSum("density", "phiGrav", curTime,  local_flag);
      }
#endif

      E_tot_old = rho_E_old - 0.5 * rho_phi_old;
      E_tot_new = rho_E_new - 0.5 * rho_phi_new;

      Array<Real> foo(2);
      foo[0] = E_tot_old;
      foo[1] = E_tot_new;

      ParallelDescriptor::ReduceRealSum(foo.dataPtr(), 2);

      E_tot_old = foo[0];
      E_tot_new = foo[1];

      if (E_tot_new > 0.0 && E_tot_new < E_tot_old) {

	jobDoneStatus = 1;

	set_job_status(&jobDoneStatus);

      }

    } else if (problem == 1) {

      // We can work out the stopping time using the formula
      // t_freefall = rotational_period / (4 * sqrt(2)).
      // We'll stop 90% of the way there because that's about
      // when the stars start coming into contact, and the
      // assumption of spherically symmetric stars breaks down.

      Real stopping_time = 0.90 * rotational_period / (4.0 * std::sqrt(2));

      if (time >= stopping_time) {

	jobDoneStatus = 1;

	set_job_status(&jobDoneStatus);

      }

    }

    // Is the job done? If so, signal this to BoxLib.

    get_job_status(&jobDoneStatus);

    if (jobDoneStatus == 1)
      stopJob();

}
#endif



//
// This function computes the center-of-mass locations (and velocities of the center of masses)
// of the primary and secondary white dwarfs.
// First, we predict where the center of mass locations will be
// for these stars based on Keplerian circular orbits.
// Then, we do a location-weighted sum that considers, for each respective star,
// only the zones that are closer to one orbit location than the other.
//

void
Castro::wdCOM (Real time, Real& mass_p, Real& mass_s, Real* com_p, Real* com_s, Real* vel_p, Real* vel_s, bool local)
{
    BL_PROFILE("Castro::wdCOM()");

    BL_ASSERT(level == 0);

    Real com_p_x = 0.0;
    Real com_p_y = 0.0;
    Real com_p_z = 0.0;
    Real com_s_x = 0.0;
    Real com_s_y = 0.0;
    Real com_s_z = 0.0;
    Real vel_p_x = 0.0;
    Real vel_p_y = 0.0;
    Real vel_p_z = 0.0;
    Real vel_s_x = 0.0;
    Real vel_s_y = 0.0;
    Real vel_s_z = 0.0;
    Real mp      = 0.0;
    Real ms      = 0.0;

    for (int lev = 0; lev <= parent->finestLevel(); lev++) {

      set_amr_info(lev, -1, -1, -1.0, -1.0);

      Castro& c_lev = getLevel(lev);

      const Real* dx = c_lev.geom.CellSize();

      // Density and momenta

      MultiFab* mfrho  = c_lev.derive("density",time,0);
      MultiFab* mfxmom = c_lev.derive("xmom",time,0);
      MultiFab* mfymom = c_lev.derive("ymom",time,0);
      MultiFab* mfzmom = c_lev.derive("zmom",time,0);

      // Effective potentials of the primary and secondary

      MultiFab* mfphip = c_lev.derive("phiEffPM_P", time, 0);
      MultiFab* mfphis = c_lev.derive("phiEffPM_S", time, 0);

      BL_ASSERT(mfrho  != 0);
      BL_ASSERT(mfxmom != 0);
      BL_ASSERT(mfymom != 0);
      BL_ASSERT(mfzmom != 0);
      BL_ASSERT(mfphip != 0);
      BL_ASSERT(mfphis != 0);

      if (lev < parent->finestLevel())
      {
          const MultiFab* mask = c_lev.getLevel(lev+1).build_fine_mask();

	  MultiFab::Multiply(*mfrho,  *mask, 0, 0, 1, 0);
	  MultiFab::Multiply(*mfxmom, *mask, 0, 0, 1, 0);
	  MultiFab::Multiply(*mfymom, *mask, 0, 0, 1, 0);
	  MultiFab::Multiply(*mfzmom, *mask, 0, 0, 1, 0);
	  MultiFab::Multiply(*mfphip, *mask, 0, 0, 1, 0);
	  MultiFab::Multiply(*mfphis, *mask, 0, 0, 1, 0);
      }

#ifdef _OPENMP
#pragma omp parallel reduction(+:com_p_x,com_p_y,com_p_z,com_s_x,com_s_y,com_s_z) \
                     reduction(+:vel_p_x,vel_p_y,vel_p_z,vel_s_x,vel_s_y,vel_s_z) \
                     reduction(+:mp, ms)
#endif
      for (MFIter mfi(*mfrho,true); mfi.isValid(); ++mfi)
	{
	  FArrayBox& fabrho  = (*mfrho )[mfi];
	  FArrayBox& fabxmom = (*mfxmom)[mfi];
	  FArrayBox& fabymom = (*mfymom)[mfi];
	  FArrayBox& fabzmom = (*mfzmom)[mfi];
	  FArrayBox& fabphip = (*mfphip)[mfi];
	  FArrayBox& fabphis = (*mfphis)[mfi];
	  FArrayBox& vol     = c_lev.volume[mfi];

	  const Box& box  = mfi.tilebox();
	  const int* lo   = box.loVect();
	  const int* hi   = box.hiVect();

	  wdcom(BL_TO_FORTRAN_3D(fabrho),
		BL_TO_FORTRAN_3D(fabxmom),
		BL_TO_FORTRAN_3D(fabymom),
		BL_TO_FORTRAN_3D(fabzmom),
		BL_TO_FORTRAN_3D(fabphip),
		BL_TO_FORTRAN_3D(fabphis),
		BL_TO_FORTRAN_3D(vol),
		ARLIM_3D(lo),ARLIM_3D(hi),
		ZFILL(dx),&time,
		&com_p_x, &com_p_y, &com_p_z,
		&com_s_x, &com_s_y, &com_s_z,
		&vel_p_x, &vel_p_y, &vel_p_z,
		&vel_s_x, &vel_s_y, &vel_s_z,
		&mp, &ms);
	}

      delete mfrho;
      delete mfxmom;
      delete mfymom;
      delete mfzmom;
      delete mfphip;
      delete mfphis;

    }

    com_p[0] = com_p_x;
    com_p[1] = com_p_y;
    com_p[2] = com_p_z;
    com_s[0] = com_s_x;
    com_s[1] = com_s_y;
    com_s[2] = com_s_z;
    vel_p[0] = vel_p_x;
    vel_p[1] = vel_p_y;
    vel_p[2] = vel_p_z;
    vel_s[0] = vel_s_x;
    vel_s[1] = vel_s_y;
    vel_s[2] = vel_s_z;
    mass_p   = mp;
    mass_s   = ms;

    if (!local) {

      const int nfoo = 14;
      Real foo[nfoo] = { 0.0 };

      foo[0] = mass_p;
      foo[1] = mass_s;

      for (int i = 0; i <=2; i++) {
	foo[i+2 ] = com_p[i];
	foo[i+5 ] = com_s[i];
	foo[i+8 ] = vel_p[i];
	foo[i+11] = vel_s[i];
      }

      ParallelDescriptor::ReduceRealSum(foo, nfoo);

      mass_p = foo[0];
      mass_s = foo[1];

      for (int i = 0; i <=2; i++) {
	com_p[i] = foo[i+2 ];
	com_s[i] = foo[i+5 ];
	vel_p[i] = foo[i+8 ];
	vel_s[i] = foo[i+11];
      }

    }

    set_amr_info(level, -1, -1, -1.0, -1.0);

}



// This function uses the known center of mass of the two white dwarfs,
// and given a density cutoff, computes the total volume of all zones
// whose density is greater or equal to that density cutoff.
// We also impose a distance requirement so that we only look 
// at zones that are within twice the original radius of the white dwarf.

void Castro::volInBoundary (Real time, Real& vol_p, Real& vol_s, Real rho_cutoff, bool local)
{
    BL_PROFILE("Castro::volInBoundary()");

    BL_ASSERT(level == 0);

    vol_p = 0.0;
    vol_s = 0.0;

    for (int lev = 0; lev <= parent->finestLevel(); lev++) {

      set_amr_info(lev, -1, -1, -1.0, -1.0);

      Castro& c_lev = getLevel(lev);

      const Real* dx   = c_lev.geom.CellSize();
      MultiFab*   mf   = c_lev.derive("density",time,0);

      // Effective potentials of the primary and secondary

      MultiFab* mfphip = c_lev.derive("phiEffPM_P", time, 0);
      MultiFab* mfphis = c_lev.derive("phiEffPM_S", time, 0);

      BL_ASSERT(mf != 0);
      BL_ASSERT(mfphip != 0);
      BL_ASSERT(mfphis != 0);

      if (lev < parent->finestLevel())
      {
	  const MultiFab* mask = c_lev.getLevel(lev+1).build_fine_mask();
	  MultiFab::Multiply(*mf, *mask, 0, 0, 1, 0);
	  MultiFab::Multiply(*mfphip, *mask, 0, 0, 1, 0);
	  MultiFab::Multiply(*mfphis, *mask, 0, 0, 1, 0);
      }

      Real vp = 0.0;
      Real vs = 0.0;

#ifdef _OPENMP
#pragma omp parallel reduction(+:vp,vs)
#endif
      for (MFIter mfi(*mf,true); mfi.isValid(); ++mfi)
      {
          FArrayBox& fab     = (*mf)[mfi];
	  FArrayBox& fabphip = (*mfphip)[mfi];
	  FArrayBox& fabphis = (*mfphis)[mfi];
	  FArrayBox& vol     = c_lev.volume[mfi];

	  Real sp = 0.0;
	  Real ss = 0.0;

	  const Box& box  = mfi.tilebox();
	  const int* lo   = box.loVect();
	  const int* hi   = box.hiVect();

	  ca_volumeindensityboundary(BL_TO_FORTRAN_3D(fab),
		                     BL_TO_FORTRAN_3D(fabphip),
				     BL_TO_FORTRAN_3D(fabphis),
				     BL_TO_FORTRAN_3D(vol),
				     ARLIM_3D(lo),ARLIM_3D(hi),
				     ZFILL(dx),&sp,&ss,&rho_cutoff);
	  vp += sp;
	  vs += ss;
      }

      delete mf;
      delete mfphis;
      delete mfphip;

      vol_p += vp;
      vol_s += vs;

    }

    if (!local) {

      const int nfoo = 2;
      Real foo[nfoo] = {vol_p, vol_s};

      ParallelDescriptor::ReduceRealSum(foo, nfoo);

      vol_p = foo[0];
      vol_s = foo[1];

    }

    set_amr_info(level, -1, -1, -1.0, -1.0);

}



//
// Calculate the gravitational wave signal.
//

void
Castro::gwstrain (Real time,
		  Real& h_plus_1, Real& h_cross_1,
		  Real& h_plus_2, Real& h_cross_2,
		  Real& h_plus_3, Real& h_cross_3,
		  bool local) {

    BL_PROFILE("Castro::gwstrain()");

    const Real* dx       = geom.CellSize();

    MultiFab* mfrho   = derive("density",time,0);
    MultiFab* mfxmom  = derive("xmom",time,0);
    MultiFab* mfymom  = derive("ymom",time,0);
    MultiFab* mfzmom  = derive("zmom",time,0);
    MultiFab* mfgravx = derive("grav_x",time,0);
    MultiFab* mfgravy = derive("grav_y",time,0);
    MultiFab* mfgravz = derive("grav_z",time,0);

    BL_ASSERT(mfrho   != 0);
    BL_ASSERT(mfxmom  != 0);
    BL_ASSERT(mfymom  != 0);
    BL_ASSERT(mfzmom  != 0);
    BL_ASSERT(mfgravx != 0);
    BL_ASSERT(mfgravy != 0);
    BL_ASSERT(mfgravz != 0);

    if (level < parent->finestLevel())
    {
	const MultiFab* mask = getLevel(level+1).build_fine_mask();

	MultiFab::Multiply(*mfrho,   *mask, 0, 0, 1, 0);
	MultiFab::Multiply(*mfxmom,  *mask, 0, 0, 1, 0);
	MultiFab::Multiply(*mfymom,  *mask, 0, 0, 1, 0);
	MultiFab::Multiply(*mfzmom,  *mask, 0, 0, 1, 0);
	MultiFab::Multiply(*mfgravx, *mask, 0, 0, 1, 0);
	MultiFab::Multiply(*mfgravy, *mask, 0, 0, 1, 0);
	MultiFab::Multiply(*mfgravz, *mask, 0, 0, 1, 0);
    }

    // Qtt stores the second time derivative of the quadrupole moment.
    // We calculate it directly rather than computing the quadrupole moment
    // and differentiating it in time, because the latter method is less accurate
    // and requires the state at other timesteps. See, e.g., Equation 5 of 
    // Loren-Aguilar et al. 2005.

    // It is a 3x3 rank-2 tensor, but BoxLib expects IntVect() to use BL_SPACEDIM
    // dimensions, so we add a redundant third index in 3D.

    Box bx( IntVect(D_DECL(0, 0, 0)), IntVect(D_DECL(2, 2, 0)) );

    FArrayBox Qtt(bx);

    Qtt.setVal(0.0);

#ifdef _OPENMP
    int nthreads = omp_get_max_threads();
    PArray<FArrayBox> priv_Qtt(nthreads, PArrayManage);
    for (int i=0; i<nthreads; i++) {
	priv_Qtt.set(i, new FArrayBox(bx));
    }
#pragma omp parallel
#endif
    {
#ifdef _OPENMP
	int tid = omp_get_thread_num();
	priv_Qtt[tid].setVal(0.0);
#endif
	for (MFIter mfi(*mfrho,true); mfi.isValid(); ++mfi) {

	    const Box& box  = mfi.tilebox();
	    const int* lo   = box.loVect();
	    const int* hi   = box.hiVect();

	    quadrupole_tensor_double_dot(BL_TO_FORTRAN_3D((*mfrho)[mfi]),
					 BL_TO_FORTRAN_3D((*mfxmom)[mfi]),
					 BL_TO_FORTRAN_3D((*mfymom)[mfi]),
					 BL_TO_FORTRAN_3D((*mfzmom)[mfi]),
					 BL_TO_FORTRAN_3D((*mfgravx)[mfi]),
					 BL_TO_FORTRAN_3D((*mfgravy)[mfi]),
					 BL_TO_FORTRAN_3D((*mfgravz)[mfi]),
					 BL_TO_FORTRAN_3D(volume[mfi]),
					 ARLIM_3D(lo),ARLIM_3D(hi),ZFILL(dx),&time,
#ifdef _OPENMP
					 priv_Qtt[tid].dataPtr());
#else
	                                 Qtt.dataPtr());
#endif
        }
    }

    delete mfrho;
    delete mfxmom;
    delete mfymom;
    delete mfzmom;
    delete mfgravx;
    delete mfgravy;
    delete mfgravz;

    // Do an OpenMP reduction on the tensor.

#ifdef _OPENMP
        int n = bx.numPts();
	Real* p = Qtt.dataPtr();
#pragma omp barrier
#pragma omp for nowait
	for (int i=0; i<n; ++i)
	{
	    for (int it=0; it<nthreads; it++) {
		const Real* pq = priv_Qtt[it].dataPtr();
		p[i] += pq[i];
	    }
	}
#endif

    // Now, do a global reduce over all processes.

    if (!local)
      ParallelDescriptor::ReduceRealSum(Qtt.dataPtr(),bx.numPts());

    // Now that we have the second time derivative of the quadrupole 
    // tensor, we can calculate the transverse-trace gauge strain tensor.

    gw_strain_tensor(&h_plus_1, &h_cross_1,
		     &h_plus_2, &h_cross_2,
		     &h_plus_3, &h_cross_3, 
		     Qtt.dataPtr(), &time);

}



// Computes standard dot-product of two three-vectors.

Real Castro::dot_product(const Real a[], const Real b[]) {

  Real c = 0.0;

  c = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];

  return c;

}



// Computes standard cross-product of two three-vectors.

void Castro::cross_product(const Real a[], const Real b[], Real c[]) {

  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];

}


// Computes norm of a three-vector.

Real Castro::norm(const Real a[]) {

  Real n = 0.0;

  n = sqrt( dot_product(a, a) );

  return n;

}



#ifdef do_problem_post_init

void Castro::problem_post_init() {

  // Get the problem number fom Fortran.

  get_problem_number(&problem);

  // Execute the post timestep diagnostics here,
  // so that the results at t = 0 and later are smooth.

  T_global_max     = 0.0;
  rho_global_max   = 0.0;
  ts_te_global_max = 0.0;

  problem_post_timestep();

}

#endif



#ifdef do_problem_post_restart

void Castro::problem_post_restart() {

  // Get the problem number from Fortran.

  get_problem_number(&problem);

}

#endif



void Castro::writeGitHashes(std::ostream& log) {

  const char* castro_hash       = buildInfoGetGitHash(1);
  const char* boxlib_hash       = buildInfoGetGitHash(2);
  const char* microphysics_hash = buildInfoGetGitHash(3);
  const char* wdmerger_hash     = buildInfoGetBuildGitHash();

  log << "# Castro       git hash: " << castro_hash       << std::endl;
  log << "# BoxLib       git hash: " << boxlib_hash       << std::endl;
  log << "# Microphysics git hash: " << microphysics_hash << std::endl;
  log << "# wdmerger     git hash: " << wdmerger_hash     << std::endl;

}
