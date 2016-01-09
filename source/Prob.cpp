#include "Castro.H"
#include "Castro_F.H"

#include "Gravity.H"

#include "ParmParse.H"

#ifdef do_problem_post_timestep
void
Castro::problem_post_timestep()
{

    if (level != 0) return;

    int finest_level = parent->finestLevel();
    Real time = state[State_Type].curTime();
    Real dt = parent->dtLevel(0);

    if (time == 0.0) dt = 0.0; // dtLevel returns the next timestep for t = 0, so overwrite    
    
    Real mass_p       = 0.0;
    Real mass_s       = 0.0;

    Real lev_mass_p   = 0.0;
    Real lev_mass_s   = 0.0;

    Real com_p[3]     = { 0.0 };
    Real com_s[3]     = { 0.0 };

    Real lev_com_p[3] = { 0.0 };
    Real lev_com_s[3] = { 0.0 };

    Real vel_p[3]     = { 0.0 };
    Real vel_s[3]     = { 0.0 };

    Real lev_vel_p[3] = { 0.0 };
    Real lev_vel_s[3] = { 0.0 };



    // Get the current stellar data
    get_star_data(com_p, com_s, vel_p, vel_s, &mass_p, &mass_s);
    
    // Update the problem center using the system bulk velocity
    update_center(&time);

    for ( int i = 0; i < 3; i++ ) {
      com_p[i] += vel_p[i] * dt;
      com_s[i] += vel_s[i] * dt;
    }

    // Now send this first estimate of the COM to Fortran, and then re-calculate
    // a more accurate result using it as a starting point.
    
    set_star_data(com_p, com_s, vel_p, vel_s, &mass_p, &mass_s);

    mass_p = 0.0;
    mass_s = 0.0;
    
    for ( int i = 0; i < 3; i++ ) {
      com_p[i] = 0.0;
      com_s[i] = 0.0;
      vel_p[i] = 0.0;
      vel_s[i] = 0.0;
    }

    for (int lev = 0; lev <= finest_level; lev++)
    {

      set_amr_info(lev, -1, -1, -1.0, -1.0);
      
      getLevel(lev).wdCOM(time, lev_mass_p, lev_mass_s, lev_com_p, lev_com_s, lev_vel_p, lev_vel_s);

      mass_p += lev_mass_p;
      mass_s += lev_mass_s;

      for ( int i = 0; i < 3; i++ ) {
	com_p[i] += lev_com_p[i];
	com_s[i] += lev_com_s[i];
	vel_p[i] += lev_vel_p[i];
	vel_s[i] += lev_vel_s[i];
      }

    }

    set_amr_info(level, -1, -1, -1.0, -1.0);    
    
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

    // Send this updated information back to the Fortran probdata module

    set_star_data(com_p, com_s, vel_p, vel_s, &mass_p, &mass_s);

    // If we are doing an initial relaxation step, determine whether the 
    // criterion for terminating the relaxation has been satisfied.

    // First, calculate the location of the L1 Lagrange point.

    Real L1[3] = { -1.0e200 };
    Real L2[3] = { -1.0e200 };
    Real L3[3] = { -1.0e200 };

    if (mass_p > 0.0 && mass_s > 0.0)
    {

      get_lagrange_points(mass_p, mass_s, com_p, com_s, L1, L2, L3);
    
      // Now cycle through the grids and determine if the L1
      // point has reached the density threshold.
    
      int relaxation_is_done = 0;

      MultiFab& S_new = get_new_data(State_Type);
    
#ifdef _OPENMP
#pragma omp parallel reduction(+:relaxation_is_done)
#endif    
      for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
	{
	  const Box& box  = mfi.tilebox();

	  const int* lo   = box.loVect();
	  const int* hi   = box.hiVect();

	  check_relaxation(BL_TO_FORTRAN_3D(S_new[mfi]),
			   ARLIM_3D(lo),ARLIM_3D(hi),
			   L1,relaxation_is_done);
	
	}

      // If any of the grids returned that it has the L1 point and
      // the density has passed the cutoff, then disable the initial relaxation.
    
      if (relaxation_is_done > 0)
	turn_off_relaxation();

    }
    
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
Castro::wdCOM (Real time, Real& mass_p, Real& mass_s, Real* com_p, Real* com_s, Real* vel_p, Real* vel_s)
{
    BL_PROFILE("Castro::wdCOM()");

    const Real* dx       = geom.CellSize();

    MultiFab*   mfrho    = derive("density",time,0);
    MultiFab*   mfxmom   = derive("xmom",time,0);
    MultiFab*   mfymom   = derive("ymom",time,0);
    MultiFab*   mfzmom   = derive("zmom",time,0);
    
    BL_ASSERT(mfrho  != 0);
    BL_ASSERT(mfxmom != 0);
    BL_ASSERT(mfymom != 0);
    BL_ASSERT(mfzmom != 0);

    if (level < parent->finestLevel())
    {
	const MultiFab* mask = getLevel(level+1).build_fine_mask();

	MultiFab::Multiply(*mfrho,  *mask, 0, 0, 1, 0);
	MultiFab::Multiply(*mfxmom, *mask, 0, 0, 1, 0);
	MultiFab::Multiply(*mfymom, *mask, 0, 0, 1, 0);
	MultiFab::Multiply(*mfzmom, *mask, 0, 0, 1, 0);
    }

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
	
        const Box& box  = mfi.tilebox();
        const int* lo   = box.loVect();
        const int* hi   = box.hiVect();

	wdcom(BL_TO_FORTRAN_3D(fabrho),
	      BL_TO_FORTRAN_3D(fabxmom),
	      BL_TO_FORTRAN_3D(fabymom),
	      BL_TO_FORTRAN_3D(fabzmom),
	      BL_TO_FORTRAN_3D(volume[mfi]),
	      ARLIM_3D(lo),ARLIM_3D(hi),ZFILL(dx),&time,
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

    ParallelDescriptor::ReduceRealSum(com_p,3);
    ParallelDescriptor::ReduceRealSum(com_s,3);
    ParallelDescriptor::ReduceRealSum(vel_p,3);
    ParallelDescriptor::ReduceRealSum(vel_s,3);

    ParallelDescriptor::ReduceRealSum(mass_p);
    ParallelDescriptor::ReduceRealSum(mass_s);

}



// This function uses the known center of mass of the two white dwarfs,
// and given a density cutoff, computes the total volume of all zones
// whose density is greater or equal to that density cutoff.
// We also impose a distance requirement so that we only look 
// at zones that are within twice the original radius of the white dwarf.

void Castro::volInBoundary (Real               time,
			    Real&              vol_p,
			    Real&              vol_s,
                            Real               rho_cutoff)
{
    BL_PROFILE("Castro::volInBoundary()");

    const Real* dx      = geom.CellSize();
    MultiFab*   mf      = derive("density",time,0);

    BL_ASSERT(mf != 0);

    if (level < parent->finestLevel())
    {
	const MultiFab* mask = getLevel(level+1).build_fine_mask();
	MultiFab::Multiply(*mf, *mask, 0, 0, 1, 0);
    }

    Real vp = 0.0;
    Real vs = 0.0;

#ifdef _OPENMP
#pragma omp parallel reduction(+:vp,vs)
#endif    
    for (MFIter mfi(*mf,true); mfi.isValid(); ++mfi)
    {
        FArrayBox& fab = (*mf)[mfi];

	Real sp = 0.0;
	Real ss = 0.0;
        const Box& box  = mfi.tilebox();
        const int* lo   = box.loVect();
        const int* hi   = box.hiVect();

	ca_volumeindensityboundary(BL_TO_FORTRAN_3D(fab),
				   BL_TO_FORTRAN_3D(volume[mfi]),
				   ARLIM_3D(lo),ARLIM_3D(hi),
				   ZFILL(dx),&sp,&ss,&rho_cutoff);
        vp += sp;
	vs += ss;
    }

    delete mf;

    vol_p = vp;
    vol_s = vs;

    ParallelDescriptor::ReduceRealSum(vol_p);
    ParallelDescriptor::ReduceRealSum(vol_s);

}



//
// Calculate the gravitational wave signal.
//

void
Castro::gwstrain (Real time,
		  Real& h_plus_1, Real& h_cross_1,
		  Real& h_plus_2, Real& h_cross_2,
		  Real& h_plus_3, Real& h_cross_3) {

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



#ifdef GRAVITY
#ifdef ROTATION
#ifdef do_problem_post_init

void Castro::problem_post_init() {

    if (level > 0)
        return;

    // Determine whether we are doing a relaxation step (from the probin file)

    int do_scf_initial_models = 0;

    get_do_scf_initial_models(do_scf_initial_models);

    if (do_scf_initial_models) {

      if (gravity->NoComposite() == 1) {
	std::cerr << "Construction of SCF initial models requires the use of multilevel gravity solves. Set gravity.no_composite = 0." << std::endl;
	BoxLib::Error();
      }

      scf_relaxation();

    }

}

#endif
#endif
#endif
