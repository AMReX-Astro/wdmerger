#include "Castro.H"
#include "Castro_F.H"

#include "Gravity.H"

#include "ParmParse.H"

int Castro::do_relax = 0;

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

	BL_FORT_PROC_CALL(WDCOM,wdcom)
            (BL_TO_FORTRAN(fabrho),
	     BL_TO_FORTRAN(fabxmom),
	     BL_TO_FORTRAN(fabymom),
	     BL_TO_FORTRAN(fabzmom),
	     lo,hi,dx,&time,
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

	BL_FORT_PROC_CALL(CA_VOLUMEINDENSITYBOUNDARY,ca_volumeindensityboundary)
	                  (BL_TO_FORTRAN(fab),lo,hi,dx,&sp,&ss,&rho_cutoff);
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
Castro::gwstrain (Real time, Real& h_plus_rot, Real& h_cross_rot,
		  Real& h_plus_star, Real& h_cross_star,
		  Real& h_plus_motion, Real& h_cross_motion) {

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

    // It is a 3x3 rank-2 tensor, but BoxLib expects our boxes to be the same 
    // dimensionality as the problem, so we add a redundant third index.

    Box bx( IntVect(0, 0, 0), IntVect(2, 2, 0) );

    FArrayBox Qtt(bx);

    Qtt.setVal(0.0);

#ifdef _OPENMP
    int nthreads = omp_get_max_threads();
    PArray<FArrayBox> priv_Qtt(nthreads);
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

	    BL_FORT_PROC_CALL(QUADRUPOLE_TENSOR_DOUBLE_DOT,quadrupole_tensor_double_dot)
	        (BL_TO_FORTRAN((*mfrho)[mfi]),
		 BL_TO_FORTRAN((*mfxmom)[mfi]),
		 BL_TO_FORTRAN((*mfymom)[mfi]),
		 BL_TO_FORTRAN((*mfzmom)[mfi]),
		 BL_TO_FORTRAN((*mfgravx)[mfi]),
		 BL_TO_FORTRAN((*mfgravy)[mfi]),
		 BL_TO_FORTRAN((*mfgravz)[mfi]),
		 lo,hi,dx,&time,
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

    BL_FORT_PROC_CALL(GW_STRAIN_TENSOR,gw_strain_tensor)
        (&h_plus_rot, &h_cross_rot,
	 &h_plus_star, &h_cross_star,
	 &h_plus_motion, &h_plus_motion, 
	 Qtt.dataPtr(), &time);

}



#ifdef GRAVITY
#ifdef ROTATION
#ifdef do_problem_post_init

void Castro::problem_post_init() {

    if (level > 0)
        return;

    ParmParse pp("castro");

    pp.query("do_relax", do_relax);

    if (do_relax) {

      if (gravity->NoComposite() == 1) {
	std::cerr << "Initial relaxation requires the use of multilevel gravity solves. Set gravity.no_composite = 0." << std::endl;
	BoxLib::Error();
      }

      int relax_type = 0;

      BL_FORT_PROC_CALL(GET_RELAX_TYPE,get_relax_type)(relax_type);

      if (relax_type == 1)
	scf_relaxation();

    }

}

#endif
#endif
#endif
