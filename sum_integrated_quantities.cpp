#include <iomanip>

#include <Castro.H>
#include <Castro_F.H>
#include <Gravity.H>

void
Castro::sum_integrated_quantities ()
{
    int finest_level  = parent->finestLevel();
    Real time         = state[State_Type].curTime();
    Real mass         = 0.0;
    Real xmom         = 0.0;
    Real ymom         = 0.0;
    Real zmom         = 0.0;
    Real rho_E        = 0.0;
    Real rho_e        = 0.0;
    Real rho_phi      = 0.0;

    Real gravitational_energy = 0.0; 
    Real kinetic_energy       = 0.0; 
    Real internal_energy      = 0.0; 
    Real rotational_energy    = 0.0; 
    Real total_energy         = 0.0; 

    Real angular_momentum[3]  = {0.0, 0.0, 0.0};
    Real moment_of_inertia[3] = {0.0, 0.0, 0.0};
    Real m_r_squared[3]       = {0.0, 0.0, 0.0};

    Real omega[3] = {0.0, 0.0, 2.0*3.1415926*rotational_frequency};

    Real mass_left    = 0.0;
    Real mass_right   = 0.0;

    Real com_xloc     = 0.0;
    Real com_xloc_l   = 0.0;
    Real com_xloc_r   = 0.0;
    Real com_xvel     = 0.0;
    Real com_xvel_l   = 0.0;
    Real com_xvel_r   = 0.0; 

    Real com_yloc     = 0.0;
    Real com_yloc_l   = 0.0;
    Real com_yloc_r   = 0.0;
    Real com_yvel     = 0.0;
    Real com_yvel_l   = 0.0;
    Real com_yvel_r   = 0.0;

    Real com_zloc     = 0.0;
    Real com_zloc_l   = 0.0;
    Real com_zloc_r   = 0.0;
    Real com_zvel     = 0.0;
    Real com_zvel_l   = 0.0;
    Real com_zvel_r   = 0.0;

    int endOfStep     =   1;

    int datawidth     =  14;
    int dataprecision =   6;

    const int* eos_ptr = &endOfStep;

    BL_FORT_PROC_CALL(GET_EOS,get_eos)(eos_ptr);    
    
    for (int lev = 0; lev <= finest_level; lev++)
    {

      // Get the current level from Castro

      Castro& ca_lev = getLevel(lev);

      // Calculate total mass, momentum and energy of system.

      mass     += ca_lev.volWgtSum("density", time);

      xmom     += ca_lev.volWgtSum("xmom", time);
      ymom     += ca_lev.volWgtSum("ymom", time);
      zmom     += ca_lev.volWgtSum("zmom", time);

      rho_E    += ca_lev.volWgtSum("rho_E", time);
      rho_e    += ca_lev.volWgtSum("rho_e", time);
      
      rho_phi  += ca_lev.volProductSum("density", "phi", time);

      for ( int i = 0; i <= 2; i++ )
        m_r_squared[i] += ca_lev.locSquaredSum("density", time, i);
      
      for ( int i = 0; i <= 2; i++ ) {
        
        int index1 = (i+1) % 3; 
        int index2 = (i+2) % 3;
        std::string name1; 
        std::string name2;
        switch (i) {
          case 0 :
            name1 = "ymom"; name2 = "zmom"; break;
          case 1 :
            name1 = "zmom"; name2 = "xmom"; break;
          case 2 :
            name1 = "xmom"; name2 = "ymom"; break;
        }

        moment_of_inertia[i] += m_r_squared[index1] + m_r_squared[index2];

        angular_momentum[i]  += ca_lev.locWgtSum(name2, time, index1) - ca_lev.locWgtSum(name1, time, index2);
       
        // Add rotational source term to angular momentum and energy

        if ( do_rotation )
          for ( int i = 0; i <= 2; i++ ) {
            angular_momentum[i] += moment_of_inertia[i] * omega[i];
            rotational_energy   += (1.0/2.0) * moment_of_inertia[i] * omega[i] * omega[i];
          }
      }
  
 
      // Calculate center of mass quantities.

      mass_left   += ca_lev.volWgtSumOneSide("density", time, 0, 0);
      mass_right  += ca_lev.volWgtSumOneSide("density", time, 1, 0);

      com_xloc    += ca_lev.locWgtSum("density", time, 0);
      com_xloc_l  += ca_lev.locWgtSumOneSide("density", time, 0, 0, 0);
      com_xloc_r  += ca_lev.locWgtSumOneSide("density", time, 0, 1, 0);
      com_xvel_l  += ca_lev.volWgtSumOneSide("xmom",    time, 0, 0);
      com_xvel_r  += ca_lev.volWgtSumOneSide("xmom",    time, 1, 0);

      com_yloc    += ca_lev.locWgtSum("density", time, 1);
      com_yloc_l  += ca_lev.locWgtSumOneSide("density", time, 1, 0, 0);
      com_yloc_r  += ca_lev.locWgtSumOneSide("density", time, 1, 1, 0);
      com_yvel_l  += ca_lev.volWgtSumOneSide("ymom",    time, 0, 0);
      com_yvel_r  += ca_lev.volWgtSumOneSide("ymom",    time, 1, 0);

      com_zloc    += ca_lev.locWgtSum("density", time, 2);
      com_zloc_l  += ca_lev.locWgtSumOneSide("density", time, 2, 0, 0);
      com_zloc_r  += ca_lev.locWgtSumOneSide("density", time, 2, 1, 0);
      com_zvel_l  += ca_lev.volWgtSumOneSide("zmom",    time, 0, 0);
      com_zvel_r  += ca_lev.volWgtSumOneSide("zmom",    time, 1, 0);
    }
        
    // Complete calculations for COM quantities

    com_xloc     = com_xloc / mass;
    com_xloc_l   = com_xloc_l / mass_left;
    com_xloc_r   = com_xloc_r / mass_right;
    com_xvel_l   = com_xvel_l / mass_left;
    com_yvel_r   = com_xvel_r / mass_right;
    com_xvel     = xmom / mass;

    com_yloc     = com_yloc / mass;
    com_yloc_l   = com_yloc_l / mass_left;
    com_yloc_r   = com_yloc_r / mass_right;
    com_zvel_l   = com_yvel_l / mass_left;
    com_zvel_r   = com_yvel_r / mass_right;
    com_yvel     = ymom / mass;

    com_zloc     = com_zloc / mass;
    com_zloc_l   = com_zloc_l / mass_left;
    com_zloc_r   = com_zloc_r / mass_right;
    com_zvel_l   = com_zvel_l / mass_left;
    com_zvel_r   = com_zvel_r / mass_right;
    com_zvel     = zmom / mass;

    const Real* ml = &mass_left;
    const Real* mr = &mass_right;
    const Real* cxl = &com_xloc_l;
    const Real* cxr = &com_xloc_r;
    const Real* cyl = &com_yloc_l;
    const Real* cyr = &com_yloc_r;
    const Real* czl = &com_zloc_l;
    const Real* czr = &com_zloc_r;

    BL_FORT_PROC_CALL(COM_SAVE,com_save)
      (ml, mr, cxl, cxr, cyl, cyr, czl, czr);

    // Complete calculations for energy

    gravitational_energy = (-1.0/2.0) * rho_phi; // avoids double counting; CASTRO uses positive phi
    internal_energy = rho_e;
    kinetic_energy = rho_E - rho_e;
    total_energy = gravitational_energy + internal_energy + kinetic_energy + rotational_energy; 
    
    // Write data out to the log.

    std::ostream& data_log1 = parent->DataLog(0);

    if ( ParallelDescriptor::IOProcessor() && data_log1.good() )
    {

      // Write header row

      if (time == 0.0) {
        data_log1 << std::setw(datawidth) << "#    TIME    ";
        data_log1 << std::setw(datawidth) << " MASS        ";
        data_log1 << std::setw(datawidth) << " XMOM        ";
        data_log1 << std::setw(datawidth) << " YMOM        ";
        data_log1 << std::setw(datawidth) << " ZMOM        ";
        data_log1 << std::setw(datawidth) << " KIN. ENERGY ";
        data_log1 << std::setw(datawidth) << " GRAV. ENERGY";
        data_log1 << std::setw(datawidth) << " INT. ENERGY ";
        data_log1 << std::setw(datawidth) << " ROT. ENERGY ";
        data_log1 << std::setw(datawidth) << " TOTAL ENERGY";
        data_log1 << std::setw(datawidth) << " ANG. MOM. X ";
        data_log1 << std::setw(datawidth) << " ANG. MOM. Y ";
        data_log1 << std::setw(datawidth) << " ANG. MOM. Z ";
        data_log1 << std::setw(datawidth) << "  X COM       ";
        data_log1 << std::setw(datawidth) << "  Y COM       ";
        data_log1 << std::setw(datawidth) << "  Z COM       ";
        data_log1 << std::setw(datawidth) << "  LEFT MASS   ";
        data_log1 << std::setw(datawidth) << "  RIGHT MASS  ";
        data_log1 << std::setw(datawidth) << "  LEFT X COM  ";
        data_log1 << std::setw(datawidth) << "  RIGHT X COM ";
        data_log1 << std::setw(datawidth) << "  LEFT Y COM  ";
        data_log1 << std::setw(datawidth) << "  RIGHT Y COM ";
        data_log1 << std::setw(datawidth) << "  LEFT Z COM  ";
        data_log1 << std::setw(datawidth) << "  RIGHT Z COM ";
        data_log1 << std::setw(datawidth) << "  X VEL       ";
        data_log1 << std::setw(datawidth) << "  Y VEL       ";
        data_log1 << std::setw(datawidth) << "  Z VEL       ";
        data_log1 << std::setw(datawidth) << "  LEFT X VEL  ";
        data_log1 << std::setw(datawidth) << "  RIGHT X VEL ";
        data_log1 << std::setw(datawidth) << "  LEFT Y VEL  ";
        data_log1 << std::setw(datawidth) << "  RIGHT Y VEL ";
        data_log1 << std::setw(datawidth) << "  LEFT Z VEL  ";
        data_log1 << std::setw(datawidth) << "  RIGHT Z VEL ";
         
        data_log1 << std::endl;
      }

      // Write data for the present time

      if ( endOfStep || time == 0.0 ) {
	data_log1 << std::fixed;

	data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << time;

	data_log1 << std::scientific;

	data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << mass;
	data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << xmom;
	data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << ymom;
	data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << zmom;
	data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << kinetic_energy;
	data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << gravitational_energy;
        data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << internal_energy;
        data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << rotational_energy;
        data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << total_energy;
        data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << angular_momentum[0];
        data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << angular_momentum[1]; 
        data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << angular_momentum[2];
	data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << com_xloc;
	data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << com_yloc;
	data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << com_zloc;
	data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << mass_left;
	data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << mass_right;
	data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << com_xloc_l;
	data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << com_xloc_r;
	data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << com_yloc_l;
	data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << com_yloc_r;
	data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << com_zloc_l;
	data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << com_zloc_r;
	data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << com_xvel;
	data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << com_yvel;
	data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << com_zvel;
	data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << com_xvel_l;
	data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << com_xvel_r;
	data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << com_yvel_l;
	data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << com_yvel_r;
	data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << com_zvel_l;
	data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << com_zvel_r;
	
	data_log1 << std::endl;
      }

    }
}

Real
Castro::volWgtSum (const std::string& name,
                   Real               time)
{
    Real        sum     = 0.0;
    const Real* dx      = geom.CellSize();
    MultiFab*   mf      = derive(name,time,0);

    BL_ASSERT(mf != 0);

    BoxArray baf;

    if (level < parent->finestLevel())
    {
        baf = parent->boxArray(level+1);
        baf.coarsen(fine_ratio);
    }

    for (MFIter mfi(*mf); mfi.isValid(); ++mfi)
    {
        FArrayBox& fab = (*mf)[mfi];

        if (level < parent->finestLevel())
        {
            std::vector< std::pair<int,Box> > isects = baf.intersections(grids[mfi.index()]);

            for (int ii = 0; ii < isects.size(); ii++)
            {
                fab.setVal(0,isects[ii].second,0,fab.nComp());
            }
        }
        Real s = 0.0;
        const Box& box  = mfi.validbox();
        const int* lo   = box.loVect();
        const int* hi   = box.hiVect();
#if(BL_SPACEDIM < 3) 
        const Real* rad = radius[mfi.index()].dataPtr();
        int irlo        = lo[0]-radius_grow;
        int irhi        = hi[0]+radius_grow;
#endif

        //
        // Note that this routine will do a volume weighted sum of
        // whatever quantity is passed in, not strictly the "mass".
        //
#if(BL_SPACEDIM == 1) 
	BL_FORT_PROC_CALL(CA_SUMMASS,ca_summass)
            (BL_TO_FORTRAN(fab),lo,hi,dx,&s,rad,irlo,irhi);
#elif(BL_SPACEDIM == 2)
	BL_FORT_PROC_CALL(CA_SUMMASS,ca_summass)
            (BL_TO_FORTRAN(fab),lo,hi,dx,&s,rad,irlo,irhi);
#elif(BL_SPACEDIM == 3)
	BL_FORT_PROC_CALL(CA_SUMMASS,ca_summass)
            (BL_TO_FORTRAN(fab),lo,hi,dx,&s);
#endif
        sum += s;
    }

    delete mf;

    ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}

Real
Castro::volWgtSquaredSum (const std::string& name,
                          Real               time)
{
    Real        sum     = 0.0;
    const Real* dx      = geom.CellSize();
    MultiFab*   mf      = derive(name,time,0);

    BL_ASSERT(mf != 0);

    BoxArray baf;

    if (level < parent->finestLevel())
    {
        baf = parent->boxArray(level+1);
        baf.coarsen(fine_ratio);
    }

    for (MFIter mfi(*mf); mfi.isValid(); ++mfi)
    {
        FArrayBox& fab = (*mf)[mfi];

        if (level < parent->finestLevel())
        {
            std::vector< std::pair<int,Box> > isects = baf.intersections(grids[mfi.index()]);

            for (int ii = 0; ii < isects.size(); ii++)
            {
                fab.setVal(0,isects[ii].second,0,fab.nComp());
            }
        }
        Real s = 0.0;
        const Box& box  = mfi.validbox();
        const int* lo   = box.loVect();
        const int* hi   = box.hiVect();
#if(BL_SPACEDIM < 3) 
        const Real* rad = radius[mfi.index()].dataPtr();
        int irlo        = lo[0]-radius_grow;
        int irhi        = hi[0]+radius_grow;
#endif

        //
        // Note that this routine will do a volume weighted sum of
        // whatever quantity is passed in, not strictly the "mass".
        //
#if(BL_SPACEDIM == 1) 
	BL_FORT_PROC_CALL(CA_SUMSQUARED,ca_summass)
            (BL_TO_FORTRAN(fab),lo,hi,dx,&s,rad,irlo,irhi);
#elif(BL_SPACEDIM == 2)
	BL_FORT_PROC_CALL(CA_SUMSQUARED,ca_summass)
            (BL_TO_FORTRAN(fab),lo,hi,dx,&s,rad,irlo,irhi);
#elif(BL_SPACEDIM == 3)
	BL_FORT_PROC_CALL(CA_SUMSQUARED,ca_summass)
            (BL_TO_FORTRAN(fab),lo,hi,dx,&s);
#endif
        sum += s;
    }

    delete mf;

    ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}

Real
Castro::locWgtSum (const std::string& name,
                   Real               time,
                   int                idir)
{
    Real        sum     = 0.0;
    const Real* dx      = geom.CellSize();
    MultiFab*   mf      = derive(name,time,0);

    BL_ASSERT(mf != 0);

    BoxArray baf;

    if (level < parent->finestLevel())
    {
        baf = parent->boxArray(level+1);
        baf.coarsen(fine_ratio);
    }

    for (MFIter mfi(*mf); mfi.isValid(); ++mfi)
    {
        FArrayBox& fab = (*mf)[mfi];

        if (level < parent->finestLevel())
        {
            std::vector< std::pair<int,Box> > isects = baf.intersections(grids[mfi.index()]);

            for (int ii = 0; ii < isects.size(); ii++)
            {
                fab.setVal(0,isects[ii].second,0,fab.nComp());
            }
        }
        Real s = 0.0;
        const Box& box  = mfi.validbox();
        const int* lo   = box.loVect();
        const int* hi   = box.hiVect();
#if (BL_SPACEDIM < 3)
        const Real* rad = radius[mfi.index()].dataPtr();
+        int irlo        = lo[0]-radius_grow;
        int irhi        = hi[0]+radius_grow;
#endif

        //
        // Note that this routine will do a volume weighted sum of
        // whatever quantity is passed in, not strictly the "mass".
        //
#if (BL_SPACEDIM == 1) 
	BL_FORT_PROC_CALL(CA_SUMLOCMASS,ca_sumlocmass)
            (BL_TO_FORTRAN(fab),lo,hi,geom.ProbLo(),dx,&s,rad,irlo,irhi,idir);
#elif (BL_SPACEDIM == 2)
        int geom_flag = Geometry::IsRZ() ? 1 : 0;
        if (idir == 0 && geom_flag == 1) {
            s = 0.0;
        } else {
	   BL_FORT_PROC_CALL(CA_SUMLOCMASS,ca_sumlocmass)
               (BL_TO_FORTRAN(fab),lo,hi,geom.ProbLo(),dx,&s,rad,irlo,irhi,idir);
        }
#else
	BL_FORT_PROC_CALL(CA_SUMLOCMASS,ca_sumlocmass)
            (BL_TO_FORTRAN(fab),lo,hi,geom.ProbLo(),dx,&s,idir);
#endif
        sum += s;
    }

    delete mf;

    ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}

Real
Castro::volWgtSumOneSide (const std::string& name,
                          Real               time, 
                          int                side,
                          int                bdir)
{
    // This function is a clone of volWgtSum except it computes the result only on half of the domain.
    // The lower half corresponds to side == 0 and the upper half corresponds to side == 1.
    // The argument bdir gives the direction along which to bisect.
    // ONLY WORKS IN THREE DIMENSIONS.

    Real        sum     = 0.0;
    const Real* dx      = geom.CellSize();
    MultiFab*   mf      = derive(name,time,0);
    const int* domlo    = geom.Domain().loVect(); 
    const int* domhi    = geom.Domain().hiVect();

    BL_ASSERT(mf != 0);

    BoxArray baf;

    if (level < parent->finestLevel())
    {
        baf = parent->boxArray(level+1);
        baf.coarsen(fine_ratio);
    }

    for (MFIter mfi(*mf); mfi.isValid(); ++mfi)
    {
        FArrayBox& fab = (*mf)[mfi];

        if (level < parent->finestLevel())
        {
            std::vector< std::pair<int,Box> > isects = baf.intersections(grids[mfi.index()]);

            for (int ii = 0; ii < isects.size(); ii++)
            {
                fab.setVal(0,isects[ii].second,0,fab.nComp());
            }
        }
        Real s = 0.0;
        const Box& box  = mfi.validbox();
        const int* lo   = box.loVect();
        const int* hi   = box.hiVect();

        int hiLeft[3]       = { *hi, *(hi+1), *(hi+2) };
        hiLeft[bdir]        = *(domhi+bdir) / 2;
        int loRight[3]      = { *lo, *(lo+1), *(lo+2) };
        loRight[bdir]       = *(domhi+bdir) / 2 + 1;

        const int* hiLeftPtr  = hiLeft;
        const int* loRightPtr = loRight;

        //
        // Note that this routine will do a volume weighted sum of
        // whatever quantity is passed in, not strictly the "mass".
        //
        
        if ( side == 0 && *(lo + bdir) <= *(hiLeftPtr + bdir) ) {
          if ( *(hi + bdir) <= *(hiLeftPtr + bdir) )
            BL_FORT_PROC_CALL(CA_SUMMASS,ca_summass)
              (BL_TO_FORTRAN(fab),lo,hi,dx,&s);
          else
            BL_FORT_PROC_CALL(CA_SUMMASS,ca_summass)
              (BL_TO_FORTRAN(fab),lo,hiLeftPtr,dx,&s);
	}  
        else if ( side == 1 && *(hi + bdir) >= *(loRightPtr + bdir) ) {
          if ( *(lo + bdir) >= *(loRightPtr + bdir) )
            BL_FORT_PROC_CALL(CA_SUMMASS,ca_summass)
              (BL_TO_FORTRAN(fab),lo,hi,dx,&s);
          else
            BL_FORT_PROC_CALL(CA_SUMMASS,ca_summass)
              (BL_TO_FORTRAN(fab),loRightPtr,hi,dx,&s);
	}
        
        sum += s;
		
    }

    delete mf;

    ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}

Real
Castro::locWgtSumOneSide (const std::string& name,
                          Real               time,
                          int                idir, 
                          int                side,
                          int                bdir)
{
  // This function is a clone of locWgtSum except that it only sums over one half of the domain.
  // The lower half corresponds to side == 0, and the upper half corresponds to side == 1.
  // The argument idir (x == 0, y == 1, z == 2) gives the direction to location weight by,
  // and the argument bdir gives the direction along which to bisect.
  // ONLY WORKS IN THREE DIMENSIONS.

    Real sum            = 0.0;
    const Real* dx      = geom.CellSize();
    MultiFab*   mf      = derive(name,time,0); 
    const int* domlo = geom.Domain().loVect(); 
    const int* domhi = geom.Domain().hiVect(); 

    BL_ASSERT(mf != 0);

    BoxArray baf;

    if (level < parent->finestLevel())
    {
        baf = parent->boxArray(level+1);
        baf.coarsen(fine_ratio);
    }

    for (MFIter mfi(*mf); mfi.isValid(); ++mfi)
    {
        FArrayBox& fab = (*mf)[mfi];

        if (level < parent->finestLevel())
        {
            std::vector< std::pair<int,Box> > isects = baf.intersections(grids[mfi.index()]);

            for (int ii = 0; ii < isects.size(); ii++)
            {
                fab.setVal(0,isects[ii].second,0,fab.nComp());
            }
        }
        Real s = 0.0;
        const Box& box  = mfi.validbox();
        const int* lo   = box.loVect();
         const int* hi   = box.hiVect();

        int hiLeft[3]       = { *hi, *(hi+1), *(hi+2) };
        hiLeft[bdir]        = *(domhi+bdir) / 2;
        int loRight[3]      = { *lo, *(lo+1), *(lo+2) };
        loRight[bdir]       = (*(domhi+bdir) / 2) + 1;

        const int* hiLeftPtr  = hiLeft;
        const int* loRightPtr = loRight;
                
        //
        // Note that this routine will do a volume weighted sum of
        // whatever quantity is passed in, not strictly the "mass".
        // 
        
        if ( (side == 0) && (*(lo + bdir) <= *(hiLeftPtr + bdir)) ) {
          if ( *(hi + bdir) <= *(hiLeftPtr + bdir) )
            BL_FORT_PROC_CALL(CA_SUMLOCMASS,ca_sumlocmass)
              (BL_TO_FORTRAN(fab),lo,hi,geom.ProbLo(),dx,&s,idir);
          else
            BL_FORT_PROC_CALL(CA_SUMLOCMASS,ca_sumlocmass)
              (BL_TO_FORTRAN(fab),lo,hiLeftPtr,geom.ProbLo(),dx,&s,idir);
	}  
        else if ( (side == 1) && (*(hi + bdir) >= *(loRightPtr + bdir)) ) {
          if ( *(lo + bdir) >= *(loRightPtr + bdir) )
            BL_FORT_PROC_CALL(CA_SUMLOCMASS,ca_sumlocmass)
              (BL_TO_FORTRAN(fab),lo,hi,geom.ProbLo(),dx,&s,idir);
          else
            BL_FORT_PROC_CALL(CA_SUMLOCMASS,ca_sumlocmass)
              (BL_TO_FORTRAN(fab),loRightPtr,hi,geom.ProbLo(),dx,&s,idir);
	}
  
        sum += s;
        
    }

    delete mf;

    ParallelDescriptor::ReduceRealSum(sum);

    return sum;

}

Real
Castro::volProductSum (const std::string& name1, 
                       const std::string& name2,
                       Real time)
{
    Real        sum     = 0.0;
    const Real* dx      = geom.CellSize();
    MultiFab*   mf1;
    MultiFab*   mf2;

    if ( name1 == "phi" )
      mf1 = gravity->get_phi_curr(level);
    else
      mf1 = derive(name1,time,0);
    
    if ( name2 == "phi" )
      mf2 = gravity->get_phi_curr(level);
    else
      mf2 = derive(name2,time,0);

    BL_ASSERT(mf1 != 0);
    BL_ASSERT(mf2 != 0);

    BoxArray baf;

    if (level < parent->finestLevel())
    {
        baf = parent->boxArray(level+1);
        baf.coarsen(fine_ratio);
    }

    for (MFIter mfi(*mf1); mfi.isValid(); ++mfi)
    {
        FArrayBox& fab1 = (*mf1)[mfi];
        FArrayBox& fab2 = (*mf2)[mfi];

        if (level < parent->finestLevel())
        {
            std::vector< std::pair<int,Box> > isects = baf.intersections(grids[mfi.index()]);

            for (int ii = 0; ii < isects.size(); ii++)
            {
                fab1.setVal(0,isects[ii].second,0,fab1.nComp());
                fab2.setVal(0,isects[ii].second,0,fab1.nComp());
            }
        }
        Real s = 0.0;
        const Box& box  = mfi.validbox();
        const int* lo   = box.loVect();
        const int* hi   = box.hiVect();

	BL_FORT_PROC_CALL(CA_SUMPRODUCT,ca_sumproduct)
	  (BL_TO_FORTRAN(fab1),BL_TO_FORTRAN(fab2),lo,hi,dx,&s);
        
        sum += s;
    }

    if ( name1 != "phi" )
      delete mf1;
    if ( name2 != "phi" )
      delete mf2;

    ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}

Real
Castro::locSquaredSum (const std::string& name,
                       Real               time,
                       int                idir)
{
    Real        sum     = 0.0;
    const Real* dx      = geom.CellSize();
    MultiFab*   mf      = derive(name,time,0);

    BL_ASSERT(mf != 0);

    BoxArray baf;

    if (level < parent->finestLevel())
    {
        baf = parent->boxArray(level+1);
        baf.coarsen(fine_ratio);
    }

    for (MFIter mfi(*mf); mfi.isValid(); ++mfi)
    {
        FArrayBox& fab = (*mf)[mfi];

        if (level < parent->finestLevel())
        {
            std::vector< std::pair<int,Box> > isects = baf.intersections(grids[mfi.index()]);

            for (int ii = 0; ii < isects.size(); ii++)
            {
                fab.setVal(0,isects[ii].second,0,fab.nComp());
            }
        }
        Real s = 0.0;
        const Box& box  = mfi.validbox();
        const int* lo   = box.loVect();
        const int* hi   = box.hiVect();

	BL_FORT_PROC_CALL(CA_SUMLOCSQUAREDMASS,ca_sumlocsquaredmass)
            (BL_TO_FORTRAN(fab),lo,hi,geom.ProbLo(),dx,&s,idir);

        sum += s;
    }

    delete mf;

    ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}
