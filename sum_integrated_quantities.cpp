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
    Real momentum[3]  = {0.0, 0.0, 0.0};
    Real rho_E        = 0.0;
    Real rho_e        = 0.0;
    Real rho_phi      = 0.0;

    Real gravitational_energy = 0.0; 
    Real kinetic_energy       = 0.0; 
    Real internal_energy      = 0.0; 
    Real total_energy         = 0.0; 

    Real angular_momentum[3]  = {0.0, 0.0, 0.0};
    Real moment_of_inertia[3] = {0.0, 0.0, 0.0};
    Real m_r_squared[3]       = {0.0, 0.0, 0.0};

    Real omega[3]     = {0.0, 0.0, 2.0*3.1415926*rotational_frequency};
    Real delta_L[3]   = {0.0, 0.0, 0.0};

    Real mass_left    = 0.0;
    Real mass_right   = 0.0;

    Real com[3]       = {0.0, 0.0, 0.0};
    Real com_l[3]     = {0.0, 0.0, 0.0};
    Real com_r[3]     = {0.0, 0.0, 0.0};
    Real delta_com[3] = {0.0, 0.0, 0.0};
 
    Real com_vel[3]   = {0.0, 0.0, 0.0};
    Real com_vel_l[3] = {0.0, 0.0, 0.0};
    Real com_vel_r[3] = {0.0, 0.0, 0.0};

    std::string name1; 
    std::string name2;

    int datawidth     =  14;
    int dataprecision =   6;
    
    for (int lev = 0; lev <= finest_level; lev++)
    {

      // Get the current level from Castro

      Castro& ca_lev = getLevel(lev);

      // Calculate center of mass quantities.

      mass_left    += ca_lev.volWgtSumOneSide("density", time, 0, 0);
      mass_right   += ca_lev.volWgtSumOneSide("density", time, 1, 0);

      for ( int i = 0; i <= 2; i++ ) {
        switch ( i ) {
	  case 0 : 
            name1 = "xmom"; break;
          case 1 :
            name1 = "ymom"; break;
          case 2 :
            name1 = "zmom"; break;
	}   

        delta_com[i]  = ca_lev.locWgtSum("density", time, i);
        com[i]       += delta_com[i];

        com_l[i]     += ca_lev.locWgtSumOneSide("density", time, i, 0, 0);
        com_r[i]     += ca_lev.locWgtSumOneSide("density", time, i, 1, 0);
        com_vel_l[i] += ca_lev.volWgtSumOneSide(name1,    time, 0, 0);
        com_vel_r[i] += ca_lev.volWgtSumOneSide(name1,    time, 1, 0);
      }

      // Calculate total mass, momentum and energy of system.

      mass     += ca_lev.volWgtSum("density", time);

      momentum[0]     += ca_lev.volWgtSum("xmom", time);
      momentum[1]     += ca_lev.volWgtSum("ymom", time);
      momentum[2]     += ca_lev.volWgtSum("zmom", time);

      rho_E    += ca_lev.volWgtSum("rho_E", time);
      rho_e    += ca_lev.volWgtSum("rho_e", time);

      if ( do_grav )      
        rho_phi  += ca_lev.volProductSum("density", "phi", time);

      // Calculate total angular momentum of system using L = r x p

      for ( int i = 0; i <= 2; i++ )
        m_r_squared[i] = ca_lev.locSquaredSum("density", time, i);
      
      for ( int i = 0; i <= 2; i++ ) {

        int index1 = (i+1) % 3; 
        int index2 = (i+2) % 3;

        switch (i) {
          case 0 :
            name1 = "ymom"; name2 = "zmom"; break;
          case 1 :
            name1 = "zmom"; name2 = "xmom"; break;
          case 2 :
            name1 = "xmom"; name2 = "ymom"; break;
        }

        moment_of_inertia[i] = m_r_squared[index1] + m_r_squared[index2];

        delta_L[i] = ca_lev.locWgtSum(name2, time, index1) - ca_lev.locWgtSum(name1, time, index2);

        angular_momentum[i]  += delta_L[i];

      }

      // Add rotation source terms

      if ( do_rotation ) {
        for ( int i = 0; i <= 2; i++ ) {
	  
          // Rotational energy == omega dot L + 0.5 * I * omega**2

          kinetic_energy     += omega[i] * delta_L[i] + (1.0/2.0) * moment_of_inertia[i] * omega[i] * omega[i];

	  // Angular momentum == (I * omega); missing a cross term which is irrelevant 
	  // since omega has only one non-zero entry.

          angular_momentum[i] += moment_of_inertia[i] * omega[i];

	  // Momentum == omega x (rho * r)

          int index1 = (i+1) % 3; 
          int index2 = (i+2) % 3;

          momentum[i] += omega[index1]*delta_com[index2] - omega[index2]*delta_com[index1];

        }
      } 
 
    }






    // Complete calculations for COM quantities

    for ( int i = 0; i <= 2; i++ ) {
      com[i]       = com[0] / mass;
      com_l[i]     = com_l[0] / mass_left;
      com_r[i]     = com_r[0] / mass_right;
      com_vel_l[i] = com_vel_l[0] / mass_left;
      com_vel_r[i] = com_vel_r[0] / mass_right;
      com_vel[i]   = momentum[i] / mass;
    } 

    const Real* ml = &mass_left;
    const Real* mr = &mass_right;
    const Real* cxl = &com_l[0];
    const Real* cxr = &com_r[0];
    const Real* cyl = &com_l[1];
    const Real* cyr = &com_r[1];
    const Real* czl = &com_l[2];
    const Real* czr = &com_r[2];

    BL_FORT_PROC_CALL(COM_SAVE,com_save)
      (ml, mr, cxl, cxr, cyl, cyr, czl, czr);







    // Complete calculations for energy

    gravitational_energy = (-1.0/2.0) * rho_phi; // avoids double counting; CASTRO uses positive phi
    internal_energy = rho_e;
    kinetic_energy += rho_E - rho_e;
    total_energy = gravitational_energy + internal_energy + kinetic_energy; 





    
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

	data_log1 << std::fixed;

	data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << time;

	data_log1 << std::scientific;

	data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << mass;
	data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << momentum[0];
	data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << momentum[1];
	data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << momentum[2];
	data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << kinetic_energy;
	data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << gravitational_energy;
        data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << internal_energy;
        data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << total_energy;
        data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << angular_momentum[0];
        data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << angular_momentum[1]; 
        data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << angular_momentum[2];
	data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << com[0];
	data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << com[1];
	data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << com[2];
	data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << mass_left;
	data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << mass_right;
	data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << com_l[0];
	data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << com_r[0];
	data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << com_l[1];
	data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << com_r[1];
	data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << com_l[2];
	data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << com_r[2];
	data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << com_vel[0];
	data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << com_vel[1];
	data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << com_vel[2];
	data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << com_vel_l[0];
	data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << com_vel_r[0];
	data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << com_vel_l[1];
	data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << com_vel_r[1];
	data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << com_vel_l[2];
	data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << com_vel_r[2];
	
	data_log1 << std::endl;
        
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
