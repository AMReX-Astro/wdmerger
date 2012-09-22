#include <iomanip>

#include <Castro.H>
#include <Castro_F.H>

void
Castro::sum_integrated_quantities ()
{
    int finest_level = parent->finestLevel();
    Real time        = state[State_Type].curTime();
    Real mass        = 0.0;
    Real xmom        = 0.0;
    Real ymom        = 0.0;
    Real zmom        = 0.0;
    Real rho_E       = 0.0;

    Real com_xloc     = 0.0;
    Real mass_left    = 0.0;
    Real mass_right   = 0.0;
    Real com_xloc_l   = 0.0;
    Real com_xloc_r   = 0.0;
    Real com_xvel     = 0.0;

    Real com_yloc_l   = 0.0;
    Real com_yloc_r   = 0.0;
    Real com_yloc     = 0.0;
    Real com_yvel     = 0.0;

    Real com_zloc_l   = 0.0;
    Real com_zloc_r   = 0.0;
    Real com_zloc     = 0.0;
    Real com_zvel     = 0.0;

    int datawidth     = 14;
    int dataprecision = 6;

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

      // If castro.show_center_of_mass == 1 in the inputs file, calculate center of mass quantities. 

      if (show_center_of_mass) {
        com_xloc      += ca_lev.locWgtSum("density", time, 0);
        mass_left     += ca_lev.volWgtSumOneSide("density", time, 0, 0);
        mass_right    += ca_lev.volWgtSumOneSide("density", time, 1, 0);
        com_xloc_l    += ca_lev.locWgtSumOneSide("density", time, 0, 0, 0);
        com_xloc_r    += ca_lev.locWgtSumOneSide("density", time, 0, 1, 0);

        com_yloc      += ca_lev.locWgtSum("density", time, 1);
        com_yloc_l    += ca_lev.locWgtSumOneSide("density", time, 1, 0, 0);
        com_yloc_r    += ca_lev.locWgtSumOneSide("density", time, 1, 1, 0);

        com_zloc      += ca_lev.locWgtSum("density", time, 2);
        com_zloc_l    += ca_lev.locWgtSumOneSide("density", time, 2, 0, 0);
        com_zloc_r    += ca_lev.locWgtSumOneSide("density", time, 2, 1, 0);
      }
        
    }

    // Complete calculations for COM quantities

    if (show_center_of_mass) {
      com_xloc = com_xloc / mass;
      com_xloc_l = com_xloc_l / mass_left;
      com_xloc_r = com_xloc_r / mass_right;
      com_xvel = xmom / mass;

      com_yloc = com_yloc / mass;
      com_yloc_l = com_yloc_l / mass_left;
      com_yloc_r = com_yloc_r / mass_right;
      com_yvel = ymom / mass;

      com_zloc = com_zloc / mass;
      com_zloc_l = com_zloc_l / mass_left;
      com_zloc_r = com_zloc_r / mass_right;
      com_zvel = zmom / mass;
    }

    // Write data out to the log. Check if data log exists first.


    if ( ParallelDescriptor::IOProcessor() )
    {
      std::ostream& data_log1 = parent->DataLog(0);
      if ( data_log1.good() )
      {
      // Write header row

      if (time == 0.0) {
        data_log1 << std::setw(datawidth) << "     TIME    ";
        data_log1 << std::setw(datawidth) << " MASS        ";
        data_log1 << std::setw(datawidth) << " XMOM        ";
        data_log1 << std::setw(datawidth) << " YMOM        ";
        data_log1 << std::setw(datawidth) << " ZMOM        ";
        data_log1 << std::setw(datawidth) << " RHO*E       ";
        if (show_center_of_mass) {
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
        } 
        data_log1 << std::endl;
      }

      // Write data for the present time

      data_log1 << std::fixed;

      data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << time;

      data_log1 << std::scientific;

      data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << mass;
      data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << xmom;
      data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << ymom;
      data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << zmom;
      data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << rho_E;
      if (show_center_of_mass) {
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
      }
      data_log1 << std::endl;
      }

    }
        
    // If castro.v == 1, then write out the fundamental quantities to stdout
    // If castro.show_center_of_mass == 1, then also write out the center of mass quantities

    if (verbose > 0 && ParallelDescriptor::IOProcessor())
    {
        std::cout << '\n';
        std::cout << "TIME= " << time << " MASS        = "   << mass  << '\n';
        std::cout << "TIME= " << time << " XMOM        = "   << xmom     << '\n';
        std::cout << "TIME= " << time << " YMOM        = "   << ymom     << '\n';
        std::cout << "TIME= " << time << " ZMOM        = "   << zmom     << '\n';
        std::cout << "TIME= " << time << " RHO*E       = "   << rho_E     << '\n';
      
        if (show_center_of_mass) {
           std::cout << "TIME= " << time << " CENTER OF MASS X-LOC       = " << com_xloc   << '\n';
	   std::cout << "TIME= " << time << " LEFT  CENTER OF MASS X-LOC = " << com_xloc_l << '\n';
	   std::cout << "TIME= " << time << " RIGHT CENTER OF MASS X-LOC = " << com_xloc_r << '\n';
           std::cout << "TIME= " << time << " CENTER OF MASS X-VEL       = " << com_xvel   << '\n';

           std::cout << "TIME= " << time << " CENTER OF MASS Y-LOC       = " << com_yloc   << '\n';
	   std::cout << "TIME= " << time << " LEFT  CENTER OF MASS Y-LOC = " << com_yloc_l << '\n';
	   std::cout << "TIME= " << time << " RIGHT CENTER OF MASS Y-LOC = " << com_yloc_r << '\n';
           std::cout << "TIME= " << time << " CENTER OF MASS Y-VEL       = " << com_yvel   << '\n';

           std::cout << "TIME= " << time << " CENTER OF MASS Z-LOC       = " << com_zloc   << '\n';
	   std::cout << "TIME= " << time << " LEFT  CENTER OF MASS Z-LOC = " << com_zloc_l << '\n';
	   std::cout << "TIME= " << time << " RIGHT CENTER OF MASS Z-LOC = " << com_zloc_r << '\n';
           std::cout << "TIME= " << time << " CENTER OF MASS Z-VEL       = " << com_zvel   << '\n';

        }
	std::cout<<'\n';
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
        int irlo        = lo[0]-radius_grow;
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




