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
    Real rho_E       = 0.0;
    Real com_xloc     = 0.0;
    Real mass_lower_x = 0.0;
    Real mass_upper_x = 0.0;
    Real com_xloc_l   = 0.0;
    Real com_xloc_u   = 0.0;
    Real com_xvel     = 0.0;
#if (BL_SPACEDIM>=2)
    Real ymom        = 0.0;
    Real com_yloc_l   = 0.0;
    Real com_yloc_u   = 0.0;
    Real com_yloc    = 0.0;
    Real com_yvel    = 0.0;
#endif
#if (BL_SPACEDIM==3)
    Real zmom        = 0.0;
    Real com_zloc_l   = 0.0;
    Real com_zloc_u   = 0.0;
    Real com_zloc    = 0.0;
    Real com_zvel    = 0.0;
#endif

    for (int lev = 0; lev <= finest_level; lev++)
    {
        Castro& ca_lev = getLevel(lev);

        mass     += ca_lev.volWgtSum("density", time);
        xmom     += ca_lev.volWgtSum("xmom", time);
#if (BL_SPACEDIM == 2)
       if (Geometry::IsRZ()) 
          xmom = 0.;
#endif

       if (show_center_of_mass) {
	 com_xloc      += ca_lev.locWgtSum("density", time, 0);
         mass_lower_x  += ca_lev.volWgtSumOneSide("density", time, 0, 0);
         mass_upper_x  += ca_lev.volWgtSumOneSide("density", time, 1, 0);
         com_xloc_l    += ca_lev.locWgtSumOneSide("density", time, 0, 0, 0);
         com_xloc_u    += ca_lev.locWgtSumOneSide("density", time, 0, 1, 0);
       }
#if (BL_SPACEDIM>=2)
       ymom     += ca_lev.volWgtSum("ymom", time);
       if (show_center_of_mass) {
         com_yloc      += ca_lev.locWgtSum("density", time, 1);
         com_yloc_l    += ca_lev.locWgtSumOneSide("density", time, 1, 0, 0);
         com_yloc_u    += ca_lev.locWgtSumOneSide("density", time, 1, 1, 0);
       }
#endif
#if (BL_SPACEDIM==3)
       zmom     += ca_lev.volWgtSum("zmom", time);
       if (show_center_of_mass) {
         com_zloc      += ca_lev.locWgtSum("density", time, 2);
         com_zloc_l    += ca_lev.locWgtSumOneSide("density", time, 2, 0, 0);
         com_zloc_u    += ca_lev.locWgtSumOneSide("density", time, 2, 1, 0);
       }
#endif
    }

    if (verbose > 0 && ParallelDescriptor::IOProcessor())
    {
        std::cout << '\n';
        std::cout << "TIME= " << time << " MASS        = "   << mass  << '\n';
        std::cout << "TIME= " << time << " XMOM        = "   << xmom     << '\n';
#if (BL_SPACEDIM>=2)
        std::cout << "TIME= " << time << " YMOM        = "   << ymom     << '\n';
#endif
#if (BL_SPACEDIM==3)
        std::cout << "TIME= " << time << " ZMOM        = "   << zmom     << '\n';
#endif
        std::cout << "TIME= " << time << " RHO*E       = "   << rho_E     << '\n';

        
        if (show_center_of_mass) {
           com_xloc = com_xloc / mass;
           com_xloc_l = com_xloc_l / mass_lower_x;
           com_xloc_u = com_xloc_u / mass_upper_x;
           com_xvel = xmom / mass;
           std::cout << "TIME= " << time << " CENTER OF MASS X-LOC = " << com_xloc  << '\n';
	   std::cout << "TIME= " << time << " LOWER CENTER OF MASS X-LOC = " << com_xloc_l << '\n';
	   std::cout << "TIME= " << time << " UPPER CENTER OF MASS X-LOC = " << com_xloc_u << '\n';
           std::cout << "TIME= " << time << " CENTER OF MASS X-VEL = " << com_xvel  << '\n';
#if (BL_SPACEDIM>=2)
           com_yloc = com_yloc / mass;
           com_yloc_l = com_yloc_l / mass_lower_x;
           com_yloc_u = com_yloc_u / mass_upper_x;
           com_yvel = ymom / mass;
           std::cout << "TIME= " << time << " CENTER OF MASS Y-LOC = " << com_yloc  << '\n';
	   std::cout << "TIME= " << time << " LOWER CENTER OF MASS Y-LOC = " << com_yloc_l << '\n';
	   std::cout << "TIME= " << time << " UPPER CENTER OF MASS Y-LOC = " << com_yloc_u << '\n';
           std::cout << "TIME= " << time << " CENTER OF MASS Y-VEL = " << com_yvel  << '\n';
#endif
#if (BL_SPACEDIM==3)
           com_zloc = com_zloc / mass;
           com_zloc_l = com_zloc_l / mass_lower_x;
           com_zloc_u = com_zloc_u / mass_upper_x;
           com_zvel = zmom / mass;
           std::cout << "TIME= " << time << " CENTER OF MASS Z-LOC = " << com_zloc  << '\n';
	   std::cout << "TIME= " << time << " LOWER CENTER OF MASS Z-LOC = " << com_zloc_l << '\n';
	   std::cout << "TIME= " << time << " UPPER CENTER OF MASS Z-LOC = " << com_zloc_u << '\n';
           std::cout << "TIME= " << time << " CENTER OF MASS Z-VEL = " << com_zvel  << '\n';
#endif
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




