/* Implementations of functions in Problem.H go here */

#include "Castro.H"
#include "Castro_F.H"
#include "Problem_F.H"

#include <iomanip>

Real Castro::volInBoundary (const std::string& name,
                            Real               time,
	                    Real               rho_cutoff)
{
    BL_PROFILE("Castro::volInBoundary()");

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

	BL_FORT_PROC_CALL(CA_VOLUMEINDENSITYBOUNDARY,ca_volumeindensityboundary)
	             	  (BL_TO_FORTRAN(fab),lo,hi,dx,&s,&rho_cutoff);
        sum += s;
    }

    delete mf;

    ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}

