#include <iomanip>

#include <Castro.H>
#include <Castro_F.H>
#include <Gravity.H>

void
Castro::center_of_mass ()
{
    int finest_level  = parent->finestLevel();
    Real time         = state[State_Type].curTime();
    Real mass         = 0.0;

    Real mass_left    = 0.0;
    Real mass_right   = 0.0;

    Real com[3]       = {0.0, 0.0, 0.0};
    Real com_l[3]     = {0.0, 0.0, 0.0};
    Real com_r[3]     = {0.0, 0.0, 0.0};
    Real delta_com[3] = {0.0, 0.0, 0.0};
    
    for (int lev = 0; lev <= finest_level; lev++)
    {
      // Get the current level from Castro
      Castro& ca_lev = getLevel(lev);

      // Calculate center of mass quantities.
      mass_left    += ca_lev.volWgtSumOneSide("density", time, 0, 0);
      mass_right   += ca_lev.volWgtSumOneSide("density", time, 1, 0);

      for ( int i = 0; i < BL_SPACEDIM; i++ ) {
        delta_com[i]  = ca_lev.locWgtSum("density", time, i);
        com[i]       += delta_com[i];

        com_l[i]     += ca_lev.locWgtSumOneSide("density", time, i, 0, 0);
        com_r[i]     += ca_lev.locWgtSumOneSide("density", time, i, 1, 0);
      }

      // Calculate total mass of system.
      mass     += ca_lev.volWgtSum("density", time);
    }

    // Complete calculations for COM quantities

    for ( int i = 0; i < BL_SPACEDIM; i++ ) {
      com[i]       = com[0] / mass;
      com_l[i]     = com_l[0] / mass_left;
      com_r[i]     = com_r[0] / mass_right;
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
}
