#include <iomanip>

#include <Castro.H>
#include <Castro_F.H>
#include <Geometry.H>

void
Castro::sum_integrated_quantities ()
{
    int finest_level  = parent->finestLevel();
    Real time         = state[State_Type].curTime();
    Real mass         = 0.0;
    Real momentum[3]  = { 0.0 };
    Real rho_E        = 0.0;
    Real rho_e        = 0.0;
    Real rho_phi      = 0.0;

    Real gravitational_energy = 0.0; 
    Real kinetic_energy       = 0.0; 
    Real internal_energy      = 0.0; 
    Real total_energy         = 0.0;

    Real angular_momentum[3]     = { 0.0 };
    Real L_grid[3]    = { 0.0 };

    Real com[3]       = { 0.0 };
    Real delta_com[3] = { 0.0 };
 
    Real com_vel[3]   = { 0.0 };

    Real vol0 = 0.0;
    Real vol1 = 0.0;
    Real vol2 = 0.0;
    Real vol3 = 0.0;
    Real vol4 = 0.0;
    Real vol5 = 0.0;
    Real vol6 = 0.0;

    int index1, index2;
    std::string name1, name2;

    int datawidth     =  23;
    int dataprecision =  15;

    Real pi = 3.1459265358979;

    for (int lev = 0; lev <= finest_level; lev++)
    {

      // Get the current level from Castro

      Castro& ca_lev = getLevel(lev);

      // Calculate center of mass quantities.

      for ( int i = 0; i <= BL_SPACEDIM-1; i++ ) {

        delta_com[i]  = ca_lev.locWgtSum("density", time, i);
        com[i]       += delta_com[i];

      }

      // Calculate total mass, momentum and energy of system.

      mass     += ca_lev.volWgtSum("density", time);

      momentum[0]     += ca_lev.volWgtSum("xmom", time);
      momentum[1]     += ca_lev.volWgtSum("ymom", time);
      momentum[2]     += ca_lev.volWgtSum("zmom", time);

      rho_E    += ca_lev.volWgtSum("rho_E", time);
      rho_e    += ca_lev.volWgtSum("rho_e", time);

#ifdef GRAVITY
      if ( do_grav ) {
        rho_phi  += ca_lev.volProductSum("density", "phi", time);
      }
#endif

      // Calculate total angular momentum on the grid using L = r x p

      for ( int i = 0; i <= 2; i++ ) {

        index1 = (i+1) % 3; 
        index2 = (i+2) % 3;

        switch (i) {
          case 0 :
            name1 = "ymom"; name2 = "zmom"; break;
          case 1 :
            name1 = "zmom"; name2 = "xmom"; break;
          case 2 :
            name1 = "xmom"; name2 = "ymom"; break;
        }

        L_grid[i] = ca_lev.locWgtSum(name2, time, index1) - ca_lev.locWgtSum(name1, time, index2);

        angular_momentum[i]  += L_grid[i];

      }

      // Calculate volume within this cutoff

      vol6 += ca_lev.volInBoundary("density", time, 1.0e6);
      vol5 += ca_lev.volInBoundary("density", time, 1.0e5);
      vol4 += ca_lev.volInBoundary("density", time, 1.0e4);
      vol3 += ca_lev.volInBoundary("density", time, 1.0e3);
      vol2 += ca_lev.volInBoundary("density", time, 1.0e2);
      vol1 += ca_lev.volInBoundary("density", time, 1.0e1);
      vol0 += ca_lev.volInBoundary("density", time, 1.0e0);
      
    }



    // Complete calculations for COM quantities

    Real center = 0.0;

    for ( int i = 0; i <= 2; i++ ) {
      center = 0.5*(Geometry::ProbLo(i) + Geometry::ProbHi(i));
      com[i]       = com[i] / mass + center;
      com_vel[i]   = momentum[i] / mass;
    } 


    // Complete calculations for energy

    gravitational_energy = (-1.0/2.0) * rho_phi; // avoids double counting; CASTRO uses positive phi
    internal_energy = rho_e;
    kinetic_energy += rho_E - rho_e;
    total_energy = gravitational_energy + internal_energy + kinetic_energy; 

    // Complete calculation for bounded volume

    Real rad0 = std::pow(vol0 * 3.0 / 4.0 / pi, 1.0/3.0);
    Real rad1 = std::pow(vol1 * 3.0 / 4.0 / pi, 1.0/3.0);
    Real rad2 = std::pow(vol2 * 3.0 / 4.0 / pi, 1.0/3.0);
    Real rad3 = std::pow(vol3 * 3.0 / 4.0 / pi, 1.0/3.0);
    Real rad4 = std::pow(vol4 * 3.0 / 4.0 / pi, 1.0/3.0);
    Real rad5 = std::pow(vol5 * 3.0 / 4.0 / pi, 1.0/3.0);
    Real rad6 = std::pow(vol6 * 3.0 / 4.0 / pi, 1.0/3.0);
    
    // Write data out to the log.

    if ( ParallelDescriptor::IOProcessor() )
    {

      // The data log is only defined on the IO processor
      // for parallel runs, so the stream should only be opened inside.

      std::ostream& data_log1 = parent->DataLog(0);
      if ( data_log1.good() ) {

        // Write header row

        if (time == 0.0) {
          data_log1 << std::setw(datawidth) << "#     TIME             ";
          data_log1 << std::setw(datawidth) << "  MASS                 ";
          data_log1 << std::setw(datawidth) << "  XMOM                 ";
          data_log1 << std::setw(datawidth) << "  YMOM                 ";
          data_log1 << std::setw(datawidth) << "  ZMOM                 ";
          data_log1 << std::setw(datawidth) << "  KINETIC ENERGY       ";
          data_log1 << std::setw(datawidth) << "  POTENTIAL ENERGY     ";
          data_log1 << std::setw(datawidth) << "  INTERNAL ENERGY      ";
          data_log1 << std::setw(datawidth) << "  TOTAL ENERGY         ";
          data_log1 << std::setw(datawidth) << "  ANGULAR MOMENTUM X   ";
          data_log1 << std::setw(datawidth) << "  ANGULAR MOMENTUM Y   ";
          data_log1 << std::setw(datawidth) << "  ANGULAR MOMENTUM Z   ";
          data_log1 << std::setw(datawidth) << "  X COM                ";
          data_log1 << std::setw(datawidth) << "  Y COM                ";
          data_log1 << std::setw(datawidth) << "  Z COM                ";
          data_log1 << std::setw(datawidth) << "  X VELOCITY           ";
          data_log1 << std::setw(datawidth) << "  Y VELOCITY           ";
          data_log1 << std::setw(datawidth) << "  Z VELOCITY           ";
	  data_log1 << std::setw(datawidth) << "  1E6 BOUNDARY RAD.    ";
	  data_log1 << std::setw(datawidth) << "  1E5 BOUNDARY RAD.    ";
	  data_log1 << std::setw(datawidth) << "  1E4 BOUNDARY RAD.    ";
	  data_log1 << std::setw(datawidth) << "  1E3 BOUNDARY RAD.    ";
	  data_log1 << std::setw(datawidth) << "  1E2 BOUNDARY RAD.    ";
	  data_log1 << std::setw(datawidth) << "  1E1 BOUNDARY RAD.    ";
	  data_log1 << std::setw(datawidth) << "  1E0 BOUNDARY RAD.    ";
         
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
	  data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << com_vel[0];
	  data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << com_vel[1];
	  data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << com_vel[2];
	  data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << rad6;
	  data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << rad5;
	  data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << rad4;
	  data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << rad3;
	  data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << rad2;
	  data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << rad1;
	  data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << rad0;
	
	  data_log1 << std::endl;
        
      }
    }
}

