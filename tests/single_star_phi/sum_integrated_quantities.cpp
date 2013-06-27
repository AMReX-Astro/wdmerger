#include <iomanip>

#include <Castro.H>
#include <Castro_F.H>
#include <Gravity.H>
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
    Real moment_of_inertia[3][3] = { 0.0 };
    Real m_r_squared[3]          = { 0.0 };

    Real omega[3]     = { 0.0, 0.0, 2.0*3.14159265358979*rotational_frequency };
    Real L_grid[3]    = { 0.0 };

    Real com[3]       = { 0.0 };
    Real delta_com[3] = { 0.0 };
 
    Real com_vel[3]   = { 0.0 };

    std::string name1; 
    std::string name2;

    int index1;
    int index2;

    int datawidth     =  23;
    int dataprecision =  15;

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
#if (BL_SPACEDIM == 3)
      momentum[2]     += ca_lev.volWgtSum("zmom", time);
#endif

      rho_E    += ca_lev.volWgtSum("rho_E", time);
      rho_e    += ca_lev.volWgtSum("rho_e", time);

#ifdef GRAVITY
      if ( do_grav && gravity->get_gravity_type() == "PoissonGrav") {
        rho_phi  += ca_lev.volProductSum("density", "phi", time);
      }
#endif

      // Calculate total angular momentum on the grid using L = r x p

#if (BL_SPACEDIM == 3)      
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
#elif (BL_SPACEDIM == 2)
	L_grid[2] = ca_lev.locWgtSum("ymom", time, 0) - ca_lev.locWgtSum("xmom", time, 1);
#endif

      // Add rotation source terms

      if ( do_rotation ) {

        // Construct (symmetric) moment of inertia tensor

	for ( int i = 0; i <= BL_SPACEDIM-1; i++ ) {
          m_r_squared[i] = ca_lev.locWgtSum2D("density", time, i, i);
        }

        for ( int i = 0; i <= 2; i++ ) {
	  for ( int j = 0; j <= 2; j++ ) {
            if ( i <= j ) {
              if ( i != j )  {
                if ( ( i < BL_SPACEDIM ) && ( j < BL_SPACEDIM ) ) // Protect against computing z direction sum in 2D
                  moment_of_inertia[i][j] = -ca_lev.locWgtSum2D("density", time, i, j);
              }
              else
                moment_of_inertia[i][j] = m_r_squared[(i+1)%3] + m_r_squared[(i+2)%3];
            }
	    else
              moment_of_inertia[i][j] = moment_of_inertia[j][i];
          }
        }

        for ( int i = 0; i <= 2; i++ ) {

	  // Momentum source from motion IN rotating frame == omega x (rho * r)

          momentum[i] += omega[(i+1)%3]*delta_com[(i+2)%3] - omega[(i+2)%3]*delta_com[(i+1)%3];

          // Rotational energy from motion IN rotating frame == omega dot L_grid

          kinetic_energy     += omega[i] * L_grid[i];

	  // Now add quantities due to motion OF rotating frame

	  for ( int j = 0; j <=2; j++ ) {
            angular_momentum[i] += moment_of_inertia[i][j] * omega[j];

            kinetic_energy += (1.0/2.0) * omega[i] * moment_of_inertia[i][j] * omega[j];
          }

        }
      } 
 
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
	
	  data_log1 << std::endl;
        
      }
    }
}

