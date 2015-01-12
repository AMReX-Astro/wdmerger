#include <cmath>

#include <iomanip>

#include <Castro.H>
#include <Castro_F.H>
#include <Geometry.H>

#include "buildInfo.H"

void
Castro::sum_integrated_quantities ()
{
    int finest_level  = parent->finestLevel();
    Real time         = state[State_Type].curTime();
    Real dt           = parent->dtLevel(0);

    if (time == 0.0) dt = 0.0; // dtLevel returns the next timestep for t = 0, so overwrite

    int step          = parent->levelSteps(0);
    Real mass         = 0.0;
    Real momentum[3]  = { 0.0 };
    Real rho_E        = 0.0;
    Real rho_e        = 0.0;
    Real rho_phi      = 0.0;

    Real gravitational_energy = 0.0; 
    Real kinetic_energy       = 0.0; 
    Real rotational_energy    = 0.0;
    Real internal_energy      = 0.0; 
    Real total_energy         = 0.0;

    Real angular_momentum[3]     = { 0.0 };
    Real moment_of_inertia[3][3] = { 0.0 };
    Real m_r_squared[3]          = { 0.0 };
#ifdef ROTATION
    Real omega[3]     = { 0.0, 0.0, 2.0 * M_PI / rotational_period };
    Real rot_kin_eng    = 0.0;
    Real rot_mom[3] = { 0.0 };
    Real rot_ang_mom[3] = { 0.0 };
#endif

    Real total_E_grid = 0.0;
    Real mom_grid[3]  = { 0.0 };
    Real L_grid[3]    = { 0.0 };

    Real mass_left    = 0.0;
    Real mass_right   = 0.0;

    Real com[3]       = { 0.0 };
    Real com_l[3]     = { 0.0 };
    Real com_r[3]     = { 0.0 };
    Real delta_com[3] = { 0.0 };
 
    Real com_vel[3]   = { 0.0 };
    Real com_vel_l[3] = { 0.0 };
    Real com_vel_r[3] = { 0.0 };

    std::string name1; 
    std::string name2;

    int index1;
    int index2;

    int datawidth     =  24;
    int dataprecision =  16;

    for (int lev = 0; lev <= finest_level; lev++)
    {

      // Get the current level from Castro

      Castro& ca_lev = getLevel(lev);

      // Calculate center of mass quantities.

      mass_left    += ca_lev.volWgtSumOneSide("density", time, 0, 0);
      mass_right   += ca_lev.volWgtSumOneSide("density", time, 1, 0);

      for ( int i = 0; i <= BL_SPACEDIM-1; i++ ) {
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

      mom_grid[0]     += ca_lev.volWgtSum("xmom", time);
      mom_grid[1]     += ca_lev.volWgtSum("ymom", time);
      mom_grid[2]     += ca_lev.volWgtSum("zmom", time);

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

      // Add rotation source terms
#ifdef ROTATION
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

          rot_mom[i] += omega[(i+1)%3]*delta_com[(i+2)%3] - omega[(i+2)%3]*delta_com[(i+1)%3];

          // Rotational energy from motion IN rotating frame == omega dot L_grid

          rot_kin_eng += omega[i] * L_grid[i];

	  // Now add quantities due to motion OF rotating frame

	  for ( int j = 0; j <=2; j++ ) {
            rot_ang_mom[i] += moment_of_inertia[i][j] * omega[j];

            rot_kin_eng += (1.0/2.0) * omega[i] * moment_of_inertia[i][j] * omega[j];
          }

        }
      } 
#endif
    }



    // Complete calculations for COM quantities

    Real center = 0.0;

    for ( int i = 0; i <= 2; i++ ) {
      center = 0.5*(Geometry::ProbLo(i) + Geometry::ProbHi(i));
      com[i]       = com[i] / mass + center;
      com_l[i]     = com_l[i] / mass_left + center; 
      com_r[i]     = com_r[i] / mass_right + center;
      com_vel_l[i] = com_vel_l[i] / mass_left;
      com_vel_r[i] = com_vel_r[i] / mass_right;
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

    // Complete calculations for energy and momenta

    gravitational_energy = (-1.0/2.0) * rho_phi; // avoids double counting; CASTRO uses positive phi
    internal_energy = rho_e;
    kinetic_energy = rho_E - rho_e;
    total_E_grid = gravitational_energy + internal_energy + kinetic_energy;
    total_energy = total_E_grid;

    for (int i = 0; i <= 2; i++) {
        momentum[i] = mom_grid[i];
        angular_momentum[i] = L_grid[i];
    }

#ifdef ROTATION

    rotational_energy = rot_kin_eng;
    total_energy += rotational_energy;
    for (int i = 0; i <= 2; i++) {
        momentum[i] += rot_mom[i];
        angular_momentum[i] += rot_ang_mom[i];
    }

#endif

    // Write data out to the log.

    if ( ParallelDescriptor::IOProcessor() )
    {

      // The data log is only defined on the IO processor
      // for parallel runs, so the stream should only be opened inside.

      std::ostream& data_log1 = parent->DataLog(0);
      if ( data_log1.good() ) {

        // Write header row

        if (time == 0.0) {

	  // Output the git commit hashes used to build the executable.

          const char* castro_hash = buildInfoGetGitHash(1);
	  const char* boxlib_hash = buildInfoGetGitHash(2);

	  data_log1 << "# Castro git hash: " << castro_hash << "\n";
	  data_log1 << "# BoxLib git hash: " << boxlib_hash << "\n";

          data_log1 << std::setw(12)        << " #  TIMESTEP";
          data_log1 << std::setw(datawidth) << "     TIME              ";
	  data_log1 << std::setw(datawidth) << "     DT                ";
          data_log1 << std::setw(datawidth) << " TOTAL ENERGY          ";
	  data_log1 << std::setw(datawidth) << " TOTAL E GRID          ";
          data_log1 << std::setw(datawidth) << " KIN. ENERGY           ";
	  data_log1 << std::setw(datawidth) << " ROT. ENERGY           ";
          data_log1 << std::setw(datawidth) << " GRAV. ENERGY          ";
          data_log1 << std::setw(datawidth) << " INT. ENERGY           ";
          data_log1 << std::setw(datawidth) << " XMOM                  ";
          data_log1 << std::setw(datawidth) << " YMOM                  ";
          data_log1 << std::setw(datawidth) << " ZMOM                  ";
	  data_log1 << std::setw(datawidth) << " XMOM GRID             ";
	  data_log1 << std::setw(datawidth) << " YMOM GRID             ";
	  data_log1 << std::setw(datawidth) << " ZMOM GRID             ";
	  data_log1 << std::setw(datawidth) << " XMOM ROT.             ";
	  data_log1 << std::setw(datawidth) << " XMOM ROT.             ";
	  data_log1 << std::setw(datawidth) << " XMOM ROT.             ";
          data_log1 << std::setw(datawidth) << " ANG. MOM. X           ";
          data_log1 << std::setw(datawidth) << " ANG. MOM. Y           ";
          data_log1 << std::setw(datawidth) << " ANG. MOM. Z           ";
          data_log1 << std::setw(datawidth) << " ANG. MOM. X GRID      ";
          data_log1 << std::setw(datawidth) << " ANG. MOM. Y GRID      ";
          data_log1 << std::setw(datawidth) << " ANG. MOM. Z GRID      ";
          data_log1 << std::setw(datawidth) << " ANG. MOM. X ROT.      ";
          data_log1 << std::setw(datawidth) << " ANG. MOM. Y ROT.      ";
          data_log1 << std::setw(datawidth) << " ANG. MOM. Z ROT.      ";
          data_log1 << std::setw(datawidth) << " MASS                  ";
          data_log1 << std::setw(datawidth) << " LEFT MASS             ";
          data_log1 << std::setw(datawidth) << " RIGHT MASS            ";
          data_log1 << std::setw(datawidth) << " X COM                 ";
          data_log1 << std::setw(datawidth) << " Y COM                 ";
          data_log1 << std::setw(datawidth) << " Z COM                 ";
          data_log1 << std::setw(datawidth) << " X COM VEL             ";
          data_log1 << std::setw(datawidth) << " Y COM VEL             ";
          data_log1 << std::setw(datawidth) << " Z COM VEL             ";
          data_log1 << std::setw(datawidth) << " LEFT X COM            ";
          data_log1 << std::setw(datawidth) << " RIGHT X COM           ";
          data_log1 << std::setw(datawidth) << " LEFT Y COM            ";
          data_log1 << std::setw(datawidth) << " RIGHT Y COM           ";
          data_log1 << std::setw(datawidth) << " LEFT Z COM            ";
          data_log1 << std::setw(datawidth) << " RIGHT Z COM           ";
          data_log1 << std::setw(datawidth) << " LEFT X VEL            ";
          data_log1 << std::setw(datawidth) << " RIGHT X VEL           ";
          data_log1 << std::setw(datawidth) << " LEFT Y VEL            ";
          data_log1 << std::setw(datawidth) << " RIGHT Y VEL           ";
          data_log1 << std::setw(datawidth) << " LEFT Z VEL            ";
          data_log1 << std::setw(datawidth) << " RIGHT Z VEL           ";
         
          data_log1 << std::endl;
        }

        // Write data for the present time

	  data_log1 << std::fixed;

	  data_log1 << std::setw(12)                                            << step;
	  data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << time;
	  data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << dt;

	  data_log1 << std::scientific;

          data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << total_energy;
          data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << total_E_grid;
	  data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << kinetic_energy;
          data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << rotational_energy;	  
	  data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << gravitational_energy;
          data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << internal_energy;
	  data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << momentum[0];
	  data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << momentum[1];
	  data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << momentum[2];
	  data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << mom_grid[0];
	  data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << mom_grid[1];
	  data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << mom_grid[2];
	  data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << rot_mom[0];
	  data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << rot_mom[1];
	  data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << rot_mom[2];
          data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << angular_momentum[0];
          data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << angular_momentum[1]; 
          data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << angular_momentum[2];
          data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << L_grid[0];
          data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << L_grid[1]; 
          data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << L_grid[2];
          data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << rot_ang_mom[0];
          data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << rot_ang_mom[1]; 
          data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << rot_ang_mom[2];
	  data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << mass;
	  data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << mass_left;
	  data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << mass_right;
	  data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << com[0];
	  data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << com[1];
	  data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << com[2];
	  data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << com_vel[0];
	  data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << com_vel[1];
	  data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << com_vel[2];
	  data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << com_l[0];
	  data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << com_r[0];
	  data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << com_l[1];
	  data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << com_r[1];
	  data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << com_l[2];
	  data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << com_r[2];
	  data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << com_vel_l[0];
	  data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << com_vel_r[0];
	  data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << com_vel_l[1];
	  data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << com_vel_r[1];
	  data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << com_vel_l[2];
	  data_log1 << std::setw(datawidth) << std::setprecision(dataprecision) << com_vel_r[2];
	
	  data_log1 << std::endl;
        
      }
    }
}

