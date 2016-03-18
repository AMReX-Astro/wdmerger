#include <cmath>

#include <iomanip>

#include <vector>

#include <Castro.H>
#include <Castro_F.H>
#include <Geometry.H>

#include <Problem.H>
#include <Problem_F.H>

#include <Gravity.H>
#include <Gravity_F.H>

#include "buildInfo.H"

void
Castro::sum_integrated_quantities ()
{
    int finest_level  = parent->finestLevel();
    Real time         = state[State_Type].curTime();
    Real dt           = parent->dtLevel(0);

    if (time == 0.0) dt = 0.0; // dtLevel returns the next timestep for t = 0, so overwrite

    int timestep = parent->levelSteps(0);

    Real mass                 = 0.0;
    Real momentum[3]          = { 0.0 };
    Real angular_momentum[3]  = { 0.0 };
    Real hybrid_momentum[3]   = { 0.0 };
    Real rho_E                = 0.0;
    Real rho_e                = 0.0;
    Real rho_K                = 0.0;
    Real rho_phi              = 0.0;
    Real rho_phirot           = 0.0;

    // Total energy on the grid, including decomposition
    // into the various components.

    Real gravitational_energy = 0.0;
    Real kinetic_energy       = 0.0;
    Real gas_energy           = 0.0;
    Real rotational_energy    = 0.0;
    Real internal_energy      = 0.0;
    Real total_energy         = 0.0;
    Real total_E_grid         = 0.0;

    // Rotation frequency.

    Real omega[3] = { 0.0 };

    get_omega_vec(omega, time);

    // Maximum temperature on the grid.

    Real T_max = 0.0;

    // Maximum density on the grid.

    Real rho_max = 0.0;

    // Maximum t_sound / t_enuc on finest level.

    Real ts_te_max = 0.0;

    // Center of mass of the system.

    Real com[3]       = { 0.0 };
    Real com_vel[3]   = { 0.0 };

    // Stellar masses.

    Real mass_p       = 0.0;
    Real mass_s       = 0.0;

    // Distance between the WDs.

    Real wd_dist[3] = { 0.0 };
    Real wd_dist_init[3] = { 0.0 };

    Real separation = 0.0;
    Real angle = 0.0;

    // Stellar centers of mass and velocities.

    Real com_p_mag = 0.0;
    Real com_s_mag = 0.0;

    Real com_p[3]     = { 0.0 };
    Real com_s[3]     = { 0.0 };

    Real vel_p_mag = 0.0;
    Real vel_s_mag = 0.0;

    Real vel_p[3] = { 0.0 };
    Real vel_s[3] = { 0.0 };

    Real vel_p_rad = 0.0;
    Real vel_s_rad = 0.0;

    Real vel_p_phi = 0.0;
    Real vel_s_phi = 0.0;

    // Gravitational free-fall timescale of the stars.
    
    Real t_ff_p = 0.0;
    Real t_ff_s = 0.0;

    // Gravitational wave amplitudes.
    
    Real h_plus_1  = 0.0;
    Real h_cross_1 = 0.0;

    Real h_plus_2  = 0.0;
    Real h_cross_2 = 0.0;

    Real h_plus_3  = 0.0;
    Real h_cross_3 = 0.0;

    // Number of species.
    
    int NumSpec;
    get_num_spec(&NumSpec);    

    // Species names and total masses on the domain.

    Real M_solar = 1.9884e33;
    
    Real species_mass[NumSpec];
    std::vector<std::string> species_names(NumSpec);
    
    std::string name1; 
    std::string name2;

    int index1;
    int index2;

    int dataprecision = 16; // Number of digits after the decimal point, for float data

    int datwidth      = 25; // Floating point data in scientific notation
    int fixwidth      = 25; // Floating point data not in scientific notation
    int intwidth      = 12; // Integer data

    int axis_1;
    int axis_2;
    int axis_3;

    // Determine various coordinate axes
    get_axes(axis_1, axis_2, axis_3);

    wd_dist_init[axis_1 - 1] = 1.0;

    // Determine the names of the species in the simulation.    

    for (int i = 0; i < NumSpec; i++) {
      species_names[i] = desc_lst[State_Type].name(FirstSpec+i);
      species_names[i] = species_names[i].substr(4,std::string::npos);
      species_mass[i]  = 0.0;	
    }

    for (int lev = 0; lev <= finest_level; lev++)
    {

      // Update the local level we're on.
      
      set_amr_info(lev, -1, -1, -1.0, -1.0);
      
      // Get the current level from Castro

      Castro& ca_lev = getLevel(lev);

      for ( int i = 0; i < 3; i++ ) {
        com[i] += ca_lev.locWgtSum("density", time, i);
      }

      // Calculate total mass, momentum, angular momentum, and energy of system.

      mass += ca_lev.volWgtSum("density", time);

      momentum[0] += ca_lev.volWgtSum("inertial_momentum_x", time);
      momentum[1] += ca_lev.volWgtSum("inertial_momentum_y", time);
      momentum[2] += ca_lev.volWgtSum("inertial_momentum_z", time);

      angular_momentum[0] += ca_lev.volWgtSum("inertial_angular_momentum_x", time);
      angular_momentum[1] += ca_lev.volWgtSum("inertial_angular_momentum_y", time);
      angular_momentum[2] += ca_lev.volWgtSum("inertial_angular_momentum_z", time);

#ifdef HYBRID_MOMENTUM
      hybrid_momentum[0] += ca_lev.volWgtSum("rmom", time);
      hybrid_momentum[1] += ca_lev.volWgtSum("lmom", time);
      hybrid_momentum[2] += ca_lev.volWgtSum("pmom", time);
#endif

      rho_E += ca_lev.volWgtSum("rho_E", time);
      rho_K += ca_lev.volWgtSum("kineng",time);
      rho_e += ca_lev.volWgtSum("rho_e", time);

#ifdef GRAVITY
      if (do_grav)
        rho_phi += ca_lev.volProductSum("density", "phiGrav", time);
#endif

#ifdef ROTATION
      if (do_rotation)
	rho_phirot += ca_lev.volProductSum("density", "phiRot", time);
#endif            
      
      // Gravitational wave signal. This is designed to add to these quantities so we can send them directly.
      ca_lev.gwstrain(time, h_plus_1, h_cross_1, h_plus_2, h_cross_2, h_plus_3, h_cross_3);

      // Integrated mass of all species on the domain.      
      for (int i = 0; i < NumSpec; i++)
	species_mass[i] += ca_lev.volWgtSum("rho_" + species_names[i], time) / M_solar;

      MultiFab& S_new = ca_lev.get_new_data(State_Type);

      // Extrema

      T_max = std::max(T_max, S_new.max(Temp));
      rho_max = std::max(rho_max, S_new.max(Density));

      if (lev == finest_level) {

        MultiFab* ts_te_MF = ca_lev.derive("t_sound_t_enuc", time, 0);
	ts_te_max = std::max(ts_te_max, ts_te_MF->max(0));
	delete ts_te_MF;

      }
    }

    // Return to the original level.
    
    set_amr_info(level, -1, -1, -1.0, -1.0);    
    
    // Complete calculations for energy and momenta

    gravitational_energy = -rho_phi; // CASTRO uses positive phi
    if (gravity->get_gravity_type() == "PoissonGrav")
      gravitational_energy *= 0.5; // avoids double counting
    internal_energy = rho_e;
    kinetic_energy = rho_K;
    gas_energy = rho_E;
    rotational_energy = rho_phirot;
    total_E_grid = gravitational_energy + rho_E;
    total_energy = total_E_grid + rotational_energy;
    
    // Complete calculations for center of mass quantities

    for ( int i = 0; i < 3; i++ ) {

      com[i]       = com[i] / mass;
      com_vel[i]   = momentum[i] / mass;

    }

    get_star_data(com_p, com_s, vel_p, vel_s, &mass_p, &mass_s, &t_ff_p, &t_ff_s);

    com_p_mag += std::pow( std::pow(com_p[0],2) + std::pow(com_p[1],2) + std::pow(com_p[2],2), 0.5 );
    com_s_mag += std::pow( std::pow(com_s[0],2) + std::pow(com_s[1],2) + std::pow(com_s[2],2), 0.5 );
    vel_p_mag += std::pow( std::pow(vel_p[0],2) + std::pow(vel_p[1],2) + std::pow(vel_p[2],2), 0.5 );
    vel_s_mag += std::pow( std::pow(vel_s[0],2) + std::pow(vel_s[1],2) + std::pow(vel_s[2],2), 0.5 );

#if (BL_SPACEDIM == 3)
    if (mass_p > 0.0) {
      vel_p_rad = (com_p[axis_1 - 1] / com_p_mag) * vel_p[axis_1 - 1] + (com_p[axis_2 - 1] / com_p_mag) * vel_p[axis_2 - 1];
      vel_p_phi = (com_p[axis_1 - 1] / com_p_mag) * vel_p[axis_2 - 1] - (com_p[axis_2 - 1] / com_p_mag) * vel_p[axis_1 - 1];
    }

    if (mass_s > 0.0) {
      vel_s_rad = (com_s[axis_1 - 1] / com_s_mag) * vel_s[axis_1 - 1] + (com_s[axis_2 - 1] / com_s_mag) * vel_s[axis_2 - 1];
      vel_s_phi = (com_s[axis_1 - 1] / com_s_mag) * vel_s[axis_2 - 1] - (com_s[axis_2 - 1] / com_s_mag) * vel_s[axis_1 - 1];
    }
#else
    if (mass_p > 0.0) {
      vel_p_rad = vel_p[axis_1 - 1];
      vel_p_phi = vel_p[axis_3 - 1];
    }

    if (mass_s > 0.0) {
      vel_s_rad = vel_s[axis_1 - 1];
      vel_s_phi = vel_s[axis_3 - 1];
    }
#endif

    if (mass_p > 0.0 && mass_s > 0.0) {
      
      // Calculate the distance between the primary and secondary.

      for ( int i = 0; i < 3; i++ ) 
	wd_dist[i] = com_s[i] - com_p[i];
    
      separation = norm(wd_dist);

      // Calculate the angle between the initial stellar axis and
      // the line currently joining the two stars. Note that this
      // neglects any motion in the plane perpendicular to the initial orbit.

      angle = atan2( wd_dist[axis_2 - 1] - wd_dist_init[axis_2 - 1],
                     wd_dist[axis_1 - 1] - wd_dist_init[axis_1 - 1] ) * 180.0 / M_PI;

      // Now let's transform from [-180, 180] to [0, 360].
      
      if (angle < 0.0) angle += 360.0;
      
    }

    // Do remaining reductions

    ParallelDescriptor::ReduceRealMax(T_max);
    ParallelDescriptor::ReduceRealMax(rho_max);
    ParallelDescriptor::ReduceRealMax(ts_te_max);

    // Write data out to the log.

    if ( ParallelDescriptor::IOProcessor() )
    {

      // The data logs are only defined on the IO processor
      // for parallel runs, so the stream should only be opened inside.

      if (parent->NumDataLogs() > 0) {

	 std::ostream& grid_log = parent->DataLog(0);

	 if ( grid_log.good() ) {

	   // Write header row

	   if (time == 0.0) {

	     // Output the git commit hashes used to build the executable.

	     const char* castro_hash   = buildInfoGetGitHash(1);
	     const char* boxlib_hash   = buildInfoGetGitHash(2);
	     const char* wdmerger_hash = buildInfoGetBuildGitHash();

	     grid_log << "# Castro   git hash: " << castro_hash   << std::endl;
	     grid_log << "# BoxLib   git hash: " << boxlib_hash   << std::endl;
	     grid_log << "# wdmerger git hash: " << wdmerger_hash << std::endl;

	     grid_log << std::setw(intwidth) << "#   TIMESTEP";
	     grid_log << std::setw(fixwidth) << "                     TIME";
	     grid_log << std::setw(datwidth) << "             TOTAL ENERGY";
	     grid_log << std::setw(datwidth) << "             TOTAL E GRID";
	     grid_log << std::setw(datwidth) << "               GAS ENERGY";
	     grid_log << std::setw(datwidth) << "              KIN. ENERGY";
	     grid_log << std::setw(datwidth) << "              ROT. ENERGY";
	     grid_log << std::setw(datwidth) << "             GRAV. ENERGY";
	     grid_log << std::setw(datwidth) << "              INT. ENERGY";
#if (BL_SPACEDIM == 3)
	     grid_log << std::setw(datwidth) << "                     XMOM";
	     grid_log << std::setw(datwidth) << "                     YMOM";
	     grid_log << std::setw(datwidth) << "                     ZMOM";
#ifdef HYBRID_MOMENTUM
	     grid_log << std::setw(datwidth) << "              HYB. MOM. R";
	     grid_log << std::setw(datwidth) << "              HYB. MOM. L";
	     grid_log << std::setw(datwidth) << "              HYB. MOM. P";
#endif
	     grid_log << std::setw(datwidth) << "              ANG. MOM. X";
	     grid_log << std::setw(datwidth) << "              ANG. MOM. Y";
	     grid_log << std::setw(datwidth) << "              ANG. MOM. Z";
#else
	     grid_log << std::setw(datwidth) << "                     RMOM";
	     grid_log << std::setw(datwidth) << "                     ZMOM";
	     grid_log << std::setw(datwidth) << "              ANG. MOM. R";
	     grid_log << std::setw(datwidth) << "              ANG. MOM. Z";
#endif
	     grid_log << std::setw(datwidth) << "                     MASS";
#if (BL_SPACEDIM == 3 )
	     grid_log << std::setw(datwidth) << "                    X COM";
	     grid_log << std::setw(datwidth) << "                    Y COM";
	     grid_log << std::setw(datwidth) << "                    Z COM";
	     grid_log << std::setw(datwidth) << "                X COM VEL";
	     grid_log << std::setw(datwidth) << "                Y COM VEL";
	     grid_log << std::setw(datwidth) << "                Z COM VEL";
#else
	     grid_log << std::setw(datwidth) << "                R COM    ";
	     grid_log << std::setw(datwidth) << "                Z COM    ";
	     grid_log << std::setw(datwidth) << "                R COM VEL";
	     grid_log << std::setw(datwidth) << "                Z COM VEL";
#endif
	     grid_log << std::setw(datwidth) << "                    T MAX";
	     grid_log << std::setw(datwidth) << "                  RHO MAX";
	     grid_log << std::setw(datwidth) << "            T_S / T_E MAX";
	     grid_log << std::setw(datwidth) << "             h_+ (axis 1)";
	     grid_log << std::setw(datwidth) << "             h_x (axis 1)";
	     grid_log << std::setw(datwidth) << "             h_+ (axis 2)";
	     grid_log << std::setw(datwidth) << "             h_x (axis 2)";
	     grid_log << std::setw(datwidth) << "             h_+ (axis 3)";
	     grid_log << std::setw(datwidth) << "             h_x (axis 3)";

	     grid_log << std::endl;
	   }

	   // Write data for the present time

	   grid_log << std::fixed;

	   grid_log << std::setw(intwidth)                                     << timestep;
	   grid_log << std::setw(fixwidth) << std::setprecision(dataprecision) << time;

	   grid_log << std::scientific;

	   grid_log << std::setw(datwidth) << std::setprecision(dataprecision) << total_energy;
	   grid_log << std::setw(datwidth) << std::setprecision(dataprecision) << total_E_grid;
	   grid_log << std::setw(datwidth) << std::setprecision(dataprecision) << gas_energy;
	   grid_log << std::setw(datwidth) << std::setprecision(dataprecision) << kinetic_energy;
	   grid_log << std::setw(datwidth) << std::setprecision(dataprecision) << rotational_energy;
	   grid_log << std::setw(datwidth) << std::setprecision(dataprecision) << gravitational_energy;
	   grid_log << std::setw(datwidth) << std::setprecision(dataprecision) << internal_energy;
	   grid_log << std::setw(datwidth) << std::setprecision(dataprecision) << momentum[0];
	   grid_log << std::setw(datwidth) << std::setprecision(dataprecision) << momentum[1];
#if (BL_SPACEDIM == 3)
	   grid_log << std::setw(datwidth) << std::setprecision(dataprecision) << momentum[2];
#endif
#ifdef HYBRID_MOMENTUM
	   grid_log << std::setw(datwidth) << std::setprecision(dataprecision) << hybrid_momentum[0];
	   grid_log << std::setw(datwidth) << std::setprecision(dataprecision) << hybrid_momentum[1];
	   grid_log << std::setw(datwidth) << std::setprecision(dataprecision) << hybrid_momentum[2];
#endif
	   grid_log << std::setw(datwidth) << std::setprecision(dataprecision) << angular_momentum[0];
	   grid_log << std::setw(datwidth) << std::setprecision(dataprecision) << angular_momentum[1];
#if (BL_SPACEDIM == 3)
	   grid_log << std::setw(datwidth) << std::setprecision(dataprecision) << angular_momentum[2];
#endif
	   grid_log << std::setw(datwidth) << std::setprecision(dataprecision) << mass;
	   grid_log << std::setw(datwidth) << std::setprecision(dataprecision) << com[0];
	   grid_log << std::setw(datwidth) << std::setprecision(dataprecision) << com[1];
#if (BL_SPACEDIM == 3)
	   grid_log << std::setw(datwidth) << std::setprecision(dataprecision) << com[2];
#endif
	   grid_log << std::setw(datwidth) << std::setprecision(dataprecision) << com_vel[0];
	   grid_log << std::setw(datwidth) << std::setprecision(dataprecision) << com_vel[1];
#if (BL_SPACEDIM == 3)
	   grid_log << std::setw(datwidth) << std::setprecision(dataprecision) << com_vel[2];
#endif
	   grid_log << std::setw(datwidth) << std::setprecision(dataprecision) << T_max;
	   grid_log << std::setw(datwidth) << std::setprecision(dataprecision) << rho_max;
	   grid_log << std::setw(datwidth) << std::setprecision(dataprecision) << ts_te_max;
	   grid_log << std::setw(datwidth) << std::setprecision(dataprecision) << h_plus_1;
	   grid_log << std::setw(datwidth) << std::setprecision(dataprecision) << h_cross_1;
	   grid_log << std::setw(datwidth) << std::setprecision(dataprecision) << h_plus_2;
	   grid_log << std::setw(datwidth) << std::setprecision(dataprecision) << h_cross_2;
	   grid_log << std::setw(datwidth) << std::setprecision(dataprecision) << h_plus_3;
	   grid_log << std::setw(datwidth) << std::setprecision(dataprecision) << h_cross_3;

	   grid_log << std::endl;
	 }
      }

      if (parent->NumDataLogs() > 1) {

	 std::ostream& star_log = parent->DataLog(1);

	 if ( star_log.good() ) {

	   if (time == 0.0) {

	     // Output the git commit hashes used to build the executable.

	     const char* castro_hash   = buildInfoGetGitHash(1);
	     const char* boxlib_hash   = buildInfoGetGitHash(2);
	     const char* wdmerger_hash = buildInfoGetBuildGitHash();

	     star_log << "# Castro   git hash: " << castro_hash   << std::endl;
	     star_log << "# BoxLib   git hash: " << boxlib_hash   << std::endl;
	     star_log << "# wdmerger git hash: " << wdmerger_hash << std::endl;

	     star_log << std::setw(intwidth) << "#   TIMESTEP";
	     star_log << std::setw(fixwidth) << "                     TIME";

	     star_log << std::setw(datwidth) << "              WD DISTANCE";
	     star_log << std::setw(fixwidth) << "                 WD ANGLE";
	     star_log << std::setw(datwidth) << "          PRIMARY MAG COM";
#if (BL_SPACEDIM == 3)
	     star_log << std::setw(datwidth) << "            PRIMARY X COM";
	     star_log << std::setw(datwidth) << "            PRIMARY Y COM";
	     star_log << std::setw(datwidth) << "            PRIMARY Z COM";
#else
	     star_log << std::setw(datwidth) << "            PRIMARY R COM";
	     star_log << std::setw(datwidth) << "            PRIMARY Z COM";
#endif
	     star_log << std::setw(datwidth) << "          PRIMARY MAG VEL";
	     star_log << std::setw(datwidth) << "          PRIMARY RAD VEL";
	     star_log << std::setw(datwidth) << "          PRIMARY ANG VEL";
#if (BL_SPACEDIM == 3)
	     star_log << std::setw(datwidth) << "            PRIMARY X VEL";
	     star_log << std::setw(datwidth) << "            PRIMARY Y VEL";
	     star_log << std::setw(datwidth) << "            PRIMARY Z VEL";
#else
	     star_log << std::setw(datwidth) << "            PRIMARY R VEL";
	     star_log << std::setw(datwidth) << "            PRIMARY Z VEL";
#endif
	     star_log << std::setw(datwidth) << "             PRIMARY MASS";
	     star_log << std::setw(datwidth) << "       PRIMARY T_FREEFALL";
	     for (int i = 0; i <= 6; ++i)
	       star_log << "       PRIMARY 1E" << i << " RADIUS";

	     star_log << std::setw(datwidth) << "        SECONDARY MAG COM";
#if (BL_SPACEDIM == 3)
	     star_log << std::setw(datwidth) << "          SECONDARY X COM";
	     star_log << std::setw(datwidth) << "          SECONDARY Y COM";
	     star_log << std::setw(datwidth) << "          SECONDARY Z COM";
#else
	     star_log << std::setw(datwidth) << "          SECONDARY R COM";
	     star_log << std::setw(datwidth) << "          SECONDARY Z COM";
#endif
	     star_log << std::setw(datwidth) << "        SECONDARY MAG VEL";
	     star_log << std::setw(datwidth) << "        SECONDARY RAD VEL";
	     star_log << std::setw(datwidth) << "        SECONDARY ANG VEL";
#if (BL_SPACEDIM == 3)
	     star_log << std::setw(datwidth) << "          SECONDARY X VEL";
	     star_log << std::setw(datwidth) << "          SECONDARY Y VEL";
	     star_log << std::setw(datwidth) << "          SECONDARY Z VEL";
#else
	     star_log << std::setw(datwidth) << "          SECONDARY R VEL";
	     star_log << std::setw(datwidth) << "          SECONDARY Z VEL";
#endif
	     star_log << std::setw(datwidth) << "           SECONDARY MASS";
	     star_log << std::setw(datwidth) << "     SECONDARY T_FREEFALL";
	     for (int i = 0; i <= 6; ++i)
	       star_log << "     SECONDARY 1E" << i << " RADIUS";

	     star_log << std::endl;

	   }

	   star_log << std::fixed;

	   star_log << std::setw(intwidth)                                     << timestep;
	   star_log << std::setw(fixwidth) << std::setprecision(dataprecision) << time;

	   star_log << std::scientific;
	   star_log << std::setw(datwidth) << std::setprecision(dataprecision) << separation;

	   star_log << std::fixed;
	   star_log << std::setw(fixwidth) << std::setprecision(dataprecision) << angle;

	   star_log << std::scientific;

	   star_log << std::setw(datwidth) << std::setprecision(dataprecision) << com_p_mag;
	   star_log << std::setw(datwidth) << std::setprecision(dataprecision) << com_p[0];
	   star_log << std::setw(datwidth) << std::setprecision(dataprecision) << com_p[1];
#if (BL_SPACEDIM == 3)
	   star_log << std::setw(datwidth) << std::setprecision(dataprecision) << com_p[2];
#endif
	   star_log << std::setw(datwidth) << std::setprecision(dataprecision) << vel_p_mag;
	   star_log << std::setw(datwidth) << std::setprecision(dataprecision) << vel_p_rad;
	   star_log << std::setw(datwidth) << std::setprecision(dataprecision) << vel_p_phi;
	   star_log << std::setw(datwidth) << std::setprecision(dataprecision) << vel_p[0];
	   star_log << std::setw(datwidth) << std::setprecision(dataprecision) << vel_p[1];
#if (BL_SPACEDIM == 3)
	   star_log << std::setw(datwidth) << std::setprecision(dataprecision) << vel_p[2];
#endif
	   star_log << std::setw(datwidth) << std::setprecision(dataprecision) << mass_p;
	   star_log << std::setw(datwidth) << std::setprecision(dataprecision) << t_ff_p;
	   for (int i = 0; i <= 6; ++i)
	       star_log << std::setw(datwidth) << std::setprecision(dataprecision) << rad_p[i];

	   star_log << std::setw(datwidth) << std::setprecision(dataprecision) << com_s_mag;
	   star_log << std::setw(datwidth) << std::setprecision(dataprecision) << com_s[0];
	   star_log << std::setw(datwidth) << std::setprecision(dataprecision) << com_s[1];
#if (BL_SPACEDIM == 3)
	   star_log << std::setw(datwidth) << std::setprecision(dataprecision) << com_s[2];
#endif
	   star_log << std::setw(datwidth) << std::setprecision(dataprecision) << vel_s_mag;
	   star_log << std::setw(datwidth) << std::setprecision(dataprecision) << vel_s_rad;
	   star_log << std::setw(datwidth) << std::setprecision(dataprecision) << vel_s_phi;
	   star_log << std::setw(datwidth) << std::setprecision(dataprecision) << vel_s[0];
	   star_log << std::setw(datwidth) << std::setprecision(dataprecision) << vel_s[1];
#if (BL_SPACEDIM == 3)
	   star_log << std::setw(datwidth) << std::setprecision(dataprecision) << vel_s[2];
#endif
	   star_log << std::setw(datwidth) << std::setprecision(dataprecision) << mass_s;
	   star_log << std::setw(datwidth) << std::setprecision(dataprecision) << t_ff_s;
	   for (int i = 0; i <= 6; ++i)
	       star_log << std::setw(datwidth) << std::setprecision(dataprecision) << rad_s[i];

	   star_log << std::endl;

	 }
      }

      if (parent->NumDataLogs() > 2) {

	 std::ostream& species_log = parent->DataLog(2);

	 if ( species_log.good() ) {

	   if (time == 0.0) {

	     // Output the git commit hashes used to build the executable.

	     const char* castro_hash   = buildInfoGetGitHash(1);
	     const char* boxlib_hash   = buildInfoGetGitHash(2);
	     const char* wdmerger_hash = buildInfoGetBuildGitHash();

	     species_log << "# Castro   git hash: " << castro_hash   << std::endl;
	     species_log << "# BoxLib   git hash: " << boxlib_hash   << std::endl;
	     species_log << "# wdmerger git hash: " << wdmerger_hash << std::endl;

	     species_log << std::setw(intwidth) << "#   TIMESTEP";
	     species_log << std::setw(fixwidth) << "                  TIME";

	     // We need to be careful here since the species names have differing numbers of characters

	     for (int i = 0; i < NumSpec; i++) {
	       std::string outString  = "";
	       std::string massString = "Mass ";
	       std::string specString = species_names[i];
               while (outString.length() + specString.length() + massString.length() < datwidth) outString += " ";
	       outString += massString;
	       outString += specString;
	       species_log << std::setw(datwidth) << outString;
	     }

	     species_log << std::endl;

	   }

	   species_log << std::fixed;

	   species_log << std::setw(intwidth)                                     << timestep;
	   species_log << std::setw(fixwidth) << std::setprecision(dataprecision) << time;

	   species_log << std::scientific;

	   for (int i = 0; i < NumSpec; i++)
	     species_log << std::setw(datwidth) << std::setprecision(dataprecision) << species_mass[i];

	   species_log << std::endl;

	 }
      }

      if (parent->NumDataLogs() > 3) {

	 std::ostream& amr_log = parent->DataLog(3);

	 if ( amr_log.good() ) {

	   if (time == 0.0) {

	     // Output the git commit hashes used to build the executable.

	     const char* castro_hash   = buildInfoGetGitHash(1);
	     const char* boxlib_hash   = buildInfoGetGitHash(2);
	     const char* wdmerger_hash = buildInfoGetBuildGitHash();

	     amr_log << "# Castro   git hash: " << castro_hash   << std::endl;
	     amr_log << "# BoxLib   git hash: " << boxlib_hash   << std::endl;
	     amr_log << "# wdmerger git hash: " << wdmerger_hash << std::endl;

	     amr_log << std::setw(intwidth) << "#   TIMESTEP";
	     amr_log << std::setw(fixwidth) << "                     TIME";
	     amr_log << std::setw(fixwidth) << "                       DT";
	     amr_log << std::setw(intwidth) << "  FINEST LEV";

	     amr_log << std::endl;

	   }

	   amr_log << std::fixed;

	   amr_log << std::setw(intwidth)                                     << timestep;
	   amr_log << std::setw(fixwidth) << std::setprecision(dataprecision) << time;
	   amr_log << std::setw(fixwidth) << std::setprecision(dataprecision) << dt;
	   amr_log << std::setw(intwidth)                                     << parent->finestLevel();

	   amr_log << std::endl;

	 }

      }

      if (parent->NumDataLogs() > 4 && level == 0) {

	 std::ostream& boundary_log = parent->DataLog(4);

	 if ( boundary_log.good() ) {

	   if (time == 0.0) {

	     // Output the git commit hashes used to build the executable.

	     const char* castro_hash   = buildInfoGetGitHash(1);
	     const char* boxlib_hash   = buildInfoGetGitHash(2);
	     const char* wdmerger_hash = buildInfoGetBuildGitHash();

	     boundary_log << "# Castro   git hash: " << castro_hash   << std::endl;
	     boundary_log << "# BoxLib   git hash: " << boxlib_hash   << std::endl;
	     boundary_log << "# wdmerger git hash: " << wdmerger_hash << std::endl;

	     boundary_log << std::setw(intwidth) << "#   TIMESTEP";
	     boundary_log << std::setw(fixwidth) << "                     TIME";
	     boundary_log << std::setw(datwidth) << "                MASS LOST";
	     boundary_log << std::setw(datwidth) << "                XMOM LOST";
	     boundary_log << std::setw(datwidth) << "                YMOM LOST";
	     boundary_log << std::setw(datwidth) << "                ZMOM LOST";
	     boundary_log << std::setw(datwidth) << "                EDEN LOST";
	     boundary_log << std::setw(datwidth) << "         ANG. MOM. X LOST";
	     boundary_log << std::setw(datwidth) << "         ANG. MOM. Y LOST";
	     boundary_log << std::setw(datwidth) << "         ANG. MOM. Z LOST";

	     boundary_log << std::endl;

	   }

	   boundary_log << std::fixed;

	   boundary_log << std::setw(intwidth)                                     << timestep;
	   boundary_log << std::setw(fixwidth) << std::setprecision(dataprecision) << time;

	   boundary_log << std::scientific;

	   for (int i = 0; i < n_lost; i++)
	     boundary_log << std::setw(datwidth) << std::setprecision(dataprecision) << material_lost_through_boundary_cumulative[i];

	   boundary_log << std::endl;

	 }

      }

    }
}

