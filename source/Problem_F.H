/* problem-specific stuff goes here */

BL_FORT_PROC_DECL(PROBLEM_CHECKPOINT,problem_checkpoint)(int * int_dir_name, int * len);

BL_FORT_PROC_DECL(PROBLEM_RESTART,problem_restart)(int * int_dir_name, int * len);

BL_FORT_PROC_DECL(WDCOM,wdcom)
    (BL_FORT_FAB_ARG(rho), 
     BL_FORT_FAB_ARG(xmom), 
     BL_FORT_FAB_ARG(ymom), 
     BL_FORT_FAB_ARG(zmom),
     const int lo[], const int hi[],
     const Real dx[], const Real* time,
     Real* com_p_x, Real* com_p_y, Real* com_p_z,
     Real* com_s_x, Real* com_s_y, Real* com_s_z,
     Real* vel_p_x, Real* vel_p_y, Real* vel_p_z,
     Real* vel_s_x, Real* vel_s_y, Real* vel_s_z,
     Real* m_p, Real* m_s);

BL_FORT_PROC_DECL(CA_VOLUMEINDENSITYBOUNDARY,ca_volumeindensityboundary)
   (BL_FORT_FAB_ARG(rho), const int lo[], const int hi[], const Real dx[],  
    Real* volp, Real* vols, Real* rho_cutoff);

BL_FORT_PROC_DECL(GET_SINGLE_STAR,get_single_star)
   (int& flag);

BL_FORT_PROC_DECL(SET_STAR_LOCATIONS,set_star_locations)
   (const Real com_p[], const Real com_s[], 
    const Real vel_p[], const Real vel_s[], 
    const Real* mass_p, const Real* mass_s);

BL_FORT_PROC_DECL(GET_STAR_LOCATIONS,get_star_locations)
   (Real com_p[], Real com_s[], 
    Real vel_p[], Real vel_s[], 
    Real* mass_p, Real* mass_s);
