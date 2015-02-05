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
     Real com_p[], Real com_s[],
     Real com_vel_p[], Real com_vel_s[],
     Real* mass_p, Real* mass_s);
