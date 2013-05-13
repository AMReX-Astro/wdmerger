/* problem-specific stuff goes here */

BL_FORT_PROC_DECL(COM_SAVE,com_save)
    ( const Real* mass_p, const Real* mass_right, 
      const Real* cxp, const Real* cxs,
      const Real* cyp, const Real* cys,
      const Real* czp, const Real* czs
    );


BL_FORT_PROC_DECL(PROBLEM_CHECKPOINT,problem_checkpoint)(int * int_dir_name, int * len);

BL_FORT_PROC_DECL(PROBLEM_RESTART,problem_restart)(int * int_dir_name, int * len);



