/* problem-specific stuff goes here */

BL_FORT_PROC_DECL(COM_SAVE,com_save)
    ( const Real* mass_left, const Real* mass_right, 
      const Real* com_xloc_l, const Real* com_xloc_r,
      const Real* com_yloc_l, const Real* com_yloc_r,
      const Real* com_zloc_l, const Real* com_zloc_r
    );

BL_FORT_PROC_DECL(GET_EOS,get_eos)
    (const int* eos);

BL_FORT_PROC_DECL(SET_EOS,set_eos)
    (const int* eos);
 
BL_FORT_PROC_DECL(CA_SUMPRODUCT,ca_sumproduct)
    ( BL_FORT_FAB_ARG(f1), BL_FORT_FAB_ARG(f2),
      const int lo[], const int hi[], const Real dx[],
      Real* s
    );

BL_FORT_PROC_DECL(CA_SUMLOCSQUAREDMASS,ca_sumlocsquaredmass)
    (BL_FORT_FAB_ARG(rho), const int lo[], const int hi[],
     const Real* problo, const Real dx[], Real* mass, const int& idir); 
