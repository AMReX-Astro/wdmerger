
BL_FORT_PROC_DECL(SCF_SETUP_RELAXATION,scf_setup_relaxation)
   (const Real dx[], const Real problo[], const Real probhi[]);

BL_FORT_PROC_DECL(SCF_GET_COEFF_INFO,scf_get_coeff_info)
(int* loc_A, int* loc_B, int* loc_C,
    Real* c_A, Real* c_B, Real* c_C);

BL_FORT_PROC_DECL(SCF_GET_OMEGASQ,scf_get_omegasq)
   (const int lo[], const int hi[], const int domlo[], const int domhi[],
    BL_FORT_FAB_ARG(state),
    BL_FORT_FAB_ARG(phi),
    const Real dx[], const Real problo[], const Real probhi[], Real* omegasq);

BL_FORT_PROC_DECL(SCF_GET_BERNOULLI_CONST,scf_get_bernoulli_const)
   (const int lo[], const int hi[], const int domlo[], const int domhi[],
    BL_FORT_FAB_ARG(state),
    BL_FORT_FAB_ARG(phi),
    const Real dx[], const Real problo[], const Real probhi[], 
    Real* bernoulli_1, Real* bernoulli_2);

BL_FORT_PROC_DECL(SCF_CONSTRUCT_ENTHALPY,scf_construct_enthalpy)
   (const int lo[], const int hi[], const int domlo[], const int domhi[],
    BL_FORT_FAB_ARG(state),
    BL_FORT_FAB_ARG(phi),
    BL_FORT_FAB_ARG(enthalpy),
    const Real dx[], const Real problo[], const Real probhi[], 
    Real* bernoulli_1, Real* bernoulli_2, Real* h_max_1, Real* h_max_2);

BL_FORT_PROC_DECL(SCF_UPDATE_DENSITY,scf_update_density)
   (const int lo[], const int hi[], const int domlo[], const int domhi[],
    BL_FORT_FAB_ARG(state),
    BL_FORT_FAB_ARG(phi),
    BL_FORT_FAB_ARG(enthalpy),
    const Real dx[], const Real problo[], const Real probhi[], 
    Real* h_max_1, Real* h_max_2,
    Real* kin_eng, Real* pot_eng, Real* int_eng,
    Real* left_mass, Real* right_mass,
    Real* delta_rho, Real* l2_norm_resid, Real* l2_norm_source);

BL_FORT_PROC_DECL(SCF_CHECK_CONVERGENCE,scf_check_convergence)
   (Real* kin_eng, Real* pot_eng, Real* int_eng,
    Real* left_mass, Real* right_mass,
    Real* delta_rho, Real* l2_norm,
    int* is_relaxed, int* num_iterations);
