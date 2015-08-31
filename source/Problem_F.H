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

BL_FORT_PROC_DECL(QUADRUPOLE_TENSOR_DOUBLE_DOT,quadrupole_tensor_double_dot)
   (BL_FORT_FAB_ARG(rho), 
    BL_FORT_FAB_ARG(xmom), 
    BL_FORT_FAB_ARG(ymom), 
    BL_FORT_FAB_ARG(zmom),
    BL_FORT_FAB_ARG(gx),
    BL_FORT_FAB_ARG(gy),
    BL_FORT_FAB_ARG(gz),
    const int lo[], const int hi[],
    const Real dx[], const Real* time,
    Real* Qtt);

BL_FORT_PROC_DECL(GW_STRAIN_TENSOR,gw_strain_tensor)
   (Real* h_plus_rot, Real* h_cross_rot, 
    Real* h_plus_star, Real* h_cross_star,
    Real* h_plus_motion, Real* h_cross_motion,
    Real* Qtt, const Real* time);

BL_FORT_PROC_DECL(SET_PERIOD,set_period)
   (const Real* period);

BL_FORT_PROC_DECL(UPDATE_CENTER,update_center)
   (const Real* time);

BL_FORT_PROC_DECL(GET_RELAX_TYPE,get_relax_type)
   (int& relax_type);

#include <scf_relax_F.H>
