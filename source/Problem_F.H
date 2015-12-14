/* problem-specific stuff goes here */

extern "C"
{

  void problem_checkpoint(int * int_dir_name, int * len);
  
  void problem_restart(int * int_dir_name, int * len);

  void wdcom
    (BL_FORT_FAB_ARG_3D(rho), 
     BL_FORT_FAB_ARG_3D(xmom), 
     BL_FORT_FAB_ARG_3D(ymom), 
     BL_FORT_FAB_ARG_3D(zmom),
     BL_FORT_FAB_ARG_3D(vol),
     const int* lo, const int* hi,
     const Real* dx, const Real* time,
     Real* com_p_x, Real* com_p_y, Real* com_p_z,
     Real* com_s_x, Real* com_s_y, Real* com_s_z,
     Real* vel_p_x, Real* vel_p_y, Real* vel_p_z,
     Real* vel_s_x, Real* vel_s_y, Real* vel_s_z,
     Real* m_p, Real* m_s);

  void ca_volumeindensityboundary
    (BL_FORT_FAB_ARG_3D(rho), BL_FORT_FAB_ARG_3D(vol),
     const int* lo, const int* hi, const Real* dx,  
     Real* volp, Real* vols, Real* rho_cutoff);

  void get_single_star(int& flag);

  void set_star_data
    (const Real com_p[], const Real com_s[], 
     const Real vel_p[], const Real vel_s[], 
     const Real* mass_p, const Real* mass_s);

  void get_star_data
    (Real com_p[], Real com_s[], 
     Real vel_p[], Real vel_s[], 
     Real* mass_p, Real* mass_s);

  void quadrupole_tensor_double_dot
    (BL_FORT_FAB_ARG_3D(rho), 
     BL_FORT_FAB_ARG_3D(xmom), 
     BL_FORT_FAB_ARG_3D(ymom), 
     BL_FORT_FAB_ARG_3D(zmom),
     BL_FORT_FAB_ARG_3D(gx),
     BL_FORT_FAB_ARG_3D(gy),
     BL_FORT_FAB_ARG_3D(gz),
     BL_FORT_FAB_ARG_3D(vol),
     const int* lo, const int* hi,
     const Real* dx, const Real* time,
     Real* Qtt);

  void gw_strain_tensor
    (Real* h_plus_1, Real* h_cross_1, 
     Real* h_plus_2, Real* h_cross_2,
     Real* h_plus_3, Real* h_cross_3,
     Real* Qtt, const Real* time);

  void set_period(const Real* period);

  void update_center(const Real* time);

  void get_omega_vec(const Real omega[], const Real& time);

  void get_axes(int& axis_1, int& axis_2, int& axis_3);

  void get_do_scf_initial_models(int& do_scf_initial_models);

}

#include <scf_relax_F.H>
