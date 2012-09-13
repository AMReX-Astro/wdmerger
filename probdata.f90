module probdata_module

  use network

  ! refinement criteria
  double precision, save ::    denerr,   dengrad
  double precision, save ::    enterr,   entgrad
  double precision, save ::    velerr,   velgrad
  double precision, save ::   temperr,  tempgrad
  double precision, save ::  presserr, pressgrad
  double precision, save ::    raderr,   radgrad
  integer         , save ::  max_denerr_lev,   max_dengrad_lev
  integer         , save ::  max_enterr_lev,   max_entgrad_lev
  integer         , save ::  max_velerr_lev,   max_velgrad_lev
  integer         , save ::  max_temperr_lev,  max_tempgrad_lev
  integer         , save ::  max_presserr_lev, max_pressgrad_lev
  integer         , save ::  max_raderr_lev,   max_radgrad_lev
  
  ! model info
  character (len=80), save :: model_A_name
  character (len=80), save :: model_B_name
  
  double precision, allocatable, save :: model_A_r(:), model_A_state(:,:)
  double precision, allocatable, save :: model_B_r(:), model_B_state(:,:)
  
  integer, save :: npts_model_A, npts_model_B
  
  double precision, save :: mass_A, mass_B
  double precision, save :: radius_A, radius_B
  
  double precision, save :: dens_ambient, temp_ambient
  double precision, save :: xn_ambient(nspec)
  double precision, save :: eint_ambient, pres_ambient

  ! grid info
  double precision, save ::  center(3)

  ! binary properties
  double precision, save :: period
  double precision, save :: a_A, a_B
  
  double precision, save :: x_cen_A, y_cen_A, z_cen_A
  double precision, save :: x_cen_B, y_cen_B, z_cen_B
  
end module probdata_module
