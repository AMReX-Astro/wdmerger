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
  double precision, save ::  starBuffer, boundaryBuffer
  
  ! model info
  character (len=80), save :: model_P_name
  character (len=80), save :: model_S_name
  
  double precision, allocatable, save :: model_P_r(:), model_P_state(:,:)
  double precision, allocatable, save :: model_S_r(:), model_S_state(:,:)
  
  integer, save :: npts_model_P, npts_model_S
  
  double precision, save :: mass_P_initial, mass_S_initial
  double precision, save :: radius_P_initial, radius_S_initial
  
  double precision, save :: dens_ambient, temp_ambient
  double precision, save :: xn_ambient(nspec)
  double precision, save :: eint_ambient, pres_ambient

  logical, save :: interp_temp

  integer, save :: nsub

  ! inertial reference frame flag
  logical, save :: inertial

  ! damping
  logical, save    :: damping
  double precision :: damping_alpha

  ! grid info
  double precision, save ::  center(3)

  ! binary properties
  double precision, save :: period
  double precision, save :: a_P_initial, a_S_initial
  
  double precision, dimension(3), save :: center_P_initial, center_S_initial

  ! tagging

  
  
end module probdata_module
