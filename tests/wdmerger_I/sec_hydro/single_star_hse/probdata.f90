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
  character (len=80), save :: model_name
  
  double precision, allocatable, save :: model_star_r(:), model_star_state(:,:)
  
  integer, save :: npts_model_star
  
  double precision, save :: total_mass
  double precision, save :: radius
  
  double precision, save :: dens_ambient, temp_ambient
  double precision, save :: xn_ambient(nspec)
  double precision, save :: eint_ambient, pres_ambient

  ! grid info
  double precision, save ::  center(3)
  
end module probdata_module
