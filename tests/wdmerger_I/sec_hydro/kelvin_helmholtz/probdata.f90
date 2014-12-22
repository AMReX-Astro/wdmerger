module probdata_module

  use network

  ! refinement criteria
  double precision ::    denerr,   dengrad
  double precision ::    enterr,   entgrad
  double precision ::    velerr,   velgrad
  double precision ::   temperr,  tempgrad
  double precision ::  presserr, pressgrad
  double precision ::    raderr,   radgrad
  integer          ::  max_denerr_lev,   max_dengrad_lev
  integer          ::  max_enterr_lev,   max_entgrad_lev
  integer          ::  max_velerr_lev,   max_velgrad_lev
  integer          ::  max_temperr_lev,  max_tempgrad_lev
  integer          ::  max_presserr_lev, max_pressgrad_lev
  integer          ::  max_raderr_lev,   max_radgrad_lev

  ! grid info
  double precision ::  center(3)

  ! problem setup data
  double precision :: rho1, rho2, pressure

  ! problem number
  integer :: problem

  ! uniform flow speed
  double precision :: bulk_velocity

end module probdata_module
