module probdata_module

  use network
  use eos_module

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
  character (len=80) :: model_P_name
  character (len=80) :: model_S_name
  
  double precision, allocatable :: model_P_r(:), model_P_state(:,:)
  double precision, allocatable :: model_S_r(:), model_S_state(:,:)
  
  double precision :: mass_P, mass_S
  double precision :: radius_P_initial, radius_S_initial

  type (eos_t) :: ambient_state
  
  logical :: interp_temp

  integer :: npts_model

  integer :: nsub

  ! inertial reference frame flag
  logical :: inertial

  ! damping
  logical          :: damping
  double precision :: damping_alpha

  ! relaxation
  logical          :: do_relax
  double precision :: relax_tau

  ! grid info
  double precision, save ::  center(3)

  ! binary properties
  double precision :: a_P_initial, a_S_initial
  
  double precision, dimension(3) :: center_P_initial, center_S_initial

  integer :: star_axis

  ! tagging

  
  
end module probdata_module
