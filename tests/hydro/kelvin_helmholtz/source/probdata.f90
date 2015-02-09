module probdata_module

  use network

  ! problem setup data
  double precision :: rho1, rho2, pressure

  ! problem number
  integer :: problem

  ! uniform flow speed
  double precision :: bulk_velocity

end module probdata_module
