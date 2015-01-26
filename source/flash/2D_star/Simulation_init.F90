!!***
!!
!! NAME
!!
!!  Simulation_init
!!
!! SYNOPSIS
!!
!!  Simulation_init()
!!
!! DESCRIPTION
!!  Initialize private data for the 3-stage flame test setup
!!
!! ARGUMENTS
!!
!!
!!
!!***

subroutine Simulation_init()

  use Simulation_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface, only : Logfile_stampMessage
  use Grid_interface, only : Grid_getGeometry
  use Logfile_interface, only : Logfile_stampMessage
  use Driver_interface, only : Driver_abortFlash
  use PhysicalConstants_interface, ONLY: PhysicalConstants_get

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"

  !--------------------------------------------------------
  !  initialize runtime parameters
  !--------------------------------------------------------

  call RuntimeParameters_get('x_velocity', sim_x_velocity)
  call RuntimeParameters_get('y_velocity', sim_y_velocity)
  call RuntimeParameters_get('star_radius', sim_star_radius)

  call RuntimeParameters_get('dens_star', sim_densStar)
  call RuntimeParameters_get('temp_star', sim_tempStar)  

  call RuntimeParameters_get('dens_fluff', sim_densFluff)
  call RuntimeParameters_get('temp_fluff', sim_tempFluff)

  call Driver_getMype(MESH_COMM, sim_meshMe)
  
end subroutine Simulation_init
