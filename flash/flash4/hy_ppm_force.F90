!!****if* source/physics/Hydro/HydroMain/split/PPM/hy_ppm_force
!!
!! NAME
!!
!!  hy_ppm_force
!!
!!
!! SYNOPSIS
!!
!!  call hy_ppm_force(integer(IN) :: numCells,
!!                    integer(IN) :: numIntCells, 
!!                    integer(IN) :: j, 
!!                    integer(IN) :: k, 
!!                    integer(IN) :: igeom, 
!!                    real(IN)    :: coord(numCells), 
!!                    real(IN)    :: radialCoord(numCells), 
!!                    real(IN)    :: thirdCoord(numCells), 
!!                    real(IN)    :: u(numCells), 
!!                    real(IN)    :: ut(numCells), 
!!                    real(IN)    :: utt(numCells), 
!!                    real(OUT)   :: fict(numCells))
!!
!! DESCRIPTION
!!
!!  compute fictitious forces (centrifugal and Coriolis) for cylindrical 
!!  and spherical coordinates
!!
!!  When a particle moves in a non-Cartesian coordinate system, these forces
!!  act to keep it moving in a straight line, in the absence of any other
!!  external forces.  They are computed by setting the Cartesian accelerations
!!  equal to 0, and transforming coordinate systems to whatever geometry 
!!  we are presently in, and solving for the accelerations in that system.
!!
!!  We note that in angular coordinates, the velocities and accelerations
!!  have physical units of cm/s and cm/s**2 respectively, instead of just
!!  1/s and 1/s**2.  
!!
!!
!! ARGUMENTS
!!
!!  numCells --    The maximum number of zones in any sweep direction
!!
!!  numIntCells -- The number of zones in the current sweep direction
!!
!!  j --         the zone index for the first transverse coord
!!  k --         the zone index for the second transverse coord
!!
!!               note, j and k are always in x, y, z, order, so if it
!!               is the x sweep, j refers to the y direction and k
!!               refers to the z direction.  If it is the y sweep, j
!!               refers to the x direction and k refers to the z
!!               direction, and for the z sweep, j refers to the x
!!               direction and k refers to the y direction
!!
!!
!!  igeom --     The geometry of the current coordinate direction
!!
!!  coord --     The zone coordinates in the current coordinate direction
!!
!!  radialCoord - The radial coordinate, no matter whether it's primary, second, or third.
!!               Actually, this is currently always the IAXIS coordinate.
!!
!!  thirdCoord - The second transverse coordinate to the sweep direction
!!               for an x sweep: z coordinates; for a y sweep: z coord; for a z sweep: ycoord
!!
!!  u, ut, utt-- The velocity in the current direction and the transverse
!!               directions
!!
!!  fict --      The resulting fictitious force
!!
!!
!! NOTES
!!
!!  This routine should be generalized so that any hydro solver can use it.
!!  
!!***



subroutine hy_ppm_force(sweepDir, numCells, numIntCells, j, k,  igeom, coord, &
             radialCoord, thirdCoord, u, ut, utt, fict)


  use Simulation_data, ONLY : sim_binaryPeriod
  use Hydro_data, ONLY : hy_geometry
  use Driver_interface, ONLY : Driver_abortFlash

!!$  use Hydro_data, ONLY : ilo,ihi,jlo,jhi,klo,khi,&
!!$                         ilo_gc,ihi_gc,jlo_gc,jhi_gc,klo_gc,khi_gc,&
!!$                         numCells

 
  implicit none

#include "constants.h"

  integer,INTENT(in) :: j, k, numCells, numIntCells, igeom, sweepDir

  real, DIMENSION(numCells),INTENT(in) :: coord, u, ut, utt
  real, DIMENSION(numCells),INTENT(in) :: radialCoord, thirdCoord
  real, DIMENSION(numCells),INTENT(out) :: fict

  integer :: i, numIntCells8

  real :: omega

  omega = 2.0 * PI / sim_binaryPeriod

  numIntCells8 = numIntCells + 8

  do i = 1, numIntCells8
    select case(sweepDir)
    case(1)
      fict(i) =  2.0e0 * omega * ut(i) + omega * omega * coord(i)
    case(2)
      fict(i) = -2.0e0 * omega * ut(i) + omega * omega * coord(i)
    case(3)
      fict(i) =  0.0e0
    end select
  enddo
  
  fict(i) = 0.0e0

  return
end subroutine hy_ppm_force



