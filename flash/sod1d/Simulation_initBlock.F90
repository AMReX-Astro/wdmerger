!!****if* source/Simulation/SimulationMain/Sod/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer(IN) :: blockID) 
!!                       
!!
!!
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.  This version sets up the Sod shock-tube
!!  problem.
!!
!!  Reference: Sod, G. A., 1978, J. Comp. Phys., 27, 1
!!
!! 
!! ARGUMENTS
!!
!!  blockID -           the number of the block to update
!!
!! PARAMETERS
!!
!!  sim_rhoLeft    Density in the left part of the grid
!!  sim_rhoRight   Density in the right part of the grid
!!  sim_TLeft      Temperature in the left part of the grid
!!  sim_TRight     Temperature in the right part of the grid
!!  sim_uLeft      fluid velocity in the left part of the grid
!!  sim_uRight     fluid velocity in the right part of the grid
!!  sim_xangle     Angle made by diaphragm normal w/x-axis (deg)
!!  sim_ yangle    Angle made by diaphragm normal w/y-axis (deg)
!!  sim_posnR      Point of intersection between the shock plane and the x-axis
!!
!!
!!***

subroutine Simulation_initBlock(blockID)

#include "constants.h"
#include "Flash.h"
#include "Eos.h"

  use Simulation_data, ONLY: sim_rhoLeft,  sim_TLeft, sim_uLeft, sim_rhoRight, sim_TRight, sim_uRight, &
     &  sim_smallX, sim_gamma, sim_smallP, xmin, xmax
     
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
    Grid_getCellCoords, Grid_putPointData
  use Eos_interface, ONLY : Eos


  implicit none

  ! compute the maximum length of a vector in each coordinate direction 
  ! (including guardcells)
  
  integer, intent(in) :: blockID
  

  integer :: i, j, k, n
  integer :: iMax, jMax, kMax
  
  real, dimension(SPECIES_BEGIN:SPECIES_END) ::  massFraction

  real :: xx, yy,  zz, xxL, xxR
  
  real :: center
  

  real,allocatable, dimension(:) ::xCenter,xLeft,xRight,yCoord,zCoord

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer :: sizeX,sizeY,sizeZ
  integer, dimension(MDIM) :: axis

  real, dimension(EOS_NUM) :: state

  real :: rhoZone, velxZone, velyZone, velzZone, tempZone, & 
       eintZone, enerZone, ekinZone
  
  logical :: gcell = .true.   
  
  ! get the integer index information for the current block
  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)
  
  sizeX = blkLimitsGC(HIGH,IAXIS)
  allocate(xLeft(sizeX))
  allocate(xRight(sizeX))
  allocate(xCenter(sizeX))
  allocate(yCoord(sizeY))
  allocate(zCoord(sizeZ))
  xCenter = 0.0
  xLeft = 0.0
  xRight = 0.0
  yCoord = 0.0
  zCoord = 0.0

  call Grid_getCellCoords(IAXIS, blockId, LEFT_EDGE, gcell, xLeft, sizeX)
  call Grid_getCellCoords(IAXIS, blockId, CENTER, gcell, xCenter, sizeX)
  call Grid_getCellCoords(IAXIS, blockId, RIGHT_EDGE, gcell, xRight, sizeX)

!------------------------------------------------------------------------------

! Loop over cells in the block.  For each, compute the physical position of 
! its left and right edge and its center as well as its physical width.  
! Then decide which side of the initial discontinuity it is on and initialize 
! the hydro variables appropriately.

  center = (xmin + xmax) / 2.0d0
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           
           ! get the cell center, left, and right positions in x
           xx  = xCenter(i)
           
           xxL = xLeft(i)
           xxR = xRight(i)
          
           
           ! initialize cells to the left of the initial shock.
           if (xxR <= center) then
              tempZone = sim_TLeft

              rhoZone = sim_rhoLeft
              velxZone = sim_uLeft 

              ! initialize cells which straddle the shock.  Treat them as though 1/2 of 
              ! the cell lay to the left and 1/2 lay to the right.
           elseif ((xxL < center) .and. (xxR > center)) then              
              tempZone = 0.5 * (sim_TLeft+sim_TRight)
            
              rhoZone = 0.5 * (sim_rhoLeft+sim_rhoRight)
              velxZone = 0.5 *(sim_uLeft+sim_uRight)
              
              ! initialize cells to the right of the initial shock.
           else              
              tempZone = sim_TRight
              
              rhoZone = sim_rhoRight
              velxZone = sim_uRight 
              
           endif

           state(:)        = 0.0d0
           state(EOS_DENS) = rhoZone
           state(EOS_TEMP) = tempZone
           state(EOS_ABAR) = 12.0d0
           state(EOS_ZBAR) = 6.0d0

           axis(IAXIS) = i

           !put in value of default species
           if (NSPECIES > 0) then
              call Grid_putPointData(blockID, CENTER, SPECIES_BEGIN, EXTERIOR, &
                   axis, 1.0e0-(NSPECIES-1)*sim_smallX)


              !if there is only 1 species, this loop will not execute
              do n = SPECIES_BEGIN+1,SPECIES_END
                 call Grid_putPointData(blockID, CENTER, n, EXTERIOR, &
                      axis, sim_smallX)
              enddo
           end if

           ! Compute the gas energy and set the gamma-values needed for the equation of 
           ! state.
           state(EOS_EKIN) = 0.5 * (velxZone**2)
 
           call Eos(MODE_DENS_TEMP,1,state,[1.0d0])

           enerZone = state(EOS_EINT) + state(EOS_EKIN)
           
           ! store the variables in the current zone via Grid put methods
           ! data is put stored one cell at a time with these calls to Grid_putData           

           call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, state(EOS_DENS))
           call Grid_putPointData(blockId, CENTER, TEMP_VAR, EXTERIOR, axis, state(EOS_TEMP))
           call Grid_putPointData(blockId, CENTER, PRES_VAR, EXTERIOR, axis, state(EOS_PRES))
           call Grid_putPointData(blockId, CENTER, VELX_VAR, EXTERIOR, axis, velxZone)

#ifdef ENER_VAR
           call Grid_putPointData(blockId, CENTER, ENER_VAR, EXTERIOR, axis, enerZone)   
#endif
#ifdef EINT_VAR
           call Grid_putPointData(blockId, CENTER, EINT_VAR, EXTERIOR, axis, state(EOS_EINT))   
#endif
#ifdef GAME_VAR          
           call Grid_putPointData(blockId, CENTER, GAME_VAR, EXTERIOR, axis, state(EOS_PRES)/(state(EOS_DENS)*state(EOS_EINT))+1.0)
#endif
#ifdef GAMC_VAR
           call Grid_putPointData(blockId, CENTER, GAMC_VAR, EXTERIOR, axis, state(EOS_GAMC))
#endif

        enddo

!! Cleanup!  Must deallocate arrays

  deallocate(xLeft)
  deallocate(xRight)
  deallocate(xCenter)
  deallocate(yCoord)
  deallocate(zCoord)

 
  return
end subroutine Simulation_initBlock
