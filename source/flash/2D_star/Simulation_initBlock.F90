!!****if* source/Simulation/SimulationMain/Flame3StageNoise/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer, intent(IN) :: blockid)
!!
!! DESCRIPTION
!!  Initialization of white dwarf and flame for Type Ia simulation
!!  This initialization is intended for modelling DDTs
!!  The (1d) hydrostatic star was read from a file
!!
!! ARGUMENTS
!!
!!   blockid : ID of block in current processor
!!
!!
!!
!!***

subroutine Simulation_initBlock(blockID)
  
  use Simulation_data
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
    Grid_getBlkBoundBox, Grid_getDeltas, Grid_putPointData, &
    Grid_getCellCoords
  use Eos_interface, ONLY : Eos

  implicit none
#include "Flash.h"
#include "constants.h"
#include "Eos.h"

  
  integer, intent(in) :: blockID

  integer :: i, j, k

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer, dimension(MDIM) :: cell
  integer :: isizeGC, jsizeGC, ksizeGC
  real, allocatable, dimension(:) :: iCoords, jCoords, kCoords
  real, dimension(MDIM) :: deltas

  real, dimension(EOS_NUM) :: state
  integer :: n

  real, dimension(SPECIES_BEGIN:SPECIES_END) :: massFraction

  real :: x, y, r

  real :: game

!==============================================================================

  ! get essential info about this block - index limits and cell coordinates
  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

  isizeGC = blkLimitsGC(HIGH,IAXIS)
  allocate(iCoords(isizeGC))
  jsizeGC = blkLimitsGC(HIGH,JAXIS)
  allocate(jCoords(jsizeGC))
  ksizeGC = blkLimitsGC(HIGH,KAXIS)
  allocate(kCoords(ksizeGC))
  call Grid_getCellCoords(IAXIS,blockID,CENTER,.true.,iCoords,isizeGC)
  call Grid_getCellCoords(JAXIS,blockID,CENTER,.true.,jCoords,jsizeGC)
  call Grid_getCellCoords(KAXIS,blockID,CENTER,.true.,kCoords,ksizeGC)

  call Grid_getDeltas(blockID, deltas)

  !-----------------------------------------------
  ! loop over all zones and init
  !-----------------------------------------------
  do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
     y = jCoords(j)
     do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
        x = iCoords(i)

        !-----------------------------------------------
        !  determine whether or not we are in the star
        !-----------------------------------------------

        r = (x**2 + y**2)**0.5d0

        if ( r .le. sim_star_radius ) then

          state(EOS_DENS) = sim_densStar
          state(EOS_TEMP) = sim_tempStar

        ! Otherwise, we're in the ambient medium

        else

          state(EOS_DENS) = sim_densFluff
          state(EOS_TEMP) = sim_tempFluff

        endif

        ! Fill the mass fraction array with only one species.

        massFraction(:) = 0.d0
        massFraction(SPECIES_BEGIN) = 1.d0

        call Eos(MODE_DENS_TEMP, 1, state, massFraction)

        GAME = state(EOS_PRES)/(state(EOS_DENS)*state(EOS_EINT))+1.0

        !-----------------------------------------------
        !  Now store all this info on the grid
        !-----------------------------------------------
        cell(IAXIS) = i
        cell(JAXIS) = j
        call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, cell, state(EOS_DENS))
        call Grid_putPointData(blockId, CENTER, TEMP_VAR, EXTERIOR, cell, state(EOS_TEMP))

        call Grid_putPointData(blockId, CENTER, VELX_VAR, EXTERIOR, cell, sim_x_velocity)
        call Grid_putPointData(blockId, CENTER, VELY_VAR, EXTERIOR, cell, sim_y_velocity)

        call Grid_putPointData(blockId, CENTER, ENER_VAR, EXTERIOR, cell, state(EOS_EINT) + &
             0.5 * (sim_x_velocity**2 + sim_y_velocity**2) )
        call Grid_putPointData(blockId, CENTER, EINT_VAR, EXTERIOR, cell, state(EOS_EINT))
        call Grid_putPointData(blockId, CENTER, PRES_VAR, EXTERIOR, cell, state(EOS_PRES))
        call Grid_putPointData(blockId, CENTER, GAMC_VAR, EXTERIOR, cell, state(EOS_GAMC))
        call Grid_putPointData(blockId, CENTER, GAME_VAR, EXTERIOR, cell, GAME)

        do n = SPECIES_BEGIN, SPECIES_END
          call Grid_putPointData(blockId, CENTER, n, EXTERIOR, cell, massFraction(n))
        enddo

     enddo
  enddo
  
  deallocate(iCoords)
  deallocate(jCoords)

  return
  
end subroutine Simulation_initBlock






