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
  use sim_local_interface, ONLY : sim_interpolate1dWd
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
  real :: radius_s, radius_p, ign_dist, idistsq, min_distsq

  real, dimension(SPECIES_BEGIN:SPECIES_END) :: massFraction

  real :: velx, vely, velz

  real :: xc12initial, xne22initial

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
  do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

           !-----------------------------------------------
           !  determine state of material at this radius if unburned
           !  from external 1-d hydrostatic model
           !-----------------------------------------------

           ! calculate distances from centers of two stars
           ! assume primary is centered at (-sim_wdp_a,0,0)
           ! assume secondary is centered at (sim_wds_a,0,0)

           radius_p = ( (iCoords(i) + sim_wdp_a)**2 + jCoords(j)**2 + kCoords(k)**2 )**0.5d0
           radius_s = ( (iCoords(i) - sim_wds_a)**2 + jCoords(j)**2 + kCoords(k)**2 )**0.5d0

           ! First set everything to fluff.
           call sim_interpolate1dWd(radius_p, deltas(IAXIS), state(EOS_DENS), state(EOS_TEMP), xc12initial, &
                                     xne22initial, 0)

           ! Then if within the two stars, interpolate and set values.
           call sim_interpolate1dWd(radius_p, deltas(IAXIS), state(EOS_DENS), state(EOS_TEMP), xc12initial, &
                                     xne22initial, 1) ! Primary
           call sim_interpolate1dWd(radius_s, deltas(IAXIS), state(EOS_DENS), state(EOS_TEMP), xc12initial, &
                                     xne22initial, 2) ! Secondary

           ! Now calculate the rotational velocity from the initial period.

           velx = jCoords(j)*(-2.0*PI/sim_binaryPeriod)
           vely = iCoords(i)*( 2.0*PI/sim_binaryPeriod)
           velz = sim_smallu

           ! Fill the mass fraction arrays.

           massFraction(:) = sim_smallx
           massFraction(C12_SPEC)  = max(xc12initial, sim_smallx)
           massFraction(NE22_SPEC) = max(xne22initial, sim_smallx)
           massFraction(O16_SPEC)  = max(1.0d0 - xc12initial - xne22initial, sim_smallx)

           call Eos(MODE_DENS_TEMP, 1, state, massFraction)

           GAME = state(EOS_PRES)/(state(EOS_DENS)*state(EOS_EINT))+1.0

           !-----------------------------------------------
           !  Now store all this info on the grid
           !-----------------------------------------------
           cell(IAXIS) = i
           cell(JAXIS) = j
           cell(KAXIS) = k
           call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, cell, state(EOS_DENS))
           call Grid_putPointData(blockId, CENTER, TEMP_VAR, EXTERIOR, cell, state(EOS_TEMP))

           call Grid_putPointData(blockId, CENTER, VELX_VAR, EXTERIOR, cell, velx)
           call Grid_putPointData(blockId, CENTER, VELY_VAR, EXTERIOR, cell, vely)
           call Grid_putPointData(blockId, CENTER, VELZ_VAR, EXTERIOR, cell, velz)

           call Grid_putPointData(blockId, CENTER, ENER_VAR, EXTERIOR, cell, state(EOS_EINT) + 0.5 * (velx**2 + vely**2 + velz**2) )
           call Grid_putPointData(blockId, CENTER, EINT_VAR, EXTERIOR, cell, state(EOS_EINT))
           call Grid_putPointData(blockId, CENTER, PRES_VAR, EXTERIOR, cell, state(EOS_PRES))
           call Grid_putPointData(blockId, CENTER, GAMC_VAR, EXTERIOR, cell, state(EOS_GAMC))
           call Grid_putPointData(blockId, CENTER, GAME_VAR, EXTERIOR, cell, GAME)

           do n = SPECIES_BEGIN, SPECIES_END
             call Grid_putPointData(blockId, CENTER, n, EXTERIOR, cell, massFraction(n))
           enddo

        enddo
     enddo
  enddo
  
  deallocate(iCoords)
  deallocate(jCoords)
  deallocate(kCoords)

  return
  
end subroutine Simulation_initBlock






