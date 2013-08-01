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

  real, dimension(EOS_NUM) :: state, unburned_state, nse_state
  integer :: n
  real :: radius_s, radius_p, ign_dist, idistsq, min_distsq

! rotation variables
  real :: velx, vely

  real :: P_l, P_lm1, hold, flame_radius, costheta
  real :: theta
  integer :: l
  real :: fsurf_distance
  real :: flam, xc12initial, xne22initial, ye, dyi_qn, dqbar_qn
  real :: yi
  real :: qbar_nse
  real :: ye_f, ye_a, yi_f, yi_a, qbar_f, qbar_a
  real :: enuc

  real :: cgsMeVperAmu = 9.6485e17

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

        !first set everything to fluff
           call sim_interpolate1dWd(radius_p, deltas(IAXIS), unburned_state(EOS_DENS), unburned_state(EOS_TEMP), xc12initial, &
                                     xne22initial, 0) !primary

        !then if within the two stars, interpolate and set values.
           call sim_interpolate1dWd(radius_p, deltas(IAXIS), unburned_state(EOS_DENS), unburned_state(EOS_TEMP), xc12initial, &
                                     xne22initial, 1) !primary
           call sim_interpolate1dWd(radius_s, deltas(IAXIS), unburned_state(EOS_DENS), unburned_state(EOS_TEMP), xc12initial, &
                                     xne22initial, 2) !secondary

           unburned_state(EOS_ABAR) = 1.0d0 / ( (1.0d0 / 12.0d0) + (1.0d0 / 16.0d0) )
           unburned_state(EOS_ZBAR) = 1.0d0 / ( (1.0d0 / 6.0d0) + (1.0d0 / 8.0d0) )

!           print *, unburned_state(EOS_ABAR), unburned_state(EOS_ZBAR)
!   now the velocity from the initial period.
!
           velx = jCoords(j)*(-2.0*PI/sim_binaryPeriod)
           vely = iCoords(i)*(2.0*PI/sim_binaryPeriod)

              state(:) = unburned_state(:)
              call Eos(MODE_DENS_TEMP, 1, state, [0.5d0, 0.5d0])

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
           call Grid_putPointData(blockId, CENTER, VELZ_VAR, EXTERIOR, cell, 0.0)

           !  usually I would just call the EOS, but we happen to have all this data
           !  so we'll just put it in.
           call Grid_putPointData(blockId, CENTER, ENER_VAR, EXTERIOR, cell, state(EOS_EINT))
           call Grid_putPointData(blockId, CENTER, EINT_VAR, EXTERIOR, cell, state(EOS_EINT))
           call Grid_putPointData(blockId, CENTER, PRES_VAR, EXTERIOR, cell, state(EOS_PRES))
           call Grid_putPointData(blockId, CENTER, GAMC_VAR, EXTERIOR, cell, state(EOS_GAMC))
           call Grid_putPointData(blockId, CENTER, GAME_VAR, EXTERIOR, cell, &
                                       state(EOS_PRES)/(state(EOS_DENS)*state(EOS_EINT))+1.0)
        enddo
     enddo
  enddo
  
  deallocate(iCoords)
  deallocate(jCoords)
  deallocate(kCoords)

  return
  
end subroutine Simulation_initBlock






