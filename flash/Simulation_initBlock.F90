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
  use Flame_interface, ONLY : Flame_getProfile, Flame_rhJump, Flame_rhJumpReactive
  use bn_paraInterface, ONLY : bn_paraFuelAshProperties
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

!acc            write(6,*)'after interpolating, density =', unburned_state(EOS_DENS)

           call bn_paraFuelAshProperties(xc12initial, xne22initial, ye_f, ye_a, yi_f, yi_a, qbar_f, qbar_a)

!acc           write(6,*)'after interpolating, dens temp xc12 xne22 ='
!acc           print *, unburned_state(EOS_DENS), unburned_state(EOS_TEMP)
!acc           print *, xc12initial, xne22initial

           unburned_state(EOS_ABAR) = 1.0/yi_f
           unburned_state(EOS_ZBAR) = ye_f/yi_f

!           print *, unburned_state(EOS_ABAR), unburned_state(EOS_ZBAR)

!   now the velocity from the initial period.
!
           velx = jCoords(j)*(-2.0*PI/sim_binaryPeriod)
           vely = iCoords(i)*(2.0*PI/sim_binaryPeriod)

           if (.not. sim_ignite) then
              ! no burned material, only unburned
              state(:) = unburned_state(:)
              call Eos(MODE_DENS_TEMP, 1, state)
              flam = 0.0
              ye = state(EOS_ZBAR)/state(EOS_ABAR)
              dyi_qn = 0.0
              dqbar_qn = 0.0
              enuc = 0.0
           else
              !-----------------------------------------------
              ! initialize, including a burned region
              !-----------------------------------------------

              ! find distance from flame surface (positive is in front of flame)

              ! find the distance to the closest ignition point specified
              ! will default to a spherical region around it
              min_distsq = HUGE(1.0)
              do n = 1, sim_ign_numpnts
                 idistsq = (iCoords(i) - sim_ignX(n))**2
                 if (NDIM >= 2) idistsq = idistsq + (jCoords(j) - sim_ignY(n))**2
                 if (NDIM == 3) idistsq = idistsq + (kCoords(k) - sim_ignZ(n))**2
                 if (idistsq < min_distsq) then
                    min_distsq = idistsq
                    flame_radius = sim_ignR(n)
                 endif
              enddo
              ign_dist = sqrt(min_distsq)

              if (sim_ignMpole .and. sim_ignSin)  &
                 call Driver_abortFlash("Simulation_initBlock: multipole and sinusoidal ignition are exclusive")

              if (sim_ignMpole) then
                 ! only CYLINDRICAL is implemented right now (check is in Simulation_init)
                 costheta = (jCoords(j)-sim_ignY(1))/ign_dist
                 ! Legendre polynomial
                 P_lm1 = 1.0
                 P_l = costheta
                 do l = 1, sim_ignMPoleMaxL-1
                    if ( l >= sim_ignMPoleMinL ) then
                       flame_radius = flame_radius + mp_A(l)*P_l
                    endif
                    ! calculate next polynomial
                    hold = P_l
                    P_l = ( (2*l+1)*costheta*P_l - l*P_lm1 )/real(l+1)
                    P_lm1 = hold
                 enddo
                 flame_radius = flame_radius + mp_A(sim_ignMpoleMaxL)*P_l
              endif

              if (sim_ignSin) then
                 ! sinusoidal (in theta) border between burned and unburned region
                 ! sim_ignitionSinN = number of periods between theta = 0 and pi
                    ! this is not a good place for this warning (reported for every block)
                    !... need to move
                    !if (sim_meshMe==MASTER_PE) then
                    !   write(6,*) "warning: sinusoidal ignition reduces to spherical in 1 dimension"
                    !endif
                 if (NDIM == 2) then
                    theta = acos( (jCoords(j)-sim_ignY(1))/ign_dist)
                    flame_radius = flame_radius - sim_ignSinA*cos(2*theta*sim_ignSinN)
                 else if (NDIM == 3) then
                    theta = acos( (kCoords(k)-sim_ignZ(1))/ign_dist)
                    flame_radius = flame_radius - sim_ignSinA*cos(2*theta*sim_ignSinN)
                 endif
              endif


              fsurf_distance = ign_dist - flame_radius

              ! determine local state in this zone
              ! assume deltas are equal
              if ( (fsurf_distance-0.5*deltas(IAXIS)) > 1.5*sim_laminarWidth ) then
                 ! whole cell unburned material
                 state(:) = unburned_state(:)
                 call Eos(MODE_DENS_TEMP, 1, state)
                 flam = 0.0
                 ye = state(EOS_ZBAR)/state(EOS_ABAR)
                 dyi_qn = 0.0
                 dqbar_qn = 0.0
                 enuc = 0.0
              else if ( (fsurf_distance+0.5*deltas(IAXIS)) < -1.5*sim_laminarWidth ) then
                 ! fully burned to NSE
                 call Flame_rhJumpReactive(unburned_state, qbar_f, state, dqbar_qn, MODE_DENS_TEMP)
                 flam = 1.0
                 dyi_qn   = 1.0/state(EOS_ABAR)
                 ye    = dyi_qn*state(EOS_ZBAR)
                 enuc = 0.0
              else
                 ! partially burned
                 ! at least one cell will fall here (necessary to get initial refinement right)
                 call Flame_getProfile(fsurf_distance, flam)
                 ! calculate properties of NSE final state
                 call Flame_rhJumpReactive(unburned_state, qbar_f, nse_state, qbar_nse, MODE_DENS_TEMP)

                 ! calculate propertise for partially burned material
                 ! note, in fact ye_f and ye_a should be equal
                 yi = yi_f*(1.0-flam) + (1.0/nse_state(EOS_ABAR))*flam
                 ye = ye_f*(1.0-flam) + (nse_state(EOS_ZBAR)/nse_state(EOS_ABAR))*flam
                 state(:) = unburned_state(:)
                 state(EOS_ABAR) = 1.0/yi
                 state(EOS_ZBAR) = ye/yi

                 ! put this in pressure equilibrium with unburned material
                 call Flame_rhJump(unburned_state, state, flam*(qbar_nse-qbar_f)*cgsMeVperAmu, 0.0, MODE_DENS_TEMP)

                 dyi_qn   = flam * 1.0/nse_state(EOS_ABAR)
                 dqbar_qn = flam * qbar_nse

                 ! to trigger refinement
                 enuc = 1.1*sim_refNogenEnucThresh
              endif

           endif ! sim_ignite
           

           !-----------------------------------------------
           !  Now store all this info on the grid
           !-----------------------------------------------
           cell(IAXIS) = i
           cell(JAXIS) = j
           cell(KAXIS) = k
           call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, cell, state(EOS_DENS))
           call Grid_putPointData(blockId, CENTER, TEMP_VAR, EXTERIOR, cell, state(EOS_TEMP))

           call Grid_putPointData(blockId, CENTER, FLAM_MSCALAR, EXTERIOR, cell, flam)

           call Grid_putPointData(blockId, CENTER, CI_MSCALAR, EXTERIOR, cell, xc12initial)
           call Grid_putPointData(blockId, CENTER, NEI_MSCALAR, EXTERIOR, cell, xne22initial)

           call Grid_putPointData(blockId, CENTER, PHFA_MSCALAR, EXTERIOR, cell, flam)
           call Grid_putPointData(blockId, CENTER, PHAQ_MSCALAR, EXTERIOR, cell, flam)
           call Grid_putPointData(blockId, CENTER, PHQN_MSCALAR, EXTERIOR, cell, flam)

           call Grid_putPointData(blockId, CENTER, YE_MSCALAR, EXTERIOR, cell, ye)
           call Grid_putPointData(blockId, CENTER, DYQN_MSCALAR, EXTERIOR, cell, dyi_qn)
           call Grid_putPointData(blockId, CENTER, DQQN_MSCALAR, EXTERIOR, cell, dqbar_qn)

           call Grid_putPointData(blockId, CENTER, ENUC_VAR, EXTERIOR, cell, enuc)

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






