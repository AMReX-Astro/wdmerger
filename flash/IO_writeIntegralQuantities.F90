!!
!! Replacement of default routine from
!!  source/IO/IOMain/IO_writeIntegralQuantities.F90
!! see that file for API spec and notes.
!!
!!
!! This includes various quantities important for type Ia simulations.
!! As a result, it makes various assumptions about the geometry when
!! calculating the energy terms.
!!
!! One slightly unusual addition is that the cumulative neutrino losses
!! are also summed and saved as an IO scaler here as well.  Each processor
!! keeps a running total of energy lost to neutrinos since the last time
!! this routine was called (generally once every two-sweep step) and and
!! these contributions are summed here.
!!
!! original routine: FLASH 3.0
!! modifications: Dean M. Townsley 2009
!!


subroutine IO_writeIntegralQuantities (isFirst, simTime)

  use IO_data, ONLY : io_restart, io_statsFileName, io_globalComm, io_globalMe
  use IO_interface, only : IO_setScalar
  use Grid_interface, ONLY : Grid_getListOfBlocks, &
    Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_getSingleCellVol, &
    Grid_releaseBlkPtr
  use Burn_data, only : c_l, N_A, m_p, m_n, m_e, bn_neutLoss, bn_neutLossThisProcStep
  use bn_paraInterface, only : bn_paraFuelAshProperties

  implicit none

#include "Flash_mpi.h"
#include "constants.h"
#include "Flash.h"
  
  real, intent(in) :: simTime

  integer, intent(in) :: isFirst

  integer :: lb, count
  
  integer :: funit = 99
  integer :: error
  
  character (len=MAX_STRING_LENGTH), save :: fname 
  
  integer :: blockList(MAXBLOCKS)

  integer :: blkLimits(HIGH, MDIM), blkLimitsGC(HIGH, MDIM)

  integer, parameter ::  nGlobalSum = 24          ! Number of globally-summed quantities
  real :: gsum(nGlobalSum) !Global summed quantities
  real :: lsum(nGlobalSum) !Global summed quantities
  real :: lmaxdens, gmaxdens
  real :: lminflamdens, gminflamdens

  integer :: i, j, k
  real :: dvol             !, del(MDIM)
  real, DIMENSION(:,:,:,:), POINTER :: solnData

  integer :: point(MDIM)
  integer :: ioStat

  real :: qbar, phi_fa, phi_aq, phi_qn, dens, dm, ekin
  real :: qbar_f, qbar_a, ye_f, ye_a, yi_f, yi_a

  real :: cgsMeVperGram   = 9.6485e17

  real :: totalNeutLossThisStep

  real :: yen, xni56

  ! Sum quantities over all locally held leaf-node blocks.
  gsum  = 0.
  lsum = 0.
  lmaxdens = 0.0
  lminflamdens = 1e10
  
  call Grid_getListOfBlocks(LEAF, blockList, count)
  
  do lb = 1, count
     !get the index limits of the block
     call Grid_getBlkIndexLimits(blockList(lb), blkLimits, blkLimitsGC)

     ! get a pointer to the current block of data
     call Grid_getBlkPtr(blockList(lb), solnData)

     ! Sum contributions from the indicated blkLimits of cells.
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              
              point(IAXIS) = i
              point(JAXIS) = j
              point(KAXIS) = k

              !! Get the cell volume for a single cell
              call Grid_getSingleCellVol(blockList(lb), EXTERIOR, point, dvol)

              dm = solnData(DENS_VAR,i,j,k)*dvol
     
              ! mass   
              lsum(1) = lsum(1) + dm

              ! momentum
              lsum(2) = lsum(2) + solnData(VELX_VAR,i,j,k)*dm
#ifdef VELY_VAR      
              lsum(3) = lsum(3) + solnData(VELY_VAR,i,j,k)*dm
#endif
#ifdef VELZ_VAR      
              lsum(4) = lsum(4) + solnData(VELZ_VAR,i,j,k)*dm
#endif

              ! thermal + kinetic energy
              lsum(5) = lsum(5) + solnData(ENER_VAR,i,j,k)*dm

              ! kinetic energy
              ekin = 0.5*solnData(VELX_VAR,i,j,k)**2
#ifdef VELY_VAR
              ekin = ekin + 0.5*solnData(VELY_VAR,i,j,k)**2
#endif
#ifdef VELZ_VAR
              ekin = ekin + 0.5*solnData(VELZ_VAR,i,j,k)**2
#endif
              lsum(6) = lsum(6) + ekin*dm

#ifdef EINT_VAR
              ! internal energy
              lsum(7) = lsum(7) + solnData(EINT_VAR,i,j,k)*dm
#endif

              ! gravitational binding energy
              lsum(8) = lsum(8) + 0.5*solnData(GPOT_VAR,i,j,k)*dm

              phi_fa = solnData(PHFA_MSCALAR,i,j,k)
              phi_aq = solnData(PHAQ_MSCALAR,i,j,k)
              phi_qn = solnData(PHQN_MSCALAR,i,j,k)

              ! nuclear energy
              ! defined as the mass-energy difference between the number
              ! of baryons given by dm/amu in the form of unbound, symmetric matter
              ! and the stuff we have here
              call bn_paraFuelAshProperties( solnData(CI_MSCALAR,i,j,k), solnData(NEI_MSCALAR,i,j,k), &
                                             ye_f, ye_a, yi_f, yi_a, qbar_f, qbar_a )
              qbar = (1.0-phi_fa)*qbar_f + (phi_fa-phi_aq)*qbar_a + solnData(DQQN_MSCALAR,i,j,k)
              lsum(9) = lsum(9) + ((solnData(YE_MSCALAR,i,j,k)-0.5)*(m_p+m_e-m_n)*c_l*c_l*N_A-qbar*cgsMeVperGram)*dm

              ! burned masses
              lsum(10) = lsum(10) + phi_fa*dm
              lsum(11) = lsum(11) + phi_aq*dm
              lsum(12) = lsum(12) + phi_qn*dm
              lsum(13) = lsum(13) + solnData(FLAM_MSCALAR,i,j,k)*dm

              ! estimate Ni56 by assuming "other" continuent has Ye of equal parts Ni58 and Fe54
              if (phi_qn > 1e-6) then
                 yen = (solnData(YE_MSCALAR,i,j,k)-(1.0-phi_qn)*ye_a)/phi_qn
                 xni56 = min(1.0, max( 0.0, (yen-0.48212)/(0.5-0.48212) ))
                 lsum(14) = lsum(14) + phi_qn*xni56*dm
              endif

              ! density bins
              dens = solnData(DENS_VAR,i,j,k)
              if (dens>1.e9) lsum(15) = lsum(15) + dm
              if (dens>3.e8) lsum(16) = lsum(16) + dm
              if (dens>1.e8) lsum(17) = lsum(17) + dm
              if (dens>3.e7) lsum(18) = lsum(18) + dm
              if (dens>2.e7) lsum(19) = lsum(19) + dm
              if (dens>1.5e7) lsum(20) = lsum(20) + dm
              if (dens>1.e7) lsum(21) = lsum(21) + dm
              if (dens>7.e6) lsum(22) = lsum(22) + dm
              if (dens>3.e6) lsum(23) = lsum(23) + dm
              if (dens>1.e6) lsum(24) = lsum(24) + dm

              ! max density
              if (dens > lmaxdens) lmaxdens = dens

              ! min density of material starting to burn in flame (for DDT)
              if ( solnData(FLAM_MSCALAR,i,j,k) > 0.001 .and. &
                   solnData(FLAM_MSCALAR,i,j,k) < 0.01 .and. dens < lminflamdens) &
                   lminflamdens = dens

           enddo
        enddo
     enddo
     call Grid_releaseBlkPtr(blockList(lb), solnData)

  enddo

  call MPI_AllReduce(bn_neutLossThisProcStep, totalNeutLossThisStep, 1, FLASH_REAL, MPI_SUM, io_globalComm, error)
  bn_neutLossThisProcStep = 0.0
  bn_neutLoss = bn_neutLoss + totalNeutLossThisStep
  call IO_setScalar("bn_neutLoss", bn_neutLoss)
  
  ! Now the MASTER_PE sums the local contributions from all of
  ! the processors and writes the total to a file.
  
  call MPI_Reduce (lsum, gsum, nGlobalSum, FLASH_REAL, MPI_SUM, & 
       &                MASTER_PE, io_globalComm, error)
  call MPI_Reduce (lmaxdens, gmaxdens, 1, FLASH_REAL, MPI_MAX, & 
       &                MASTER_PE, io_globalComm, error)
  call MPI_Reduce (lminflamdens, gminflamdens, 1, FLASH_REAL, MPI_MIN, & 
       &                MASTER_PE, io_globalComm, error)
  
  

  if (io_globalMe == MASTER_PE) then
     
     ! create the file from scratch if it is a not a restart simulation, 
     ! otherwise append to the end of the file
     
     !No mater what, we are opening the file. Check to see if already there
     ioStat = 0
     open(funit, file=trim(io_statsFileName), position='APPEND', status='OLD', iostat=ioStat)
     if(ioStat .NE. 0) then
        !print *, 'FILE FOUND'
        open(funit, file=trim(io_statsFileName), position='APPEND')
     endif
     
     if (isFirst .EQ. 1 .AND. (.NOT. io_restart .or. ioStat .NE. 0)) then
        
        write (funit, 10)               &
             '#time                     ', &
             'mass                      ', &
             'x-momentum                ', &
             'y-momentum                ', & 
             'z-momentum                ', &
             'E_internal+kinetic        ', &
             'E_kinetic (from vel)      ', &
             'E_internal                ', &
             'E_grav                    ', &
             'E_nuc                     ', &
             'E_neutloss                ', &
             'mass burned               ', &
             'mass burned to NSQE       ', &
             'mass burned to NSE        ', &
             'mass burned by flame      ', &
             'estimated Ni56 mass       ', &
             'mass with dens > 1e9      ', &
             'mass with dens > 3e8      ', &
             'mass with dens > 1e8      ', &
             'mass with dens > 3e7      ', &
             'mass with dens > 2e7      ', &
             'mass with dens > 1.5e7    ', &
             'mass with dens > 1e7      ', &
             'mass with dens > 7e6      ', &
             'mass with dens > 3e6      ', &
             'mass with dens > 1e6      ', &
             'maximum density           ', &
             'minimum flame density     '


        
        
10         format (2x,50(a25, :, 1X))

     else if(isFirst .EQ. 1) then
        write (funit, 11) 
11      format('# simulation restarted')
     endif
     
     
     ! Write the global sums to the file.
     write (funit, 12) simtime, gsum(1:9), bn_neutLoss, gsum(10:nGlobalSum), gmaxdens, gminflamdens
12   format (1x, 50(es25.18, :, 1x))
 
     close (funit)          ! Close the file.
     
  endif
  
  call MPI_Barrier (io_globalComm, error)
  
  !=============================================================================
  
  return
end subroutine IO_writeIntegralQuantities



