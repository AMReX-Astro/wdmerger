!!****if* source/Simulation/SimulationMain/Flame3StageNoise_cleanup/Simulation_init
!!
!! NAME
!!
!!  Simulation_initAfterDomain
!!
!! SYNOPSIS
!!
!!  Simulation_initAfterDomain(integer,intent(IN)  :: mype)
!!
!! DESCRIPTION
!!   Perform initialization which comes after domain has been filled star.
!!   This is to put in a convection field.
!!
!!   Loop over all the blocks, check if each is within the convection field.
!!   If it is, read in the relevant piece of the field and map it in.
!!   Don't want to read the whole convection field in because it would induce
!!   a big overhead of about 500 MB.
!!
!! ARGUMENTS
!!
!!   mype : my Processor Number
!!
!!
!!
!!***

subroutine Simulation_initAfterDomain()

  use Simulation_data
  use sim_local_interface, ONLY : sim_interpolate1dWd
  use Logfile_interface, only : Logfile_stampMessage
  use Grid_interface, only : Grid_getListOfBlocks, Grid_getBlkIndexLimits, &
                             Grid_getCellCoords, Grid_getBlkPtr, &
                             Grid_releaseBlkPtr
  use IO_data, only : io_useCollectiveHDF5

  implicit none

#include "constants.h"
#include "Flash.h"
#include "mpif.h"

  !----------------------------------
  integer :: numLeafBlks, maxLeafBlks, blk, blockID, ierr
  integer :: leafBlkIDs(MAXBLOCKS)
  integer :: blkLimits(2,MDIM), blkLimitsGC(2,MDIM)
  real, pointer, dimension(:,:,:,:)                    :: blkData

  integer :: isize, jsize, ksize, i, j, k
  real, allocatable, dimension(:) :: iCoords, jCoords, kCoords
  real, allocatable, dimension(:,:,:) :: velx, vely, velz

  integer :: ext_i, ext_j, ext_k
  real :: radius, dens, temp, xc12initial, xne22initial, f
  logical :: use_collective

  if ( .not. sim_read_turbfield ) return

  use_collective = io_useCollectiveHDF5

  call Grid_getListOfBlocks(LEAF, leafBlkIDs, numLeafBlks)


  call Logfile_stampMessage("[Simulation_initAfterDomain] opening turbulence field file")
  ! call into c routines to open file and such
  call sim_turb_field_setup( sim_turbfield_bbox(IAXIS,LOW),  &
                             sim_turbfield_bbox(IAXIS,HIGH), &
                             sim_turbfield_bbox(JAXIS,LOW),  &
                             sim_turbfield_bbox(JAXIS,HIGH), &
                             sim_turbfield_bbox(KAXIS,LOW),  &
                             sim_turbfield_bbox(KAXIS,HIGH), &
                             sim_smooth_level, sim_turbfield_filename, &
                             use_collective)
  ! use_collective is set to false if TURB_PARALLELIO is not defined
  if ( use_collective ) then
     call MPI_Allreduce(numLeafBlks, maxLeafBlks, 1, MPI_Integer, &
                        MPI_MAX, MPI_COMM_WORLD, ierr)
  else
     maxLeafBlks = numLeafBlks
  endif

  do blk = 1, numLeafBlks
     blockID = leafBlkIDs(blk)

     ! get essential info about this block - index limits and cell coordinates
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

     isize = blkLimits(HIGH,IAXIS)-blkLimits(LOW,IAXIS)+1
     jsize = blkLimits(HIGH,JAXIS)-blkLimits(LOW,JAXIS)+1
     ksize = blkLimits(HIGH,KAXIS)-blkLimits(LOW,KAXIS)+1

     allocate(iCoords(isize))
     allocate(jCoords(jsize))
     allocate(kCoords(ksize))
     call Grid_getCellCoords(IAXIS,blockID,CENTER,.false.,iCoords,isize)
     call Grid_getCellCoords(JAXIS,blockID,CENTER,.false.,jCoords,jsize)
     call Grid_getCellCoords(KAXIS,blockID,CENTER,.false.,kCoords,ksize)

  !   print *, "starting block with corner x,y,z", iCoords(1), jCoords(1), kCoords(1)

     allocate(velx(isize,jsize,ksize))
     allocate(vely(isize,jsize,ksize))
     allocate(velz(isize,jsize,ksize))

     call sim_turb_field_get_vel(isize, jsize, ksize, velx, vely, velz, iCoords, jCoords, kCoords)

     ! now write to the grid

     call Grid_getBlkPtr(blockID, blkData, CENTER)

     do k = 1, ksize
        ext_k = blkLimits(LOW,KAXIS) + k - 1
        do j = 1, jsize
           ext_j = blkLimits(LOW,JAXIS) + j - 1
           do i = 1, isize
              ext_i = blkLimits(LOW,IAXIS) + i - 1

              ! velocity field magnitude is parameterized by temperature of the background star
              radius = sqrt( iCoords(i)**2 + jCoords(j)**2 + kCoords(k)**2 )
!acc              call sim_interpolate1dWd(radius, iCoords(2)-iCoords(1), dens, temp, xc12initial, xne22initial)

              if (temp < 1.0001*sim_vrms_T0) then
                 f=0.0
              else
                 f = sim_vrms_reduced * ( (temp-0.9999*sim_vrms_T0)/sim_vrms_T0*2 )**sim_vrms_alpha
              endif

              blkData(VELX_VAR,ext_i,ext_j,ext_k) = f*velx(i,j,k)
              blkData(VELY_VAR,ext_i,ext_j,ext_k) = f*vely(i,j,k)
              blkData(VELZ_VAR,ext_i,ext_j,ext_k) = f*velz(i,j,k)

              blkData(ENER_VAR,ext_i,ext_j,ext_k) = blkData(EINT_VAR,ext_i,ext_j,ext_k) &
                                 + 0.5*( blkData(VELX_VAR,ext_i,ext_j,ext_k)**2 &
                                        +blkData(VELY_VAR,ext_i,ext_j,ext_k)**2 &
                                        +blkData(VELZ_VAR,ext_i,ext_j,ext_k)**2 )
              
           enddo
        enddo
     enddo

     call Grid_releaseBlkPtr(blockID, blkData)

     deallocate(velx)
     deallocate(vely)
     deallocate(velz)

     deallocate(iCoords)
     deallocate(jCoords)
     deallocate(kCoords)

  enddo ! end loop over leaf blocks

  do blk = numLeafBlks + 1, maxLeafBlks

     call sim_turb_field_null_read()

  enddo

  call sim_turb_field_teardown
  call Logfile_stampMessage("[Simulation_initAfterDomain] turbulence field file closed")

end subroutine Simulation_initAfterDomain
