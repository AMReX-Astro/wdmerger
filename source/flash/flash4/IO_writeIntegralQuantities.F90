!!****if* source/IO/IOMain/IO_writeIntegralQuantities
!!
!!
!!  NAME
!!    IO_writeIntegralQuantities
!!
!!  SYNOPSIS
!!    call IO_writeIntegralQuantities(integer(in) :: isFirst,
!!                                    real(in)    :: simTime)
!!
!!  DESCRIPTION
!!
!!   Compute the values of integral quantities (eg. total energy)
!!   and write them to an ASCII file.  If this is the initial step,
!!   create the file and write a header to it before writing the data.
!!
!!   Presently, this supports 1, 2, and 3-d Cartesian geometry and 2-d
!!   cylindrical geometry (r,z).  More geometries can be added by
!!   modifying the volume of each zone (dvol).
!!
!!   Users should modify this routine if they want to store any
!!   quantities other than default values in the flash.dat.  Make sure
!!   to modify the nGlobalSum parameter to match the number of
!!   quantities written.  Also make sure to modify the header to match
!!   the names of quantities with those calculated in the lsum and
!!   gsum arrays.
!!  
!!  ARGUMENTS
!!    
!!   isFirst - if 1 then write header info plus data, otherwise just write data
!!   simTime - simulation time
!!
!!
!!***

!!REORDER(4):solnData

subroutine IO_writeIntegralQuantities ( isFirst, simTime)

  use Simulation_data, ONLY : sim_binaryPeriod, sim_inertial
  use IO_data, ONLY : io_restart, io_statsFileName, io_globalComm
  use Grid_interface, ONLY : Grid_getListOfBlocks, &
    Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_getSingleCellVol, &
    Grid_releaseBlkPtr, Grid_getCellCoords

   use IO_data, ONLY : io_globalMe
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
  integer :: isizeGC, jsizeGC, ksizeGC

  integer, parameter ::  nGlobalSum = 22          ! Number of globally-summed quantities
  real :: gsum(nGlobalSum) ! Global summed quantities
  real :: lsum(nGlobalSum) ! Locally summed quantities

  integer :: i, j, k
  real :: dvol, dm
  real, DIMENSION(:,:,:,:), POINTER :: solnData
  real, DIMENSION(:), ALLOCATABLE :: iCoords, jCoords, kCoords
  real :: x, y, z, vx, vy, vz

  integer :: point(MDIM)
  integer :: ioStat
  integer :: massIdx = 1
  integer :: xmomIdx = 2, ymomIdx = 3, zmomIdx = 4
  integer :: kinEngIdx = 5, gravEngIdx = 6, intEngIdx = 7, totEngIdx = 8
  integer :: angMomXIdx = 9, angMomYIdx = 10, angMomZIdx = 11
  integer :: xcomIdx = 12, ycomIdx = 13, zcomIdx = 14
  integer :: leftMassIdx = 15, rightMassIdx = 16
  integer :: xcomlIdx = 17, ycomlIdx = 18, zcomlIdx = 19
  integer :: xcomrIdx = 20, ycomrIdx = 21, zcomrIdx = 22

  real    :: moment_of_inertia(0:2, 0:2)
  real    :: loc(0:2), vel(0:2), omega(0:2), L_grid(0:2), delta_com(0:2)
  integer :: m, n, idx1, idx2

  omega(0) = 0.0e0
  omega(1) = 0.0e0
  omega(2) = 2.0 * PI / sim_binaryPeriod

  ! Sum quantities over all locally held leaf-node blocks.
  gsum = 0.
  lsum = 0.
  
  call Grid_getListOfBlocks(LEAF, blockList, count)
  
  do lb = 1, count
     !get the index limits of the block
     call Grid_getBlkIndexLimits(blockList(lb), blkLimits, blkLimitsGC)

     ! Get the coordinates

     isizeGC = blkLimitsGC(HIGH,IAXIS)
     allocate(iCoords(isizeGC))
     jsizeGC = blkLimitsGC(HIGH,JAXIS)
     allocate(jCoords(jsizeGC))
     ksizeGC = blkLimitsGC(HIGH,KAXIS)
     allocate(kCoords(ksizeGC))
     call Grid_getCellCoords(IAXIS,blockList(lb),CENTER,.true.,iCoords,isizeGC)
     call Grid_getCellCoords(JAXIS,blockList(lb),CENTER,.true.,jCoords,jsizeGC)
     call Grid_getCellCoords(KAXIS,blockList(lb),CENTER,.true.,kCoords,ksizeGC)

     ! get a pointer to the current block of data
     call Grid_getBlkPtr(blockList(lb), solnData)

     ! Sum contributions from the indicated blkLimits of cells.
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        loc(2) = kCoords(k)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           loc(1) = jCoords(j)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              loc(0) = iCoords(i) 
              point(IAXIS) = i
              point(JAXIS) = j
              point(KAXIS) = k

!! Get the cell volume and mass for a single cell
              call Grid_getSingleCellVol(blockList(lb), EXTERIOR, point, dvol)
              dm = solnData(DENS_VAR,i,j,k) * dvol

              ! mass   
              lsum(massIdx) = lsum(massIdx) + dm

              if ( loc(0) > 0.0e0 ) then
                 lsum(leftMassIdx) = lsum(leftMassIdx) + dm
              elseif ( loc(0) < 0.0e0 ) then
                 lsum(rightMassIdx) = lsum(rightMassIdx) + dm
              endif
    
              ! momentum
              vel(0) = solnData(VELX_VAR,i,j,k)
              vel(1) = solnData(VELY_VAR,i,j,k)
              vel(2) = solnData(VELZ_VAR,i,j,k)

              lsum(xmomIdx) = lsum(xmomIdx) + vel(0) * dm

              lsum(ymomIdx) = lsum(ymomIdx) + vel(1) * dm
            
              lsum(zmomIdx) = lsum(zmomIdx) + vel(2) * dm

              ! kinetic energy
              lsum(kinEngIdx) = lsum(kinEngIdx) + 0.5 * dm * ( vel(0)**2 + vel(1)**2 + vel(2)**2 )

#ifdef GPOT_VAR
              ! gravitational energy
              lsum(gravEngIdx) = lsum(gravEngIdx) + 0.5 * solnData(GPOT_VAR,i,j,k) * dm
#endif
              ! internal energy
              lsum(intEngIdx) = lsum(intEngIdx) + dm * solnData(EINT_VAR,i,j,k)


              ! center of mass

              delta_com(0) = dm * loc(0)
              delta_com(1) = dm * loc(1)
              delta_com(2) = dm * loc(2)

              do m = 0, 2

                lsum(xcomIdx+m) = lsum(xcomIdx+m) + delta_com(m)
             
                if ( loc(0) > 0.0e0 ) then
                   lsum(xcomlIdx+2*m) = lsum(xcomlIdx+2*m) + delta_com(m)
                elseif ( loc(0) < 0.0e0 ) then
                   lsum(xcomrIdx+2*m) = lsum(xcomrIdx+2*m) + delta_com(m)
                endif

              enddo

              ! Angular momentum

              do m = 0, 2

                 idx1 = MOD(m+1,3)
                 idx2 = MOD(m+2,3)

                 L_grid(m) = dm * (vel(idx2) * loc(idx1) - vel(idx1) * loc(idx2))

                 lsum(angMomXIdx+m) = lsum(angMomXIdx+m) + L_grid(m)

              enddo

              ! Now add quantities if in rotating reference frame. First, motion in the frame.

              if ( .not. sim_inertial ) then

               do m = 0, 2

                  idx1 = MOD(m+1,3)
                  idx2 = MOD(m+2,3)

                  ! Momentum source = omega x (m * r)
 
                  lsum(xmomIdx+m) = lsum(xmomIdx+m) + (omega(idx1)*delta_com(idx2) + omega(idx2)*delta_com(idx1))

                  ! Kinetic energy source = omega . L

                  lsum(kinEngIdx) = lsum(kinEngIdx) + omega(m) * L_grid(m)

               enddo

               ! Now we'll add quantities from motion of the frame itself. It helps to calculate a 
               ! moment of inertia tensor for the current cell first.

               do m = 0, 2
                  do n = 0, 2

                     if ( m <= n ) then
                        if ( m /= n ) then
                           moment_of_inertia(m,n) = -dm * loc(m) * loc(n)
                        else
                           moment_of_inertia(m,n) =  dm * ( loc(MOD(m+1,3))**2 + loc(MOD(m+2,3))**2 )
                        endif
                     else
                        moment_of_inertia(m,n) = moment_of_inertia(n,m)
                     endif

                  enddo
               enddo

               do m = 0, 2
                  do n = 0, 2

                    lsum(angMomXIdx+m) = lsum(angMomXIdx+m) + moment_of_inertia(m,n) * omega(n)
                    lsum(kinEngIdx)  = lsum(kinEngIdx)  + 0.5e0 * omega(m) * moment_of_inertia(m,n) * omega(n)

                  enddo
               enddo

             endif
           
           enddo
        enddo
     enddo

     call Grid_releaseBlkPtr(blockList(lb), solnData)

     deallocate(iCoords)
     deallocate(jCoords)
     deallocate(kCoords)

  enddo
  
  
  ! Now the MASTER_PE sums the local contributions from all of
  ! the processors and writes the total to a file.
  
  call MPI_Reduce (lsum, gsum, nGlobalSum, FLASH_REAL, MPI_SUM, & 
       &                MASTER_PE, io_globalComm, error)
  

  ! Calculate total energy from kinetic, internal and gravitational

  gsum(totEngIdx) = gsum(gravEngIdx) + gsum(intEngIdx) + gsum(kinEngIdx)

  ! Divide by masses for COMs

  gsum(xcomIdx:zcomIdx) = gsum(xcomIdx:zcomIdx) / gsum(massIdx)

  gsum(xcomlIdx:zcomlIdx:2) = gsum(xcomlIdx:zcomlIdx:2) / gsum(leftMassIdx)

  gsum(xcomrIdx:zcomrIdx:2) = gsum(xcomrIdx:zcomrIdx:2) / gsum(rightMassIdx)

  if (io_globalMe  == MASTER_PE) then
     
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
             '# TIME             ', &
             'MASS               ', &
             'XMOM               ', &
             'YMOM               ', & 
             'ZMOM               ', &
             'KIN. ENERGY        ', &
             'GRAV. ENERGY       ', &
             'INT. ENERGY        ', &
             'TOTAL ENERGY       ', &
             'ANG. MOM. X        ', &
             'ANG. MOM. Y        ', &
             'ANG. MOM. Z        ', &
             'X COM              ', &
             'Y COM              ', &
             'Z COM              ', &
             'LEFT MASS          ', &
             'RIGHT MASS         ', &
             'LEFT X COM         ', &
             'RIGHT X COM        ', &
             'LEFT Y COM         ', &
             'RIGHT Y COM        ', &
             'LEFT Z COM         ', &
             'RIGHT Z COM        '
       
        
10         format (2x,50(a13, :, 1X))

     else if(isFirst .EQ. 1) then
        write (funit, 11) 
11      format('# simulation restarted')
     endif
     
     
     write (funit, 12) simtime, gsum      ! Write the global sums to the file.
12   format (1x, 50(es13.6, :, 1x))
 
     close (funit)          ! Close the file.
     
  endif
  
  call MPI_Barrier (io_globalComm, error)
  
  !=============================================================================
  
  return
end subroutine IO_writeIntegralQuantities



