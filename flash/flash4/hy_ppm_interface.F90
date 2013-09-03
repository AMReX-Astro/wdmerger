!!****ih* source/physics/Hydro/HydroMain/split/PPM/hy_ppm_interface
!!
!! NAME
!!
!!  hy_ppm_interface
!!
!! SYNOPSIS
!!
!!  use hy_ppm_interface
!!
!! DESCRIPTION
!!
!!  Interface module for internal use within the split PPM Hydro implementation.
!!
!!***

! solnData depends on the ordering on unk
!!REORDER(4): solnData, tempFlx

Module hy_ppm_interface 
  implicit none

#include "constants.h"
#include "Flash.h"
  
  interface
     subroutine hy_ppm_block ( hy_meshMe,blockID,sweepDir, dt, dtOld, &
          blkLimits,blkLimitsGC,bcs, &
          numCells,numguard, &
          primaryCoord     , &
          primaryLeftCoord , &
          primaryRghtCoord , &
          primaryDx        , &
          secondCoord      , &
          thirdCoord       , &
          radialCoord      , &
          ugrid            , &
          tempArea, tempGrav1d_o, tempGrav1d, &
          tempDtDx, tempFict, tempAreaLeft,   &
          tempFlx,           & 
          shock, solnData )

       implicit none

       !! ------------
       !! ---- ARGUMENTS
       integer, intent(IN) :: hy_meshMe, blockID
       integer, intent(IN) :: sweepDir
       real,    intent(IN) :: dt, dtOld
       integer, intent(IN) :: numCells,numguard
       integer, intent(IN),dimension(2,MDIM) :: blkLimitsGC,blkLimits,bcs
       real,    pointer :: solnData(:,:,:,:) 

#ifdef FIXEDBLOCKSIZE
       real, intent(OUT), DIMENSION(GRID_ILO_GC:GRID_IHI_GC,   &
            GRID_JLO_GC:GRID_JHI_GC,          &
            GRID_KLO_GC:GRID_KHI_GC) ::       &
            tempArea,       &
            tempGrav1d_o,   &
            tempGrav1d,     &
            tempDtDx,       &
            tempFict,       &
            tempAreaLeft

       real, intent(IN), DIMENSION(GRID_ILO_GC:GRID_IHI_GC,    &
            GRID_JLO_GC:GRID_JHI_GC,          &
            GRID_KLO_GC:GRID_KHI_GC) :: &
            shock

       real, intent(OUT), DIMENSION(NFLUXES,                   &
            GRID_ILO_GC:GRID_IHI_GC,     &
            GRID_JLO_GC:GRID_JHI_GC,     &
            GRID_KLO_GC:GRID_KHI_GC) ::  &
            tempFlx

       real,intent(IN), DIMENSION(MAXCELLS) :: primaryCoord ,  &
            primaryLeftCoord , &
            primaryRghtCoord , &
            primaryDx        , &
            secondCoord      , &
            thirdCoord       , &
            radialCoord     , &
            ugrid
#else
       real, intent(OUT), DIMENSION(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),&
            blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),          &
            blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)) ::       &
            tempArea,       &
            tempGrav1d_o,   &
            tempGrav1d,     &
            tempDtDx,       &
            tempFict,       &
            tempAreaLeft
       real, intent(IN), DIMENSION(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),&
            blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),&
            blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)) :: shock

       real, intent(OUT), DIMENSION(NFLUXES,                    &
            blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),     &
            blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),     &
            blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)) ::  &
            tempFlx  
       real, intent(IN), DIMENSION(numCells) :: primaryCoord ,  &
            primaryLeftCoord , &
            primaryRghtCoord , &
            primaryDx        , &
            secondCoord      , &
            thirdCoord       , &
            radialCoord      , &
            ugrid
#endif
     end subroutine hy_ppm_block
  end interface

  
  
  interface
     subroutine hy_ppm_updateSoln(rangeSwitch, &
          xyzswp, dt,                          &
          blkLimits,blkLimitsGC,numCells,      &
          tempArea, tempGrav1d_o, tempGrav1d,  &
          tempDtDx, tempFict,                  &
          tempFlx,  solnData )
       
       implicit none
       integer, intent(IN) :: rangeSwitch
       integer, intent(IN) :: xyzswp
       real,    intent(IN) :: dt
       integer, intent(IN) :: numCells
       integer, intent(IN),dimension(2,MDIM)::blkLimitsGC,blkLimits
#ifdef FIXEDBLOCKSIZE
       real, intent(IN), DIMENSION(GRID_ILO_GC:GRID_IHI_GC, &
            GRID_JLO_GC:GRID_JHI_GC, &
            GRID_KLO_GC:GRID_KHI_GC  ) :: &
            tempArea, tempGrav1d_o,  &
            tempGrav1d, &
            tempDtDx, tempFict
       real, intent(IN), DIMENSION(NFLUXES,GRID_ILO_GC:GRID_IHI_GC, &
            GRID_JLO_GC:GRID_JHI_GC, &
            GRID_KLO_GC:GRID_KHI_GC) :: tempFlx
#else
       real, intent(IN), DIMENSION(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
            blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
            blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)  ) :: &
            tempArea, tempGrav1d_o, &
            tempGrav1d, &
            tempDtDx, tempFict
       real, intent(IN), DIMENSION(NFLUXES,blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
            blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
            blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)) :: tempFlx
       
#endif
       real, pointer :: solnData(:,:,:,:) 
       
     end subroutine hy_ppm_updateSoln
  end interface
  
  
  interface
     subroutine hy_ppm_getTemporaryData&
          (axis, blockID, size, area, dtdx, grav, ngrav, fict, areaLeft)
       integer, intent(in) :: axis, blockID
       integer, intent(in), dimension(3) :: size
       real, intent(out), dimension(size(1),size(2),size(3)) :: &
            area,dtdx,grav,ngrav,fict,areaLeft
     end subroutine hy_ppm_getTemporaryData
  end interface

  interface
     subroutine hy_ppm_putTemporaryData&
          (axis, blockID, size, area, dtdx, grav, ngrav, fict, areaLeft)
       integer, intent(in) :: axis, blockID
       integer, intent(in), dimension(3) :: size
       real, intent(in), dimension(size(1),size(2),size(3)) :: &
            area,dtdx,grav,ngrav,fict,areaLeft
     end subroutine hy_ppm_putTemporaryData
  end interface

  interface
     subroutine hy_ppm_sweep ( blockCount, blockList, &
                               timeEndAdv, dt, dtOld, &
                               sweepDir )


       !! ----------------------
       !! ---- INPUT ARGUMENTS
       integer, INTENT(IN) :: blockCount
       real,    INTENT(IN) :: timeEndAdv, dt, dtOld
       integer, INTENT(IN) :: sweepDir
       integer, INTENT(IN), dimension(blockCount) :: blockList
     end subroutine hy_ppm_sweep
  end interface


  interface
     subroutine hy_ppm_force(xyzswp, numCells, numIntCells, j, k,  igeom, coord, &
          radialCoord, thirdCoord, u, ut, utt, fict)
       implicit none
       integer,INTENT(in) :: j, k, numCells, numIntCells, igeom, xyzswp
       real, DIMENSION(numCells),INTENT(in) :: coord, u, ut, utt
       real, DIMENSION(numCells),INTENT(in) :: radialCoord, thirdCoord
       real, DIMENSION(numCells),INTENT(out) :: fict
     end subroutine hy_ppm_force
  end interface

  interface
     subroutine hy_ppm_geom (numIntCells, numCells,j, k, xyzswp, igeom, &
          areal, arear, area, dx, dvol, xl, xr, radialCoord, thirdCoord)
       implicit none
       integer, INTENT(in)              :: j, k, xyzswp, igeom, numIntCells, numCells
       real, DIMENSION(numCells),INTENT(in)    :: dx,xl, xr
       real, DIMENSION(numCells),INTENT(in) :: radialCoord, thirdCoord
       real, DIMENSION(numCells),INTENT(out) :: areal, arear, area, dvol
     end subroutine hy_ppm_geom
  end interface

  interface
     subroutine hy_ppm_completeGeomFactors (numIntCells4, numCells, igeom, &
          dx, r, dvol, cvol, area, areaLeft)
       implicit none
       integer, INTENT(in)                   :: igeom, numIntCells4, numCells
       real, DIMENSION(numCells),INTENT(in)  :: dx
       real, INTENT(in)                      :: r !        the radial coordinate.
       real, DIMENSION(numCells),intent(out) :: area, dvol
       real, DIMENSION(numCells),intent(in)  :: areaLeft, cvol
     end subroutine hy_ppm_completeGeomFactors
  end interface

end Module hy_ppm_interface
