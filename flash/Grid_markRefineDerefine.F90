!!
!! Dean M. Townsley 2012
!!
!! Replacement of default refinement marking from
!!  source/Grid/Grid_markRefineDerefine.F90
!! see that file for API spec and notes.
!!
!! The main difference between the default refinement marking and this
!! implementation is the usage of several different maximum refinement levels
!! which depend on the local physical variables (density, energy
!! generation rate).  This requires that the parent-child
!! reconciliation (i.e. a child should not derefine if its parent is marked
!! for refinement) be performed after all refinement checks have been
!! performed.  This is somewhat complex because parents may be on different
!! processors than their children, so refinement information must be
!! communicated.
!!
!! Also in this version, specialized for turbulent initial conditions,
!! adds a central region of elevated refinement over a specified area
!! with a specified refinement level
!!
!! The procedure is thus:
!!   1. Use standard routines to mark based on gradients
!!   2. Impose additional limits on refinement level
!!      2.1. impose specified max refinement in fluff
!!         where fluff is determined by a density threshold
!!         regions with maxdens < threshold are derefined
!!      2.2. impose specified max refinement outside energy-generating regions
!!         regions are determined by an enuc threshold
!!         note that the margin is on the other side of the
!!         threshold compared to the Fluff threshold above
!!         i.e. regions with enuc > threshold can fully refine
!!      2.3. boost refinement if necessary in central region
!!      2.4. standard refinement min and max
!!   3. Do parent-child consistency
!!
!!***

subroutine Grid_markRefineDerefine()

  use Grid_data, ONLY : gr_refine_cutoff, gr_derefine_cutoff,&
                        gr_refine_filter,&
                        gr_numRefineVars,gr_refine_var, &
                        gr_maxRefine, gr_meshComm, gr_meshMe
  use tree, ONLY : newchild, refine, derefine, stay, lrefine, &
                   lrefine_min, &
                   nodetype, parent, child, nchild
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_getListOfBlocks,&
       Grid_releaseBlkPtr, Grid_getBlkIndexLimits, Grid_fillGuardCells, &
       Grid_getBlkBoundBox
  use Simulation_data, ONLY : sim_refFluffDensThresh, sim_refFluffMargin, &
       sim_refFluffLevel, sim_refNogenEnucThresh, sim_refNogenMargin, &
       sim_refNogenLevel, sim_refNogenFldtThresh, &
       sim_refCentRegionDist, sim_refCentRegionLevel

  implicit none

#include "Flash_mpi.h"
#include "constants.h"
#include "Flash.h"

  real :: ref_cut,deref_cut,ref_filter
  integer       :: l,i,iref

  logical :: gcMask(NUNK_VARS)
  
  integer, dimension(2,MDIM) :: blkLimits
  integer, dimension(2,MDIM) :: blkLimitsGC
  integer ,dimension(MAXBLOCKS) :: blkList
  integer :: blkCount, iblk, j, k, error
  real, dimension(2,MDIM) :: boundBox

  real, dimension(:,:,:,:), pointer :: solnData

  real :: maxdens, maxenuc, maxfldt
  real :: minxd, minyd, minzd

  ! for reconciling parents and children
  logical, dimension(MAXBLOCKS) :: refine_parent
  integer :: nsend, nrecv, ierr
  ! need enough receive slots for all parents (up to 1 each block)
  integer, dimension(MAXBLOCKS) :: recvreq
  integer, dimension(MPI_STATUS_SIZE,MAXBLOCKS) :: recvstat
  ! need enough send slots for all children (up to 2**NDIM each block)
  integer, dimension((K1D+1)*(K2D+1)*(K3D+1)*MAXBLOCKS) :: sendreq
  integer, dimension(MPI_STATUS_SIZE,(K1D+1)*(K2D+1)*(K3D+1)*MAXBLOCKS) :: sendstat

  newchild(:) = .FALSE.
  refine(:)   = .FALSE.
  derefine(:) = .FALSE.
  stay(:)     = .FALSE.

  !----------------------------------------------------
  !  1. Use standard routines to mark based on gradients
  !----------------------------------------------------

  gcMask = .FALSE.
  do i = 1,gr_numRefineVars
     iref = gr_refine_var(i)
     if (iref > 0) gcMask(iref) = .TRUE.
  end do

  ! refinement max/min levels are based on these variables
  gcMask(DENS_VAR) = .true.
  gcMask(ENUC_VAR) = .true.
  gcMask(FLDT_VAR) = .true.

  call Grid_fillGuardCells(CENTER, ALLDIR, eosMode=MODE_DENS_EI, &
       doEos=.true., maskSize=NUNK_VARS, mask=gcMask, makeMaskConsistent=.true.)


  ! default is to derefine unless a criteria requires resolution
  ! note only mark active blocks, otherwise nothing works
  call Grid_getListOfBlocks(ACTIVE_BLKS, blkList,blkCount)
  do iblk = 1, blkCount
     derefine(blkList(iblk)) = .true.
  end do
  do l = 1,gr_numRefineVars
     iref = gr_refine_var(l)
     ref_cut = gr_refine_cutoff(l)
     deref_cut = gr_derefine_cutoff(l)
     ref_filter = gr_refine_filter(l)
     call gr_markRefineDerefine(iref,ref_cut,deref_cut,ref_filter)
  end do

#ifdef FLASH_GRID_PARAMESH2
  ! Make sure lrefine_min and lrefine_max are obeyed - KW
  if (gr_numRefineVars .LE. 0) then
     call gr_markRefineDerefine(-1, 0.0, 0.0, 0.0)
  end if
#endif

  !----------------------------------------------------
  !  2. Impose additional limits on refinement level
  !----------------------------------------------------
  call Grid_getListOfBlocks(ACTIVE_BLKS, blkList,blkCount)
  do iblk = 1, blkCount
     ! get info about block
     call Grid_getBlkindexLimits(blkList(iblk), blkLimits, blkLimitsGC)
     call Grid_getBlkPtr(blkList(iblk),solnData,CENTER)
     call Grid_getBlkBoundBox( blkList(iblk), boundBox)

     ! find maximum density and enuc in block + guardcells
     maxdens = 0.0
     maxenuc = 0.0
     maxfldt = 0.0
     do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
        do j = blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)
           do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)
              maxdens = max(solnData(DENS_VAR,i,j,k),maxdens)
              maxenuc = max(abs(solnData(ENUC_VAR,i,j,k)),maxenuc)
              maxfldt = max(abs(solnData(FLDT_VAR,i,j,k)),maxfldt)
           enddo
        enddo
     enddo

     call Grid_releaseBlkPtr(blkList(iblk),solnData)

     !----------------------------------------------------
     !  2.1. impose specified max refinement in fluff
     !     where fluff is determined by a density threshold
     !  regions with maxdens < threshold are derefined
     !----------------------------------------------------

     ! below threshold force derefinement to sim_refFluffLevel
     if ( maxdens < sim_refFluffDensThresh ) then
        if ( lrefine(blkList(iblk)) > sim_refFluffLevel ) then
           refine(blkList(iblk)) = .false.
           derefine(blkList(iblk)) = .true.
        else if ( lrefine(blkList(iblk)) == sim_refFluffLevel ) then
           refine(blkList(iblk)) = .false.
        endif
     ! above threshold but within margin prevent refinement past FluffLevel
     ! but don't force derefinement
     else if ( maxdens < sim_refFluffDensThresh*(1.0+sim_refFluffMargin) &
               .and. lrefine(blkList(iblk)) >= sim_refFluffLevel ) then
        refine(blkList(iblk)) = .false.
     endif
     ! else allow refinement (lrefine_max is enforced below)

     !----------------------------------------------------
     !  2.2. impose specified max refinement outside energy-generating regions
     !     regions are determined by an enuc threshold
     !  note that the margin is on the other side of the
     !  threshold compared to the Fluff threshold above
     !   i.e. regions with enuc > threshold can fully refine
     !----------------------------------------------------

     ! below threshold and margin force derefinement to NegenLevel
     if ( maxenuc < sim_refNogenEnucThresh*(1.0-sim_refNogenMargin) .and. &
          maxfldt < sim_refNogenFldtThresh*(1.0-sim_refNogenMargin) ) then
        if ( lrefine(blkList(iblk)) > sim_refNogenLevel ) then
           refine(blkList(iblk)) = .false.
           derefine(blkList(iblk)) = .true.
        else if ( lrefine(blkList(iblk)) == sim_refNogenLevel ) then
           refine(blkList(iblk)) = .false.
        endif
     ! below threshold but within margin prevent refinement past NogenLevel
     ! but don't force derefinement
     else if ( maxenuc < sim_refNogenEnucThresh .and. maxfldt < sim_refNogenFldtThresh &
               .and. lrefine(blkList(iblk)) >= sim_refNogenLevel ) then
        refine(blkList(iblk)) = .false.
     endif
     ! else allow refinement (lrefine_max is enforced below)

     !----------------------------------------------------
     !  2.3. allow a region of defined (minimum) refinement level
     !       nominally for initial turbulence field
     !       All blocks with cells closer than sim_refCentRegionDist to the axis
     !       will set to be no lower than sim_refCentRegionLevel
     !----------------------------------------------------

     ! nearest distances to axis planes
     minxd = min( abs(boundBox(1,IAXIS)), abs(boundBox(2,IAXIS)) )
     minyd = min( abs(boundBox(1,JAXIS)), abs(boundBox(2,JAXIS)) )
     minzd = min( abs(boundBox(1,KAXIS)), abs(boundBox(2,KAXIS)) )
     if ( minxd <= sim_refCentRegionDist .and. &
          minyd <= sim_refCentRegionDist .and. &
          minzd <= sim_refCentRegionDist ) then
        if ( lrefine(blkList(iblk)) < sim_refCentRegionLevel ) then
           refine(blkList(iblk)) = .true.
           derefine(blkList(iblk)) = .false.
        else if ( lrefine(blkList(iblk)) == sim_refCentRegionLevel ) then
           derefine(blkList(iblk)) = .false.
        endif
     endif

     !----------------------------------------------------
     !  2.4. standard refinement min and max
     !----------------------------------------------------

     if (lrefine(blkList(iblk)) < lrefine_min) then
        refine(blkList(iblk)) = .true.
        derefine(blkList(iblk)) = .false.
     else if (lrefine(blkList(iblk)) == lrefine_min) then
        derefine(blkList(iblk)) = .false.
     endif
     if (lrefine(blkList(iblk)) > gr_maxRefine) then
        refine(blkList(iblk)) = .false.
        derefine(blkList(iblk)) = .true.
     else if (lrefine(blkList(iblk)) == gr_maxRefine) then
        refine(blkList(iblk)) = .false.
     endif
     
  enddo

  !---------------------------------------------------------------
  ! 3. Do parent-child consistency
  !
  ! for children that are marked derefine, check if parent is 
  ! marked refine and if so unmark derefine
  !---------------------------------------------------------------
  call Grid_getListOfBlocks(ALL_BLKS, blkList,blkCount)
  refine_parent(:) = .false.
  ! open (async) message recieve if parent is off-procssor
  ! otherwise fill directly
  !    message id is child block number on local processor
  nrecv = 0
  do iblk = 1, blkCount
     i = blkList(iblk)
     if (parent(1,i) > 0) then
        if (parent(2,i)/=gr_meshMe) then
           nrecv = nrecv+1
           call MPI_IRecv(refine_parent(i), 1, MPI_LOGICAL, parent(2,i), &
                                    i, gr_meshComm, recvreq(nrecv), ierr)
        else
           refine_parent(i) = refine(parent(1,i))
        endif
     endif
  end do
  ! parents send refine flag to each off-processor child
  nsend = 0
  do iblk = 1, blkCount
     i = blkList(iblk)
     do j = 1,nchild
        if (child(1,j,i) > 0) then
           if (child(2,j,i) /= gr_meshMe) then
              nsend = nsend + 1
              call MPI_ISend(refine(i), 1, MPI_LOGICAL, child(2,j,i), &
                                  child(1,j,i), gr_meshComm, sendreq(nsend), ierr)
           endif
        endif
     enddo
  enddo
  ! wait to recieve all parent info
  if (nrecv > 0) then
     call MPI_Waitall(nrecv,recvreq,recvstat,ierr)
  endif
  ! now reconcile, deferring to parent marked for refine
  do iblk = 1, blkCount
     i = blkList(iblk)
     if (nodetype(i) == LEAF .and. derefine(i) .and. refine_parent(i) ) then
        derefine(i) = .false.
     endif
  enddo
  ! wait until last to ask to have refine() buffer back under our control since
  ! we are not modifying it
  if (nsend > 0) then
     call MPI_Waitall(nsend,sendreq,sendstat,ierr)
  endif


  return
end subroutine Grid_markRefineDerefine
