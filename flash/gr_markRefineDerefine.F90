!! Replacement of single-variable refinement marking routine from
!!  source/Grid/GridMain/paramesh/paramesh4/gr_markRefineDerefine.F90
!! see that file for main functional description
!!
!! This version is modified so that the parent-child reconciliation is
!! performed outside this subroutine (see notes at end).  Thus each
!! block is checked and marked if refinement is required based on the
!! input parameters (see main documentation in above subroutine for their
!! meaning).  The refinement test, based on a combination of derivatives
!! is the same as the default routine.
!!
!! Originial code: FLASH 3.0
!! Modifications: Dean M. Townsley 2009
!! Updates for flash4.0 (removed mpirank argument): Dean Townsley 2012
!!

subroutine gr_markRefineDerefine(iref,refine_cutoff,derefine_cutoff,refine_filter)

  use paramesh_dimensions
  use physicaldata, ONLY : gcell_on_cc, unk, unk1, no_permanent_guardcells
!  use workspace
  use tree
  use paramesh_interfaces, ONLY : amr_guardcell
  use Grid_data, ONLY: gr_geometry, gr_oneBlock, gr_meshMe

  use paramesh_interfaces, ONLY: amr_1blk_guardcell

  implicit none

  include 'mpif.h'
#include "Flash.h"
#include "constants.h"  

  integer, intent(IN) :: iref
  real, intent(IN) :: refine_cutoff, derefine_cutoff, refine_filter
  integer, parameter :: SQNDIM = NDIM*NDIM
  
  real delx,dely,delz
  real dely_f, delz_f
  real delu(mdim,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,kl_bnd1:ku_bnd1)
  real delua(mdim,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,kl_bnd1:ku_bnd1)

  real delu2(SQNDIM), delu3(SQNDIM), delu4(SQNDIM)

  real num,denom,error(maxblocks),error_par(maxblocks), errort
  real J_sq,xx,dd
  real error_max
  real vx,vy,vz,bv,bx,by,bz,bx1,by1,bz1,rhbi,rhi,bk,et
  real px,py,pz,rho,b2,v2,gf,vv,gr,dx,vt,press
  real times,time_exe

  integer lb,i,j,k,nref
  integer ineigh,ineigh_proc,ierr
  integer kstart,kend,jstart,jend,istart,iend
  integer nsend,nrecv
  integer reqr(maxblocks),reqs(maxblocks)
!
!       for the message buffering
!
  integer stags(maxblocks_tr)
  real    svals(maxblocks_tr)
  integer sprocs(maxblocks_tr)
  integer rtags(maxblocks_tr)
  real    rvals(maxblocks_tr)
  integer rprocs(maxblocks_tr)
  
  logical, parameter :: msgbuffer = .false.
  integer :: ii, jj, kk, ip, jp, kp
  integer :: ndim2
  integer :: statr(MPI_STATUS_SIZE,maxblocks)
  integer :: stats(MPI_STATUS_SIZE,maxblocks)

  real, pointer :: solnData(:,:,:)
  logical :: gcell_on_cc_backup(NUNK_VARS)
  integer :: idest, iopt, nlayers, icoord
  logical :: lcc, lfc, lec, lnc, l_srl_only, ldiag

!==============================================================================

  if (no_permanent_guardcells) gcell_on_cc_backup = gcell_on_cc

  if (.not. no_permanent_guardcells) then
! A non-directional guardcell fill for CENTER (and also EOS calls for
! all block cells, including guardcells, if any refinement variables
! refine_var_# require this to be current) must have been performed
! when this routine is invoked. Moreover, there must not be any
! intervening calls that would modify the solution data in unk (at
! least, for the variables to be used for refinement criteria).
!
! If that is not true, in the simplest case uncommenting the following line
! should be uncommented here:
!!$  call Grid_fillGuardCells(CENTER,ALLDIR)
!
! It would be better to have the caller of this routine do what's necessary,
! or, if the call has to be here, to call Grid_fillGuardCells with the
! optional mask argument indicating that only guardcells fore the refinement
! variable(s) are needed.

! The following is unused code - it copies only the inner cells to work.
!!$  do lb = 1,lnblocks
!!$     do k = NGUARD*K3D+1,NGUARD*K3D+NZB
!!$        do j = NGUARD*K2D+1,NGUARD*K2D+NYB
!!$           do i = NGUARD+1,NGUARD+NXB
!!$              work(i,j,k,lb,1) = unk(iref,i,j,k,lb)
!!$           end do
!!$        end do
!!$     end do
!!$  end do

! We are using more cell layers, including guardcells, from unk.

     
!     work(:,:,:,:,1)=unk(iref,:,:,:,:)

  end if
  !==============================================================================

  ndim2 = ndim*ndim


#define XCOORD(I) (gr_oneBlock(lb)%firstAxisCoords(CENTER,I))
#define YCOORD(I) (gr_oneBlock(lb)%secondAxisCoords(CENTER,I))

  do lb = 1,lnblocks


     error(lb) = 0.

     if (nodetype(lb).eq.1.or.nodetype(lb).eq.2) then

        if (no_permanent_guardcells) then
           !just get the iref var:
           gcell_on_cc = .false.
           gcell_on_cc(iref) = .true.

           idest = 1
           iopt = 1
           nlayers = NGUARD
           lcc = .true.
           lfc = .false.
           lec = .false.
           lnc = .false.
           l_srl_only = .false.
           icoord = 0
           ldiag = .true.

           call amr_1blk_guardcell(gr_meshMe,iopt,nlayers,lb,gr_meshMe, &
                lcc,lfc,lec,lnc, &
                l_srl_only,icoord,ldiag)

           solnData => unk1(iref,:,:,:,idest)

        else
           solnData => unk(iref,:,:,:,lb)
        end if


        delx = 0.5e0*float(NXB)/bsize(1,lb)
        
#if N_DIM >= 2
        dely = 0.5e0*float(NYB)/bsize(2,lb)
        dely_f = dely
#endif
        
#if N_DIM == 3
        delz = 0.5e0*float(NZB)/bsize(3,lb)
        delz_f = delz
#endif
        
        ! Compute first derivatives
        
        do k = 1+K3D,NZB+(K3D*((2*NGUARD)-1))
           do j = 1+K2D,NYB+(K2D*((2*NGUARD)-1))
              do i = 2,NXB+(2*NGUARD)-1
                 
                 if (gr_geometry == SPHERICAL) &
                      delx = 1.0/(XCOORD(i+1) - XCOORD(i-1))

                 ! d/dx
                 delu(1,i,j,k) = solnData(i+1,j,k) - solnData(i-1,j,k)
                 delu(1,i,j,k) = delu(1,i,j,k)*delx
                 
                 delua(1,i,j,k) = abs(solnData(i+1,j,k)) + &
                      abs(solnData(i-1,j,k))
                 delua(1,i,j,k) = delua(1,i,j,k)*delx
                 
#if N_DIM >= 2
                 ! d/dy
                 if ((gr_geometry == SPHERICAL) .or. (gr_geometry == POLAR)) then
                    dely_f = dely/XCOORD(i)
                 end if

                 delu(2,i,j,k) = solnData(i,j+1,k) - solnData(i,j-1,k)
                 delu(2,i,j,k) = delu(2,i,j,k)*dely_f
                 
                 delua(2,i,j,k) = abs(solnData(i,j+1,k)) + &
                      abs(solnData(i,j-1,k))
                 delua(2,i,j,k) = delua(2,i,j,k)*dely_f
#endif
                 
#if N_DIM == 3
                 ! d/dz
                 if (gr_geometry == SPHERICAL) then
                    delz_f = delz/(  XCOORD(i) &
                              &    * sin(YCOORD(j))  )
                 else if (gr_geometry == CYLINDRICAL) then
                    delz_f = delz/XCOORD(i)
                 end if
                 delu(3,i,j,k) = solnData(i,j,k+1) -  solnData(i,j,k-1)
                 delu(3,i,j,k) = delu(3,i,j,k)*delz_f
                 
                 delua(3,i,j,k) = abs(solnData(i,j,k+1)) + &
                      abs(solnData(i,j,k-1))
                 delua(3,i,j,k) = delua(3,i,j,k)*delz_f
#endif
                 
              end do
           end do
        end do
        
        ! Compute second derivatives
        
        ! Test if at a block boundary
        
        ! Two guardcells
        kstart = 2*K3D+1
        kend   = NZB+(K3D*((2*NGUARD)-2))
        jstart = 2*K2D+1
        jend   = NYB+(K2D*((2*NGUARD)-2))
        istart = 3
        iend   = NXB+(2*NGUARD)-2
        ! One guardcell
        !            kstart = 2*K3D+1+K3D
        !            kend   = NZB+(K3D*((2*NGUARD)-2))-K3D
        !            jstart = 2*K2D+1+K2D
        !            jend   = NYB+(K2D*((2*NGUARD)-2))-K2D
        !            istart = NGUARD
        !            iend   = NXB+(2*NGUARD)-3
        ! No guardcells
        !            kstart = K3D*NGUARD+1
        !            kend   = NZB+K3D*NGUARD
        !            jstart = K2D*NGUARD+1
        !            jend   = NYB+K2D*NGUARD
        !            istart = NGUARD+1
        !            iend   = NXB+NGUARD
        
        if (neigh(1,1,lb).le.-20) istart = NGUARD+1
        if (neigh(1,2,lb).le.-20) iend   = NGUARD+NXB
        
#if N_DIM >= 2
        if (neigh(1,3,lb).le.-20) jstart = NGUARD*K2D+1
        if (neigh(1,4,lb).le.-20) jend   = NGUARD*K2D+NYB
#endif
#if N_DIM == 3
        if (neigh(1,5,lb).le.-20) kstart = NGUARD*K3D+1
        if (neigh(1,6,lb).le.-20) kend   = NGUARD*K3D+NZB
#endif
        
        do k = kstart,kend
           do j = jstart,jend
              do i = istart,iend
                 
                 if (gr_geometry == SPHERICAL) &
                      delx = 1.0/(XCOORD(i+1) - XCOORD(i-1))

                 ! d/dxdx
                 delu2(1) = delu(1,i+1,j,k) - delu(1,i-1,j,k)
                 delu2(1) = delu2(1)*delx
                 
                 delu3(1) = abs(delu(1,i+1,j,k)) + abs(delu(1,i-1,j,k))
                 delu3(1) = delu3(1)*delx
                 
                 delu4(1) = delua(1,i+1,j,k) + delua(1,i-1,j,k)
                 delu4(1) = delu4(1)*delx
                 
#if N_DIM >= 2
                 if ((gr_geometry == SPHERICAL) .or. (gr_geometry == POLAR)) then
                    dely_f = dely/XCOORD(i)
                 end if

                 ! d/dydx
                 delu2(2) = delu(1,i,j+1,k) - delu(1,i,j-1,k)
                 delu2(2) = delu2(2)*dely_f
                 
                 delu3(2) = abs(delu(1,i,j+1,k)) + abs(delu(1,i,j-1,k))
                 delu3(2) = delu3(2)*dely_f
                 
                 delu4(2) = delua(1,i,j+1,k) + delua(1,i,j-1,k)
                 delu4(2) = delu4(2)*dely_f
                 
                 ! d/dxdy
                 delu2(3) = delu(2,i+1,j,k) - delu(2,i-1,j,k)
                 delu2(3) = delu2(3)*delx
                 
                 delu3(3) = abs(delu(2,i+1,j,k)) + abs(delu(2,i-1,j,k))
                 delu3(3) = delu3(3)*delx
                 
                 delu4(3) = delua(2,i+1,j,k) + delua(2,i-1,j,k)
                 delu4(3) = delu4(3)*delx
                 
                 ! d/dydy
                 delu2(4) = delu(2,i,j+1,k) - delu(2,i,j-1,k)
                 delu2(4) = delu2(4)*dely_f
                 
                 delu3(4) = abs(delu(2,i,j+1,k)) +  &
                      &                          abs(delu(2,i,j-1,k))
                 delu3(4) = delu3(4)*dely_f
                 
                 delu4(4) = delua(2,i,j+1,k) + delua(2,i,j-1,k)
                 delu4(4) = delu4(4)*dely_f
#endif
                 
#if N_DIM == 3
                 if (gr_geometry == SPHERICAL) then
                    delz_f = delz/(  XCOORD(i) &
                              &    * sin(YCOORD(j))  )
                 else if (gr_geometry == CYLINDRICAL) then
                    delz_f = delz/XCOORD(i)
                 end if

                 ! d/dzdx
                 delu2(5) = delu(1,i,j,k+1) - delu(1,i,j,k-1)
                 delu2(5) = delu2(5)*delz_f
                 
                 delu3(5) = abs(delu(1,i,j,k+1)) + abs(delu(1,i,j,k-1))
                 delu3(5) = delu3(5)*delz_f
                 
                 delu4(5) = delua(1,i,j,k+1) + delua(1,i,j,k-1)
                 delu4(5) = delu4(5)*delz_f
                 
                 ! d/dzdy
                 delu2(6) = delu(2,i,j,k+1) - delu(2,i,j,k-1)
                 delu2(6) = delu2(6)*delz_f
                 
                 delu3(6) = abs(delu(2,i,j,k+1)) + abs(delu(2,i,j,k-1))
                 delu3(6) = delu3(6)*delz_f
                 
                 delu4(6) = delua(2,i,j,k+1) + delua(2,i,j,k-1)
                 delu4(6) = delu4(6)*delz_f
                 
                 ! d/dxdz
                 delu2(7) = delu(3,i+1,j,k) - delu(3,i-1,j,k)
                 delu2(7) = delu2(7)*delx
                 
                 delu3(7) = abs(delu(3,i+1,j,k)) + abs(delu(3,i-1,j,k))
                 delu3(7) = delu3(7)*delx
                 
                 delu4(7) = delua(3,i+1,j,k) + delua(3,i-1,j,k)
                 delu4(7) = delu4(7)*delx
                 
                 ! d/dydz
                 delu2(8) = delu(3,i,j+1,k) - delu(3,i,j-1,k)
                 delu2(8) = delu2(8)*dely_f
                 
                 delu3(8) = abs(delu(3,i,j+1,k)) + abs(delu(3,i,j-1,k))
                 delu3(8) = delu3(8)*dely_f
                 
                 delu4(8) = delua(3,i,j+1,k) + delua(3,i,j-1,k)
                 delu4(8) = delu4(8)*dely_f
                 
                 ! d/dzdz
                 delu2(9) = delu(3,i,j,k+1) - delu(3,i,j,k-1)
                 delu2(9) = delu2(9)*delz_f
                 
                 delu3(9) = abs(delu(3,i,j,k+1)) + abs(delu(3,i,j,k-1))
                 delu3(9) = delu3(9)*delz_f
                 
                 delu4(9) = delua(3,i,j,k+1) + delua(3,i,j,k-1)
                 delu4(9) = delu4(9)*delz_f
#endif

! compute the error
                 num = 0.
                 denom = 0.
                 
                 do kk = 1, ndim2
                    num = num + delu2(kk)**2
                    denom = denom + (delu3(kk) + &
                         (refine_filter*delu4(kk)+1.e-20))**2
                 end do
                    
                    ! mz -- compare the square of the error
                 error(lb) = max(error(lb),num/denom)
                    
              end do
           end do
        end do
           
           ! store the maximum error for the current block
        error(lb) = sqrt(error(lb))
           
     end if
        
  end do
     
  
! MARK FOR REFINEMENT OR DEREFINEMENT

! Note that max refine/derefine and reconciliation between
!  parents and children must be done separately after all
!  refinement tests have been applied

  do lb = 1,lnblocks

     if (nodetype(lb).eq.1 .or. nodetype(lb).eq.2) then

        ! default state is derefine
        ! we must indicate if derefinement is forbidden
        ! or if refinement is required
        ! this is effectively or'ed over all criteria
        ! any one can prevent derefinement or request refinement

        ! test for refinement
        if (error(lb) .gt. refine_cutoff) then
           derefine(lb) = .false.
           refine(lb) = .true.
        else if ( error(lb).gt.derefine_cutoff ) then
           derefine(lb) = .false.
        endif
        
     end if
     
  end do

  !restore to the state when we came in
  if (no_permanent_guardcells) gcell_on_cc = gcell_on_cc_backup
  !=========================================================================
  return
end subroutine gr_markRefineDerefine














