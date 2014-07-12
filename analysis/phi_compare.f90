! Calculates the error associated with various approaches
! for the gravity boundary conditions.
!
!
program fwdbcs

  use f2kcli
  use bl_space
  use bl_error_module
  use bl_constants_module
  use bl_IO_module
  use plotfile_module

  implicit none

  type(plotfile) pfComp
  integer :: unit
  integer :: i, j, ii, jj, kk
  real(kind=dp_t) :: c(0:1,0:2)
  real(kind=dp_t) :: num1, num2, den1, den2
  integer :: rr

  real(kind=dp_t) :: r_axis, xctr, yctr, zctr, v_rotate, magvel

  real(kind=dp_t) :: dx(MAX_SPACEDIM)

  real(kind=dp_t), pointer :: pComp(:,:,:,:)

  integer :: phi_comp, true_comp

  real(kind=dp_t) :: residNorm, trueNorm

  logical, allocatable :: imask(:,:,:)
  integer :: lo(MAX_SPACEDIM), hi(MAX_SPACEDIM)
  integer :: flo(MAX_SPACEDIM), fhi(MAX_SPACEDIM)

  real(kind=dp_t) :: xlo(3), x, y, z
  real(kind=dp_t) :: phi, true, dphi, L0, L2

  integer :: narg, farg
  character(len=256) :: trueName, compName

  real(kind=dp_t) :: time


  unit = unit_new()

  narg = command_argument_count()

  call get_command_argument(1, value = compName)

  call build(pfComp, trim(compName), unit)

  time = pfComp%tm

  ! Figure out the variable index for phi.

  phi_comp  = plotfile_var_index(pfComp, "phiGrav")
  true_comp = plotfile_var_index(pfComp, "adv_0")
 
  ! Get dx for the coarse level. We need this for the L2 norm.

  dx = plotfile_get_dx(pfComp, 1)

  ! Get the index bounds for the coarse level.

  flo = lwb(plotfile_get_pd_box(pfComp, 1))
  fhi = upb(plotfile_get_pd_box(pfComp, 1))

  residNorm = 0.0d0
  trueNorm = 0.0d0

  xlo(:) = -1.6

  L0 = 0.0

  do j = 1, nboxes(pfComp, 1)
        
    call fab_bind(pfComp, 1, j)

    lo = lwb(get_box(pfComp, 1, j))
    hi = upb(get_box(pfComp, 1, j))

    ! get a pointer to the current patch
    pComp => dataptr(pfComp, 1, j)

    ! loop over all of the zones in the patch.
    do kk = lbound(pComp,dim=3), ubound(pComp,dim=3)
      z = xlo(3) + (dble(kk)+0.5) * dx(3)
      do jj = lbound(pComp,dim=2), ubound(pComp,dim=2)
        y = xlo(2) + (dble(jj)+0.5) * dx(2)
        do ii = lbound(pComp,dim=1), ubound(pComp,dim=1)
          x = xlo(1) + (dble(ii)+0.5) * dx(1)

!          if (kk .eq. flo(3) .or. kk .eq. fhi(3) .or. &
!              jj .eq. flo(2) .or. jj .eq. fhi(2) .or. &
!              ii .eq. flo(1) .or. ii .eq. fhi(1)) then

          phi  = pComp(ii,jj,kk,phi_comp)
          true = pComp(ii,jj,kk,true_comp)

          dphi = phi - true

          residNorm = residNorm + ( dphi )**2.0d0
          trueNorm  = trueNorm  + ( true )**2.0d0

          if (abs(dphi/true) .gt. L0) L0 = abs(dphi/true)

!          endif  

        enddo
      enddo
    enddo

    call fab_unbind(pfComp, 1, j)
                
  end do

  residNorm = ( dx(1) * dx(2) * dx(3) * residNorm )**0.5d0
  trueNorm  = ( dx(1) * dx(2) * dx(3) * trueNorm  )**0.5d0

  L2 = residNorm / trueNorm

100 format(1x, 5(g20.10))

  !write (*, 100) pf%tm, kinetic_energy, internal_energy, -potential_energy, &
  !     kinetic_energy + internal_energy - potential_energy, max_mach

  print *, "Rel. L2 Norm = ", L2
  print *, "Inf norm = ", L0

  call destroy(pfComp)

end program fwdbcs
