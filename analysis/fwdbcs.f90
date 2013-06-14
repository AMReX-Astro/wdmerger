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

  type(plotfile) pfTrue, pfComp
  integer :: unit
  integer :: i, j, ii, jj, kk
  real(kind=dp_t) :: xx, yy, zz
  integer :: rr, r1

  real(kind=dp_t) :: r_axis, xctr, yctr, zctr, v_rotate, magvel

  real(kind=dp_t) :: dx(MAX_SPACEDIM)

  real(kind=dp_t), pointer :: pTrue(:,:,:,:), pComp(:,:,:,:)

  real(kind=dp_t) :: kinetic_energy, potential_energy, internal_energy, mach_number, max_mach, rad_star,xyz_rad

  integer :: phiGrav_comp

  real(kind=dp_t) :: residNorm, trueNorm

  logical, allocatable :: imask(:,:,:)
  integer :: lo(MAX_SPACEDIM), hi(MAX_SPACEDIM)
  integer :: flo(MAX_SPACEDIM), fhi(MAX_SPACEDIM)

  integer :: narg, farg
  character(len=256) :: trueName, compName

  real(kind=dp_t) :: time


  unit = unit_new()

  narg = command_argument_count()

  call get_command_argument(1, value = trueName)
  call get_command_argument(2, value = compName)

  call build(pfTrue, trim(trueName), unit)
  call build(pfComp, trim(compName), unit)

  time = pfTrue%tm

  if ( pfComp%tm /= time ) then
    print *, "Times between plotfiles do not match. Exiting."
    stop
  endif


  ! Figure out the variable index for phi.

  phiGrav_comp = plotfile_var_index(pfTrue, "phiGrav")
  

  ! Get dx for the coarse level. We need this for the L2 norm.

  dx = plotfile_get_dx(pfTrue, 1)

  ! Get the index bounds for the coarse level.

  flo = lwb(plotfile_get_pd_box(pfTrue, 1))
  fhi = upb(plotfile_get_pd_box(pfTrue, 1))

  residNorm = 0.0d0
  trueNorm = 0.0d0

  do j = 1, nboxes(pfTrue, 1)
        
    call fab_bind(pfTrue, 1, j)
    call fab_bind(pfComp, 1, j)

    lo = lwb(get_box(pfTrue, 1, j))
    hi = upb(get_box(pfTrue, 1, j))


    ! get a pointer to the current patch
    pTrue => dataptr(pfTrue, 1, j)
    pComp => dataptr(pfComp, 1, j)

    ! loop over all of the zones in the patch.
    do kk = lbound(pTrue,dim=3), ubound(pTrue,dim=3)

      do jj = lbound(pTrue,dim=2), ubound(pTrue,dim=2)

        do ii = lbound(pTrue,dim=1), ubound(pTrue,dim=1)

          residNorm = residNorm + ( pTrue(ii,jj,kk,phiGrav_comp) - pComp(ii,jj,kk,phiGrav_comp) )**2.0d0
          trueNorm  = trueNorm  + ( pTrue(ii,jj,kk,phiGrav_comp) )**2.0d0
  
        enddo
      enddo
    enddo

    call fab_unbind(pfTrue, 1, j)
    call fab_unbind(pfComp, 1, j)
                
  end do

  residNorm = ( residNorm )**0.5d0
  trueNorm  = ( trueNorm  )**0.5d0

100 format(1x, 5(g20.10))

  !write (*, 100) pf%tm, kinetic_energy, internal_energy, -potential_energy, &
  !     kinetic_energy + internal_energy - potential_energy, max_mach

  print *, "Norm = ", residNorm / trueNorm

  call destroy(pfTrue)
  call destroy(pfComp)

end program fwdbcs
