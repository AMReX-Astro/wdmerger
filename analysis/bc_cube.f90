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
  integer :: i, j, k, l, m, n, ii, jj, kk
  real(kind=dp_t) :: c(0:1,0:2)
  real(kind=dp_t) :: num1, num2, den1, den2
  integer :: rr

  real(kind=dp_t) :: r_axis, xctr, yctr, zctr, v_rotate, magvel

  real(kind=dp_t) :: dx(MAX_SPACEDIM)

  real(kind=dp_t), pointer :: pComp(:,:,:,:)

  integer :: phiGrav_comp

  real(kind=dp_t) :: residNorm, trueNorm

  logical, allocatable :: imask(:,:,:)
  integer :: lo(MAX_SPACEDIM), hi(MAX_SPACEDIM)
  integer :: flo(MAX_SPACEDIM), fhi(MAX_SPACEDIM)

  real(kind=dp_t) :: xlo(3), x, y, z
  real(kind=dp_t) :: phi, dphi, L0
  real(kind=dp_t) :: Gconst = 6.67428e-8_dp_t

  integer :: narg, farg
  character(len=256) :: trueName, compName

  real(kind=dp_t) :: time


  unit = unit_new()

  narg = command_argument_count()

  call get_command_argument(1, value = compName)

  call build(pfComp, trim(compName), unit)

  time = pfComp%tm

  ! Figure out the variable index for phi.

  phiGrav_comp = plotfile_var_index(pfComp, "phiGrav")
 
  ! Get dx for the coarse level. We need this for the L2 norm.

  dx = plotfile_get_dx(pfComp, 1)

  ! Get the index bounds for the coarse level.

  flo = lwb(plotfile_get_pd_box(pfComp, 1))
  fhi = upb(plotfile_get_pd_box(pfComp, 1))

  residNorm = 0.0d0
  trueNorm = 0.0d0

  xlo(:) = -1.6

  L0 = 0.0

  do m = 1, nboxes(pfComp, 1)
        
    call fab_bind(pfComp, 1, m)

    lo = lwb(get_box(pfComp, 1, m))
    hi = upb(get_box(pfComp, 1, m))

    ! get a pointer to the current patch
    pComp => dataptr(pfComp, 1, m)

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
 
          c(0,2) = -0.5 - z
          c(1,2) =  0.5 - z
          c(0,1) = -0.5 - y
          c(1,1) =  0.5 - y
          c(0,0) = -0.5 - x
          c(1,0) =  0.5 - x

          phi = 0.0

          do i = 0, 1
             do j = 0, 1
                 do l = 0, 2

                      num1 = ( (c(i,l)**2 + c(j,mod(l+1,3))**2 + c(1,mod(l+2,3))**2)**(0.5) + c(1,mod(l+2,3)) )**3
                      num2 = ( (c(i,l)**2 + c(j,mod(l+1,3))**2 + c(0,mod(l+2,3))**2)**(0.5) - c(0,mod(l+2,3)) )
                      den1 = ( (c(i,l)**2 + c(j,mod(l+1,3))**2 + c(1,mod(l+2,3))**2)**(0.5) - c(1,mod(l+2,3)) )
                      den2 = ( (c(i,l)**2 + c(j,mod(l+1,3))**2 + c(0,mod(l+2,3))**2)**(0.5) + c(0,mod(l+2,3)) )**3

                      phi = phi + 0.5 * (-1)**(i+j) * ( c(i,l) * c(j,mod(l+1,3)) * &
                                         log( num1 * num2 / (den1 * den2) ) )

                enddo
             enddo
          enddo

          do i = 0, 1
             do j = 0, 1
                do k = 0, 1
                   do l = 0, 2
                          phi = phi + (-1)**(i+j+k+1) * (c(i,l)**2 * &
                                      atan2( c(i,l) * c(k,mod(l+2,3)), c(i,l)**2 + c(j,mod(l+1,3))**2 + &
                                             c(j,mod(l+1,3))*(c(i,l)**2 + c(j,mod(l+1,3))**2 + c(k,mod(l+2,3))**2)**(0.5) ) )
                   enddo
                enddo
             enddo
          enddo

          phi = phi * 0.5 * Gconst

!          print *, ii, jj, kk, x, y, z, dx, phi, pComp(ii,jj,kk,phiGrav_comp)

          dphi = pComp(ii,jj,kk,phiGrav_comp) - phi

          residNorm = residNorm + ( dphi )**2.0d0
          trueNorm  = trueNorm  + (  phi )**2.0d0

          if (abs(dphi/phi) .gt. L0) L0 = abs(dphi/phi)

!          endif  

        enddo
      enddo
    enddo

    call fab_unbind(pfComp, 1, m)
                
  end do

  residNorm = ( dx(1) * dx(2) * dx(3) * residNorm )**0.5d0
  trueNorm  = ( dx(1) * dx(2) * dx(3) * trueNorm  )**0.5d0

100 format(1x, 5(g20.10))

  !write (*, 100) pf%tm, kinetic_energy, internal_energy, -potential_energy, &
  !     kinetic_energy + internal_energy - potential_energy, max_mach

  print *, "Rel. L2 Norm = ", residNorm / trueNorm
  print *, "Inf norm = ", L0

  call destroy(pfComp)

end program fwdbcs
