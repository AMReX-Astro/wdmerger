! Analysis routines for the 3-d WD merger Castro problem.
! These are based on fsedov3d_sph.f90
!
!
program fwdmerger

  use f2kcli
  use bl_space
  use bl_error_module
  use bl_constants_module
  use bl_IO_module
  use plotfile_module

  implicit none

  type(plotfile) pf
  integer :: unit
  integer :: i, j, ii, jj, kk
  real(kind=dp_t) :: xx, yy, zz
  integer :: rr, r1

  real(kind=dp_t) :: r_axis, xctr, yctr, zctr, v_rotate, magvel

  real(kind=dp_t) :: dx(MAX_SPACEDIM)

  real(kind=dp_t), pointer :: p(:,:,:,:)

  real(kind=dp_t) :: kinetic_energy, potential_energy, internal_energy, mach_number, max_mach, rad_star,xyz_rad

  integer :: dens_comp, rhoe_comp, u_comp, v_comp, w_comp, phi_comp, mach_comp

  real(kind=dp_t) :: f_rotate

  logical, allocatable :: imask(:,:,:)
  integer :: lo(MAX_SPACEDIM), hi(MAX_SPACEDIM)
  integer :: flo(MAX_SPACEDIM), fhi(MAX_SPACEDIM)

  integer :: narg, farg
  character(len=256) :: fname

  real(kind=dp_t) :: time


  unit = unit_new()


  ! set the defaults
  f_rotate = 0.005

  narg = command_argument_count()

  farg = 1
  do while ( farg <= narg )
     call get_command_argument(farg, value = fname)

     select case (fname)

     case ('--f_rotate')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname,*) f_rotate

     case default
        exit

     end select
     farg = farg + 1
  end do

  call get_command_argument(farg, value = fname)

  call build(pf, trim(fname), unit)

  time = pf%tm


  ! figure out the variable indices
  
  ! density
  dens_comp = plotfile_var_index(pf, "density")
  
  ! velocity
  u_comp = plotfile_var_index(pf, "x_velocity")
  v_comp = plotfile_var_index(pf, "y_velocity")
  w_comp = plotfile_var_index(pf, "z_velocity")

  ! internal energy
  rhoe_comp = plotfile_var_index(pf, "rho_e")

  ! gravitational potential
  phi_comp = plotfile_var_index(pf, "phiGrav")
  
  !Mach Number
  mach_comp = plotfile_var_index(pf, "MachNumber")



  if ( dens_comp < 0 .or. &
       u_comp < 0 .or. v_comp < 0 .or. w_comp < 0 .or. &
       rhoe_comp < 0 .or. phi_comp < 0) then
     call bl_error("ERROR: varaible(s) not defined")
  endif


  ! get dx for the coarse level.  
  dx = plotfile_get_dx(pf, 1)


  ! get the index bounds for the finest level.  Note, lo and hi are
  ! ZERO based indicies
  flo = lwb(plotfile_get_pd_box(pf, pf%flevel))
  fhi = upb(plotfile_get_pd_box(pf, pf%flevel))


  ! the default for the center of the star will be the geometric center
  ! of the domain
  xctr = HALF*(pf%phi(1) + pf%plo(1))
  yctr = HALF*(pf%phi(2) + pf%plo(2))
  zctr = HALF*(pf%phi(3) + pf%plo(3))


  ! imask will be set to false if we've already output the data.
  ! Note, imask is defined in terms of the finest level.  As we loop
  ! over levels, we will compare to the finest level index space to
  ! determine if we've already output here
  allocate(imask(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3)))



  !----------------------------------------------------------------------------
  ! compute the angle averaged quantities
  !----------------------------------------------------------------------------

  kinetic_energy = ZERO
  internal_energy = ZERO
  potential_energy = ZERO
  mach_number = ZERO
  max_mach = ZERO
  !Radius of the star. Model radius for secondary is 7.60156250e8, Model radius for primary is 9.69531250e8
  rad_star = 7.60156250e8

  imask(:,:,:) = .true.


  ! loop over the data, starting at the finest grid, and if we haven't
  ! already stored data in that grid location (according to imask),
  ! store it.  


  ! r1 is the factor between the current level grid spacing and the
  ! FINEST level
  r1  = 1

  do i = pf%flevel, 1, -1

     ! rr is the factor between the COARSEST level grid spacing and
     ! the current level
     rr = product(pf%refrat(1:i-1,1))

     do j = 1, nboxes(pf, i)
        
        ! read in the data 1 patch at a time -- read in all the variables
        call fab_bind(pf, i, j)

        lo = lwb(get_box(pf, i, j))
        hi = upb(get_box(pf, i, j))


        ! get a pointer to the current patch
        p => dataptr(pf, i, j)

        
        ! loop over all of the zones in the patch.  Here, we convert
        ! the cell-centered indices at the current level into the
        ! corresponding RANGE on the finest level, and test if we've
        ! stored data in any of those locations.  If we haven't then
        ! we store this level's data and mark that range as filled.
        do kk = lbound(p,dim=3), ubound(p,dim=3)
           zz = (kk + HALF)*dx(3)/rr

           do jj = lbound(p,dim=2), ubound(p,dim=2)
              yy = (jj + HALF)*dx(2)/rr

              do ii = lbound(p,dim=1), ubound(p,dim=1)
                 xx = (ii + HALF)*dx(1)/rr

                 if ( any(imask(ii*r1:(ii+1)*r1-1, &
                                jj*r1:(jj+1)*r1-1, &
                                kk*r1:(kk+1)*r1-1) ) ) then


                    ! weight the zone's data by its size

                    ! note, for p(:,:,:,n), n refers to index of the
                    ! variable as found via plotfile_var_index

                    r_axis = sqrt((xx-xctr)**2 + (yy-yctr)**2)

                    v_rotate = 2.0_dp_t*M_PI*f_rotate*r_axis

                    magvel = (p(ii,jj,kk,u_comp) - v_rotate*(yy-yctr)/r_axis)**2 + &
                             (p(ii,jj,kk,v_comp) + v_rotate*(xx-xctr)/r_axis)**2 + &
                             (p(ii,jj,kk,w_comp))**2
                    magvel = sqrt(magvel)

                    kinetic_energy = kinetic_energy + &
                         0.5_dp_t*p(ii,jj,kk,dens_comp) * magvel**2 * &
                            (dx(1)/rr)*(dx(2)/rr)*(dx(3)/rr)

                    potential_energy = potential_energy + &
                         0.5_dp_t*p(ii,jj,kk,dens_comp) * p(ii,jj,kk,phi_comp) * &
                            (dx(1)/rr)*(dx(2)/rr)*(dx(3)/rr)

                    internal_energy = internal_energy + &
                         p(ii,jj,kk,rhoe_comp) * &
                         (dx(1)/rr)*(dx(2)/rr)*(dx(3)/rr)

                    
                    mach_number = p(ii,jj,kk,mach_comp)
                    xyz_rad = sqrt((xx-xctr)**2 + (yy-yctr)**2 + (zz-zctr)**2)
                    if (mach_number > max_mach) then
                       if (xyz_rad <= rad_star) then
                          max_mach = mach_number
                       endif
                    endif

                    imask(ii*r1:(ii+1)*r1-1, &
                          jj*r1:(jj+1)*r1-1, &
                          kk*r1:(kk+1)*r1-1) = .false.
                    
                 end if

              end do
           enddo
        enddo

        call fab_unbind(pf, i, j)
                
     end do

     ! adjust r1 for the next lowest level
     if ( i /= 1 ) r1 = r1*pf%refrat(i-1,1)
  end do


  ! normalize

  ! we already normalized the kinetic energy by multiplying 
  ! each zone by dV

100 format(1x, 5(g20.10))

  write (*, 100) pf%tm, kinetic_energy, internal_energy, -potential_energy, &
       kinetic_energy + internal_energy - potential_energy, max_mach


  call destroy(pf)

end program fwdmerger
