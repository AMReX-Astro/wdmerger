! This model contains routines to:
! Generate a 1D isothermal white dwarf in hydrostatic equilibrium
! Interpolate a 1D WD model onto a 3D Cartesian grid

module initial_model_module

use bl_types
use bl_constants_module
use bl_error_module, only: bl_error
use eos_module, only: eos_input_rt, eos
use eos_type_module, only: eos_t
use network, only: nspec
use model_parser_module, only: itemp_model, idens_model, ipres_model, ispec_model
use fundamental_constants_module, only: Gconst, M_solar
use interpolate_module, only: interpolate

contains

  subroutine init_1d(model_r, model_hse, nx, dx, hse_tol, mass_tol, &
                     radius, mass, central_density, temp_core, xn_core, ambient_state)

    implicit none

    ! Arguments

    integer,          intent(in   ) :: nx
    double precision, intent(in   ) :: dx
    double precision, intent(in   ) :: hse_tol, mass_tol
    double precision, intent(in   ) :: temp_core, xn_core(nspec)
    double precision, intent(inout) :: radius
    double precision, intent(inout) :: model_hse(nx,3+nspec)
    double precision, intent(inout) :: model_r(nx)
    double precision, intent(inout) :: mass, central_density

    type (eos_t), intent(in) :: ambient_state

    ! Local variables

    double precision :: temp_base, delta

    double precision :: xzn_hse(nx), xznl(nx), xznr(nx), M_enclosed(nx)

    integer :: i, n

    double precision :: rho_c, rho_c_old, mass_wd, mass_wd_old, drho_c

    integer, parameter :: nvar = 3 + nspec

    double precision :: dens_zone, temp_zone, pres_zone, entropy
    double precision :: dpd, dpt, dsd, dst

    double precision :: p_want, drho, dtemp, delx
    double precision :: entropy_base

    double precision :: g_zone

    integer, parameter :: MAX_ITER = 250

    integer :: iter, iter_mass

    integer :: icutoff

    logical :: converged_hse, fluff, mass_converged

    double precision, dimension(nspec) :: xn

    double precision :: smallx, smallt

    double precision :: M_tot

    double precision :: max_hse_error, dpdr, rhog

    integer :: narg

    type (eos_t) :: eos_state

    ! First, make sure we haven't overspecified with both the mass and the central density,
    ! and make sure we've specified at least one.

    if (mass > ZERO .and. central_density > ZERO) then
       call bl_error('Error: Cannot specify both mass and central density in the initial model generator.')
    else if (mass < ZERO .and. central_density < ZERO) then
       call bl_error('Error: Must specify either mass or central density in the initial model generator.')
    endif

    ! If we are specifying the mass, then we don't know what WD central density 
    ! will give the desired total mass, 
    ! so we need to do a secant iteration over central density.
    ! rho_c_old is the 'old' guess for
    ! the central density and rho_c is the current guess.  After 2
    ! loops, we can start estimating the density required to yield our
    ! desired mass.

    ! If instead we are specifying the central density, then we only need to do a 
    ! single HSE integration.

    if (mass > ZERO) then

       ! Convert the WD mass into solar masses
       M_tot = mass * M_solar

       rho_c_old = -ONE
       rho_c     = 1.d7     ! A reasonable starting guess for moderate-mass WDs

    elseif (central_density > ZERO) then

       rho_c_old = central_density
       rho_c     = central_density

    endif


    !----------------------------------------------------------------------------
    ! Create a 1-D uniform grid.
    !----------------------------------------------------------------------------

    do i = 1, nx
       xznl(i) = (dble(i) - ONE)*dx
       xznr(i) = (dble(i))*dx
       model_r(i) = HALF*(xznl(i) + xznr(i))
    enddo

    mass_converged = .false.


    do iter_mass = 1, MAX_ITER

       fluff = .false.

       ! We start at the center of the WD and integrate outward.  Initialize
       ! the central conditions.
       eos_state%T     = temp_core
       eos_state%rho   = rho_c
       eos_state%xn(:) = xn_core(:)

       ! (T, rho) -> (p, s)    
       call eos(eos_input_rt, eos_state, .false.)

       ! Make the initial guess be completely uniform
       model_hse(:,idens_model) = eos_state%rho
       model_hse(:,itemp_model) = eos_state%T
       model_hse(:,ipres_model) = eos_state%p

       do i = 1, nspec
          model_hse(:,ispec_model-1+i) = eos_state%xn(i)
       enddo

       ! Keep track of the mass enclosed below the current zone
       M_enclosed(1) = FOUR3RD*M_PI*(xznr(1)**3 - xznl(1)**3)*model_hse(1,idens_model)


       !-------------------------------------------------------------------------
       ! HSE + entropy solve
       !-------------------------------------------------------------------------
       do i = 2, nx

          delx = model_r(i) - model_r(i-1)

          ! As the initial guess for the density, use the previous zone
          dens_zone = model_hse(i-1,idens_model)

          temp_zone = temp_core
          xn(:) = xn_core(:)

          g_zone = -Gconst*M_enclosed(i-1)/(xznl(i)*xznl(i))


          !----------------------------------------------------------------------
          ! Iteration loop
          !----------------------------------------------------------------------

          ! Start off the Newton loop by assuming that the zone has not converged
          converged_hse = .FALSE.

          if (.not. fluff) then

             do iter = 1, MAX_ITER

                ! The core is isothermal, so we just need to constrain
                ! the density and pressure to agree with the EOS and HSE.

                ! We difference HSE about the interface between the current
                ! zone and the one just inside.
                p_want = model_hse(i-1,ipres_model) + &
                     delx*0.5*(dens_zone + model_hse(i-1,idens_model))*g_zone

                eos_state%T     = temp_zone
                eos_state%rho   = dens_zone
                eos_state%xn(:) = xn(:)

                ! (T, rho) -> (p, s)
                call eos(eos_input_rt, eos_state, .false.)

                entropy = eos_state%s
                pres_zone = eos_state%p

                dpd = eos_state%dpdr

                drho = (p_want - pres_zone)/(dpd - 0.5*delx*g_zone)

                dens_zone = max(0.9*dens_zone, &
                     min(dens_zone + drho, 1.1*dens_zone))

                if (abs(drho) < hse_tol*dens_zone) then
                   converged_hse = .TRUE.
                   exit
                endif

                if (dens_zone < ambient_state % rho) then

                   icutoff = i
                   dens_zone = ambient_state % rho
                   temp_zone = ambient_state % T
                   converged_hse = .TRUE.
                   fluff = .TRUE.
                   exit

                endif

             enddo

             if (.NOT. converged_hse) then

                print *, 'Error zone', i, ' did not converge in init_1d'
                print *, dens_zone, temp_zone
                print *, p_want
                print *, drho
                call bl_error('Error: HSE non-convergence')

             endif

          else
             dens_zone = ambient_state % rho
             temp_zone = ambient_state % T
          endif


          ! Call the EOS one more time for this zone and then go on
          ! to the next.
          eos_state%T     = temp_zone
          eos_state%rho   = dens_zone
          eos_state%xn(:) = xn(:)

          ! (T, rho) -> (p, s)    
          call eos(eos_input_rt, eos_state, .false.)

          pres_zone = eos_state%p


          ! Update the thermodynamics in this zone.
          model_hse(i,idens_model) = dens_zone
          model_hse(i,itemp_model) = temp_zone
          model_hse(i,ipres_model) = pres_zone

          model_hse(i,ispec_model:ispec_model-1+nspec) = xn(:)

          M_enclosed(i) = M_enclosed(i-1) + &
               FOUR3RD*M_PI*(xznr(i) - xznl(i))* &
               (xznr(i)**2 +xznl(i)*xznr(i) + xznl(i)**2)*model_hse(i,idens_model)

       enddo  ! End loop over zones


       mass_wd = FOUR3RD*M_PI*(xznr(1)**3 - xznl(1)**3)*model_hse(1,idens_model)

       do i = 2, icutoff
          mass_wd = mass_wd + &
               FOUR3RD*M_PI*(xznr(i) - xznl(i))* &
               (xznr(i)**2 +xznl(i)*xznr(i) + xznl(i)**2)*model_hse(i,idens_model)
       enddo

       if (rho_c_old < ZERO) then
          ! Not enough iterations yet -- store the old central density and
          ! mass and pick a new value.
          rho_c_old = rho_c
          mass_wd_old = mass_wd

          rho_c = HALF*rho_c_old

       else
          ! Check if we have converged. If we specified the central density, we can just stop here.
          if ( abs(mass_wd - M_tot)/M_tot < mass_tol .or. central_density > ZERO ) then
             mass_converged = .true.
             exit
          endif

          ! Do a secant iteration:
          ! M_tot = M(rho_c) + dM/drho |_rho_c x drho + ...        
          drho_c = (M_tot - mass_wd)/ &
               ( (mass_wd  - mass_wd_old)/(rho_c - rho_c_old) )

          rho_c_old = rho_c
          mass_wd_old = mass_wd

          rho_c = min(1.1d0*rho_c_old, &
                      max((rho_c + drho_c), 0.9d0*rho_c_old))

       endif     

    enddo  ! End mass constraint loop

    if (.not. mass_converged) then
       print *, 'ERROR: WD mass did not converge'
       call bl_error("ERROR: mass did not converge")
    endif

    if (central_density < ZERO) then
       central_density = model_hse(1,idens_model)
    endif

    radius = model_r(icutoff)
    mass = mass_wd / M_solar

  end subroutine init_1d

  ! Takes a one-dimensional stellar model and interpolates it to a point in
  ! 3D Cartesian grid space. Optionally does a sub-grid-scale interpolation
  ! if nsub > 1 (set in the probin file).

  subroutine interpolate_3d_from_1d(model_r, model_state, npts_model, loc, dx, state, nsub_in)

    implicit none

    integer,          intent(in) :: npts_model

    double precision, intent(in) :: model_r(npts_model)
    double precision, intent(in) :: model_state(npts_model,3+nspec)
    double precision, intent(in) :: loc(3), dx(3)

    type (eos_t) :: state

    integer, optional, intent(in) :: nsub_in

    integer :: i, j, k, n
    integer :: nsub
    double precision :: x, y, z, r

    if (present(nsub_in)) then
       nsub = nsub_in
    else
       nsub = 1
    endif
    
    state % rho = ZERO 
    state % p   = ZERO
    state % T   = ZERO
    state % xn  = ZERO

    ! We perform a sub-grid-scale interpolation, where 
    ! nsub determines the number of intervals we split the zone into.
    ! If nsub = 1, we simply interpolate using the cell-center location.
    ! As an example, if nsub = 3, then the sampled locations will be
    ! k = 0 --> z = loc(3) - dx(3) / 3   (1/6 of way from left edge of zone)
    ! k = 1 --> z = loc(3)               (halfway between left and right edge)
    ! k = 2 --> z = loc(3) + dx(3) / 3   (1/6 of way from right edge of zone)

    do k = 0, nsub-1
       z = loc(3) + dble(k + HALF * (1 - nsub)) * dx(3) / nsub

       do j = 0, nsub-1
          y = loc(2) + dble(j + HALF * (1 - nsub)) * dx(2) / nsub

          do i = 0, nsub-1
             x = loc(1) + dble(i + HALF * (1 - nsub)) * dx(1) / nsub

             r = (x**2 + y**2 + z**2)**(0.5)

             state % rho = state % rho + interpolate(r, npts_model, model_r, model_state(:,idens_model))
             state % T   = state % T   + interpolate(r, npts_model, model_r, model_state(:,itemp_model))

             do n = 1, nspec
                state % xn(n) = state % xn(n) + interpolate(r, npts_model, model_r, model_state(:,ispec_model-1+n))
             enddo

          enddo
       enddo
    enddo

    ! Now normalize by the number of intervals, and complete the thermodynamics.

    state % rho = state % rho / (nsub**3)
    state % T   = state % T   / (nsub**3)
    state % xn  = state % xn  / (nsub**3)

    call eos(eos_input_rt, state)
                    
  end subroutine interpolate_3d_from_1d
  
end module initial_model_module
