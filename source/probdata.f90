module probdata_module

  use network, only: nspec, network_species_index
  use eos_type_module, only: eos_t
  use eos_module, only: eos_input_rt, eos
  use bl_constants_module, only: ZERO, THIRD, HALF, ONE, TWO, THREE, M_PI, FOUR
  use fundamental_constants_module, only: Gconst, M_solar, AU
  use initial_model_module, only: initial_model

  ! Initial stellar properties
  ! Note that the envelope mass is included within the total mass of the star

  double precision, save :: mass_P = ONE
  double precision, save :: mass_S = ONE
  double precision, save :: central_density_P = -ONE
  double precision, save :: central_density_S = -ONE
  double precision, save :: stellar_temp = 1.0d7
  double precision, save :: primary_envelope_mass, secondary_envelope_mass
  double precision, save :: primary_envelope_comp(nspec), secondary_envelope_comp(nspec)



  ! Ambient medium

  double precision, save :: ambient_density = 1.0d-4
  double precision, save :: ambient_temp = 1.0d7
  double precision, save :: ambient_comp(nspec)



  ! Smallest allowed mass fraction

  double precision, save :: smallx = 1.0d-12

  ! Smallest allowed velocity on the grid

  double precision, save :: smallu = 1.0d-12



  ! Parameters for nterpolation from 1D model to 3D model:

  ! Number of sub-grid-scale zones to use

  integer, save :: nsub = 1

  ! Default to interpolation that preserves temperature; otherwise, use pressure

  logical, save :: interp_temp = .true.



  ! Whether or not to give stars an initial orbital velocity
  ! consistent with their Keplerian orbit speed.

  logical, save :: no_orbital_kick = .false.



  ! Collision parameters

  ! Whether or not to give the stars an initial velocity
  ! consistent with the free-fall speed.

  logical, save :: collision = .false.

  ! For a collision, number of (secondary) WD radii to 
  ! separate the WDs by.

  double precision, save :: collision_separation = 4.0d0

  ! For a collision, the impact parameter measured in
  ! units of the primary's initial radius.

  double precision, save :: collision_impact_parameter = 0.0d0



  ! Binary orbit properties

  double precision, save :: r_P_initial, r_S_initial, a_P_initial, a_S_initial, a  
  double precision, save :: v_P_r, v_S_r, v_P_phi, v_S_phi
  double precision, save :: center_P_initial(3), center_S_initial(3)
  double precision, save :: orbital_eccentricity = 0.0d0
  double precision, save :: orbital_angle = 0.0d0



  ! Axis is in orbital plane; we measure angle with respect to this axis. Normally the x axis.

  integer, save :: axis_1 = 1

  ! Perpendicular axis in the orbital plane. Normally the y axis.

  integer, save :: axis_2 = 2

  ! Perpendicular to both other axies. Normally the z axis and also the rotation axis.

  integer, save :: axis_3 = 3

  ! Location of the physical center of the problem, as a fraction of domain size

  double precision, save :: center_fracx = HALF
  double precision, save :: center_fracy = HALF
  double precision, save :: center_fracz = HALF

  ! Bulk system motion

  double precision, save :: bulk_velx = ZERO
  double precision, save :: bulk_vely = ZERO
  double precision, save :: bulk_velz = ZERO

  ! Whether we're doing an initialization or a restart

  integer, save :: init

  ! Are we doing a single star simulation?

  logical, save :: single_star

  ! Should we override the domain boundary conditions with
  ! ambient material?

  logical, save :: fill_ambient_bc = .false.



  ! 1D initial models

  type (initial_model) :: model_P, model_S

  ! For the grid spacing for our model, we'll use 
  ! 6.25 km. No simulation we do is likely to have a resolution
  ! higher than that inside the stars (it represents
  ! three jumps by a factor of four compared to our 
  ! normal coarse grid resolution). By using 4096
  ! grid points, the size of the 1D domain will be 2.56e9 cm,
  ! which is larger than any reasonable mass white dwarf.

  double precision, save :: initial_model_dx = 6.25d5
  integer, save          :: initial_model_npts = 4096

  ! initial_model_mass_tol is tolerance used for getting the total WD mass 
  ! equal to the desired mass. It can be reasonably small, since there
  ! will always be a central density value that can give the desired
  ! WD mass on the grid we use.

  double precision, save :: initial_model_mass_tol = 1.d-6

  ! hse_tol is the tolerance used when iterating over a zone to force
  ! it into HSE by adjusting the current density (and possibly
  ! temperature).  hse_tol should be very small (~ 1.e-10).

  double precision, save :: initial_model_hse_tol = 1.d-10



  ! Composition properties of initial models.
  ! We follow the prescription of Dan et al. 2012 for determining
  ! the composition of the WDs. In this approach, below a certain 
  ! mass there are pure He WDs; just above that are hybrid WDs
  ! with pure CO cores and a He mantle; above that are pure CO WDs
  ! with slightly more oxygen than carbon; and above that are 
  ! ONeMg WDs. All masses are in solar masses.

  double precision, save :: max_he_wd_mass = 0.45d0
  double precision, save :: max_hybrid_wd_mass = 0.6d0
  double precision, save :: hybrid_wd_he_shell_mass = 0.1d0
  double precision, save :: max_co_wd_mass = 1.05d0
  double precision, save :: co_wd_he_shell_mass = 0.0d0

  double precision, save :: hybrid_wd_c_frac = 0.50d0
  double precision, save :: hybrid_wd_o_frac = 0.50d0

  double precision, save :: co_wd_c_frac = 0.40d0
  double precision, save :: co_wd_o_frac = 0.60d0

  double precision, save :: onemg_wd_o_frac  = 0.60d0
  double precision, save :: onemg_wd_ne_frac = 0.35d0
  double precision, save :: onemg_wd_mg_frac = 0.05d0


  ! Tagging criteria

  double precision, save :: max_tagging_radius = 0.75d0
  double precision, save :: stellar_density_threshold = 1.0d0
  double precision, save :: temperature_tagging_threshold = 5.0d8



  ! Stores the center of mass location of the stars throughout the run

  double precision, save :: com_P(3), com_S(3)
  double precision, save :: vel_P(3), vel_S(3)

  ! Stores the effective Roche radii

  double precision, save :: roche_rad_P, roche_rad_S



  ! Relaxation parameters

  logical, save          :: do_initial_relaxation = .false.
  double precision, save :: relaxation_timescale = 0.001
  double precision, save :: relaxation_density_cutoff = 1.0d0

  ! Distance (in kpc) used for calculation of the gravitational wave amplitude
  ! (this wil be calculated along all three coordinate axes).

  double precision, save :: gw_dist = 10.0d0

contains

  ! This routine calls all of the other subroutines at the beginning
  ! of a simulation to fill in the basic problem data.

  subroutine initialize_problem(init_in)

    use bl_error_module, only: bl_error
    use prob_params_module, only: dim

    implicit none

    integer :: init_in

    ! Safety check: we can't run this problem in one dimension.
    if (dim .eq. 1) then
       call bl_error("Cannot run wdmerger problem in one dimension. Exiting.")
    endif

    init = init_in

    ! Read in the namelist to set problem parameters.

    call read_namelist

    ! Establish binary parameters and create initial models.

    call binary_setup

    ! Set small_pres and small_ener.

    call set_small

  end subroutine initialize_problem



  ! This routine reads in the namelist

  subroutine read_namelist

    use meth_params_module
    use prob_params_module, only: dim, coord_type
    use problem_io_module, only: probin

    implicit none

    integer :: untin

    namelist /fortin/ &
         mass_P, mass_S, &
         central_density_P, central_density_S, &
         nsub, &
         no_orbital_kick, &
         collision, &
         collision_separation, &
         collision_impact_parameter, &
         interp_temp, &
         do_initial_relaxation, &
         relaxation_timescale, &
         relaxation_density_cutoff, &
         ambient_density, &
         stellar_temp, ambient_temp, &
         max_he_wd_mass, &
         max_hybrid_wd_mass, hybrid_wd_he_shell_mass, &
         max_co_wd_mass, &
         co_wd_he_shell_mass, &
         hybrid_wd_c_frac, hybrid_wd_o_frac, &
         co_wd_c_frac, co_wd_o_frac, &
         onemg_wd_o_frac, onemg_wd_ne_frac, onemg_wd_mg_frac, &
         orbital_eccentricity, orbital_angle, &
         axis_1, axis_2, axis_3, &
         max_tagging_radius, stellar_density_threshold, &
         temperature_tagging_threshold, &
         bulk_velx, bulk_vely, bulk_velz, &
         smallx, smallu, &
         center_fracx, center_fracy, center_fracz, &
         initial_model_dx, &
         initial_model_npts, &
         initial_model_mass_tol, &
         initial_model_hse_tol, &
         gw_dist, &
         fill_ambient_bc

    ! Read namelist to override the module defaults.

    untin = 9 
    open(untin,file=probin,form='formatted',status='old')
    read(untin,fortin)
    close(unit=untin)

    ! Convert masses from solar masses to grams.

    mass_P = mass_P * M_solar
    mass_S = mass_S * M_solar

    max_he_wd_mass = max_he_wd_mass * M_solar
    max_hybrid_wd_mass = max_hybrid_wd_mass * M_solar
    max_co_wd_mass = max_co_wd_mass * M_solar

    hybrid_wd_he_shell_mass = hybrid_wd_he_shell_mass * M_solar
    co_wd_he_shell_mass = co_wd_he_shell_mass * M_solar

    if (mass_S < ZERO .and. central_density_S < ZERO) single_star = .true.

    ! Make sure that the primary mass is really the larger mass

    call ensure_primary_mass_larger

    ! We enforce that the orbital plane is (z, phi) for two-dimensional problems,
    ! with rotation about z = 0 along the radial axis.

    if (dim .eq. 2) then

       if (coord_type .ne. 1) then
          call bl_error("We only support cylindrical coordinates in two dimensions. Set coord_type == 1.")
       endif

       axis_1 = 2
       axis_2 = 3
       rot_axis = 1

    endif

    ! Make sure we have a sensible collision impact parameter.

    if (collision_impact_parameter > 1.0) then
       call bl_error("Impact parameter must be less than one in our specified units.")
    endif

    ! Don't do a collision in a rotating reference frame.

    if (collision .and. do_rotation .eq. 1) then
       call bl_error("The collision problem does not make sense in a rotating reference frame.")
    endif

    ! Make sure we have a sensible eccentricity.

    if (orbital_eccentricity >= 1.0) then
       call bl_error("Orbital eccentricity cannot be larger than one.")
    endif

    ! Make sure we have a sensible angle. Then convert it to radians.

    if (orbital_angle < 0.0 .or. orbital_angle > 360.0) then
       call bl_error("Orbital angle must be between 0 and 360 degrees.")
    endif

    orbital_angle = orbital_angle * M_PI / 180.0

  end subroutine read_namelist



  ! Calculate small_pres and small_ener

  subroutine set_small

    use meth_params_module, only: small_temp, small_pres, small_dens, small_ener

    implicit none

    type (eos_t) :: eos_state

    ! Given the inputs of small_dens and small_temp, figure out small_pres.

    eos_state % rho = small_dens
    eos_state % T   = small_temp
    eos_state % xn  = ambient_comp

    call eos(eos_input_rt, eos_state)

    small_pres = eos_state % p
    small_ener = eos_state % e

  end subroutine set_small



  ! Returns the ambient state

  subroutine get_ambient(ambient_state)

    implicit none

    type (eos_t) :: ambient_state

    ! Define ambient state, using a composition that is an 
    ! even mixture of the primary and secondary composition, 
    ! and then call the EOS to get internal energy and pressure.

    ambient_state % rho = ambient_density
    ambient_state % T   = ambient_temp
    ambient_state % xn  = ambient_comp

    call eos(eos_input_rt, ambient_state)

  end subroutine get_ambient



  ! Set up a binary simulation

  subroutine binary_setup

    use meth_params_module, only: do_rotation, rot_period
    use initial_model_module, only: initialize_model, establish_hse
    use prob_params_module, only: center, problo, probhi, dim
    use meth_params_module, only: rot_axis
    use rotation_frequency_module, only: get_omega
    use math_module, only: cross_product
    use binary_module, only: get_roche_radii
    use problem_io_module, only: ioproc

    implicit none

    double precision :: v_ff, collision_offset
    double precision :: mu
    double precision :: omega(3)

    omega = get_omega(ZERO)

    ! Set up the center variable. We want it to be at 
    ! problo + center_frac * domain_width in each direction.
    ! center_frac is 1/2 by default, so the problem
    ! would be set up exactly in the center of the domain.
    ! Note that we override this for 2D axisymmetric, as the
    ! radial coordinate must be centered at zero for the problem to make sense.

    if (dim .eq. 3) then

       center(1) = problo(1) + center_fracx * (probhi(1) - problo(1))
       center(2) = problo(2) + center_fracy * (probhi(2) - problo(2))
       center(3) = problo(3) + center_fracz * (probhi(3) - problo(3))

    else

       center(1) = problo(1)
       center(2) = problo(2) + center_fracz * (probhi(2) - problo(2))
       center(3) = ZERO

    endif

    ! Set some default values for these quantities;
    ! we'll update them soon.

    center_P_initial = center
    center_S_initial = center

    com_P = center
    com_S = center

    vel_P = ZERO
    vel_S = ZERO

    ! Allocate arrays to hold the stellar models.

    call initialize_model(model_P, initial_model_dx, initial_model_npts, initial_model_mass_tol, initial_model_hse_tol)
    call initialize_model(model_S, initial_model_dx, initial_model_npts, initial_model_mass_tol, initial_model_hse_tol)

    model_P % min_density = ambient_density
    model_S % min_density = ambient_density

    model_P % central_temp = stellar_temp
    model_S % central_temp = stellar_temp



    ! Fill in the model's physical details.
    ! If we're integrating to reach a desired mass, set the composition accordingly.
    ! If instead we're fixing the central density, then first we'll assume the composition is
    ! that of a solar mass WD as a initial guess, and get the corresponding mass. 
    ! Then we set the composition to match this preliminary mass, and we'll get a final mass later.

    if (mass_P > ZERO) then

       model_P % mass = mass_P

       call set_wd_composition(model_P)

    elseif (central_density_P > ZERO) then

       model_P % mass = M_solar

       call set_wd_composition(model_P)

       model_P % central_density = central_density_P

       call establish_hse(model_P)

       call set_wd_composition(model_P)

    else

       call bl_error("Must specify either a positive primary mass or a positive primary central density.")

    endif



    if (.not. single_star) then

       if (mass_S > ZERO) then

          model_S % mass = mass_S

          call set_wd_composition(model_S)

       elseif (central_density_S > ZERO) then

          model_S % mass = M_solar

          call set_wd_composition(model_S)

          model_S % central_density = central_density_S

          call establish_hse(model_S)

          call set_wd_composition(model_S)

       else

          call bl_error("If we are doing a binary calculation, we must specify either a ", &
                        "positive secondary mass or a positive secondary central density.")

       endif

       ambient_comp = (model_P % envelope_comp + model_S % envelope_comp) / 2

    else

       ambient_comp = model_P % envelope_comp

    endif



    roche_rad_P = ZERO
    roche_rad_S = ZERO

    ! Generate primary and secondary WD models.

    call establish_hse(model_P)

    if (ioproc .and. init == 1) then

       ! Set the color to bold green for printing to terminal in this section. See:
       ! http://stackoverflow.com/questions/6402700/coloured-terminal-output-from-fortran

       print *, ''//achar(27)//'[1;32m'

       write (*,1001) model_P % mass / M_solar, model_P % central_density, model_P % radius
       1001 format ("Generated initial model for primary WD of mass ", f4.2, &
                    " solar masses, central density ", ES8.2, " g cm**-3, and radius ", ES8.2, " cm.")
       print *, ""
    endif

    mass_P = model_P % mass
    roche_rad_P = model_P % radius

    if (.not. single_star) then

       call establish_hse(model_S)

       if (ioproc .and. init == 1) then
          write (*,1002) model_S % mass / M_solar, model_S % central_density, model_S % radius
          1002 format ("Generated initial model for secondary WD of mass ", f4.2, &
                       " solar masses, central density ", ES8.2, " g cm**-3, and radius ", ES8.2, " cm.")
          print *, ""
       endif

       mass_S = model_S % mass
       roche_rad_S = model_S % radius

       ! For a collision, set the stars up with a free-fall velocity at a specified number
       ! of secondary radii; for a circular orbit, use Kepler's third law.

       if (collision) then

          axis_2 = axis_1

          collision_separation = collision_separation * model_S % radius

          call freefall_velocity(mass_P + mass_S, collision_separation, v_ff)

          vel_P(axis_1) =  (mass_P / (mass_S + mass_P)) * v_ff
          vel_S(axis_1) = -(mass_S / (mass_S + mass_P)) * v_ff

          r_P_initial = -(mass_P / (mass_S + mass_P)) * collision_separation
          r_S_initial =  (mass_S / (mass_S + mass_P)) * collision_separation

          a = r_S_initial - r_P_initial

          center_P_initial(axis_1) = center_P_initial(axis_1) + r_P_initial
          center_S_initial(axis_1) = center_S_initial(axis_1) + r_S_initial

          ! We also permit a non-zero impact parameter b in the direction perpendicular
          ! to the motion of the stars. This is measured in units of the radius of the
          ! primary, so that b > 1 doesn't make any sense as the stars won't collide.
          ! Since the secondary's radius is greater than the primary's, measuring in the
          ! units of the primary's radius will guarantee contact.

          collision_offset = collision_impact_parameter * model_P % radius

          center_P_initial(axis_2) = center_P_initial(axis_2) - collision_offset
          center_S_initial(axis_2) = center_S_initial(axis_2) + collision_offset                  

       else

          call kepler_third_law(model_P % radius, model_P % mass, model_S % radius, model_S % mass, &
                                rot_period, orbital_eccentricity, orbital_angle, &
                                a, r_P_initial, r_S_initial, v_P_r, v_S_r, v_P_phi, v_S_phi)

          if (ioproc .and. init == 1) then
             write (*,1003) a, a / AU
             write (*,1004) r_P_initial, r_P_initial / AU
             write (*,1005) r_S_initial, r_S_initial / AU
1003         format ("Generated binary orbit of distance ", ES8.2, " cm = ", ES8.2, " AU.")
1004         format ("The primary orbits the center of mass at distance ", ES9.2, " cm = ", ES9.2, " AU.")
1005         format ("The secondary orbits the center of mass at distance ", ES9.2, " cm = ", ES9.2, " AU.")
          endif

          ! Star center positions -- we'll put them in the midplane, with the center of mass at the center of the domain.

          center_P_initial(axis_1) = center_P_initial(axis_1) + r_P_initial * cos(orbital_angle)
          center_P_initial(axis_2) = center_P_initial(axis_2) + r_P_initial * sin(orbital_angle)

          center_S_initial(axis_1) = center_S_initial(axis_1) + r_S_initial * cos(orbital_angle)
          center_S_initial(axis_2) = center_S_initial(axis_2) + r_S_initial * sin(orbital_angle)           

          ! Star velocities, from Kepler's third law.

          vel_P(axis_1) = v_P_r   * cos(orbital_angle) - v_P_phi * sin(orbital_angle)
          vel_P(axis_2) = v_P_phi * cos(orbital_angle) + v_P_r   * sin(orbital_angle)

          vel_S(axis_1) = v_S_r   * cos(orbital_angle) - v_S_phi * sin(orbital_angle)
          vel_S(axis_2) = v_S_phi * cos(orbital_angle) + v_S_r   * sin(orbital_angle)

          ! Subtract off the rigid rotation rate. Note that this is also OK
          ! for the inertial frame because there we're going to give everything
          ! a kick of rigid body rotation in ca_initdata.

          if ( .not. no_orbital_kick ) then
             vel_P = vel_P - cross_product(omega, center_P_initial)
             vel_S = vel_S - cross_product(omega, center_S_initial)
          endif

       endif

       ! Compute initial Roche radii

       call get_roche_radii(mass_S / mass_P, roche_rad_S, roche_rad_P, a)


    endif

    ! Reset the terminal color to its previous state.

    if (ioproc .and. init == 1) then
       print *, ''//achar(27)//'[0m'
    endif

    com_P = center_P_initial
    com_S = center_S_initial

    ! Safety check: make sure the stars are actually inside the computational domain.

    if ( ( center_P_initial(1) - model_P % radius .lt. problo(1) .and. dim .eq. 3 ) .or. &
         center_P_initial(1) + model_P % radius .gt. probhi(1) .or. &
         center_P_initial(2) - model_P % radius .lt. problo(2) .or. &
         center_P_initial(2) + model_P % radius .gt. probhi(2) .or. &
         ( center_P_initial(3) - model_P % radius .lt. problo(3) .and. dim .eq. 3 ) .or. &
         ( center_P_initial(3) + model_P % radius .gt. probhi(3) .and. dim .eq. 3 ) ) then
       call bl_error("Primary does not fit inside the domain.")
    endif

    if ( ( center_S_initial(1) - model_S % radius .lt. problo(1) .and. dim .eq. 3 ) .or. &
         center_S_initial(1) + model_S % radius .gt. probhi(1) .or. &
         center_S_initial(2) - model_S % radius .lt. problo(2) .or. &
         center_S_initial(2) + model_S % radius .gt. probhi(2) .or. &
         ( center_S_initial(3) - model_S % radius .lt. problo(3) .and. dim .eq. 3 ) .or. &
         ( center_S_initial(3) + model_S % radius .gt. probhi(3) .and. dim .eq. 3 ) ) then
       call bl_error("Secondary does not fit inside the domain.")
    endif


  end subroutine binary_setup



  ! Given a WD mass, set its core and envelope composition.

  subroutine set_wd_composition(model)

    type (initial_model), intent(inout) :: model

    integer :: iHe4, iC12, iO16, iNe20, iMg24

    iHe4 = network_species_index("helium-4")
    iC12 = network_species_index("carbon-12")
    iO16 = network_species_index("oxygen-16")
    iNe20 = network_species_index("neon-20")
    iMg24 = network_species_index("magnesium-24")

    if (iHe4 < 0) call bl_error("Must have He4 in the nuclear network.")
    if (iC12 < 0) call bl_error("Must have C12 in the nuclear network.")
    if (iO16 < 0) call bl_error("Must have O16 in the nuclear network.")
    if (iNe20 < 0) call bl_error("Must have Ne20 in the nuclear network.")
    if (iMg24 < 0) call bl_error("Must have Mg24 in the nuclear network.")

    model % core_comp = smallx
    model % envelope_comp = smallx

    model % envelope_mass = ZERO

    ! Here we follow the prescription of Dan et al. 2012.

    if (model % mass > ZERO .and. model % mass < max_he_wd_mass) then

       model % core_comp(iHe4) = ONE

       model % envelope_comp = model % core_comp

    else if (model % mass >= max_he_wd_mass .and. model % mass < max_hybrid_wd_mass) then

       model % core_comp(iC12) = hybrid_wd_c_frac
       model % core_comp(iO16) = hybrid_wd_o_frac

       model % envelope_mass = hybrid_wd_he_shell_mass

       if (model % envelope_mass > ZERO) then
          model % envelope_comp(iHe4) = ONE
       else
          model % envelope_comp = model % core_comp
       endif

    else if (model % mass >= max_hybrid_wd_mass .and. model % mass < max_co_wd_mass) then

       model % core_comp(iC12) = co_wd_c_frac
       model % core_comp(iO16) = co_wd_o_frac

       model % envelope_mass = co_wd_he_shell_mass

       if (model % envelope_mass > ZERO) then
          model % envelope_comp(iHe4) = ONE
       else
          model % envelope_comp = model % core_comp
       endif

    else if (model % mass > max_co_wd_mass) then

       model % core_comp(iO16)  = onemg_wd_o_frac
       model % core_comp(iNe20) = onemg_wd_ne_frac
       model % core_comp(iMg24) = onemg_wd_mg_frac

       model % envelope_comp = model % core_comp

    endif

    ! Normalize compositions so that they sum to one.

     model % core_comp = model % core_comp / sum(model % core_comp)
     model % envelope_comp = model % envelope_comp / sum(model % envelope_comp)

  end subroutine set_wd_composition



  ! Accepts the masses of two stars (in solar masses)
  ! and the orbital period of a system,
  ! and returns the semimajor axis of the orbit (in cm),
  ! as well as the distances a_1 and a_2 from the center of mass.

  subroutine kepler_third_law(radius_1, mass_1, radius_2, mass_2, period, eccentricity, phi, a, r_1, r_2, v_1r, v_2r, v_1p, v_2p)

    use prob_params_module, only: problo, probhi

    implicit none

    double precision, intent(in   ) :: mass_1, mass_2, period, eccentricity, phi, radius_1, radius_2
    double precision, intent(inout) :: a, r_1, r_2, v_1r, v_2r, v_1p, v_2p

    double precision :: length

    double precision :: mu, M ! Reduced mass, total mass
    double precision :: r     ! Position
    double precision :: v_r, v_phi ! Radial and azimuthal velocity

    ! Definitions of total and reduced mass

    M  = mass_1 + mass_2
    mu = mass_1 * mass_2 / M    

    ! First, solve for the orbit in the reduced one-body problem, where
    ! an object of mass mu orbits an object with mass M located at r = 0.
    ! For this we follow Carroll and Ostlie, Chapter 2, but many texts discuss this.
    ! Note that we use the convention that phi measures angle from aphelion,
    ! which is opposite to the convention they use.

    a = (Gconst * M * period**2 / (FOUR * M_PI**2))**THIRD ! C + O, Equation 2.37

    r = a * (ONE - eccentricity**2) / (ONE - eccentricity * cos(phi)) ! C + O, Equation 2.3

    ! To get the radial and azimuthal velocity, we take the appropriate derivatives of the above.
    ! v_r = dr / dt = dr / d(phi) * d(phi) / dt, with d(phi) / dt being derived from
    ! C + O, Equation 2.30 for the angular momentum, and the fact that L = mu * r**2 * d(phi) / dt.

    v_r   = -TWO * M_PI * a * eccentricity * sin(phi) / (period * (ONE - eccentricity**2)**HALF)
    v_phi =  TWO * M_PI * a * (ONE - eccentricity * cos(phi)) / (period * (ONE - eccentricity**2)**HALF)

    ! Now convert everything back to the binary frame, using C+O, Equation 2.23 and 2.24. This applies
    ! to the velocities as well as the positions because the factor in front of r_1 and r_2 is constant.

    r_1  = -(mu / mass_1) * r
    r_2  =  (mu / mass_2) * r

    v_1r = -(mu / mass_1) * v_r
    v_2r =  (mu / mass_2) * v_r

    v_1p = -(mu / mass_1) * v_phi
    v_2p =  (mu / mass_2) * v_phi

    ! Make sure the domain is big enough to hold stars in an orbit this size

    length = (r_2 - r_1) + radius_1 + radius_2

    if (length > (probhi(axis_1)-problo(axis_1))) then
        call bl_error("ERROR: The domain width is too small to include the binary orbit.")
    endif

     ! Make sure the stars are not touching.
     if (radius_1 + radius_2 > a) then
        call bl_error("ERROR: Stars are touching!")
     endif

  end subroutine kepler_third_law



  ! Given total mass of a binary system and the initial separation of
  ! two point particles, obtain the velocity at this separation 
  ! assuming the point masses fell in from infinity. This will
  ! be the velocity in the frame where the center of mass is stationary.

  subroutine freefall_velocity(mass, distance, vel)

    implicit none

    double precision, intent(in   ) :: mass, distance
    double precision, intent(inout) :: vel

    vel = (TWO * Gconst * mass / distance)**HALF

  end subroutine freefall_velocity



  ! This routine checks to see if the primary mass is actually larger
  ! than the secondary mass, and switches them if not.

  subroutine ensure_primary_mass_larger

    use problem_io_module, only: ioproc

    implicit none

    double precision :: temp_mass

    ! We want the primary WD to be more massive. If what we're calling
    ! the primary is less massive, switch the stars.

    if ( mass_P < mass_S ) then

      if (ioproc) then
        print *, "Primary mass is less than secondary mass; switching the stars so that the primary is more massive."
      endif

      temp_mass = mass_P
      mass_P = mass_S
      mass_S = temp_mass

    endif

  end subroutine ensure_primary_mass_larger



  ! Given a zone state, fill it with ambient material.

   subroutine fill_ambient(state, loc, time)

    use bl_constants_module, only: ZERO
    use meth_params_module, only: NVAR, URHO, UMX, UMZ, UTEMP, UEINT, UEDEN, UFS, do_rotation
    use network, only: nspec
    use rotation_frequency_module, only: get_omega
    use math_module, only: cross_product

    implicit none

    double precision :: state(NVAR)
    double precision :: loc(3), time

    type (eos_t) :: ambient_state
    double precision :: omega(3)

    omega = get_omega(time)

    call get_ambient(ambient_state)

    state(URHO) = ambient_state % rho
    state(UTEMP) = ambient_state % T
    state(UFS:UFS-1+nspec) = ambient_state % rho * ambient_state % xn(:)                 

    ! If we're in the inertial frame, give the material the rigid-body rotation speed.
    ! Otherwise set it to zero.

    if ( (do_rotation .ne. 1) .and. (.not. no_orbital_kick) ) then

       state(UMX:UMZ) = state(URHO) * cross_product(omega, loc)

    else

       state(UMX:UMZ) = ZERO

    endif

    state(UEINT) = ambient_state % rho * ambient_state % e
    state(UEDEN) = state(UEINT) + HALF * sum(state(UMX:UMZ)**2) / state(URHO)

  end subroutine fill_ambient



  ! If we are in a rotating reference frame, then rotate a vector
  ! by an amount corresponding to the time that has passed
  ! since the beginning of the simulation.

  function inertial_rotation(vec, time) result(vec_i)

    use rotation_frequency_module, only: get_omega
    use meth_params_module, only: do_rotation, rot_period, rot_period_dot

    implicit none

    double precision :: vec(3), time

    double precision :: vec_i(3)

    double precision :: omega(3), theta(3), rot_matrix(3,3)


    ! To get the angle, we integrate omega over the time of the
    ! simulation. Since the time rate of change is linear in the
    ! period, let's work that variable. At some time t the current
    ! period P is given by P = P_0 + Pdot * t. Then:
    !
    ! theta(t) = int( omega(t) * dt )
    !      theta = int( omega(t) dt )
    !            = (2 * pi / P_0) * int( dt / (1 + (dPdt / P_0) * t) )
    !
    ! if dPdt = 0, then theta = 2 * pi * t / P_0 = omega_0 * t, as expected.
    ! if dPdt > 0, then theta = (2 * pi / P_0) * (P_0 / dPdt) * ln| (dPdt / P_0) * t + 1 |
    ! Note that if dPdt << P_0, then we have ln(1 + x) = x, and we again
    ! recover the original expression as expected.

    if (do_rotation .eq. 1) then

       if (abs(rot_period_dot) > ZERO .and. time > ZERO) then
          theta = get_omega(ZERO) * (rot_period / rot_period_dot) * &
                  log( abs( (rot_period_dot / rot_period) * time + 1 ) )
       else
          theta = get_omega(ZERO) * time
       endif

       omega = get_omega(time)

    else

       omega = ZERO
       theta = ZERO

    endif

    ! This is the 3D rotation matrix for converting between reference frames.
    ! It is the composition of rotations along the x, y, and z axes. Therefore 
    ! it allows for the case where we are rotating about multiple axes. Normally 
    ! we use the right-hand convention for constructing the usual rotation matrix, 
    ! but this is the transpose of that rotation matrix to account for the fact 
    ! that we are rotating *back* to the inertial frame, rather than from the 
    ! inertial frame to the rotating frame.

    rot_matrix(1,1) =  cos(theta(2)) * cos(theta(3))
    rot_matrix(1,2) = -cos(theta(2)) * sin(theta(3))
    rot_matrix(1,3) =  sin(theta(2))
    rot_matrix(2,1) =  cos(theta(1)) * sin(theta(3)) + sin(theta(1)) * sin(theta(2)) * cos(theta(3))
    rot_matrix(2,2) =  cos(theta(1)) * cos(theta(3)) - sin(theta(1)) * sin(theta(2)) * sin(theta(3))
    rot_matrix(2,3) = -sin(theta(1)) * cos(theta(2))
    rot_matrix(3,1) =  sin(theta(1)) * sin(theta(3)) - cos(theta(1)) * sin(theta(2)) * cos(theta(3))
    rot_matrix(3,2) =  sin(theta(1)) * cos(theta(3)) + cos(theta(1)) * sin(theta(2)) * sin(theta(3))
    rot_matrix(3,3) =  cos(theta(1)) * cos(theta(2))

    vec_i = matmul(rot_matrix, vec)

  end function inertial_rotation



  ! Given a rotating frame velocity, get the inertial frame velocity.
  ! Note that this will simply return the original velocity if we're
  ! already in the inertial frame, since omega = 0.

  function inertial_velocity(loc, vel, time) result (vel_i)

    use rotation_frequency_module, only: get_omega
    use math_module, only: cross_product

    implicit none

    double precision :: loc(3), vel(3), time
    double precision :: omega(3)

    double precision :: vel_i(3)

    omega = get_omega(time)

    vel_i = vel + cross_product(omega, loc)

  end function inertial_velocity



  ! Return the locations of the stellar centers of mass

  subroutine get_star_data(P_com, S_com, P_vel, S_vel, P_mass, S_mass) bind(C)

    implicit none

    double precision, intent(inout) :: P_com(3), S_com(3)
    double precision, intent(inout) :: P_vel(3), S_vel(3)
    double precision, intent(inout) :: P_mass, S_mass

    P_com = com_P
    S_com = com_S

    P_vel = vel_P
    S_vel = vel_S

    P_mass = mass_P
    S_mass = mass_S

  end subroutine get_star_data



  ! Set the locations of the stellar centers of mass

  subroutine set_star_data(P_com, S_com, P_vel, S_vel, P_mass, S_mass) bind(C)

    use bl_constants_module, only: TENTH, ZERO
    use prob_params_module, only: center
    use binary_module, only: get_roche_radii

    implicit none

    double precision, intent(in) :: P_com(3), S_com(3)
    double precision, intent(in) :: P_vel(3), S_vel(3)
    double precision, intent(in) :: P_mass, S_mass

    double precision :: r

    r = ZERO

    if (mass_P > ZERO) then

       com_P = P_com
       vel_P = P_vel
       mass_P = P_mass

    endif

    if (mass_S > ZERO) then

       com_S = S_com
       vel_S = S_vel
       mass_S = S_mass

    endif

    if (mass_P > ZERO .and. mass_S > ZERO) then

       r = sum((com_P-com_S)**2)**(0.5)

       call get_roche_radii(mass_S/mass_P, roche_rad_S, roche_rad_P, r)

       ! Beyond a certain point, it doesn't make sense to track the stars separately
       ! anymore. We'll set the secondary to a fixed constant and keep it there
       ! if its Roche radius becomes smaller than 10% of the primary's. Also, for exactly 
       ! equal mass systems sometimes it is the primary that disrupts, perhaps
       ! just due to numerical noise, so do the same check for the primary.

       if (roche_rad_S < TENTH * roche_rad_P) then
          com_S = center
          vel_S = ZERO
          mass_S = ZERO
          roche_rad_S = ZERO
       else if (roche_rad_P < TENTH * roche_rad_S) then
          com_P = center
          vel_P = ZERO
          mass_S = ZERO
          roche_rad_P = ZERO
       endif

    endif

  end subroutine set_star_data



  ! Check whether we should stop the initial relaxation.
  ! If so, set do_initial_relaxation to false, which will effectively
  ! turn off the external source terms.

  subroutine check_relaxation(state, s_lo, s_hi, lo, hi, L1, is_done) bind(C)

    use meth_params_module, only: URHO, NVAR
    use castro_util_module, only: position_to_index

    implicit none

    integer :: lo(3), hi(3)
    integer :: s_lo(3), s_hi(3)
    integer :: is_done
    double precision :: L1(3)
    double precision :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)

    integer :: L1_idx(3)

    ! Convert the Lagrange point into a grid index.

    L1_idx = position_to_index(L1)

    ! Check if the Lagrange point is in this box.

    if ( lo(1) .le. L1_idx(1) .and. hi(1) .ge. L1_idx(1) .and. &
         lo(2) .le. L1_idx(2) .and. hi(2) .ge. L1_idx(2) .and. &
         lo(3) .le. L1_idx(3) .and. hi(3) .ge. L1_idx(3) ) then

       ! If so, check if the density in that zone is above the cutoff.

       if (state(L1_idx(1),L1_idx(2),L1_idx(3),URHO) > relaxation_density_cutoff) then

          is_done = 1

       endif

    endif

  end subroutine check_relaxation



  subroutine turn_off_relaxation() bind(C)

    use problem_io_module, only: ioproc

    implicit none

    do_initial_relaxation = .false.

    if (ioproc) then
       print *, "Initial relaxation phase terminated."
    endif

  end subroutine turn_off_relaxation



  subroutine get_axes(axis_1_in, axis_2_in, axis_3_in) bind(C)

    implicit none

    integer :: axis_1_in, axis_2_in, axis_3_in

    axis_1_in = axis_1
    axis_2_in = axis_2
    axis_3_in = axis_3

  end subroutine get_axes



  ! Return whether we're doing a single star simulation or not.

  subroutine get_single_star(flag) bind(C)

    implicit none

    integer :: flag

    flag = 0

    if (single_star) flag = 1

  end subroutine get_single_star

end module probdata_module
