module probdata_module

  use network, only: nspec, network_species_index
  use eos_type_module, only: eos_t
  use eos_module, only: eos_input_rt, eos
  use bl_constants_module, only: ZERO, THIRD, HALF, ONE, TWO, THREE, M_PI, FOUR
  use fundamental_constants_module, only: Gconst, M_solar, AU
  use initial_model_module, only: initial_model

  ! Probin file
  character (len=:), allocatable :: probin

  ! For determining if we are the I/O processor.
  integer :: ioproc
  
  ! Initial binary orbit characteristics
  ! Note that the envelope mass is included within the total mass of the star
  double precision :: mass_P, mass_S
  double precision :: central_density_P, central_density_S
  double precision :: orbital_speed_P, orbital_speed_S
  double precision :: stellar_temp
  double precision :: primary_envelope_mass, secondary_envelope_mass
  double precision :: primary_envelope_comp(nspec), secondary_envelope_comp(nspec)

  ! Ambient medium
  double precision :: ambient_density, ambient_temp, ambient_comp(nspec)

  ! Smallest allowed mass fraction
  double precision :: smallx

  ! Smallest allowed velocity on the grid
  double precision :: smallu
  
  ! Controls interpolation from 1D model to 3D model
  integer :: nsub
  logical :: interp_temp

  ! Whether or not to give stars an initial orbital velocity
  ! consistent with their Keplerian orbit speed.
  logical :: no_orbital_kick

  ! Whether or not to give the stars an initial velocity
  ! consistent with the free-fall speed.
  logical :: collision

  ! For a collision, number of (secondary) WD radii to 
  ! separate the WDs by.
  double precision :: collision_separation

  ! Damping
  logical          :: damping
  double precision :: damping_alpha

  ! Binary properties
  double precision :: a_P_initial, a_S_initial, a  
  double precision :: center_P_initial(3), center_S_initial(3)

  integer :: star_axis
  integer :: initial_motion_dir

  ! Location of the physical center of the problem, as a fraction of domain size
  double precision :: center_fracx, center_fracy, center_fracz

  ! Bulk system motion
  double precision :: bulk_velx, bulk_vely, bulk_velz

  ! Whether we're doing an initialization or a restart
  integer :: init

  ! Are we doing a single star simulation?
  logical :: single_star

  ! Should we override the domain boundary conditions with
  ! ambient material?
  logical :: fill_ambient_bc
  
  ! 1D initial models

  type (initial_model) :: model_P, model_S

  double precision :: initial_model_dx
  integer          :: initial_model_npts
  double precision :: initial_model_mass_tol
  double precision :: initial_model_hse_tol

  ! Composition properties of initial models
  double precision :: max_he_wd_mass
  double precision :: max_hybrid_co_wd_mass, hybrid_co_wd_he_shell_mass
  double precision :: max_co_wd_mass, co_wd_he_shell_mass

  ! Tagging criteria
  double precision :: maxTaggingRadius

  ! Stores the center of mass location of the stars throughout the run
  double precision :: com_P(3), com_S(3)
  double precision :: vel_P(3), vel_S(3)

  ! Stores the effective Roche radii
  double precision :: roche_rad_P, roche_rad_S

  ! Relaxation parameters
  logical          :: do_relax
  integer          :: relax_type ! 1 = SCF

  ! Input parameters for SCF relaxation
  double precision :: scf_d_A, scf_d_B, scf_d_C
  double precision :: scf_relax_tol
  
  ! Internal data for SCF relaxation
  double precision :: scf_h_max_P, scf_h_max_S
  double precision :: scf_enthalpy_min
  double precision :: scf_rpos(3,3)         ! Relative position of points A, B, and C
  double precision :: scf_d_vector(3,3)     ! Positions of points relative to system center
  double precision :: scf_c(0:1,0:1,0:1,3)  ! Interpolation coefficients for points
  integer          :: scf_rloc(3,3)         ! Indices of zones nearby to these points

  ! Sponge parameters
  double precision :: sponge_dist      ! Fraction of the domain size interior to which we should not sponge
  double precision :: sponge_width     ! Interval over which we should smooth out the sponging (in units of the domain size)
  double precision :: sponge_timescale ! Typical timescale to use in determining sponging strength

  ! Distance (in kpc) used for calculation of the gravitational wave amplitude
  ! (this wil be calculated along all three coordinate axes).
  double precision :: gw_dist

contains

  ! This routine calls all of the other subroutines at the beginning
  ! of a simulation to fill in the basic problem data.

  subroutine initialize(name, namlen, init_in)

    use bl_error_module, only: bl_error

    implicit none

    integer :: namlen, i, init_in
    integer :: name(namlen)

    init = init_in

    ! Build "probin" filename -- the name of the file containing the fortin namelist.
    allocate(character(len=namlen) :: probin)
    do i = 1, namlen
       probin(i:i) = char(name(i))
    enddo

    ! Read in the namelist to set problem parameters.

    call read_namelist

    ! Determine if we are the I/O processor, and save it to the ioproc variable.

    call get_ioproc

    ! Establish binary parameters and create initial models.

    call binary_setup

    ! Set small_pres and small_ener.

    call set_small

  end subroutine



  ! This routine reads in the namelist

  subroutine read_namelist

    use meth_params_module

    implicit none

    integer :: untin

    namelist /fortin/ &
         mass_P, mass_S, &
         central_density_P, central_density_S, &
         nsub, &
         no_orbital_kick, &
         collision, &
         collision_separation, &
         interp_temp, &
         damping, damping_alpha, &
         do_relax, &
         relax_type, &
         scf_d_A, scf_d_B, scf_d_C, &
         scf_relax_tol, &
         ambient_density, &
         stellar_temp, ambient_temp, &
         max_he_wd_mass, &
         max_hybrid_co_wd_mass, hybrid_co_wd_he_shell_mass, &
         max_co_wd_mass, &
         star_axis, initial_motion_dir, &
         maxTaggingRadius, &
         bulk_velx, bulk_vely, bulk_velz, &
         smallx, smallu, &
         center_fracx, center_fracy, center_fracz, &
         initial_model_dx, &
         initial_model_npts, &
         initial_model_mass_tol, &
         initial_model_hse_tol, &
         sponge_dist, sponge_width, sponge_timescale, &
         gw_dist, &
         fill_ambient_bc

    nsub = 1

    central_density_P = -ONE
    central_density_S = -ONE

    mass_P = 1.0d0
    mass_S = 1.0d0

    stellar_temp = 1.0d7
    ambient_temp = 1.0d7

    ! We follow the prescription of Dan et al. 2012 for determining
    ! the composition of the WDs. In this approach, below a certain 
    ! mass there are pure He WDs; just above that are hybrid WDs
    ! with pure CO cores and a He mantle; above that are pure CO WDs
    ! with slightly more oxygen than carbon; and above that are 
    ! ONeMg WDs. All masses are in solar masses.

    max_he_wd_mass = 0.45d0
    max_hybrid_co_wd_mass = 0.6d0
    hybrid_co_wd_he_shell_mass = 0.1d0
    max_co_wd_mass = 1.05d0
    co_wd_he_shell_mass = 0.0d0

    smallx = 1.d-12
    smallu = 1.d-12

    ambient_density = 1.d-4

    no_orbital_kick = .false.

    collision = .false.
    collision_separation = 4.0

    interp_temp = .false.
    damping  = .false.
    star_axis = 1
    initial_motion_dir = 2

    maxTaggingRadius = 0.75d0

    do_relax = .false.
    relax_type = 1

    scf_d_A = 1.0d9
    scf_d_B = 1.0d9
    scf_d_C = 1.8d9
    scf_relax_tol = 1.d-3

    bulk_velx = ZERO
    bulk_vely = ZERO
    bulk_velz = ZERO

    center_fracx = HALF
    center_fracy = HALF
    center_fracz = HALF

    ! For the grid spacing for our model, we'll use 
    ! 6.25 km. No simulation we do is likely to have a resolution
    ! higher than that inside the stars (it represents
    ! three jumps by a factor of four compared to our 
    ! normal coarse grid resolution). By using 2048
    ! grid points, the size of the 1D domain will be 2.56e9 cm,
    ! which is larger than any reasonable mass white dwarf.

    initial_model_dx = 6.25d5
    initial_model_npts = 4096

    ! initial_model_mass_tol is tolerance used for getting the total WD mass 
    ! equal to the desired mass. It can be reasonably small, since there
    ! will always be a central density value that can give the desired
    ! WD mass on the grid we use.

    initial_model_mass_tol = 1.d-6

    ! hse_tol is the tolerance used when iterating over a zone to force
    ! it into HSE by adjusting the current density (and possibly
    ! temperature).  hse_tol should be very small (~ 1.e-10).

    initial_model_hse_tol = 1.d-10

    sponge_dist      = 0.75d0
    sponge_width     = 0.1d0
    sponge_timescale = 0.01d0

    gw_dist = 10.0 ! kpc

    fill_ambient_bc = .false.
    
    ! Read namelist to override the defaults.

    untin = 9 
    open(untin,file=probin,form='formatted',status='old')
    read(untin,fortin)
    close(unit=untin)

    ! Convert masses from solar masses to grams.

    mass_P = mass_P * M_solar
    mass_S = mass_S * M_solar

    max_he_wd_mass = max_he_wd_mass * M_solar
    max_hybrid_co_wd_mass = max_hybrid_co_wd_mass * M_solar
    max_co_wd_mass = max_co_wd_mass * M_solar

    hybrid_co_wd_he_shell_mass = hybrid_co_wd_he_shell_mass * M_solar
    co_wd_he_shell_mass = co_wd_he_shell_mass * M_solar

    if (mass_S < ZERO .and. central_density_S < ZERO) single_star = .true.

    ! Make sure that the primary mass is really the larger mass

    call ensure_primary_mass_larger

  end subroutine read_namelist



  ! Determine if we are the I/O processor, and save it to the ioproc variable

  subroutine get_ioproc

    implicit none

    call bl_pd_is_ioproc(ioproc)

  end subroutine get_ioproc



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
    use prob_params_module, only: center, problo, probhi
    use meth_params_module, only: rot_axis

    implicit none

    type (eos_t) :: ambient_state
    double precision :: v_ff
    double precision :: mass_temp

    ! Set up the center variable. We want it to be at 
    ! problo + center_frac * domain_width in each direction.
    ! center_frac is 1/2 by default, so the problem
    ! would be set up exactly in the center of the domain.

    center(1) = problo(1) + center_fracx * (probhi(1) - problo(1))
    center(2) = problo(2) + center_fracy * (probhi(2) - problo(2))
    center(3) = problo(3) + center_fracz * (probhi(3) - problo(3))

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

    if (ioproc == 1 .and. init == 1) then
       
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

       if (ioproc == 1 .and. init == 1) then
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

           initial_motion_dir = star_axis

           collision_separation = collision_separation * model_S % radius

           call freefall_velocity(mass_P + mass_S, collision_separation, v_ff)

           vel_P(star_axis) =  (mass_P / (mass_S + mass_P)) * v_ff
           vel_S(star_axis) = -(mass_S / (mass_S + mass_P)) * v_ff

           a_P_initial = (mass_P / (mass_S + mass_P)) * collision_separation
           a_S_initial = (mass_S / (mass_S + mass_P)) * collision_separation

           a = a_P_initial + a_S_initial

       else

           ! Now, if we're not doing a relaxation step, then we want to put the WDs 
           ! where they ought to be by Kepler's third law. But if we are, then we 
           ! need a better first guess. The central location for each WD should be
           ! equal to their inner distance plus their radius.

           if (do_relax .and. relax_type .eq. 1) then

              a_P_initial = scf_d_A + model_P % radius
              a_S_initial = scf_d_B + model_S % radius

              a = a_P_initial + a_S_initial

           else

              call kepler_third_law(model_P % radius, model_P % mass, model_S % radius, model_S % mass, &
                                    rot_period, a_P_initial, a_S_initial, a, &
                                    orbital_speed_P, orbital_speed_S)

           endif

           if (ioproc == 1 .and. init == 1) then
              write (*,1003) a, a / AU
              write (*,1004) a_P_initial, a_P_initial / AU
              write (*,1005) a_S_initial, a_S_initial / AU
              1003 format ("Generated binary orbit of distance ", ES8.2, " cm = ", ES8.2, " AU.")
              1004 format ("The primary orbits the center of mass at distance ", ES8.2, " cm = ", ES8.2, " AU.")
              1005 format ("The secondary orbits the center of mass at distance ", ES8.2, " cm = ", ES8.2, " AU.")
           endif      

           ! Direction of initial motion -- it is the position axis
           ! that is not star axis and that is also not the rotation axis.

           if (rot_axis .eq. 3) then
              if (star_axis .eq. 1) initial_motion_dir = 2
              if (star_axis .eq. 2) initial_motion_dir = 1
           else if (rot_axis .eq. 2) then
              if (star_axis .eq. 1) initial_motion_dir = 3
              if (star_axis .eq. 3) initial_motion_dir = 1
           else if (rot_axis .eq. 1) then
              if (star_axis .eq. 2) initial_motion_dir = 3
              if (star_axis .eq. 3) initial_motion_dir = 2
           else
              call bl_error("Error: probdata module: invalid choice for rot_axis.")
           endif

           if ( (do_rotation .ne. 1) .and. (.not. no_orbital_kick) ) then
              vel_P(initial_motion_dir) = - orbital_speed_P
              vel_S(initial_motion_dir) =   orbital_speed_S
           endif

       endif

       ! Star center positions -- we'll put them in the midplane on the
       ! axis specified by star_axis, with the center of mass at the center of the domain.

       center_P_initial(star_axis) = center_P_initial(star_axis) - a_P_initial
       center_S_initial(star_axis) = center_S_initial(star_axis) + a_S_initial

       ! Compute initial Roche radii

       call get_roche_radii(mass_S / mass_P, roche_rad_S, roche_rad_P, a)

       
    endif

    ! Reset the terminal color to its previous state.

    if (ioproc == 1 .and. init == 1) then
       print *, ''//achar(27)//'[0m'
    endif

    com_P = center_P_initial
    com_S = center_S_initial

    ! Safety check: make sure the stars are actually inside the computational domain.

    if (center_P_initial(1) - model_P % radius .lt. problo(1) .or. &
        center_P_initial(1) + model_P % radius .gt. probhi(1) .or. &
        center_P_initial(2) - model_P % radius .lt. problo(2) .or. &
        center_P_initial(2) + model_P % radius .gt. probhi(2) .or. &
        center_P_initial(3) - model_P % radius .lt. problo(3) .or. &
        center_P_initial(3) + model_P % radius .gt. probhi(3)) then
       call bl_error("Primary does not fit inside the domain.")
    endif

    if (center_S_initial(1) - model_S % radius .lt. problo(1) .or. &
        center_S_initial(1) + model_S % radius .gt. probhi(1) .or. &
        center_S_initial(2) - model_S % radius .lt. problo(2) .or. &
        center_S_initial(2) + model_S % radius .gt. probhi(2) .or. &
        center_S_initial(3) - model_S % radius .lt. problo(3) .or. &
        center_S_initial(3) + model_S % radius .gt. probhi(3)) then
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

    model % core_comp = smallx
    model % envelope_comp = smallx

    model % envelope_mass = ZERO

    ! Here we follow the prescription of Dan et al. 2012.

    if (model % mass > ZERO .and. model % mass < max_he_wd_mass) then

       model % core_comp(iHe4) = ONE
       
       model % envelope_comp = model % core_comp

    else if (model % mass >= max_he_wd_mass .and. model % mass < max_hybrid_co_wd_mass) then

       model % core_comp(iC12) = HALF
       model % core_comp(iO16) = HALF

       model % envelope_mass = hybrid_co_wd_he_shell_mass

       if (model % envelope_mass > ZERO) then
          model % envelope_comp(iHe4) = ONE
       else
          model % envelope_comp = model % core_comp
       endif

    else if (model % mass >= max_hybrid_co_wd_mass .and. model % mass < max_co_wd_mass) then
         
       model % core_comp(iC12) = 0.4d0
       model % core_comp(iO16) = 0.6d0

       model % envelope_mass = co_wd_he_shell_mass

       if (model % envelope_mass > ZERO) then
          model % envelope_comp(iHe4) = ONE
       else
          model % envelope_comp = model % core_comp
       endif

    else if (model % mass > max_co_wd_mass) then

       model % core_comp(iO16)  = 0.60d0
       model % core_comp(iNe20) = 0.35d0
       model % core_comp(iMg24) = 0.05d0

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

  subroutine kepler_third_law(radius_1, mass_1, radius_2, mass_2, period, a_1, a_2, a, v_1, v_2)

    use prob_params_module, only: problo, probhi

    implicit none

    double precision, intent(in   ) :: mass_1, mass_2, period, radius_1, radius_2
    double precision, intent(inout) :: a, a_1, a_2
    double precision, intent(inout) :: v_1, v_2

    double precision :: length

    ! Evaluate Kepler's third law

    a = (Gconst*(mass_1 + mass_2)*period**2/(FOUR*M_PI**2))**THIRD

    a_2 = a/(ONE + mass_2/mass_1)
    a_1 = (mass_2/mass_1)*a_2

    ! Calculate the orbital speeds. This is simply equal to
    ! orbital circumference over the orbital period.

    v_1 = TWO * M_PI * a_1 / period
    v_2 = TWO * M_PI * a_2 / period

    ! Make sure the domain is big enough to hold stars in an orbit this size

    length = a + radius_1 + radius_2

    if (length > (probhi(1)-problo(1)) .or. &
        length > (probhi(2)-problo(2)) .or. &
        length > (probhi(3)-problo(3))) then
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

    implicit none

    double precision :: temp_mass

    ! We want the primary WD to be more massive. If what we're calling
    ! the primary is less massive, switch the stars.

    if ( mass_P < mass_S ) then

      if (ioproc == 1) then
        print *, "Primary mass is less than secondary mass; switching the stars so that the primary is more massive."
      endif

      temp_mass = mass_P
      mass_P = mass_S
      mass_S = temp_mass

    endif

  end subroutine ensure_primary_mass_larger

end module probdata_module
