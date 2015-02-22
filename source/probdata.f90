module probdata_module

  use network, only: nspec, network_species_index
  use eos_type_module, only: eos_t
  use eos_module, only: eos_input_rt, eos

  ! Probin file
  character (len=:), allocatable :: probin

  ! Determine if we are the I/O processor
  integer :: ioproc
  
  ! Model info
  character (len=80) :: model_P_name
  character (len=80) :: model_S_name
  
  double precision, allocatable :: model_P_r(:), model_P_state(:,:)
  double precision, allocatable :: model_S_r(:), model_S_state(:,:)

  integer :: npts_model

  ! Initial binary orbit characteristics

  double precision :: mass_P_initial, mass_S_initial
  double precision :: central_density_P, central_density_S
  double precision :: radius_P_initial, radius_S_initial
  double precision :: orbital_speed_P, orbital_speed_S
  double precision :: stellar_temp, stellar_comp(nspec)

  ! Smallest allowed mass fraction
  double precision :: smallx

  ! Smallest allowed velocity on the grid
  double precision :: smallu
  
  ! Controls interpolation from 1D model to 3D model
  integer :: nsub
  logical :: interp_temp

  ! Whether or not to give stars an initial orbital velocity
  ! consistent with their Keplerian orbit speed.
  logical :: orbital_kick

  ! Damping
  logical          :: damping
  double precision :: damping_alpha

  ! Relaxation
  logical          :: do_relax
  double precision :: relax_tau

  ! Binary properties
  double precision :: a_P_initial, a_S_initial, a
  
  double precision, dimension(3) :: center_P_initial, center_S_initial

  integer :: star_axis
  integer :: initial_motion_dir

  ! Density of ambient medium
  double precision :: ambient_density

  ! Bulk system motion
  double precision :: bulk_velx, bulk_vely, bulk_velz

  ! Whether we're doing an initialization or a restart
  integer :: init

  ! Are we doing a single star simulation?
  logical :: single_star

  ! Tagging criteria
  double precision :: maxTaggingRadius

  ! Stores the center of mass location of the stars throughout the run
  double precision :: com_P(3), com_S(3)
  double precision :: vel_P(3), vel_S(3)
  double precision :: mass_P,   mass_S

  ! Stores the effective Roche radii
  double precision :: roche_rad_P, roche_rad_S

contains

  ! This routine calls all of the other subroutines at the beginning
  ! of a simulation to fill in the basic problem data.

  subroutine initialize(name, namlen, init_in)

    use bl_constants_module, only: ZERO
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

    ! Set small_pres and small_ener.

    call set_small

    ! Determine if we are the I/O processor, and save it to the ioproc variable.

    call get_ioproc

    ! Establish binary parameters and create initial models.

    call binary_setup

  end subroutine



  ! This routine reads in the namelist

  subroutine read_namelist

    use bl_constants_module, only: ZERO, HALF, ONE
    use meth_params_module

    implicit none

    integer :: untin, i

    integer :: iC12, iO16
    double precision :: stellar_C12, stellar_O16
    double precision :: mass_p, mass_s

    namelist /fortin/ &
         mass_p, mass_s, &
         central_density_P, central_density_S, &
         nsub, &
         orbital_kick, &
         interp_temp, &
         damping, damping_alpha, &
         do_relax, relax_tau, &
         ambient_density, &
         stellar_temp, stellar_C12, stellar_O16, &
         star_axis, &
         maxTaggingRadius, &
         bulk_velx, bulk_vely, bulk_velz

    maxTaggingRadius = 0.75d0

    nsub = 1

    central_density_P = -ONE
    central_density_S = -ONE

    mass_P_initial = 1.0d0
    mass_S_initial = 1.0d0

    stellar_temp = 1.0d7
    stellar_C12  = HALF
    stellar_O16  = HALF

    smallx = 1.d-12
    smallu = 1.d-12

    ambient_density = 1.d-4

    orbital_kick = .false.
    interp_temp = .false.
    damping  = .false.
    do_relax = .false.
    star_axis = 1

    bulk_velx = ZERO
    bulk_vely = ZERO
    bulk_velz = ZERO

    ! Read namelist to override the defaults
    untin = 9 
    open(untin,file=probin,form='formatted',status='old')
    read(untin,fortin)
    close(unit=untin)

    mass_P_initial = mass_p
    mass_S_initial = mass_s

    if (mass_S_initial < ZERO) single_star = .true.

    ! Fill in stellar composition using C12 and O16

    iC12 = network_species_index("carbon-12")
    iO16 = network_species_index("oxygen-16")

    stellar_comp(1:nspec) = smallx
    stellar_comp(iC12)    = stellar_C12 - (nspec - 2) * smallx
    stellar_comp(iO16)    = stellar_O16 - (nspec - 2) * smallx

    ! Make sure that the primary mass is really the larger mass

    call ensure_primary_mass_larger

  end subroutine read_namelist



  ! Determine if we are the I/O processor, and save it to the ioproc variable

  subroutine get_ioproc

    implicit none

    ! For outputting -- determine if we are the IO processor
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
    eos_state % xn  = stellar_comp
 
    call eos(eos_input_rt, eos_state, .false.)

    small_pres = eos_state % p
    small_ener = eos_state % e

  end subroutine set_small



  ! Returns the ambient state

  subroutine get_ambient(ambient_state)

    implicit none

    type (eos_t) :: ambient_state

    ! Define ambient state and call EOS to get eint and pressure

    ambient_state % rho = ambient_density
    ambient_state % T   = stellar_temp
    ambient_state % xn  = stellar_comp

    call eos(eos_input_rt, ambient_state, .false.)

  end subroutine get_ambient



  ! Set up a binary simulation

  subroutine binary_setup

    use bl_constants_module, only: ZERO, ONE
    use fundamental_constants_module, only: M_solar
    use meth_params_module, only: do_rotation, rot_period
    use initial_model_module, only: init_1d
    use prob_params_module, only: center
    use meth_params_module, only: rot_axis

    implicit none

    double precision :: dx
    type (eos_t) :: ambient_state

    double precision :: loc_P(3), loc_S(3)

    call get_ambient(ambient_state)

    center_P_initial = center
    center_S_initial = center

    radius_P_initial = ZERO
    radius_S_initial = ZERO

    npts_model = 1024

    dx = 1.0d6

    com_P = center
    com_S = center

    vel_P = ZERO
    vel_S = ZERO

    mass_P = mass_P_initial * M_solar
    mass_S = mass_S_initial * M_solar

    ! Allocate arrays to hold the stellar models

    allocate(model_P_state(npts_model,3+nspec))
    allocate(model_S_state(npts_model,3+nspec))

    allocate(model_P_r(npts_model))
    allocate(model_S_r(npts_model))

    ! Generate primary and secondary WD

    call init_1d(model_P_r, model_P_state, npts_model, dx, radius_P_initial, mass_P_initial, &
                 central_density_P, stellar_temp, stellar_comp, ambient_state)

    if (ioproc == 1 .and. init == 1) then
        print *, "Generated initial model for primary WD of mass", mass_P_initial, &
                 ", central density", central_density_P, ", and radius", radius_P_initial
    endif

    if (.not. single_star) then

       call init_1d(model_S_r, model_S_state, npts_model, dx, radius_S_initial, mass_S_initial, &
                    central_density_S, stellar_temp, stellar_comp, ambient_state)

       if (ioproc == 1 .and. init == 1) then
          print *, "Generated initial model for secondary WD of mass", mass_S_initial, &
                   ", central density", central_density_S, ", and radius", radius_S_initial
       endif

       ! Get the orbit from Kepler's third law

       call kepler_third_law(radius_P_initial, mass_P_initial, radius_S_initial, mass_S_initial, &
                             rot_period, a_P_initial, a_S_initial, a, &
                             orbital_speed_P, orbital_speed_S)

       ! Star center positions -- we'll put them in the midplane on the
       ! axis specified by star_axis, with the center of mass at the center of the domain.
       ! If we're only doing a single star, which we know because the secondary mass is
       ! less than zero, then just leave the star in the center.

       center_P_initial(star_axis) = center_P_initial(star_axis) - a_P_initial
       center_S_initial(star_axis) = center_S_initial(star_axis) + a_S_initial

       com_P = center_P_initial
       com_S = center_S_initial

       ! Direction of initial motion -- essentially it's the position axis
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

       if ( (do_rotation .ne. 1) .and. orbital_kick ) then
          vel_P(initial_motion_dir) = - orbital_speed_P
          vel_S(initial_motion_dir) =   orbital_speed_S
       endif

    endif

  end subroutine binary_setup



  ! Accepts the masses of two stars (in solar masses)
  ! and the orbital period of a system,
  ! and returns the semimajor axis of the orbit (in cm),
  ! as well as the distances a_1 and a_2 from the center of mass.

  subroutine kepler_third_law(radius_1, mass_1, radius_2, mass_2, period, a_1, a_2, a, v_1, v_2)

    use bl_constants_module, only: ONE, TWO, FOUR, THIRD, M_PI
    use fundamental_constants_module, only: Gconst, M_solar
    use prob_params_module, only: problo, probhi

    implicit none

    double precision, intent(in   ) :: mass_1, mass_2, period, radius_1, radius_2
    double precision, intent(inout) :: a, a_1, a_2
    double precision, intent(inout) :: v_1, v_2

    double precision :: length

    ! Evaluate Kepler's third law

    a = (Gconst*(mass_1 + mass_2)*M_solar*period**2/(FOUR*M_PI**2))**THIRD

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



  ! This routine checks to see if the primary mass is actually larger
  ! than the secondary mass, and switches them if not.

  subroutine ensure_primary_mass_larger

    implicit none

    double precision :: temp_mass

    ! We want the primary WD to be more massive. If what we're calling
    ! the primary is less massive, switch the stars.

    if ( mass_P_initial < mass_S_initial ) then

      if (ioproc == 1) then
        print *, "Primary mass is less than secondary mass; switching the stars."
      endif

      temp_mass = mass_P_initial
      mass_P_initial = mass_S_initial
      mass_S_initial = temp_mass

    endif

  end subroutine ensure_primary_mass_larger
  
end module probdata_module
