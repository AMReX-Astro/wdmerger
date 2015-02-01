module probdata_module

  use network, only: nspec, network_species_index
  use eos_type_module, only: eos_t
  use eos_module, only: eos_input_rt, eos

  ! Probin file
  character (len=:), allocatable :: probin

  ! Refinement criteria
  double precision ::    denerr,   dengrad
  double precision ::    enterr,   entgrad
  double precision ::    velerr,   velgrad
  double precision ::   temperr,  tempgrad
  double precision ::  presserr, pressgrad
  double precision ::    raderr,   radgrad
  integer          ::  max_denerr_lev,   max_dengrad_lev
  integer          ::  max_enterr_lev,   max_entgrad_lev
  integer          ::  max_velerr_lev,   max_velgrad_lev
  integer          ::  max_temperr_lev,  max_tempgrad_lev
  integer          ::  max_presserr_lev, max_pressgrad_lev
  integer          ::  max_raderr_lev,   max_radgrad_lev
  double precision ::  starBuffer, boundaryBuffer

  ! Determine if we are the I/O processor
  integer :: ioproc
  
  ! Model info
  character (len=80) :: model_P_name
  character (len=80) :: model_S_name
  
  double precision, allocatable :: model_P_r(:), model_P_state(:,:)
  double precision, allocatable :: model_S_r(:), model_S_state(:,:)

  integer :: npts_model

  double precision :: mass_P, mass_S
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

  ! Inertial reference frame flag
  logical :: inertial

  ! Damping
  logical          :: damping
  double precision :: damping_alpha

  ! Relaxation
  logical          :: do_relax
  double precision :: relax_tau

  ! Grid info
  double precision :: center(3)

  ! Binary properties
  double precision :: a_P_initial, a_S_initial, a
  
  double precision, dimension(3) :: center_P_initial, center_S_initial

  integer :: star_axis

  ! Density of ambient medium
  double precision :: ambient_density

  ! Bulk system motion

  double precision :: bulk_velx, bulk_vely, bulk_velz

contains

  ! This routine calls all of the other subroutines at the beginning
  ! of a simulation to fill in the basic problem data.

  subroutine initialize(name, namlen, init)

    use bl_constants_module, only: ZERO
    use bl_error_module, only: bl_error

    implicit none

    integer :: namlen, i, init
    integer :: name(namlen)

    ! Build "probin" filename -- the name of the file containing the fortin namelist.
    allocate(character(len=namlen) :: probin)
    do i = 1, namlen
       probin(i:i) = char(name(i))
    enddo

    ! Read in the namelist to set problem parameters.

    call read_namelist

    ! Set small_pres and small_ener.

    call set_small

    ! Complete any grid related calculations like defining its center.

    call grid_data

    ! Determine if we are the I/O processor, and save it to the ioproc variable.

    call get_ioproc

    ! Establish binary parameters and create initial models.

    if (init == 1) then 
       call binary_setup
    endif

  end subroutine



  ! This routine reads in the namelist

  subroutine read_namelist

    use bl_constants_module, only: ZERO, HALF, ONE
    use meth_params_module

    implicit none

    integer :: untin, i

    integer :: iC12, iO16
    double precision :: stellar_C12, stellar_O16

    namelist /fortin/ &
         mass_p, mass_s, &
         central_density_P, central_density_S, &
         nsub, &
         inertial, &
         interp_temp, &
         damping, damping_alpha, &
         do_relax, relax_tau, &
         ambient_density, &
         stellar_temp, stellar_C12, stellar_O16, &
         denerr,     dengrad,   max_denerr_lev,   max_dengrad_lev, &
         velerr,     velgrad,   max_velerr_lev,   max_velgrad_lev, &
         presserr, pressgrad, max_presserr_lev, max_pressgrad_lev, &
         temperr,   tempgrad,  max_temperr_lev,  max_tempgrad_lev, &
         starBuffer, boundaryBuffer, star_axis, &
         bulk_velx, bulk_vely, bulk_velz

    ! Set namelist defaults
    denerr = 1.d20
    dengrad = 1.d20
    max_denerr_lev = 10
    max_dengrad_lev = 10

    presserr = 1.d20
    pressgrad = 1.d20
    max_presserr_lev = -1
    max_pressgrad_lev = -1

    velerr  = 1.d0
    velgrad = 1.d20
    max_velerr_lev = -1
    max_velgrad_lev = -1

    temperr  = 1.d0
    tempgrad = 1.d20
    max_temperr_lev = -1
    max_tempgrad_lev = -1

    starBuffer = 1.5d0
    boundaryBuffer = 0.6d0

    nsub = 1

    central_density_P = -ONE
    central_density_S = -ONE

    mass_p = 1.0
    mass_s = 1.0

    stellar_temp = 1.0d7
    stellar_C12  = HALF
    stellar_O16  = HALF

    smallx = 1.d-10
    smallu = 1.d-12

    ambient_density = 1.d-4

    inertial = .false.
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

    use meth_params_module, only: small_temp, small_pres, small_dens !, small_ener

    implicit none

    type (eos_t) :: eos_state

    ! Given the inputs of small_dens and small_temp, figure out small_pres.
 
    eos_state % rho = small_dens
    eos_state % T   = small_temp
    eos_state % xn  = stellar_comp
 
    call eos(eos_input_rt, eos_state, .false.)

    small_pres = eos_state % p
!    small_ener = eos_state % e

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



  ! Do any grid related calculations

  subroutine grid_data

    use bl_constants_module, only: HALF
    use prob_params_module, only: xmin, xmax, ymin, ymax, zmin, zmax

    implicit none

    center(1) = HALF * (xmax - xmin)
    center(2) = HALF * (ymax - ymin)
    center(3) = HALF * (zmax - zmin)

  end subroutine grid_data



  ! Set up a binary simulation

  subroutine binary_setup

    use bl_constants_module, only: ZERO, ONE
    use meth_params_module, only: rot_period
    use initial_model_module, only: init_1d

    implicit none

    double precision :: dx
    type (eos_t) :: ambient_state

    call get_ambient(ambient_state)

    center_P_initial = center
    center_S_initial = center

    radius_P_initial = ZERO
    radius_S_initial = ZERO

    npts_model = 1024

    dx = 1.0d6

    ! Allocate arrays to hold the stellar models

    allocate(model_P_state(npts_model,3+nspec))
    allocate(model_S_state(npts_model,3+nspec))

    allocate(model_P_r(npts_model))
    allocate(model_S_r(npts_model))

    ! Generate primary and secondary WD

    call init_1d(model_P_r, model_P_state, npts_model, dx, radius_P_initial, mass_P, &
                 central_density_P, stellar_temp, stellar_comp, ambient_state)

    if (ioproc == 1) then
        print *, "Generated initial model for primary WD of mass", mass_P, &
                 ", central density", central_density_P, ", and radius", radius_P_initial
    endif

    if (mass_S > ZERO) then

       call init_1d(model_S_r, model_S_state, npts_model, dx, radius_S_initial, mass_S, &
                    central_density_S, stellar_temp, stellar_comp, ambient_state)

       if (ioproc == 1) then
          print *, "Generated initial model for secondary WD of mass", mass_S, &
                   ", central density", central_density_S, ", and radius", radius_S_initial
       endif

       ! Get the orbit from Kepler's third law

       call kepler_third_law(radius_P_initial, mass_P, radius_S_initial, mass_S, &
                             rot_period, a_P_initial, a_S_initial, a, &
                             orbital_speed_P, orbital_speed_S)

       ! Star center positions -- we'll put them in the midplane on the
       ! axis specified by star_axis, with the center of mass at the center of the domain.
       ! If we're only doing a single star, which we know because the secondary mass is
       ! less than zero, then just leave the star in the center.

       center_P_initial(star_axis) = center_P_initial(star_axis) - a_P_initial
       center_S_initial(star_axis) = center_S_initial(star_axis) + a_S_initial

    endif

  end subroutine binary_setup



  ! Accepts the masses of two stars (in solar masses)
  ! and the orbital period of a system,
  ! and returns the semimajor axis of the orbit (in cm),
  ! as well as the distances a_1 and a_2 from the center of mass.

  subroutine kepler_third_law(radius_1, mass_1, radius_2, mass_2, period, a_1, a_2, a, v_1, v_2)

    use bl_constants_module, only: ONE, TWO, FOUR, THIRD, M_PI
    use fundamental_constants_module, only: Gconst, M_solar
    use prob_params_module, only: xmin, xmax, ymin, ymax, zmin, zmax

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

    if (length > (xmax-xmin) .or. &
        length > (ymax-ymin) .or. &
        length > (zmax-zmin)) then
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

    if ( mass_p < mass_s ) then

      if (ioproc == 1) then
        print *, "Primary mass is less than secondary mass; switching the stars."
      endif

      temp_mass = mass_p
      mass_p = mass_s
      mass_s = temp_mass

    endif

  end subroutine ensure_primary_mass_larger

end module probdata_module
