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

  ! Determine if we are the I/O processor
  integer :: ioproc
  
  double precision, allocatable :: model_r(:), model_state(:,:)

  integer :: npts_model

  double precision :: P_c, rho_c, mass, radius
  double precision :: polytrope_temp, polytrope_comp(nspec)

  ! Fraction by which we decrease the central pressure
  double precision :: eta, central_frac

  ! Smallest allowed mass fraction
  double precision :: smallx

  ! Smallest allowed velocity on the grid
  double precision :: smallu
  
  ! EOS state type that describes the ambient gas around the polytrope
  type (eos_t) :: ambient_state
  double precision :: ambient_density
  
  ! Controls interpolation from 1D model to 3D model
  integer :: nsub
  logical :: interp_temp

  ! Grid info
  double precision :: center(3)

contains

  ! This routine calls all of the other subroutines at the beginning
  ! of a simulation to fill in the basic problem data.

  subroutine initialize(name, namlen)

    use bl_constants_module, only: ZERO
    use bl_error_module, only: bl_error

    implicit none

    integer :: namlen, i
    integer :: name(namlen)

    ! Build "probin" filename -- the name of the file containing the fortin namelist.
    allocate(character(len=namlen) :: probin)
    do i = 1, namlen
       probin(i:i) = char(name(i))
    enddo

    ! Determine if we are the I/O processor, and save it to the ioproc variable

    call get_ioproc

    ! Read in the namelist to set problem parameters

    call read_namelist

    ! Finalize ambient state, and get small_pres

    call set_ambient_and_small

    ! Complete any grid related calculations like defining its center

    call grid_data

    ! Set up the polytropic sphere

    call polytrope_setup

  end subroutine



  ! This routine reads in the namelist

  subroutine read_namelist

    use bl_constants_module, only: ZERO, HALF
    use meth_params_module

    implicit none

    integer :: untin, i

    integer :: iC12, iO16
    double precision :: polytrope_C12, polytrope_O16

    namelist /fortin/ &
         eta, rho_c, central_frac, &
         ambient_density, &
         nsub, &
         polytrope_temp, polytrope_C12, polytrope_O16, &
         denerr,     dengrad,   max_denerr_lev,   max_dengrad_lev, &
         velerr,     velgrad,   max_velerr_lev,   max_velgrad_lev, &
         presserr, pressgrad, max_presserr_lev, max_pressgrad_lev, &
         temperr,   tempgrad,  max_temperr_lev,  max_tempgrad_lev

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

    nsub = 1

    eta = 0.93d0
    central_frac = 0.2d0

    rho_c = 1.0d7

    ambient_density = 1.0d-10

    polytrope_temp = 1.0d7
    polytrope_C12  = HALF
    polytrope_O16  = HALF

    smallx = 1.d-10

    ! Read namelist to override the defaults
    untin = 9 
    open(untin,file=probin,form='formatted',status='old')
    read(untin,fortin)
    close(unit=untin)

    ! Fill in stellar composition using C12 and O16

    iC12 = network_species_index("carbon-12")
    iO16 = network_species_index("oxygen-16")

    polytrope_comp(1:nspec) = smallx
    polytrope_comp(iC12)    = polytrope_C12 - (nspec - 2) * smallx
    polytrope_comp(iO16)    = polytrope_O16 - (nspec - 2) * smallx

  end subroutine read_namelist



  ! Determine if we are the I/O processor, and save it to the ioproc variable

  subroutine get_ioproc

    implicit none

    ! For outputting -- determine if we are the IO processor
    call bl_pd_is_ioproc(ioproc)

  end subroutine get_ioproc



  ! Define the ambient state, and calculate small_pres

  subroutine set_ambient_and_small

    use meth_params_module, only: small_temp, small_pres, small_dens, small_ener

    implicit none

    type (eos_t) :: eos_state

    smallu = 1.d-12

    ! Define ambient state and call EOS to get eint and pressure

    ambient_state % rho = ambient_density * rho_c
    ambient_state % T   = polytrope_temp
    ambient_state % xn  = polytrope_comp

    call eos(eos_input_rt, ambient_state, .false.)

    ! Given the inputs of small_dens and small_temp, figure out small_pres.
 
    eos_state % rho = small_dens
    eos_state % T   = small_temp
    eos_state % xn  = polytrope_comp
 
    call eos(eos_input_rt, eos_state, .false.)

    small_pres = eos_state % p
    small_ener = eos_state % e

  end subroutine set_ambient_and_small



  ! Do any grid related calculations

  subroutine grid_data

    use bl_constants_module, only: HALF
    use prob_params_module, only: xmin, xmax, ymin, ymax, zmin, zmax

    implicit none

    center(1) = HALF * (xmax - xmin)
    center(2) = HALF * (ymax - ymin)
    center(3) = HALF * (zmax - zmin)

  end subroutine grid_data



  subroutine polytrope_setup

    use bl_constants_module, only: ZERO
    use initial_model_module, only: init_1d
    use model_parser_module, only: itemp_model, idens_model, ipres_model, ispec_model
    use eos_module, only: eos_set_polytrope_parameters, eos_get_polytrope_parameters

    implicit none

    double precision :: dx

    integer :: polytrope
    double precision :: gamma, K_const, mu_e

    radius = ZERO

    npts_model = 1024

    dx = 1.0d7

    mass = -1.d0

    ! Allocate arrays to hold the stellar models

    allocate(model_state(npts_model,3+nspec))
    allocate(model_r(npts_model))

    ! Generate the polytropic sphere

    call init_1d(model_r, model_state, npts_model, dx, radius, mass, rho_c, &
                 polytrope_temp, polytrope_comp, ambient_state)

    if (ioproc == 1) then
        print *, "Generated initial model for WD of mass", mass, &
                 ", central density", model_state(1,idens_model), &
                 ", central pressure", model_state(1,ipres_model), &
                 ", and radius", radius

    endif

  end subroutine polytrope_setup

end module probdata_module
