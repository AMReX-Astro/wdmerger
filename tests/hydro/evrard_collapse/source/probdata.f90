module probdata_module

  use network, only: nspec, network_species_index
  use eos_type_module, only: eos_t
  use eos_module, only: eos_input_rt, eos

  ! Probin file
  character (len=:), allocatable :: probin

  ! Determine if we are the I/O processor
  integer :: ioproc
  
  double precision, allocatable :: model_r(:), model_state(:,:)

  integer :: npts_model

  double precision :: mass, radius
  double precision :: stellar_comp(nspec)

  ! Smallest allowed mass fraction
  double precision :: smallx

  ! Smallest allowed velocity on the grid
  double precision :: smallu
  
  ! Density of ambient gas around the star
  double precision :: ambient_density
  
  ! Controls interpolation from 1D model to 3D model
  integer :: nsub

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

    call set_small

  end subroutine



  ! This routine reads in the namelist

  subroutine read_namelist

    use bl_constants_module, only: ZERO, HALF
    use meth_params_module

    implicit none

    integer :: untin, i

    integer :: iC12, iO16
    double precision :: stellar_C12, stellar_O16

    namelist /fortin/ &
         mass, radius, &
         ambient_density, &
         nsub, &
         stellar_C12, stellar_O16

    ! Set namelist defaults

    nsub = 1

    mass   = 1.0d0
    radius = 1.0d9

    ambient_density = 1.0d-4

    stellar_C12  = HALF
    stellar_O16  = HALF

    smallx = 1.d-10

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

  end subroutine read_namelist



  ! Determine if we are the I/O processor, and save it to the ioproc variable

  subroutine get_ioproc

    implicit none

    ! For outputting -- determine if we are the IO processor
    call bl_pd_is_ioproc(ioproc)

  end subroutine get_ioproc



  ! Define the ambient state, and calculate small_pres

  subroutine set_small

    use meth_params_module, only: small_temp, small_pres, small_dens, small_ener

    implicit none

    type (eos_t) :: eos_state

    smallu = 1.d-12

    ! Given the inputs of small_dens and small_temp, figure out small_pres.
 
    eos_state % rho = small_dens
    eos_state % T   = small_temp
    eos_state % xn  = stellar_comp
 
    call eos(eos_input_rt, eos_state, .false.)

    small_pres = eos_state % p
    small_ener = eos_state % e

  end subroutine set_small

end module probdata_module
