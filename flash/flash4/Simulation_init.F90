!!***
!!
!! NAME
!!
!!  Simulation_init
!!
!! SYNOPSIS
!!
!!  Simulation_init()
!!
!! DESCRIPTION
!!  Initialize private data for the 3-stage flame test setup
!!
!! ARGUMENTS
!!
!!
!!
!!***

subroutine Simulation_init()

  use Simulation_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface, only : Logfile_stampMessage
  use Grid_interface, only : Grid_getGeometry
  use Logfile_interface, only : Logfile_stampMessage
  use Driver_interface, only : Driver_abortFlash
  use PhysicalConstants_interface, ONLY: PhysicalConstants_get

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"

  character(len=256) :: initialWDpFileName
  character(len=256) :: initialWDsFileName
  integer :: istat, i
  real :: radius, r_l, r_r

  logical :: found_radius

  real, parameter :: fourthirds = 4./3.
  real, parameter :: third = 1./3.
  real, save :: newton ! Newton's gravity

  !--------------------------------------------------------
  !  initialize runtime parameters and some other constants
  !--------------------------------------------------------
  call RuntimeParameters_get('initialWDpFile', initialWDpFileName)
  call RuntimeParameters_get('initialWDsFile', initialWDsFileName)

  call RuntimeParameters_get('binary_period', sim_binaryPeriod)
  call RuntimeParameters_get('inertial', sim_inertial)
  call RuntimeParameters_get('single_star', sim_single_star)

  call RuntimeParameters_get('dens_fluff', sim_densFluff)
  call RuntimeParameters_get('temp_fluff', sim_tempFluff)
  call RuntimeParameters_get('xc12_fluff', sim_xc12Fluff)
  call RuntimeParameters_get('xne22_fluff', sim_xne22Fluff)

  call RuntimeParameters_get('smallx', sim_smallx)
  call RuntimeParameters_get('smallp', sim_smallp)
  call RuntimeParameters_get('smallt', sim_smallt)
  call RuntimeParameters_get('smallu', sim_smallu)
  call RuntimeParameters_get('smalle', sim_smalle)

  ! Store physical constants:
  call PhysicalConstants_get("Newton", newton)

  call Driver_getMype(MESH_COMM, sim_meshMe)

  !--------------------------------------------------------
  !  read in 1d initial wd profile
  !--------------------------------------------------------

  ! First we'll read in the file for the primary WD.

  if (sim_meshMe == MASTER_PE) then
    write(6,*) 'Reading primary WD file'
  endif

  call Logfile_stampMessage('[Simulation_init] Reading initial 1-d WD profile, primary')
  open(unit=2,file=initialWDpFileName,status='OLD',iostat=istat)
  if (istat /= 0) call Driver_abortFlash('Unable to open initial WD profile, primary')

  ! eat header
  read(2,*)
  read(2,*) sim_wdp_npnts

  allocate(sim_wdp_rad_tab(sim_wdp_npnts),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate space for initial WD")
  allocate(sim_wdp_dens_tab(sim_wdp_npnts),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate space for initial WD")
  allocate(sim_wdp_temp_tab(sim_wdp_npnts),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate space for initial WD")
  allocate(sim_wdp_c12_tab(sim_wdp_npnts),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate space for initial WD")
  allocate(sim_wdp_ne22_tab(sim_wdp_npnts),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate space for initial WD")

  r_l = 0.0
  r_r = 0.0
  radius = 0.0
  sim_wdp_mass = 0.0
  found_radius = .false.

  do i = 1, sim_wdp_npnts
     read(2,*) sim_wdp_rad_tab(i), sim_wdp_dens_tab(i), sim_wdp_temp_tab(i), sim_wdp_c12_tab(i), sim_wdp_ne22_tab(i)
  enddo

  do i = 1, sim_wdp_npnts-1

     if ( i /= 1 ) then
       r_l = 0.5 * ( sim_wdp_rad_tab(i) + sim_wdp_rad_tab(i-1) )
     endif

     r_r = 0.5 * ( sim_wdp_rad_tab(i) + sim_wdp_rad_tab(i+1) )

     sim_wdp_mass = sim_wdp_mass + (r_r - r_l)*(r_r**2 + r_r*r_l + r_l**2)*sim_wdp_dens_tab(i) 

     ! Determine the radius of the primary by finding when it comes close to the fluff density.

     if ( (sim_wdp_dens_tab(i) <= 1.10d0 * sim_densFluff) .and. (.not. found_radius) ) then
       sim_wdp_radius = sim_wdp_rad_tab(i)
       found_radius = .true.
     endif

  enddo
  close(2)

  if ( .not. found_radius ) then
    sim_wdp_radius = sim_wdp_rad_tab(sim_wdp_npnts)
  endif
  sim_wdp_mass = sim_wdp_mass*fourthirds*PI
  sim_wdp_dr_inv = 1.0 / (r_r - r_l)

  if (sim_meshMe == MASTER_PE) then
    write(6,*) 'Mass of primary =', sim_wdp_mass
  endif

  ! now secondary

  if (sim_meshMe == MASTER_PE) then
    write(6,*) 'Reading secondary WD file'
  endif

  call Logfile_stampMessage('[Simulation_init] Reading initial 1-d WD profile,secondary')
  open(unit=2,file=initialWDsFileName,status='OLD',iostat=istat)
  if (istat /= 0) call Driver_abortFlash('Unable to open initial WD profile, secondary')

  ! eat header
  read(2,*)
  read(2,*) sim_wds_npnts

  allocate(sim_wds_rad_tab(sim_wds_npnts),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate space for initial WD")
  allocate(sim_wds_dens_tab(sim_wds_npnts),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate space for initial WD")
  allocate(sim_wds_temp_tab(sim_wds_npnts),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate space for initial WD")
  allocate(sim_wds_c12_tab(sim_wds_npnts),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate space for initial WD")
  allocate(sim_wds_ne22_tab(sim_wds_npnts),STAT=istat)
    if (istat /= 0) call Driver_abortFlash("Cannot allocate space for initial WD")

  r_l = 0.0
  r_r = 0.0
  radius = 0.0
  sim_wds_mass = 0.0
  found_radius = .false.

  do i = 1, sim_wds_npnts
     read(2,*) sim_wds_rad_tab(i), sim_wds_dens_tab(i), sim_wds_temp_tab(i), sim_wds_c12_tab(i), sim_wds_ne22_tab(i)
  enddo

  do i = 1, sim_wds_npnts-1

     if ( i /= 1 ) then
       r_l = 0.5 * ( sim_wds_rad_tab(i) + sim_wds_rad_tab(i-1) )
     endif

     r_r = 0.5 * ( sim_wds_rad_tab(i) + sim_wds_rad_tab(i+1) )

     sim_wds_mass = sim_wds_mass + (r_r - r_l)*(r_r**2 + r_r*r_l + r_l**2)*sim_wds_dens_tab(i) 

     ! Determine the radius of the secondary by finding when it comes close to the fluff density.

     if ( (sim_wds_dens_tab(i) <= 1.10d0 * sim_densFluff) .and. (.not. found_radius) ) then
       sim_wds_radius = sim_wds_rad_tab(i)
       found_radius = .true.
     endif

  enddo
  close(2)

  if ( .not. found_radius ) then
    sim_wds_radius = sim_wds_rad_tab(sim_wds_npnts)
  endif

  sim_wds_mass = sim_wds_mass*fourthirds*PI
  sim_wds_dr_inv = 1.0 / (r_r - r_l)


  if (sim_meshMe == MASTER_PE) then
    write(6,*) 'Mass of secondary =', sim_wds_mass
  endif

  !----------------------------------------------------------
  ! Now solve the Kepler problem for the configurations of the stars.
  ! This is in line with how this is solved in the wdmerger problem for CASTRO.
  !----------------------------------------------------------

  ! Orbital separation a = ( G * M * P**2 / (4 * pi**2) )**3
  sim_wd_sep = (newton*(sim_wdp_mass + sim_wds_mass)*sim_binaryPeriod**2/(4.0*PI**2))**third
  
  ! Semimajor axes
  sim_wds_a = sim_wd_sep / (1.0 + sim_wds_mass/sim_wdp_mass)
  sim_wdp_a = (sim_wds_mass/sim_wdp_mass)*sim_wds_a

  if (sim_meshMe == MASTER_PE) then
    print *, 'Semimajor axis: '          , sim_wd_sep
    print *, 'Primary semimajor axis: '  , sim_wdp_a
    print *, 'Secondary semimajor axis: ', sim_wds_a
  endif

  if (.not. sim_single_star) then

    ! Make sure the stars are not touching.
    if (sim_wdp_radius + sim_wds_radius > sim_wd_sep) then
      call Driver_abortFlash('Stars are touching.')
    endif

  endif
  
  ! add check to see if grid is big enough here?



end subroutine Simulation_init
