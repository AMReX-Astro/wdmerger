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
  use Flame_interface, ONLY : Flame_rhJumpReactive, Flame_getWidth, Flame_laminarSpeed
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use bn_paraInterface, ONLY : bn_paraFuelAshProperties
  use Logfile_interface, only : Logfile_stampMessage
  use Grid_interface, only : Grid_getGeometry
  use sim_local_interface, only : sim_LCGRandomIterate
  use Logfile_interface, only : Logfile_stampMessage
  use Driver_interface, only : Driver_abortFlash
  use PhysicalConstants_interface, ONLY: PhysicalConstants_get

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"

  character(len=256) :: initialWDpFileName
  character(len=256) :: initialWDsFileName
  character(len=4096) :: ignitionFileName
  character(len=80) :: reportbuf
  integer :: istat, i
  real :: radius, r_l, r_r
  integer :: ignfileunit=2

  integer :: try, l
  logical :: accept, ignitionFile
  logical :: found_radius
  real    :: u, v, r
  real    :: deltaplus, deltaminus, costheta, P_lm1, P_l, hold

  real, parameter :: fourthirds = 4./3.
  real, parameter :: third = 1./3.
  real, save :: newton ! Newton's gravity
!  real, save :: period ! Binary Period
!fix below
!  period = 40.0

  !--------------------------------------------------------
  !  initialize runtime parameters and some other constants
  !--------------------------------------------------------
  call RuntimeParameters_get( 'initialWDpFile', initialWDpFileName)
  call RuntimeParameters_get( 'initialWDsFile', initialWDsFileName)

  call RuntimeParameters_get( 'binary_period', sim_binaryPeriod)

  call RuntimeParameters_get('dens_fluff', sim_densFluff)
  call RuntimeParameters_get('temp_fluff', sim_tempFluff)
  call RuntimeParameters_get('xc12_fluff', sim_xc12Fluff)
  call RuntimeParameters_get('xne22_fluff', sim_xne22Fluff)

  call RuntimeParameters_get('ignite', sim_ignite)
  call RuntimeParameters_get('ignition_file', ignitionFile)
  
  call RuntimeParameters_get('refFluffDensThresh', sim_refFluffDensThresh)
  call RuntimeParameters_get('refFluffMargin', sim_refFluffMargin)
  call RuntimeParameters_get('refFluffLevel', sim_refFluffLevel)

  call RuntimeParameters_get('refNogenEnucThresh', sim_refNogenEnucThresh)
  call RuntimeParameters_get('refNogenFldtThresh', sim_refNogenFldtThresh)
  call RuntimeParameters_get('refNogenMargin', sim_refNogenMargin)
  call RuntimeParameters_get('refNogenLevel', sim_refNogenLevel)

  call RuntimeParameters_get('refCentRegionDist', sim_refCentRegionDist)
  call RuntimeParameters_get('refCentRegionLevel', sim_refCentRegionLevel)
  
  ! Why is G here but not pi?
  ! Store physical constants:
  call PhysicalConstants_get("Newton", newton)

  ! only need to get width of artificial flame once
  call Flame_getWidth(sim_laminarWidth)

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

  ! Make sure the stars are not touching.
  if (sim_wdp_radius + sim_wds_radius > sim_wd_sep) then
    call Driver_abortFlash('Stars are touching.')
  endif
  
  ! add check to see if grid is big enough here?


  !----------------------------------------------------------
  ! Now get and set up information about ignition
  !----------------------------------------------------------
  if (ignitionFile) then
     call Logfile_stampMessage('[Simulation_init] Reading ignition points from file')
     call RuntimeParameters_get('ignition_file_name', ignitionFileName)
     open(unit=ignfileunit,file=ignitionFileName,status='OLD',iostat=istat)
     if (istat /= 0) call Driver_abortFlash('Unable to open ignition points file')
     ! one-line header ignored
     read(ignfileunit,*)
     read(ignfileunit,*) sim_ign_numpnts
     allocate(sim_ignX(sim_ign_numpnts),STAT=istat)
       if (istat /= 0) call Driver_abortFlash("Cannot allocate space for ignition points")
     allocate(sim_ignY(sim_ign_numpnts),STAT=istat)
       if (istat /= 0) call Driver_abortFlash("Cannot allocate space for ignition points")
     allocate(sim_ignZ(sim_ign_numpnts),STAT=istat)
       if (istat /= 0) call Driver_abortFlash("Cannot allocate space for ignition points")
     allocate(sim_ignR(sim_ign_numpnts),STAT=istat)
       if (istat /= 0) call Driver_abortFlash("Cannot allocate space for ignition points")
     do i = 1, sim_ign_numpnts
        read(ignfileunit,*) sim_ignX(i), sim_ignY(i), sim_ignZ(i), sim_ignR(i)
     enddo
     close(ignfileunit)
  else
     ! just a single ignition point, read from parameter file
     sim_ign_numpnts = 1
     allocate(sim_ignX(1),STAT=istat)
       if (istat /= 0) call Driver_abortFlash("Cannot allocate space for ignition points")
     allocate(sim_ignY(1),STAT=istat)
       if (istat /= 0) call Driver_abortFlash("Cannot allocate space for ignition points")
     allocate(sim_ignZ(1),STAT=istat)
       if (istat /= 0) call Driver_abortFlash("Cannot allocate space for ignition points")
     allocate(sim_ignR(1),STAT=istat)
       if (istat /= 0) call Driver_abortFlash("Cannot allocate space for ignition points")
     call RuntimeParameters_get('x_match', sim_ignX(1))
     call RuntimeParameters_get('y_match', sim_ignY(1))
     call RuntimeParameters_get('z_match', sim_ignZ(1))
     call RuntimeParameters_get('r_match', sim_ignR(1))
  endif
   

  call RuntimeParameters_get('ignMPole', sim_ignMpole)
  call RuntimeParameters_get('ignMPoleA', sim_ignMpoleA)
  call RuntimeParameters_get('ignMPoleMinL', sim_ignMpoleMinL)
  call RuntimeParameters_get('ignMPoleMaxL', sim_ignMpoleMaxL)
  call RuntimeParameters_get('ignMPoleSeed', sim_ignMpoleSeed)

  if (sim_ignMPole) then
     ! set up coefficients of multipole contributions
     ! check some things first to catch mistakes
     if (NDIM < 2) call Driver_abortFlash("[Simulation_init] multipole ignition not sensible in 1-dimension")
     call Grid_getGeometry(sim_geom)
     if (sim_geom/=CYLINDRICAL) call Driver_abortFlash("[Simulation_init] multipole ignition only implemented in 2d")

     allocate(mp_A(sim_ignMpoleMinL:sim_ignMpoleMaxL), stat=istat)
     if (istat/=0) call Driver_abortFlash("Unable to allocate mp_A in Simulation_init")

     try = 0
     accept = .false.
     do while ( .not. accept )

        try = try + 1

        ! generate set of random coefficients
        do i = sim_ignMpoleMinL, sim_ignMpoleMaxL
           ! generate a normally distributed number, r, via the Box-Muller method
           call sim_LCGRandomIterate(sim_ignMPoleSeed)
           ! avoid taking log of zero
           if (u==0) u = 1
           u = real(sim_ignMPoleSeed)/2147483646.0
           call sim_LCGRandomIterate(sim_ignMPoleSeed)
           v = real(sim_ignMPoleSeed)/2147483646.0
           r = sqrt(-2*log(u))*cos(2*PI*v)
          
           ! now set amplitude of this component, this will be multiplied by
           ! the legendre polynomial of order i to construct the initial flame surface
           ! normalization is chosen to match spherical harmonics ( with m=0 ) as
           ! given in Jackson.
           mp_A(i) = r*sim_ignMpoleA * sqrt((2*i+1)/4.0/PI)

           ! another initialization which turned out to be non-stardard
           ! this was used for some early testing and is just kept for reference
           !call sim_LCGRandomIterate(sim_ignMPoleSeed)
           ! random phase angle, between 0 and pi since we will not use imaginary part
           !alpha = real(sim_ignMPoleSeed)/2147483646.0 * PI
           ! coefficients have amplitudes of equal complex amplitude, with random
           ! phase.  We then use real part.  Normalization is for spherical
           ! harmonics in Jackson
           !mp_A(i) = sqrt((2*i+1)/4.0/PI)*cos(alpha)*sim_ignMpoleA*2/(sim_ignMpoleMaxL-sim_ignMpoleMinL+1)

           !print *, i,alpha/PI, mp_A(i)
        enddo

        ! now check that this is negative on + and - axis to help prevent pathologies
        deltaplus = 0.0
        costheta = 1.0
        ! Legendre polynomial
        P_lm1 = 1.0
        P_l = costheta
        do l = 1, sim_ignMPoleMaxL-1
           if ( l >= sim_ignMPoleMinL ) then
              deltaplus = deltaplus + mp_A(l)*P_l
           endif
           ! calculate next polynomial
           hold = P_l
           P_l = ( (2*l+1)*costheta*P_l - l*P_lm1 )/real(l+1)
           P_lm1 = hold
        enddo
        deltaplus = deltaplus + mp_A(sim_ignMpoleMaxL)*P_l

        deltaminus = 0.0
        costheta = -1.0
        ! Legendre polynomial
        P_lm1 = 1.0
        P_l = costheta
        do l = 1, sim_ignMPoleMaxL-1
           if ( l >= sim_ignMPoleMinL ) then
              deltaminus = deltaminus + mp_A(l)*P_l
           endif
           ! calculate next polynomial
           hold = P_l
           P_l = ( (2*l+1)*costheta*P_l - l*P_lm1 )/real(l+1)
           P_lm1 = hold
        enddo
        deltaminus = deltaminus + mp_A(sim_ignMpoleMaxL)*P_l

        if (deltaplus < 0.0 .and. deltaminus < 0.0) accept = .true.

     enddo


     ! need to report ending seed for statistically independent follow-on run
     write(reportbuf, *) '[Simulation_init] Multipole ignition coefficients initialized'
     call Logfile_stampMessage(reportbuf)
     write(reportbuf, *) '[Simulation_init] required iterations: ', try
     call Logfile_stampMessage(reportbuf)
     write(reportbuf, *) '[Simulation_init] endpoint ignMPoleSeed = ', sim_ignMPoleSeed
     call Logfile_stampMessage(reportbuf)
  endif

  call RuntimeParameters_get('ignSin', sim_ignSin)
  call RuntimeParameters_get('ignSinN', sim_ignSinN)
  call RuntimeParameters_get('ignSinA', sim_ignSinA)

  call RuntimeParameters_get('sim_smooth_level', sim_smooth_level)
  call RuntimeParameters_get('sim_vrms_center', sim_vrms_center)
  call RuntimeParameters_get('sim_vrms_Tc', sim_vrms_Tc)
  call RuntimeParameters_get('sim_vrms_T0', sim_vrms_T0)

  ! calculate factors for scaling velocity field
  sim_vrms_reduced = sim_vrms_center / 10**(0.8)
  sim_vrms_alpha = log(10**(0.8)) / log( (sim_vrms_Tc-sim_vrms_T0)/sim_vrms_T0*2 )

  call RuntimeParameters_get('sim_turbfield_filename',sim_turbfield_filename)
  call RuntimeParameters_get('sim_read_turbfield', sim_read_turbfield)
  call RuntimeParameters_get('sim_turbfield_xmin',sim_turbfield_bbox(IAXIS,LOW) )
  call RuntimeParameters_get('sim_turbfield_xmax',sim_turbfield_bbox(IAXIS,HIGH) )
  call RuntimeParameters_get('sim_turbfield_ymin',sim_turbfield_bbox(JAXIS,LOW) )
  call RuntimeParameters_get('sim_turbfield_ymax',sim_turbfield_bbox(JAXIS,HIGH) )
  call RuntimeParameters_get('sim_turbfield_zmin',sim_turbfield_bbox(KAXIS,LOW) )
  call RuntimeParameters_get('sim_turbfield_zmax',sim_turbfield_bbox(KAXIS,HIGH) )
 

end subroutine Simulation_init
