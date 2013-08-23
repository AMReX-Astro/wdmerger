module Simulation_data
#include "Flash.h"
#include "constants.h"
#include "Eos.h"

  real, allocatable,dimension(:), save :: sim_wdp_rad_tab, sim_wdp_dens_tab, sim_wdp_temp_tab, sim_wdp_c12_tab, sim_wdp_ne22_tab
  real, save :: sim_wdp_dr_inv, sim_wdp_mass, sim_wdp_radius
  integer, save :: sim_wdp_npnts

  real, allocatable, dimension(:), save :: sim_wds_rad_tab, sim_wds_dens_tab, sim_wds_temp_tab, sim_wds_c12_tab, sim_wds_ne22_tab
  real, save :: sim_wds_dr_inv, sim_wds_mass, sim_wds_radius
  integer, save :: sim_wds_npnts

  real, save :: sim_wd_sep, sim_wdp_a, sim_wds_a

  real, save :: sim_densFluff, sim_tempFluff, sim_xc12Fluff, sim_xne22Fluff

  real, save :: sim_binaryPeriod

  real, save :: sim_smallx, sim_smallt, sim_smallu, sim_smalle, sim_smallp

  logical, save :: sim_ignite

  logical, save :: sim_ignitionFile
  integer, save :: sim_ign_numpnts
  real, allocatable, dimension(:), save :: sim_ignX, sim_ignY, sim_ignZ, sim_ignR

  logical, save :: sim_ignMPole
  real, save :: sim_ignMPoleA
  integer, save :: sim_ignMpoleMinL, sim_ignMpoleMaxL, sim_ignMpoleSeed
  integer, save :: sim_geom
  real, allocatable, dimension(:), save  :: mp_A

  logical, save :: sim_ignSin
  real, save    :: sim_ignSinN, sim_ignSinA

  real,    save           :: sim_laminarWidth

  real, save :: sim_refFluffDensThresh, sim_refFluffMargin
  real, save :: sim_refNogenEnucThresh, sim_refNogenFldtThresh, sim_refNogenMargin
  integer, save :: sim_refFluffLevel, sim_refNogenLevel
  real, save :: sim_refCentRegionDist
  integer, save :: sim_refCentRegionLevel

  real, save :: sim_vrms_reduced, sim_vrms_center, sim_vrms_Tc, sim_vrms_T0, sim_vrms_alpha
  logical, save :: sim_read_turbfield
  integer, save :: sim_smooth_level
  character(len=4096), save :: sim_turbfield_filename
  real, save :: sim_turbfield_bbox(IAXIS:KAXIS,LOW:HIGH)

  integer, save :: sim_meshMe

end module Simulation_data
