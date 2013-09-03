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

  logical, save :: sim_inertial

  integer, save :: sim_meshMe

end module Simulation_data
