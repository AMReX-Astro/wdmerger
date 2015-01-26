module Simulation_data
#include "Flash.h"
#include "constants.h"
#include "Eos.h"

  real, save :: sim_star_radius

  real, save :: sim_densFluff, sim_tempFluff, sim_densStar, sim_tempStar

  real, save :: sim_x_velocity, sim_y_velocity

  integer, save :: sim_meshMe

end module Simulation_data
