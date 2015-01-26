!!****if* source/Simulation/SimulationMain/Sod/Simulation_data
!!
!! NAME
!!  Simulation_data
!!
!! SYNOPSIS
!!
!!  use Simulation_data
!!
!! DESCRIPTION
!!
!!  Store the simulation data for the Sod problem
!!
!! ARGUMENTS
!!
!!
!! PARAMETERS
!!
!!  sim_rhoLeft    Density in the left part of the grid
!!  sim_rhoRight   Density in the right part of the grid
!!  sim_TLeft      Temperature in the left part of the grid
!!  sim_TRight     Temperature in the right part of the grid
!!  sim_uLeft      fluid velocity in the left part of the grid
!!  sim_uRight     fluid velocity in the right part of the grid
!!
!!   
!!
!!***

module Simulation_data
#include "Flash.h"
  implicit none

  !! *** Runtime Parameters *** !!

  real, save :: sim_rhoLeft, sim_rhoRight, sim_TLeft, sim_TRight
  real, save :: sim_uLeft, sim_uRight, sim_xAngle, sim_yAngle, sim_posn
  real, save :: sim_gamma, sim_smallP, sim_smallX

  !! *** Variables pertaining to Simulation Setup 'Sod' *** !!
  logical, save :: sim_gCell
  real, save :: xmin, xmax

  integer, save :: sim_meshMe
end module Simulation_data


