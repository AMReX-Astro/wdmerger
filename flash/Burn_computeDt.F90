!!****f* source/physics/sourceTerms/Burn/Burn_computeDt
!!
!! NAME
!!
!!  Burn_computeDt
!!
!! SYNOPSIS
!!
!!  Burn_computeDt(integer(IN) :: blockID,
!!                 integer(IN) :: blkLimits(2,MDIM)
!!                 integer(IN) :: blkLimitsGC(2,MDIM)
!!                 real,pointer::  solnData(:,:,:,:),   
!!                 real(OUT)   :: dt_burn, 
!!                 real(OUT)   :: dt_minloc(5)) 
!!
!!
!! DESCRIPTION
!!
!!  compute a burning timestep limiter, by trying to force the energy
!!  generation from burning to be smaller than the internal energy
!!  in a zone.
!!
!!   The timestep limiter would be:
!!
!!                                      eint
!!             dt     =  enucDtFactor * -----
!!               burn                   enuc
!!
!!  enuc is energy/volume/s, so the time factor is already in there, and we
!!  are actually doing
!!
!!             
!!                                      eint
!!             dt     =  enucDtFactor * -----    * dt
!!               burn                   enuc*dt
!!
!!  enuc*dt is the amount of energy / volume deposited in a zone by burning. 
!!  eint is the internal energy / volume in that zone.  If enuc*dt is 2x
!!  eint, then we want a timestep that is half the size.  
!!
!!  enucDtFactor is a prefactor to scaling the timestep.  In general, we aim
!!  for enuc*dt < enucDtFactor * eint.  For good coupling between the hydro
!!  and the burner, enucDtFactor should be < 1.
!!
!!
!! ARGUMENTS
!!
!!  blockID       --  local block ID
!!  blkLimits     --  the indices for the interior endpoints of the block
!!  blkLimitsGC   --  the indices for endpoints including the guardcells
!!  solnData      --  the physical, solution data from grid
!!  dt_burn       --  variable to hold timestep constraint
!!  dt_minloc(5)  --  array to hold limiting zone info:  zone indices
!!                    (i,j,k), block ID, PE number
!!
!!
!! PARAMETERS
!!
!!  enucDtFactor    A parameter, such that enuc*dt < enucDtFactor * eint,
!!                  that is, the energy release from burning divided by
!!                  the internal energy in that zone is < enucDtFactor.
!!
!! SEE ALSO
!!
!!  Driver_computeDt
!!
!! NOTE
!!  
!!  This routine is not implemented yet due to conflicts with the required
!!  magnitude of enucDtFactor
!!
!!***

subroutine Burn_computeDt(blockID,                       & 
                           blkLimits,blkLimitsGC,        &
                           solnData,                     &
                           dt_burn, dt_minloc)

! this version is the stub and should therefore not do anything.
! value of dt_burn is unmodified so it is not limiting.

#include "constants.h"


  implicit none


  !! arguments
  integer, intent(IN)   :: blockID
  integer, intent(IN),dimension(2,MDIM)::blkLimits,blkLimitsGC
  real, pointer           :: solnData(:,:,:,:) 
  real, intent(INOUT)     :: dt_burn
  integer, intent(INOUT)  :: dt_minloc(5)

  return

end subroutine Burn_computeDt


