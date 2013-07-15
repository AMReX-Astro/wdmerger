
! Dean Townsley 2008
!
! This is a linear congruential generator pseudorandom number generator.  This
! particular form has a seed/state/input/output (all the same thing) that is a
! signed integer, and therefore can be implemented portably and the seed can be
! easily recorded and specified in a parameter file.
! This particular LCG is described in
!   Monte Carlo Methods in Statistical physics
!   Newman and Barkema
!   Oxford University Press, 2004
!   page 386-388

! Note that this is suitable for modest-quality random numbers, not a large
! number of high quality ones

subroutine sim_LCGRandomIterate(state)
  implicit none
  integer, intent(inout) :: state

  integer, parameter :: a = 16807
  integer, parameter :: q = 127773
  integer, parameter :: r = 2836
  integer, parameter :: m = 2147483647

  state = a*mod(state, q) - (state/q)*r

  if (state < 0) state = state + m

end subroutine
