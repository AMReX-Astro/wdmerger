module binary_module

  implicit none

contains

  ! Given the mass ratio q of two stars (assumed to be q = M_1 / M_2), 
  ! compute the effective Roche radii of the stars, normalized to unity, 
  ! using the approximate formula of Eggleton (1983). We then 
  ! scale them appropriately using the current orbital distance.
  
  subroutine get_roche_radii(mass_ratio, r_1, r_2, a) bind(C)

    use bl_constants_module, only: ONE, TWO3RD, THIRD

    implicit none

    double precision, intent(in   ) :: mass_ratio, a
    double precision, intent(inout) :: r_1, r_2

    double precision :: q
    double precision :: c1, c2

    c1 = 0.49d0
    c2 = 0.60d0

    q = mass_ratio

    r_1 = a * c1 * q**(TWO3RD) / (c2 * q**(TWO3RD) + LOG(ONE + q**(THIRD)))

    q = ONE / q

    r_2 = a * c1 * q**(TWO3RD) / (c2 * q**(TWO3RD) + LOG(ONE + q**(THIRD)))

  end subroutine get_roche_radii


  
  ! Calculate Lagrange points. In each case we give the zone index
  ! closest to it (assuming we're on the coarse grid).

  subroutine get_lagrange_points(mass_1, mass_2, com_1, com_2, L1_idx) bind(C)

    use bl_constants_module
    use prob_params_module, only: dx_level
    use amrinfo_module, only: amr_level

    implicit none

    double precision :: mass_1, mass_2
    double precision :: com_1(3), com_2(3)
    integer          :: L1_idx(3)
    
    double precision :: r2 ! Distance from L1 to secondary
    double precision :: R  ! Distance between secondary and primary

    double precision :: r2_old ! For the root find
    double precision :: tolerance = 1.0d-6
    integer          :: max_iters = 100

    integer          :: i

    ! Don't try to calculate the Lagrange point if the secondary
    ! is already gone.

    if (mass_2 < ZERO) return

    R = sqrt(sum((com_2 - com_1)**2))

    ! Do a root-find over the quintic equation for L1. The initial
    ! guess will be the case of equal mass, i.e. that the Lagrange point
    ! is midway in between the stars.

    r2 = HALF * R

    do i = 1, max_iters

       r2_old = r2

       r2 = r2_old - L1(mass_1, mass_2, r2_old, R) / dL1dr(mass_1, mass_2, r2_old, R)

       if (abs( (r2 - r2_old) / r2_old ) < tolerance) then
          exit
       endif

    enddo

    if (i > max_iters) then
       call bl_error("L1 Lagrange point root find unable to converge.")
    endif

    ! Now convert this radial distance between the two stars
    ! into a zone index. Remember that it has to be along the
    ! line joining the two stars.

    L1_idx(:) = (com_1(:) + r2 * abs(com_2(:) - com_1(:)) / R) / dx_level(:,amr_level)

  end subroutine get_lagrange_points



  function L1(M1, M2, r2, R) result(func)

    use bl_constants_module

    implicit none

    double precision :: M1, M2, r2, R
    double precision :: func

    func = M1 * r2**2 * R**5 - M2 * (R - r2)**2 * R**5 &
         - M1 * (R - r2)**2 * r2**2 * R**3 + (M1 + M2) * r2**3  * R**2 * (R - r2)**2

  end function L1



  function dL1dr(M1, M2, r2, R) result(func)

    use bl_constants_module

    implicit none

    double precision :: M1, M2, r2, R
    double precision :: func

    func = TWO * M1 * r2 * R**5 + TWO * M2 * (R - r2) * R**5 &
         + TWO * M1 * (R - r2) * r2**2 * R**3 - TWO * M1 * (R - r2) * r2 * R**3 &
         + THREE * (M1 + M2) * r2**2  * R**2 * (R - r2)**2 - TWO * (M1 + M2) * r2**3 * R**2 * (R - r2)

  end function dL1dr

end module binary_module
