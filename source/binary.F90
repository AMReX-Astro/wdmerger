module binary_module

  ! This module contains a number of routines that are designed to
  ! calculate generic properties of binary orbits.
  
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

  subroutine get_lagrange_points(mass_1, mass_2, com_1, com_2, &
                                 L1, L2, L3) bind(C)

    use bl_constants_module
    use prob_params_module, only: dx_level
    use amrinfo_module, only: amr_level

    implicit none

    double precision :: mass_1, mass_2
    double precision :: com_1(3), com_2(3)
    double precision :: L1(3), L2(3), L3(3)
    
    double precision :: r ! Distance from Lagrange point to primary
    double precision :: a ! Distance between secondary and primary

    double precision :: r1, r2

    ! Don't try to calculate the Lagrange points if the secondary
    ! is already gone.

    if (mass_2 < ZERO) return

    a = sqrt(sum((com_2 - com_1)**2))

    r1 = -sqrt(sum(com_1**2))
    r2 = sqrt(sum(com_2**2))

    ! Do a root-find over the quintic equation for L1. 

    r = ZERO

    call lagrange_iterate(r, mass_1, mass_2, r1, r2, a, r_min = r1, r_max = r2)
    
    ! Now turn this radial distance into a grid coordinate.

    L1 = r * (com_2 - com_1) / a

    ! Repeat for L2 Lagrange point.

    r = r2 + HALF * a

    call lagrange_iterate(r, mass_1, mass_2, r1, r2, a, r_min = r2)
    
    L2 = r * (com_2 - com_1) / a

    ! Repeat for L3 Lagrange point.

    r = r1 - HALF * a

    call lagrange_iterate(r, mass_1, mass_2, r1, r2, a, r_max = r1)
    
    L3 = r * (com_2 - com_1) / a

  end subroutine get_lagrange_points



  ! Iterate over the force balance equation to find a Lagrange point.
  ! r_min and r_max set the domain over which we want to find the answer,
  ! since in general the equation has multiple roots.
  ! We assume that r comes in with a valid starting guess value.
  
  subroutine lagrange_iterate(r, mass_1, mass_2, r1, r2, a, r_min, r_max)

    use bl_constants_module, only: HALF
    
    implicit none

    ! Input variables
    
    double precision :: r
    double precision :: mass_1, mass_2, r1, r2, a
    double precision, optional :: r_min, r_max
    
    ! Root-find parameters
    
    double precision :: tolerance = 1.0d-8
    integer          :: max_iters = 100    

    ! Local variables
    
    double precision :: r_old    
    integer :: i

    do i = 1, max_iters

       r_old = r

       r = r_old - fL(mass_1, mass_2, r1, r2, r_old, a) / fdLdr(mass_1, mass_2, r1, r2, r_old, a)

       ! Force this to stay > r_max and < r_max.

       if (present(r_max)) then
          if (r > r_max) then
             r = r_old + HALF * (r_max - r_old)
          endif
       endif

       if (present(r_min)) then
          if (r < r_min) then
             r = r_old + HALF * (r_min - r_old)
          endif
       endif

       ! Check to see if we've convernged.

       if (abs( (r - r_old) / r_old ) < tolerance) then
          exit
       endif

    enddo

    if (i > max_iters) then
       call bl_error("Lagrange point root find unable to converge.")
    endif    
    
  end subroutine lagrange_iterate
  


  function fL(M1, M2, r1, r2, r, a)

    use bl_constants_module

    implicit none

    double precision :: M1, M2, r1, r2, r, a
    double precision :: fL

    double precision :: g1, g2, c

    g1 = gforce(M1, r - r1)
    g2 = gforce(M2, r - r2)

    c  = cforce(M1 + M2, a, r)
    
    fL = g1 + g2 + c

  end function fL



  function fdLdr(M1, M2, r1, r2, r, a)

    use bl_constants_module

    implicit none

    double precision :: M1, M2, r1, r2, r, a
    double precision :: fdLdr

    double precision :: dg1, dg2, dc
    
    dg1 = dgforcedr(M1, r - r1)
    dg2 = dgforcedr(M2, r - r2)

    dc  = dcforcedr(M1 + M2, a, r)
    
    fdLdr = dg1 + dg2 + dc

  end function fdLdr



  function gforce(M, r)

    use bl_constants_module, only: ONE    
    use fundamental_constants_module, only: Gconst

    double precision :: M, r
    
    double precision :: gforce

    gforce = -Gconst * M / r**2 * sign(ONE,r)

  end function gforce


  function dgforcedr(M, r)

    use bl_constants_module, only: ONE, TWO
    use fundamental_constants_module, only: Gconst

    double precision :: M, r
    
    double precision :: dgforcedr

    dgforcedr = TWO * Gconst * M / r**3 * sign(ONE,r)

  end function dgforcedr


  
  function cforce(M, a, r)

    use bl_constants_module, only: ONE    
    use fundamental_constants_module, only: Gconst

    double precision :: M, a, r
    
    double precision :: cforce

    cforce = Gconst * M / a**3 * r

  end function cforce



  function dcforcedr(M, a, r)

    use bl_constants_module, only: ONE
    use fundamental_constants_module, only: Gconst

    double precision :: M, a, r    
    
    double precision :: dcforcedr

    dcforcedr = Gconst * M / a**3

  end function dcforcedr
  
end module binary_module
