module probdata_module

!     Sod variables
      double precision, save ::  p_s, u_s, rho_s, p_a, u_a, rho_a, rhoe_s, rhoe_a, width
      double precision, save :: T_s, T_a

      logical, save :: use_Tinit

!     These help specify which specific problem
      integer        , save ::  probtype,idir

end module probdata_module
