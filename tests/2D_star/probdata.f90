module probdata_module

!     These determine the refinement criteria
      double precision, save :: denerr,  dengrad
      double precision, save :: velerr,  velgrad
      double precision, save :: presserr,pressgrad
      double precision, save :: temperr,tempgrad
      double precision, save :: raderr,radgrad
      integer         , save :: max_denerr_lev   ,max_dengrad_lev
      integer         , save :: max_velerr_lev   ,max_velgrad_lev
      integer         , save :: max_presserr_lev, max_pressgrad_lev
      integer         , save :: max_temperr_lev, max_tempgrad_lev
      integer         , save :: max_raderr_lev, max_radgrad_lev

!     Sod variables
      double precision, save ::  p_s, u_s, rho_s, p_a, u_a, rho_a, rhoe_s, rhoe_a, width
      double precision, save :: T_s, T_a

      logical, save :: use_Tinit

!     These help specify which specific problem
      integer        , save ::  probtype,idir

      double precision, save ::  center(3)

      
end module probdata_module
