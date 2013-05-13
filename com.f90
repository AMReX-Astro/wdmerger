module com

  double precision,               save :: mass_p,    mass_s
  double precision, dimension(3), save :: com_loc_p, com_loc_s

  real, save :: time = 0.0

end module com

subroutine com_save(mp,ms,cxp,cxs,cyp,cys,czp,czs)

  use com

  implicit none

  double precision, intent(in) :: mp,  ms
  double precision, intent(in) :: cxp, cxs
  double precision, intent(in) :: cyp, cys
  double precision, intent(in) :: czp, czs

  mass_p       = mp
  mass_s       = ms
  com_loc_p(1) = cxp
  com_loc_s(1) = cxs
  com_loc_p(2) = cyp
  com_loc_s(2) = cys
  com_loc_p(3) = czp
  com_loc_s(3) = czs

end subroutine com_save

