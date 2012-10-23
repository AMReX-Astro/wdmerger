module com

  double precision, save :: mass_left, mass_right
  double precision, save :: com_xloc_l, com_xloc_r
  double precision, save :: com_yloc_l, com_yloc_r
  double precision, save :: com_zloc_l, com_zloc_r

  real, save :: time = 0.0

  integer, save :: endOfStep = 0

end module com

subroutine com_save(ml,mr,cxl,cxr,cyl,cyr,czl,czr)

  use com

  implicit none

  double precision :: ml, mr
  double precision :: cxl, cxr
  double precision :: cyl, cyr
  double precision :: czl, czr

  mass_left = ml
  mass_right = mr
  com_xloc_l = cxl
  com_xloc_r = cxr
  com_yloc_l = cyl
  com_yloc_r = cyr
  com_zloc_l = czl
  com_zloc_r = czr

end subroutine com_save



subroutine set_eos(eos)

  use com

  implicit none

  integer :: eos

  endOfStep = eos

end subroutine set_eos



subroutine get_eos(eos)

  use com

  implicit none

  integer :: eos

  eos = endOfStep

end subroutine get_eos
