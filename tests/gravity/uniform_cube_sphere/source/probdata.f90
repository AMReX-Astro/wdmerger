module probdata_module

  use network

  double precision, save :: dens_ambient, temp_ambient
  double precision, save :: xn_ambient(nspec)
  double precision, save :: eint_ambient, pres_ambient

  ! cube data
  double precision, save :: density, diameter
  double precision, save :: problem

  double precision, save :: ambient_temp, ambient_dens

end module probdata_module
