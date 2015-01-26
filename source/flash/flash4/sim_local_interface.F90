module sim_local_interface

  interface interpolate
     double precision function interpolate(r, npts_model, model_r, model_var)
     implicit none
     double precision, intent(in) :: r
     integer,          intent(in) :: npts_model
     double precision, intent(in) :: model_r(npts_model), model_var(npts_model)
     end function
  end interface

end module

