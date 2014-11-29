program make_star

implicit none

real, parameter :: pi = 3.14159265359
integer :: npts = 1024
integer :: nvars = 7
integer :: lun = 10
integer :: i
real    :: radius, density, temperature, pressure
real    :: mass, volume
real    :: C12 = 1.0, O16 = 0.0, He4 = 0.0, Fe56 = 0.0
real    :: wd_rad

open (unit = lun, file="uniform_WD")

write (lun,*) "# npts = ", 2*npts
write (lun,*) "# num of variables = ", nvars
write (lun,*) "# density"
write (lun,*) "# temperature"
write (lun,*) "# pressure"
write (lun,*) "# helium-4"
write (lun,*) "# carbon-12"
write (lun,*) "# oxygen-16"
write (lun,*) "# iron-56"

wd_rad = 9.0e9
mass = 0.6 * 1.99e33
volume = (4.0/3.0) * pi * (wd_rad)**3

do i = 1,2*npts

  radius = (1.0*i) * wd_rad / (1.0*npts)
  if ( i <= npts ) then
    density = mass / volume
  else
    density = 1.0e-6
  endif
  temperature = 1.0e7
  pressure = 1.0e15

  write (10,*) radius, density, temperature, pressure, He4, C12, O16, Fe56

enddo

end program
