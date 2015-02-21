! problem-specific Fortran stuff goes here

subroutine problem_checkpoint(int_dir_name, len)

  ! called by the IO processor during checkpoint

  use bl_IO_module
  use probdata_module, only: com_P, com_S, vel_P, vel_S, mass_P, mass_S

  implicit none

  integer :: len
  integer :: int_dir_name(len)
  character (len=len) :: dir

  integer :: i, un

  ! dir will be the string name of the checkpoint directory
  do i = 1, len
     dir(i:i) = char(int_dir_name(i))
  enddo

  un = unit_new()
  open (unit=un, file=trim(dir)//"/COM", status="unknown")

100 format(1x, g30.20, 1x, g30.20)

  write (un,100) mass_P, mass_S
  write (un,100) com_P(1), com_S(1)
  write (un,100) com_P(2), com_S(2)
  write (un,100) com_P(3), com_S(3)
  write (un,100) vel_P(1), vel_S(1)
  write (un,100) vel_P(2), vel_S(2)
  write (un,100) vel_P(3), vel_S(3)

  close (un)

end subroutine problem_checkpoint


subroutine problem_restart(int_dir_name, len)

  ! called by ALL processors during restart 

  use bl_IO_module
  use probdata_module, only: com_P, com_S, vel_P, vel_S, mass_P, mass_S

  implicit none

  integer :: len
  integer :: int_dir_name(len)
  character (len=len) :: dir

  integer :: i, un

  ! dir will be the string name of the checkpoint directory
  do i = 1, len
     dir(i:i) = char(int_dir_name(i))
  enddo

  un = unit_new()
  open (unit=un, file=trim(dir)//"/COM", status="old")

100 format(1x, g30.20, 1x, g30.20)

  ! here we read in, and this sets the values in the com module
  read (un,100) mass_P, mass_S
  read (un,100) com_P(1), com_S(1)
  read (un,100) com_P(2), com_S(2)
  read (un,100) com_P(3), com_S(3)
  read (un,100) vel_P(1), vel_S(1)
  read (un,100) vel_P(2), vel_S(2)
  read (un,100) vel_P(3), vel_S(3)

  close (un)

end subroutine problem_restart



! Return whether we're doing a single star simulation or not.

subroutine get_single_star(flag)

 use probdata_module, only: single_star

 implicit none

 integer :: flag

 flag = 0

 if (single_star) flag = 1

end subroutine



! Return the mass-weighted center of mass and velocity for the primary and secondary, for a given FAB.
! Since this will rely on a sum over processors, we should only add to the relevant variables
! in anticipation of a MPI reduction, and not overwrite them.
! Note that what we are doing here is to use the old center of mass location to predict the new one,
! so com_P and com_S from the probdata_module are different from com_p_x, com_p_y, etc., which are 
! going back to C++ for the full MPI reduce. We'll then update the module locations accordingly.
! The same is true for mass_S and mass_P versus m_s and m_p.

subroutine wdcom(rho,  r_l1, r_l2, r_l3, r_h1, r_h2, r_h3, &
                 xmom, x_l1, x_l2, x_l3, x_h1, x_h2, x_h3, &
                 ymom, y_l1, y_l2, y_l3, y_h1, y_h2, y_h3, &
                 zmom, z_l1, z_l2, z_l3, z_h1, z_h2, z_h3, &
                 lo, hi, dx, time, &
                 com_p_x, com_p_y, com_p_z, & 
                 com_s_x, com_s_y, com_s_z, &
                 vel_p_x, vel_p_y, vel_p_z, &
                 vel_s_x, vel_s_y, vel_s_z, &
                 m_p, m_s)

  use bl_constants_module, only: HALF, ZERO, ONE
  use prob_params_module, only: problo, center
  use probdata_module, only: mass_P, mass_S, com_P, com_S, single_star

  implicit none

  integer         , intent(in   ) :: r_l1, r_l2, r_l3, r_h1, r_h2, r_h3
  integer         , intent(in   ) :: x_l1, x_l2, x_l3, x_h1, x_h2, x_h3
  integer         , intent(in   ) :: y_l1, y_l2, y_l3, y_h1, y_h2, y_h3
  integer         , intent(in   ) :: z_l1, z_l2, z_l3, z_h1, z_h2, z_h3

  double precision, intent(in   ) :: rho(r_l1:r_h1,r_l2:r_h2,r_l3:r_h3)
  double precision, intent(in   ) :: xmom(x_l1:x_h1,x_l2:x_h2,x_l3:x_h3)
  double precision, intent(in   ) :: ymom(y_l1:y_h1,y_l2:y_h2,y_l3:y_h3)
  double precision, intent(in   ) :: zmom(z_l1:z_h1,z_l2:z_h2,z_l3:z_h3)
  integer         , intent(in   ) :: lo(3), hi(3)
  double precision, intent(in   ) :: dx(3), time
  double precision, intent(inout) :: com_p_x, com_p_y, com_p_z
  double precision, intent(inout) :: com_s_x, com_s_y, com_s_z
  double precision, intent(inout) :: vel_p_x, vel_p_y, vel_p_z
  double precision, intent(inout) :: vel_s_x, vel_s_y, vel_s_z
  double precision, intent(inout) :: m_p, m_s

  integer          :: i, j, k
  double precision :: x, y, z, r_p, r_s
  double precision :: dV, dm
  double precision :: q, roche_rad_P, roche_rad_S, wd_dist
  double precision :: vx, vy, vz, rhoInv

  ! Volume of a zone

  dV = dx(1) * dx(2) * dx(3)

  q = mass_S / mass_P

  wd_dist = ((com_P(1)-com_S(1))**2 + (com_P(2)-com_S(2))**2 + (com_P(3)-com_S(3))**2)**HALF

  ! Get normalized Roche radii, then rescale them to the actual distance between the WDs

  call get_roche_radii(q, roche_rad_S, roche_rad_P)

  print *, q, mass_S, mass_P, com_P, com_S, roche_rad_S, roche_rad_P, wd_dist

  roche_rad_P = roche_rad_P * wd_dist
  roche_rad_S = roche_rad_S * wd_dist

  ! Now, add to the COM locations and velocities depending on whether we're closer
  ! to the primary or the secondary.

  do k = lo(3), hi(3)
     z = problo(3) + (dble(k) + HALF) * dx(3)
     do j = lo(2), hi(2)
        y = problo(2) + (dble(j) + HALF) * dx(2)
        do i = lo(1), hi(1)
           x = problo(1) + (dble(i) + HALF) * dx(1)

           r_P = ( (x - com_P(1))**2 + (y - com_P(2))**2 + (z - com_P(3))**2 )**0.5
           r_S = ( (x - com_S(1))**2 + (y - com_S(2))**2 + (z - com_S(3))**2 )**0.5

           dm = rho(i,j,k) * dV
           
           rhoInv = ONE / rho(i,j,k)

           vx = xmom(i,j,k) * rhoInv
           vy = ymom(i,j,k) * rhoInv
           vz = zmom(i,j,k) * rhoInv

           if (r_P < roche_rad_P .or. single_star) then

              m_p = m_p + dm

              com_p_x = com_p_x + dm * x
              com_p_y = com_p_y + dm * y
              com_p_z = com_p_z + dm * z

              vel_p_x = vel_p_x + dm * vx
              vel_p_y = vel_p_y + dm * vy
              vel_p_z = vel_p_z + dm * vz

           else if (r_S < roche_rad_S) then

              m_s = m_s + dm

              com_s_x = com_s_x + dm * x
              com_s_y = com_s_y + dm * y
              com_s_z = com_s_z + dm * z

              vel_s_x = vel_s_x + dm * vx
              vel_s_y = vel_s_y + dm * vy
              vel_s_z = vel_s_z + dm * vz

           endif

        enddo
     enddo
  enddo

end subroutine wdcom



! This function uses the known center of mass of the two white dwarfs,
! and given a density cutoff, computes the total volume of all zones
! whose density is greater or equal to that density cutoff.
! We also impose a distance requirement so that we only look 
! at zones that are within the effective Roche radius of the white dwarf.

subroutine ca_volumeindensityboundary(rho,r_l1,r_l2,r_l3,r_h1,r_h2,r_h3,lo,hi,dx,&
                                      volp,vols,rho_cutoff)

  use bl_constants_module
  use probdata_module, only: mass_P, mass_S, com_P, com_S, single_star
  use prob_params_module, only: problo

  implicit none
  integer          :: r_l1,r_l2,r_l3,r_h1,r_h2,r_h3
  integer          :: lo(3), hi(3)
  double precision :: volp, vols, rho_cutoff, dx(3)
  double precision :: rho(r_l1:r_h1,r_l2:r_h2,r_l3:r_h3)

  integer          :: i, j, k
  double precision :: x, y, z
  double precision :: dV
  double precision :: r_P, r_S

  double precision :: roche_rad_P, roche_rad_S
  double precision :: q, wd_dist

  dV  = dx(1)*dx(2)*dx(3)
  volp = ZERO
  vols = ZERO

  q = mass_S / mass_P

  wd_dist = ((com_P(1)-com_S(1))**2 + (com_P(2)-com_S(2))**2 + (com_P(3)-com_S(3))**2)**HALF

  ! Get normalized Roche radii, then rescale them to the actual distance between the WDs

  call get_roche_radii(q, roche_rad_S, roche_rad_P)

  roche_rad_P = roche_rad_P * wd_dist
  roche_rad_S = roche_rad_S * wd_dist

  do k = lo(3), hi(3)
     z = problo(3) + (dble(k) + HALF) * dx(3)
     do j = lo(2), hi(2)
        y = problo(2) + (dble(j) + HALF) * dx(2)
        do i = lo(1), hi(1)
           x = problo(1) + (dble(i) + HALF) * dx(1)

           if (rho(i,j,k) > rho_cutoff) then

              r_P = ( (x - com_p(1))**2 + (y - com_p(2))**2 + (z - com_p(3))**2 )**0.5
              r_S = ( (x - com_s(1))**2 + (y - com_s(2))**2 + (z - com_s(3))**2 )**0.5
              
              if (r_P < roche_rad_P .or. single_star) then
                 volp = volp + dV
              elseif (r_S < roche_rad_S) then
                 vols = vols + dV
              endif
             
           endif

        enddo
     enddo
  enddo

end subroutine ca_volumeindensityboundary



! Return the locations of the stellar centers of mass

subroutine get_star_locations(P_com, S_com, P_vel, S_vel, P_mass, S_mass)

  use probdata_module, only: com_P, com_S, vel_P, vel_S, mass_P, mass_S

  implicit none

  double precision, intent (inout) :: P_com(3), S_com(3)
  double precision, intent (inout) :: P_vel(3), S_vel(3)
  double precision, intent (inout) :: P_mass, S_mass

  P_com = com_P
  S_com = com_S

  P_vel = vel_P
  S_vel = vel_S

  P_mass = mass_P
  S_mass = mass_S

end subroutine get_star_locations



! Set the locations of the stellar centers of mass

subroutine set_star_locations(P_com, S_com, P_vel, S_vel, P_mass, S_mass)

  use probdata_module, only: com_P, com_S, vel_P, vel_S, mass_P, mass_S

  implicit none

  double precision, intent (in) :: P_com(3), S_com(3)
  double precision, intent (in) :: P_vel(3), S_vel(3)
  double precision, intent (in) :: P_mass, S_mass

  com_P = P_com
  com_S = S_com

  vel_P = P_vel
  vel_S = S_vel

  mass_P = P_mass
  mass_S = S_mass

end subroutine set_star_locations



! Given the mass ratio q of two stars (assumed to be q = M_1 / M_2), 
! compute the effective Roche radii of the stars, normalized to unity, 
! using the approximate formula of Eggleton (1983).

subroutine get_roche_radii(mass_ratio, r_1, r_2)

  use bl_constants_module, only: ONE, TWO3RD, THIRD

  implicit none

  double precision, intent(in   ) :: mass_ratio
  double precision, intent(inout) :: r_1, r_2

  double precision :: q
  double precision :: c1, c2

  c1 = 0.49d0
  c2 = 0.60d0

  q = mass_ratio

  r_1 = c1 * q**(TWO3RD) / (c2 * q**(TWO3RD) + LOG(ONE + q**(THIRD)))

  q = ONE / q

  r_2 = c1 * q**(TWO3RD) / (c2 * q**(TWO3RD) + LOG(ONE + q**(THIRD)))

end subroutine get_roche_radii
