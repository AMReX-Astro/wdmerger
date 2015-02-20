! problem-specific Fortran stuff goes here

subroutine problem_checkpoint(int_dir_name, len)

  ! called by the IO processor during checkpoint

  use bl_IO_module

  implicit none

  integer :: len
  integer :: int_dir_name(len)
  character (len=len) :: dir

  integer :: i, un

  ! dir will be the string name of the checkpoint directory
  do i = 1, len
     dir(i:i) = char(int_dir_name(i))
  enddo

end subroutine problem_checkpoint


subroutine problem_restart(int_dir_name, len)

  ! called by ALL processors during restart 

  use bl_IO_module

  implicit none

  integer :: len
  integer :: int_dir_name(len)
  character (len=len) :: dir

  integer :: i, un

  ! dir will be the string name of the checkpoint directory
  do i = 1, len
     dir(i:i) = char(int_dir_name(i))
  enddo

end subroutine problem_restart



! Returrn whether we're doing a single star simulation or not.

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

subroutine wdcom(rho,  r_l1, r_l2, r_l3, r_h1, r_h2, r_h3, &
                 xmom, x_l1, x_l2, x_l3, x_h1, x_h2, x_h3, &
                 ymom, y_l1, y_l2, y_l3, y_h1, y_h2, y_h3, &
                 zmom, z_l1, z_l2, z_l3, z_h1, z_h2, z_h3, &
                 lo, hi, dx, time, &
                 com_p_x, com_p_y, com_p_z, & 
                 com_s_x, com_s_y, com_s_z, &
                 vel_p_x, vel_p_y, vel_p_z, &
                 vel_s_x, vel_s_y, vel_s_z, &
                 mass_p, mass_s)

  use bl_constants_module, only: HALF, ZERO
  use prob_params_module, only: problo, center
  use probdata_module, only: get_star_locations, mass_S_initial

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
  double precision, intent(inout) :: mass_p, mass_s

  integer          :: i, j, k
  double precision :: x, y, z, loc_p(3), loc_s(3), r_p, r_s
  double precision :: dV, dm

  ! First, get the predicted locations of the primary and secondary.

  call get_star_locations(time, loc_p, loc_s)

  ! Volume of a zone

  dV = dx(1) * dx(2) * dx(3)

  ! Now, add to the COM locations and velocities depending on whether we're closer
  ! to the primary or the secondary.

  do k = lo(3), hi(3)
     z = problo(3) + (dble(k) + HALF) * dx(3)
     do j = lo(2), hi(2)
        y = problo(2) + (dble(j) + HALF) * dx(2)
        do i = lo(1), hi(1)
           x = problo(1) + (dble(i) + HALF) * dx(1)

           r_P = ( (x - loc_p(1))**2 + (y - loc_p(2))**2 + (z - loc_p(3))**2 )**0.5
           r_S = ( (x - loc_s(1))**2 + (y - loc_s(2))**2 + (z - loc_s(3))**2 )**0.5

           dm = rho(i,j,k) * dV

           if (r_P < r_S .or. mass_S_initial < ZERO) then

              mass_p = mass_p + dm

              com_p_x = com_p_x + dm * x
              com_p_y = com_p_y + dm * y
              com_p_z = com_p_z + dm * z

              vel_p_x = vel_p_x + xmom(i,j,k) * dV
              vel_p_y = vel_p_y + ymom(i,j,k) * dV
              vel_p_z = vel_p_z + zmom(i,j,k) * dV              

           else if (r_S < r_P) then

              mass_s = mass_s + dm

              com_s_x = com_s_x + dm * x
              com_s_y = com_s_y + dm * y
              com_s_z = com_s_z + dm * z

              vel_s_x = vel_s_x + xmom(i,j,k) * dV
              vel_s_y = vel_s_y + ymom(i,j,k) * dV
              vel_s_z = vel_s_z + zmom(i,j,k) * dV

           endif

        enddo
     enddo
  enddo

end subroutine wdcom



! This function uses the known center of mass of the two white dwarfs,
! and given a density cutoff, computes the total volume of all zones
! whose density is greater or equal to that density cutoff.
! We also impose a distance requirement so that we only look 
! at zones that are within twice the original radius of the white dwarf.

subroutine ca_volumeindensityboundary(rho,r_l1,r_l2,r_l3,r_h1,r_h2,r_h3,lo,hi,dx,&
                                      com_p,com_s,volp,vols,rho_cutoff)

  use bl_constants_module
  use probdata_module, only: radius_P_initial, radius_S_initial
  use prob_params_module, only: problo

  implicit none
  integer          :: r_l1,r_l2,r_l3,r_h1,r_h2,r_h3
  integer          :: lo(3), hi(3)
  double precision :: volp, vols, rho_cutoff, dx(3)
  double precision :: rho(r_l1:r_h1,r_l2:r_h2,r_l3:r_h3)
  double precision :: com_p(3), com_s(3)

  integer          :: i, j, k
  double precision :: x, y, z
  double precision :: dV
  double precision :: r_P, r_S

  dV  = dx(1)*dx(2)*dx(3)
  volp = ZERO
  vols = ZERO

  do k = lo(3), hi(3)
     z = problo(3) + (dble(k) + HALF) * dx(3)
     do j = lo(2), hi(2)
        y = problo(2) + (dble(j) + HALF) * dx(2)
        do i = lo(1), hi(1)
           x = problo(1) + (dble(i) + HALF) * dx(1)

           if (rho(i,j,k) > rho_cutoff) then

              r_P = ( (x - com_p(1))**2 + (y - com_p(2))**2 + (z - com_p(3))**2 )**0.5
              r_S = ( (x - com_s(1))**2 + (y - com_s(2))**2 + (z - com_s(3))**2 )**0.5
              
              if (r_P < TWO * radius_P_initial) then
                 volp = volp + dV
              elseif (r_S < TWO * radius_S_initial) then
                 vols = vols + dV
              endif
             
           endif

        enddo
     enddo
  enddo

end subroutine ca_volumeindensityboundary
