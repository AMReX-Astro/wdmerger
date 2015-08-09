! problem-specific Fortran stuff goes here

subroutine problem_checkpoint(int_dir_name, len)

  ! called by the IO processor during checkpoint

  use bl_IO_module
  use probdata_module, only: com_P, com_S, vel_P, vel_S, mass_P, mass_S
  use prob_params_module, only: center

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
200 format(1x, g30.20, 1x, g30.20, 1x, g30.20)

  write (un,200) center(1), center(2), center(3)
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
  use prob_params_module, only: center

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
200 format(1x, g30.20, 1x, g30.20, 1x, g30.20)

  read (un,200) center(1), center(2), center(3)
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
  use prob_params_module, only: problo
  use probdata_module, only: mass_P, mass_S, com_P, com_S, single_star, roche_rad_P, roche_rad_S

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

  ! If the stars have merged, then the Roche radius of the secondary will be zero (or effectively close).
  ! In such a situation, we want to terminate this calculation gracefully. The simplest fix is simply by just returning zero for 
  ! everything;
  ! there is 

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

           r_P = ( (x - com_P(1))**2 + (y - com_P(2))**2 + (z - com_P(3))**2 )**0.5
           r_S = ( (x - com_S(1))**2 + (y - com_S(2))**2 + (z - com_S(3))**2 )**0.5

           dm = rho(i,j,k) * dV
           
           if (r_P < roche_rad_P .or. single_star) then

              m_p = m_p + dm

              com_p_x = com_p_x + dm * x
              com_p_y = com_p_y + dm * y
              com_p_z = com_p_z + dm * z

              vel_p_x = vel_p_x + xmom(i,j,k) * dV
              vel_p_y = vel_p_y + ymom(i,j,k) * dV
              vel_p_z = vel_p_z + zmom(i,j,k) * dV

           else if (r_S < roche_rad_S) then

              m_s = m_s + dm

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
! at zones that are within the effective Roche radius of the white dwarf.

subroutine ca_volumeindensityboundary(rho,r_l1,r_l2,r_l3,r_h1,r_h2,r_h3,lo,hi,dx,&
                                      volp,vols,rho_cutoff)

  use bl_constants_module
  use probdata_module, only: mass_P, mass_S, com_P, com_S, roche_rad_P, roche_rad_S, single_star
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
  use probdata_module, only: com_P, com_S, vel_P, vel_S, mass_P, mass_S, roche_rad_P, roche_rad_S, single_star
  use prob_params_module, only: center
  use bl_constants_module, only: TENTH, ZERO

  implicit none

  double precision, intent (in) :: P_com(3), S_com(3)
  double precision, intent (in) :: P_vel(3), S_vel(3)
  double precision, intent (in) :: P_mass, S_mass

  double precision :: r

  com_P = P_com
  vel_P = P_vel
  mass_P = P_mass

  r = ZERO

  if (.not. single_star) then

     com_S = S_com
     vel_S = S_vel
     mass_S = S_mass

     r = sum((com_P-com_S)**2)**(0.5)

     call get_roche_radii(mass_S/mass_P, roche_rad_S, roche_rad_P, r)

     ! Beyond a certain point, it doesn't make sense to track the stars separately
     ! anymore. We'll set the secondary to a fixed constant and keep it there
     ! if its Roche radius becomes smaller than 10% of the primary's.

     if (roche_rad_S .lt. TENTH * roche_rad_P) then
        com_S = center
        vel_S = ZERO
        mass_S = ZERO
        roche_rad_S = ZERO
        single_star = .true.
     endif

  endif

end subroutine set_star_locations



! Calculate the second time derivative of the quadrupole moment tensor,
! according to the formula in Equation 6.5 of Blanchet, Damour and Schafer 1990.
! It involves integrating the mass distribution and then taking the symmetric 
! trace-free part of the tensor. We can do the latter operation here since the 
! integral is a linear operator and each part of the domain contributes independently.

subroutine quadrupole_tensor_double_dot(rho,  r_l1, r_l2, r_l3, r_h1, r_h2, r_h3, &
                                        xmom, x_l1, x_l2, x_l3, x_h1, x_h2, x_h3, &
                                        ymom, y_l1, y_l2, y_l3, y_h1, y_h2, y_h3, &
                                        zmom, z_l1, z_l2, z_l3, z_h1, z_h2, z_h3, &
                                        gx,  gx_l1, gx_l2, gx_l3, gx_h1, gx_h2, gx_h3, &
                                        gy,  gy_l1, gy_l2, gy_l3, gy_h1, gy_h2, gy_h3, &
                                        gz,  gz_l1, gz_l2, gz_l3, gz_h1, gz_h2, gz_h3, &
                                        lo, hi, dx, time, &
                                        Qtt)

  use bl_constants_module, only: ZERO, THIRD, HALF, ONE, TWO
  use prob_params_module, only: center, problo
  use meth_params_module, only: do_rotation, rot_period
  use rot_sources_module, only: get_omega, cross_product
  
  implicit none

  integer         , intent(in   ) :: r_l1, r_l2, r_l3, r_h1, r_h2, r_h3
  integer         , intent(in   ) :: x_l1, x_l2, x_l3, x_h1, x_h2, x_h3
  integer         , intent(in   ) :: y_l1, y_l2, y_l3, y_h1, y_h2, y_h3
  integer         , intent(in   ) :: z_l1, z_l2, z_l3, z_h1, z_h2, z_h3
  integer         , intent(in   ) :: gx_l1, gx_l2, gx_l3, gx_h1, gx_h2, gx_h3
  integer         , intent(in   ) :: gy_l1, gy_l2, gy_l3, gy_h1, gy_h2, gy_h3
  integer         , intent(in   ) :: gz_l1, gz_l2, gz_l3, gz_h1, gz_h2, gz_h3

  double precision, intent(in   ) :: rho(r_l1:r_h1,r_l2:r_h2,r_l3:r_h3)
  double precision, intent(in   ) :: xmom(x_l1:x_h1,x_l2:x_h2,x_l3:x_h3)
  double precision, intent(in   ) :: ymom(y_l1:y_h1,y_l2:y_h2,y_l3:y_h3)
  double precision, intent(in   ) :: zmom(z_l1:z_h1,z_l2:z_h2,z_l3:z_h3)
  double precision, intent(in   ) :: gx(gx_l1:gx_h1,gx_l2:gx_h2,gx_l3:gx_h3)
  double precision, intent(in   ) :: gy(gy_l1:gy_h1,gy_l2:gy_h2,gy_l3:gy_h3)
  double precision, intent(in   ) :: gz(gz_l1:gz_h1,gz_l2:gz_h2,gz_l3:gz_h3)
  integer         , intent(in   ) :: lo(3), hi(3)
  double precision, intent(in   ) :: dx(3), time
  double precision, intent(inout) :: Qtt(3,3)

  integer          :: i, j, k, l, m, dir
  double precision :: dV
  double precision :: r(3), pos(3), vel(3), grav(3), rhoInv
  double precision :: dQtt(3,3)
  double precision :: omega(3)
  double precision :: theta(3), rot_matrix(3,3)

  dV = dx(1) * dx(2) * dx(3)
  dQtt(:,:) = ZERO

  if (do_rotation .eq. 1) then
     omega = get_omega()
  else
     omega = ZERO
  endif

  ! If we are in a rotating reference frame, then we want
  ! to rotate this by an amount corresponding to the time
  ! that has passed since the beginning of the simulation.
  ! This way it will correspond to an inertial observer
  ! at large distances. There is some clipping due to material
  ! that leaves or enters the domain boundaries, but this should
  ! average out over time, and will be negligible in magnitude.
        
  theta = omega * time

  ! This is the 3D rotation matrix for converting between reference frames.
  ! It is the composition of rotations along the x, y, and z axes. Therefore 
  ! it allows for the case where we are rotating about multiple axes. Normally 
  ! we use the right-hand convention for constructing the usual rotation matrix, 
  ! but this is the transpose of that rotation matrix to account for the fact 
  ! that we are rotating *back* to the inertial frame, rather than from the 
  ! inertial frame to the rotating frame.

  rot_matrix(1,1) =  cos(theta(2)) * cos(theta(3))
  rot_matrix(1,2) = -cos(theta(2)) * sin(theta(3))
  rot_matrix(1,3) =  sin(theta(2))
  rot_matrix(2,1) =  cos(theta(1)) * sin(theta(3)) + sin(theta(1)) * sin(theta(2)) * cos(theta(3))
  rot_matrix(2,2) =  cos(theta(1)) * cos(theta(3)) - sin(theta(1)) * sin(theta(2)) * sin(theta(3))
  rot_matrix(2,3) = -sin(theta(1)) * cos(theta(2))
  rot_matrix(3,1) =  sin(theta(1)) * sin(theta(3)) - cos(theta(1)) * sin(theta(2)) * cos(theta(3))
  rot_matrix(3,2) =  sin(theta(1)) * cos(theta(3)) + cos(theta(1)) * sin(theta(2)) * sin(theta(3))
  rot_matrix(3,3) =  cos(theta(1)) * cos(theta(2))

  do k = lo(3), hi(3)
     r(3) = problo(3) + (k + HALF) * dx(3) - center(3)
     do j = lo(2), hi(2)
        r(2) = problo(2) + (j + HALF) * dx(2) - center(2)
        do i = lo(1), hi(1)
           r(1) = problo(1) + (i + HALF) * dx(1) - center(1)

           if (rho(i,j,k) > ZERO) then
              rhoInv = ONE / rho(i,j,k)
           else
              rhoInv = ZERO
           endif

           vel(1) = xmom(i,j,k) * rhoInv
           vel(2) = ymom(i,j,k) * rhoInv
           vel(3) = zmom(i,j,k) * rhoInv

           ! Account for rotation, if there is any. These will leave 
           ! r and vel and changed, if not.

           ! Remember that this is going from the rotating frame back to the inertial frame.

           pos = matmul(rot_matrix, r)

           ! For constructing the velocity in the inertial frame, we need to 
           ! account for the fact that we have rotated the system already, so that 
           ! the r in omega x r is actually the position in the inertial frame, and 
           ! not the usual position in the rotating frame. It has to be on physical 
           ! grounds, because for binary orbits where the stars aren't moving, that 
           ! r never changes, and so the contribution from rotation would never change.
           ! But it must, since the motion vector of the stars changes in the inertial 
           ! frame depending on where we are in the orbit.

           vel = vel + cross_product(omega, pos)

           grav(1) = gx(i,j,k)
           grav(2) = gy(i,j,k)
           grav(3) = gz(i,j,k)

           ! We need to rotate the gravitational field to be consistent with the rotated position.

           grav = matmul(rot_matrix, grav)

           do m = 1, 3
              do l = 1, 3
                 dQtt(l,m) = dQtt(l,m) + TWO * rho(i,j,k) * dV * (vel(l) * vel(m) + pos(l) * grav(m))
              enddo
           enddo

        enddo
     enddo
  enddo

  ! Now take the symmetric trace-free part of the quadrupole moment.
  ! The operator is defined in Equation 6.7 of Blanchet et al. (1990):
  ! STF(A^{ij}) = 1/2 A^{ij} + 1/2 A^{ji} - 1/3 delta^{ij} sum_{k} A^{kk}.

  do l = 1, 3
     do m = 1, 3

        Qtt(l,m) = Qtt(l,m) + HALF * dQtt(l,m) + HALF * dQtt(m,l)
        Qtt(l,l) = Qtt(l,l) - THIRD * dQtt(m,m)

     enddo
  enddo

end subroutine quadrupole_tensor_double_dot



! Given the above quadrupole tensor, calculate the strain tensor.

subroutine gw_strain_tensor(h_plus_rot, h_cross_rot, h_plus_star, h_cross_star, h_plus_motion, h_cross_motion, Qtt, time)

  use bl_constants_module, only: ZERO, HALF, ONE, TWO
  use fundamental_constants_module, only: Gconst, c_light, parsec
  use probdata_module, only: gw_dist, star_axis, initial_motion_dir
  use meth_params_module, only: rot_axis

  implicit none

  double precision, intent(inout) :: h_plus_rot, h_cross_rot, h_plus_star, h_cross_star, h_plus_motion, h_cross_motion
  double precision, intent(in   ) :: Qtt(3,3)
  double precision, intent(in   ) :: time
  
  integer :: i, j, k, l, m, dir
  double precision :: h(3,3), proj(3,3,3,3), delta(3,3), n(3), r
  double precision :: dist(3)
  double precision :: omega(3)
  
  ! Standard Kronecker delta.

  delta(:,:) = ZERO

  do i = 1, 3
     delta(i,i) = ONE
  enddo

  ! Unit vector for the wave is simply the distance 
  ! vector to the observer normalized by the total distance.
  ! We are going to repeat this process by looking along 
  ! all three coordinate axes.

  do dir = 1, 3

     dist(:) = ZERO
     dist(dir) = gw_dist

     r = sqrt(sum(dist**2))

     n(:) = dist(:) / r

     h = ZERO

     ! Projection operator onto the unit vector n.

     do l = 1, 3
        do k = 1, 3
           do j = 1, 3
              do i = 1, 3
                 proj(i,j,k,l) = (delta(i,k) - n(i) * n(k)) * (delta(j,l) - n(j) * n(l)) &
                               - HALF * (delta(i,j) - n(i) * n(j)) * (delta(k,l) - n(k) * n(l))
              enddo
           enddo
        enddo
     enddo

     ! Now we can calculate the strain tensor.

     do l = 1, 3
        do k = 1, 3
           do j = 1, 3
              do i = 1, 3
                 h(i,j) = h(i,j) + proj(i,j,k,l) * Qtt(k,l)
              enddo
           enddo
        enddo
     enddo

     ! Finally multiply by the coefficients.

     r = r * parsec * 1d3 ! Convert from kpc to cm.

     h(:,:) = h(:,:) * TWO * Gconst / (c_light**4 * r)

     ! If rot_axis == 3, then h_+ = h_{11} = -h_{22} and h_x = h_{12} = h_{21}.
     ! Analogous statements hold along the other axes.

     ! We are adding here so that this calculation makes sense on multiple levels.

     if (dir .eq. rot_axis) then

        h_plus_rot = h_plus_rot + h(star_axis,star_axis)
        h_cross_rot = h_cross_rot + h(star_axis,initial_motion_dir)

     else if (dir .eq. star_axis) then

        h_plus_star = h_plus_star + h(initial_motion_dir,initial_motion_dir)
        h_cross_star = h_cross_star + h(initial_motion_dir,rot_axis)

     else if (dir .eq. initial_motion_dir) then

        h_plus_motion = h_plus_motion + h(rot_axis,rot_axis)
        h_cross_motion = h_cross_motion + h(rot_axis,star_axis)

     endif

  enddo

end subroutine gw_strain_tensor



! Given the mass ratio q of two stars (assumed to be q = M_1 / M_2), 
! compute the effective Roche radii of the stars, normalized to unity, 
! using the approximate formula of Eggleton (1983). We then 
! scale them appropriately using the current orbital distance.

subroutine get_roche_radii(mass_ratio, r_1, r_2, a)

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



subroutine setup_scf_relaxation(dx, problo, probhi)

    use bl_constants_module, only: HALF, ONE, TWO
    use prob_params_module, only: center
    use eos_module
    use probdata_module

    implicit none

    double precision :: dx(3), problo(3), probhi(3)

    double precision :: x, y, z
    integer          :: n
    integer          :: ncell(3)

    double precision :: pos_l(3)
    type (eos_t) :: eos_state

    ! Following Swesty, Wang and Calder (2000), and Motl et al. (2002), 
    ! we need to fix three points to uniquely determine an equilibrium 
    ! configuration for two unequal mass stars. We can do this by 
    ! specifying d_A, the distance from the center of mass to the 
    ! inner point of the primary; d_B, the distance from the center 
    ! of mass to the inner point of the secondary; and d_C, the 
    ! distance from the center of mass to the outer point of the secondary. 
    ! We give these widths in physical units, so in general these will 
    ! be locations on the grid that are not coincident with a corner. 
    ! If a location is at point (x, y, z), then we find the eight 
    ! zone centers that surround this point, and do a tri-linear
    ! reconstruction to estimate the relative weight of each of the eight 
    ! zone centers in determining the value at that point.

    ncell = NINT( (probhi - problo) / dx )

    do n = 1, 3
       d_vector(:,n) = center
       if (n .eq. 1) then
          d_vector(star_axis,n) = d_vector(star_axis,n) - d_A
       else if (n .eq. 2) then
          d_vector(star_axis,n) = d_vector(star_axis,n) + d_B
       else if (n .eq. 3) then
          d_vector(star_axis,n) = d_vector(star_axis,n) + d_C
       endif
    enddo

    ! Locate the zone centers that bracket each point at the 
    ! lower left corner. Note that the INT function rounds down,
    ! which is what we want here since we want the lower left.

    do n = 1, 3
       rloc(:,n) = INT( (d_vector(:,n) + dx(:)/TWO) / dx(:)) + ncell / 2 - 1
    enddo

    ! Obtain the location of these points relative to the cube surrounding them.
    ! The lower left corner is at (0,0,0) and the upper right corner is at (1,1,1).

    do n = 1, 3
       pos_l = (rloc(:,n) - ncell / 2 + HALF) * dx
       rpos(:,n) = (d_vector(:,n) - pos_l) / dx(:)
    enddo

    ! Determine the tri-linear coefficients

    do n = 1, 3
       x = rpos(1,n)
       y = rpos(2,n)
       z = rpos(3,n)
       c(0,0,0,n) = (ONE - x) * (ONE - y) * (ONE - z)
       c(1,0,0,n) = x         * (ONE - y) * (ONE - z)
       c(0,1,0,n) = (ONE - x) * y         * (ONE - z)
       c(1,1,0,n) = x         * y         * (ONE - z)
       c(0,0,1,n) = (ONE - x) * (ONE - y) * z
       c(1,0,1,n) = x         * (ONE - y) * z
       c(0,1,1,n) = (ONE - x) * y         * z
       c(1,1,1,n) = x         * y         * z
    enddo

    ! Convert the maximum densities into maximum enthalpies.

    eos_state % T   = stellar_temp
    eos_state % rho = central_density_P
    eos_state % xn  = stellar_comp

    call eos(eos_input_rt, eos_state)

    h_max_P = eos_state % h

    eos_state % T   = stellar_temp
    eos_state % rho = central_density_S
    eos_state % xn  = stellar_comp

    call eos(eos_input_rt, eos_state)
    
    h_max_S = eos_state % h

    ! Determine the lowest possible enthalpy that can be 
    ! obtained by the EOS; this is useful for protecting
    ! against trying to compute a corresponding temperature
    ! in zones where the enthalpy is just too low for convergence.

    eos_state % T   = stellar_temp
    eos_state % rho = ambient_density
    eos_state % xn  = stellar_comp

    do while (eos_state % rho < 1.d11)

      eos_state % rho = eos_state % rho * 1.1

      call eos(eos_input_rt, eos_state)

      if (eos_state % h < enthalpy_min) enthalpy_min = eos_state % h

    enddo

end subroutine setup_scf_relaxation



subroutine get_coeff_info(loc_A,loc_B,loc_C,c_A,c_B,c_C)

  use probdata_module, only: rloc, c

  implicit none

  double precision, intent(inout) :: loc_A(3), loc_B(3), loc_C(3)
  double precision, intent(inout) :: c_A(2,2,2), c_B(2,2,2), c_C(2,2,2)
  
  loc_A(:) = rloc(:,1)
  loc_B(:) = rloc(:,2)
  loc_C(:) = rloc(:,3)

  c_A(:,:,:) = c(:,:,:,1)
  c_B(:,:,:) = c(:,:,:,2)
  c_C(:,:,:) = c(:,:,:,3)

end subroutine get_coeff_info



subroutine get_omegasq(lo,hi,domlo,domhi, &
                       state,state_l1,state_l2,state_l3,state_h1,state_h2,state_h3, &
                       phi,phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3, &
                       dx,problo,probhi,omegasq)

    use bl_constants_module, only: ONE, TWO
    use meth_params_module, only: NVAR, URHO
    use prob_params_module, only: center
    use probdata_module, only: c, d_vector, rloc

    implicit none
    
    integer :: lo(3), hi(3), domlo(3), domhi(3)
    integer :: state_l1,state_h1,state_l2,state_h2,state_l3,state_h3
    integer :: phi_l1,phi_h1,phi_l2,phi_h2,phi_l3,phi_h3
    double precision :: problo(3), probhi(3), dx(3)
    double precision :: state(state_l1:state_h1,state_l2:state_h2,state_l3:state_h3,NVAR)
    double precision :: phi(phi_l1:phi_h1,phi_l2:phi_h2,phi_l3:phi_h3)
    double precision :: omegasq

    integer :: i, j, k

    phi(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = -phi(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))    

    do k = rloc(3,2), rloc(3,2)+1
       do j = rloc(2,2), rloc(2,2)+1
          do i = rloc(1,2), rloc(1,2)+1
             if (i .ge. lo(1) .and. j .ge. lo(2) .and. k .ge. lo(3) .and. &
                 i .le. hi(1) .and. j .le. hi(2) .and. k .le. hi(3)) then
                omegasq = omegasq - c(i-rloc(1,2),j-rloc(2,2),k-rloc(3,2),2) &
                        * TWO * phi(i,j,k) / (sum(d_vector(:,3)**2)-sum(d_vector(:,2)**2))
             endif
          enddo
       enddo
    enddo

    do k = rloc(3,3), rloc(3,3) + 1
       do j = rloc(2,3), rloc(2,3) + 1
          do i = rloc(1,3), rloc(1,3) + 1
             if (i .ge. lo(1) .and. j .ge. lo(2) .and. k .ge. lo(3) .and. &
                 i .le. hi(1) .and. j .le. hi(2) .and. k .le. hi(3)) then
                omegasq = omegasq + c(i-rloc(1,3),j-rloc(2,3),k-rloc(3,3),3) &
                        * TWO * phi(i,j,k) / (sum(d_vector(:,3)**2)-sum(d_vector(:,2)**2))
             endif
          enddo
       enddo
    enddo

    phi(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = -phi(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))

end subroutine get_omegasq



subroutine set_period(period)

  use meth_params_module, only: rot_period

  implicit none

  double precision :: period

  rot_period = period

end subroutine set_period



subroutine get_bernoulli_const(lo,hi,domlo,domhi, &
                               state,state_l1,state_l2,state_l3,state_h1,state_h2,state_h3, &
                               phi,phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3, &
                               dx,problo,probhi,bernoulli_1,bernoulli_2)

    use bl_constants_module, only: HALF, ONE, TWO, M_PI
    use meth_params_module, only: NVAR, URHO, rot_period
    use prob_params_module, only: center
    use probdata_module, only: c, d_vector, rloc

    implicit none
    
    integer :: lo(3), hi(3), domlo(3), domhi(3)
    integer :: state_l1,state_h1,state_l2,state_h2,state_l3,state_h3
    integer :: phi_l1,phi_h1,phi_l2,phi_h2,phi_l3,phi_h3
    double precision :: problo(3), probhi(3), dx(3)
    double precision :: state(state_l1:state_h1,state_l2:state_h2,state_l3:state_h3,NVAR)
    double precision :: phi(phi_l1:phi_h1,phi_l2:phi_h2,phi_l3:phi_h3)
    double precision :: bernoulli_1, bernoulli_2

    integer :: i, j, k
    double precision :: omega

    omega = TWO * M_PI / rot_period

    phi(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = -phi(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))

    do k = rloc(3,1), rloc(3,1)+1
       do j = rloc(2,1), rloc(2,1)+1
          do i = rloc(1,1), rloc(1,1)+1
             if (i .ge. lo(1) .and. j .ge. lo(2) .and. k .ge. lo(3) .and. &
                 i .le. hi(1) .and. j .le. hi(2) .and. k .le. hi(3)) then
                bernoulli_1 = bernoulli_1 + c(i-rloc(1,1),j-rloc(2,1),k-rloc(3,1),1) &
                            * (phi(i,j,k) - HALF * omega**2 * sum(d_vector(:,1)**2))
             endif
          enddo
       enddo
    enddo

    do k = rloc(3,2), rloc(3,2) + 1
       do j = rloc(2,2), rloc(2,2) + 1
          do i = rloc(1,2), rloc(1,2) + 1
             if (i .ge. lo(1) .and. j .ge. lo(2) .and. k .ge. lo(3) .and. &
                 i .le. hi(1) .and. j .le. hi(2) .and. k .le. hi(3)) then
                bernoulli_2 = bernoulli_2 + c(i-rloc(1,2),j-rloc(2,2),k-rloc(3,2),2) &
                            * (phi(i,j,k) - HALF * omega**2 * sum(d_vector(:,2)**2))
             endif
          enddo
       enddo
    enddo

    phi(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = -phi(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))

end subroutine get_bernoulli_const



subroutine construct_enthalpy(lo,hi,domlo,domhi, &
                              state,state_l1,state_l2,state_l3,state_h1,state_h2,state_h3, &
                              phi,phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3, &
                              enthalpy,h_l1,h_l2,h_l3,h_h1,h_h2,h_h3, &
                              dx,problo,probhi,&
                              bernoulli_1,bernoulli_2,h_max_1,h_max_2)

    use bl_constants_module, only: ZERO, HALF, ONE, TWO, M_PI
    use meth_params_module, only: NVAR, URHO, rot_period
    use prob_params_module, only: center
    use probdata_module, only: star_axis

    implicit none
    
    integer :: lo(3), hi(3), domlo(3), domhi(3)
    integer :: state_l1,state_h1,state_l2,state_h2,state_l3,state_h3
    integer :: phi_l1,phi_h1,phi_l2,phi_h2,phi_l3,phi_h3
    integer :: h_l1,h_h1,h_l2,h_h2,h_l3,h_h3
    double precision :: problo(3), probhi(3), dx(3)
    double precision :: state(state_l1:state_h1,state_l2:state_h2,state_l3:state_h3,NVAR)
    double precision :: phi(phi_l1:phi_h1,phi_l2:phi_h2,phi_l3:phi_h3)
    double precision :: enthalpy(h_l1:h_h1,h_l2:h_h2,h_l3:h_h3)
    double precision :: bernoulli_1, bernoulli_2, h_max_1,h_max_2

    integer :: i, j, k
    double precision :: r(3)
    double precision :: omega

    omega = TWO * M_PI / rot_period

    phi(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = -phi(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))

    do k = lo(3), hi(3)
       r(3) = problo(3) + (dble(k) + HALF) * dx(3) - center(3)
       do j = lo(2), hi(2)
          r(2) = problo(2) + (dble(j) + HALF) * dx(2) - center(2)
          do i = lo(1), hi(1)
             r(1) = problo(1) + (dble(i) + HALF) * dx(1) - center(1)

             if (r(star_axis) < ZERO) then
                enthalpy(i,j,k) = bernoulli_1 - phi(i,j,k) + HALF * omega**2 * sum(r**2)
                if (enthalpy(i,j,k) > h_max_1) then
                   h_max_1 = enthalpy(i,j,k)
                endif
             else
                enthalpy(i,j,k) = bernoulli_2 - phi(i,j,k) + HALF * omega**2 * sum(r**2)
                if (enthalpy(i,j,k) > h_max_2) then
                   h_max_2 = enthalpy(i,j,k)
                endif
             endif

          enddo
       enddo
    enddo

    phi(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = -phi(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))

end subroutine construct_enthalpy



subroutine update_density(lo,hi,domlo,domhi, &
                          state,state_l1,state_l2,state_l3,state_h1,state_h2,state_h3, &
                          phi,phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3, &
                          enthalpy,h_l1,h_l2,h_l3,h_h1,h_h2,h_h3, &
                          dx,problo,probhi, &
                          h_max_1,h_max_2, &
                          kin_eng, pot_eng, int_eng, &
                          left_mass, right_mass, &
                          delta_rho, l2_norm_resid, l2_norm_source)

    use bl_constants_module, only: ZERO, ONE, TWO, M_PI
    use meth_params_module, only: NVAR, URHO, UTEMP, UMX, UMY, UMZ, UEDEN, UEINT, UFS, rot_period
    use network, only: nspec
    use prob_params_module, only: center
    use probdata_module, only: star_axis, get_ambient, h_max_P, h_max_S, enthalpy_min
    use eos_module

    implicit none
    
    integer :: lo(3), hi(3), domlo(3), domhi(3)
    integer :: state_l1,state_h1,state_l2,state_h2,state_l3,state_h3
    integer :: phi_l1,phi_h1,phi_l2,phi_h2,phi_l3,phi_h3
    integer :: h_l1,h_h1,h_l2,h_h2,h_l3,h_h3
    double precision :: problo(3), probhi(3), dx(3)
    double precision :: state(state_l1:state_h1,state_l2:state_h2,state_l3:state_h3,NVAR)
    double precision :: phi(phi_l1:phi_h1,phi_l2:phi_h2,phi_l3:phi_h3)
    double precision :: enthalpy(h_l1:h_h1,h_l2:h_h2,h_l3:h_h3)
    double precision :: h_max_1,h_max_2
    double precision :: kin_eng, pot_eng, int_eng
    double precision :: left_mass, right_mass
    double precision :: delta_rho, l2_norm_resid, l2_norm_source

    integer :: i, j, k
    double precision :: r(3)
    double precision :: old_rho, drho
    double precision :: dV
    double precision :: omega
    double precision :: max_dist

    type (eos_t) :: eos_state, ambient_state

    phi(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = -phi(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))

    dV = dx(1) * dx(2) * dx(3)

    max_dist = 0.75 * max(maxval(abs(probhi-center)), maxval(abs(problo-center)))

    omega = TWO * M_PI / rot_period

    call get_ambient(ambient_state)

    do k = lo(3), hi(3)
       r(3) = problo(3) + (dble(k) + HALF) * dx(3) - center(3)
       do j = lo(2), hi(2)
          r(2) = problo(2) + (dble(j) + HALF) * dx(2) - center(2)
          do i = lo(1), hi(1)
             r(1) = problo(1) + (dble(i) + HALF) * dx(1) - center(1)

             old_rho = state(i,j,k,URHO)

             ! Rescale the enthalpy by the maximum value.

             if (r(star_axis) .lt. ZERO) then
                enthalpy(i,j,k) = h_max_P * (enthalpy(i,j,k) / h_max_1)
             else
                enthalpy(i,j,k) = h_max_S * (enthalpy(i,j,k) / h_max_2)
             endif

             ! We only want to call the EOS for zones with enthalpy > 0,
             ! but for distances far enough from the center, the rotation
             ! term can overcome the other terms and make the enthalpy 
             ! spuriously positive. So we'll only consider zones with 75%
             ! of the distance from the center.

             if (enthalpy(i,j,k) > enthalpy_min .and. sum(r**2) .lt. max_dist**2) then

                eos_state % T   = state(i,j,k,UTEMP)
                eos_state % h   = enthalpy(i,j,k)
                eos_state % xn  = state(i,j,k,UFS:UFS+nspec-1) / state(i,j,k,URHO)
                eos_state % rho = state(i,j,k,URHO) ! Initial guess for the EOS

                ! eos_state % loc = (/ i, j, k /)

                call eos(eos_input_th, eos_state)

             else

                eos_state = ambient_state

             endif

             state(i,j,k,URHO) = eos_state % rho
             state(i,j,k,UEINT) = eos_state % rho * eos_state % e
             state(i,j,k,UFS:UFS+nspec-1) = eos_state % rho * eos_state % xn

             state(i,j,k,UMX:UMZ) = ZERO

             state(i,j,k,UEDEN) = state(i,j,k,UEINT) + HALF * (state(i,j,k,UMX)**2 + &
                                  state(i,j,k,UMY)**2 + state(i,j,k,UMZ)**2) / state(i,j,k,URHO)

             ! Convergence tests and diagnostic quantities

             drho = abs( state(i,j,k,URHO) - old_rho ) / old_rho
             if (drho > delta_rho) delta_rho = drho

             l2_norm_resid = l2_norm_resid + dV * (state(i,j,k,URHO) - old_rho)**2
             l2_norm_source = l2_norm_source + dV * old_rho**2

             kin_eng = kin_eng + HALF * omega**2 * state(i,j,k,URHO) * (r(1)**2 + r(2)**2) * dV

             pot_eng = pot_eng + HALF * state(i,j,k,URHO) * phi(i,j,k) * dV

             int_eng = int_eng + eos_state % p * dV

             if (r(star_axis) < ZERO) then
               left_mass = left_mass + state(i,j,k,URHO) * dV
             else
               right_mass = right_mass + state(i,j,k,URHO) * dV
             endif

          enddo
       enddo
    enddo

    phi(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = -phi(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))

end subroutine update_density
  


  



subroutine check_convergence(kin_eng, pot_eng, int_eng, &
                             left_mass, right_mass, &
                             delta_rho, l2_norm, &
                             is_relaxed, num_iterations)

  use probdata_module, only: relax_tol
  use meth_params_module, only: rot_period
  use multifab_module
  use bl_constants_module, only: THREE
  use fundamental_constants_module, only: M_solar
  use eos_module

  implicit none

  integer :: is_relaxed
  integer :: num_iterations
  integer :: ioproc

  double precision :: kin_eng, pot_eng, int_eng
  double precision :: left_mass, right_mass
  double precision :: delta_rho, l2_norm

  double precision :: virial_error

  virial_error = abs(TWO * kin_eng + pot_eng + THREE * int_eng) / abs(pot_eng)

  if (l2_norm .lt. relax_tol) then
    is_relaxed = 1
  endif

  call bl_pd_is_ioproc(ioproc)
  if (ioproc == 1) then
    print *, ""
    print *, ""
    print *, "  Relaxation iterations completed:", num_iterations
    print *, "  Maximum change in rho:", delta_rho
    print *, "  L2 Norm of Residual (relative to old state):", l2_norm
    print *, "  Current value of rot_period:", rot_period
    print *, "  Kinetic energy:", kin_eng 
    print *, "  Potential energy:", pot_eng
    print *, "  Internal energy:", int_eng
    print *, "  Virial error:", virial_error
    print *, "  Mass (M_sun) on left side of grid:", left_mass / M_solar
    print *, "  Mass (M_sun) on right side of grid:", right_mass / M_solar
    if (is_relaxed .eq. 1) print *, "  Relaxation completed!"
    print *, ""
    print *, ""
  endif

end subroutine check_convergence



! Update the location of center using the bulk velocity. 

subroutine update_center(time)

  use probdata_module, only: bulk_velx, bulk_vely, bulk_velz, &
                             center_fracx, center_fracy, center_fracz
  use prob_params_module, only: center, problo, probhi

  implicit none

  double precision, intent(in) :: time

  ! Determine the original location of the center.

  center(1) = problo(1) + center_fracx * (probhi(1) - problo(1))
  center(2) = problo(2) + center_fracy * (probhi(2) - problo(2))
  center(3) = problo(3) + center_fracz * (probhi(3) - problo(3))

  ! Now update using the time passed since the beginning of the simulation.

  center(1) = center(1) + bulk_velx * time
  center(2) = center(2) + bulk_vely * time
  center(3) = center(3) + bulk_velz * time

end subroutine update_center
