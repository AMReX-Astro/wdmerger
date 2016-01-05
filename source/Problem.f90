! Problem-specific Fortran routines that are designed to interact with C++

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



! Return the mass-weighted center of mass and velocity for the primary and secondary, for a given FAB.
! Since this will rely on a sum over processors, we should only add to the relevant variables
! in anticipation of a MPI reduction, and not overwrite them.
! Note that what we are doing here is to use the old center of mass location to predict the new one,
! so com_P and com_S from the probdata_module are different from com_p_x, com_p_y, etc., which are 
! going back to C++ for the full MPI reduce. We'll then update the module locations accordingly.
! The same is true for mass_S and mass_P versus m_s and m_p.

subroutine wdcom(rho,  r_lo, r_hi, &
                 xmom, px_lo, px_hi, &
                 ymom, py_lo, py_hi, &
                 zmom, pz_lo, pz_hi, &
                 vol,  vo_lo, vo_hi, &
                 lo, hi, dx, time, &
                 com_p_x, com_p_y, com_p_z, & 
                 com_s_x, com_s_y, com_s_z, &
                 vel_p_x, vel_p_y, vel_p_z, &
                 vel_s_x, vel_s_y, vel_s_z, &
                 m_p, m_s) bind(C)

  use bl_constants_module, only: HALF, ZERO, ONE
  use prob_params_module, only: problo
  use probdata_module, only: mass_P, mass_S, com_P, com_S, single_star

  implicit none

  integer         , intent(in   ) :: r_lo(3), r_hi(3)
  integer         , intent(in   ) :: px_lo(3), px_hi(3)
  integer         , intent(in   ) :: py_lo(3), py_hi(3)
  integer         , intent(in   ) :: pz_lo(3), pz_hi(3)
  integer         , intent(in   ) :: vo_lo(3), vo_hi(3)
  
  double precision, intent(in   ) :: rho(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3))
  double precision, intent(in   ) :: xmom(px_lo(1):px_hi(1),px_lo(2):px_hi(2),px_lo(3):px_hi(3))
  double precision, intent(in   ) :: ymom(py_lo(1):py_hi(1),py_lo(2):py_hi(2),py_lo(3):py_hi(3))
  double precision, intent(in   ) :: zmom(pz_lo(1):pz_hi(1),pz_lo(2):pz_hi(2),pz_lo(3):pz_hi(3))  
  double precision, intent(in   ) :: vol(vo_lo(1):vo_hi(1),vo_lo(2):vo_hi(2),vo_lo(3):vo_hi(3))
  
  integer         , intent(in   ) :: lo(3), hi(3)
  double precision, intent(in   ) :: dx(3), time
  double precision, intent(inout) :: com_p_x, com_p_y, com_p_z
  double precision, intent(inout) :: com_s_x, com_s_y, com_s_z
  double precision, intent(inout) :: vel_p_x, vel_p_y, vel_p_z
  double precision, intent(inout) :: vel_s_x, vel_s_y, vel_s_z
  double precision, intent(inout) :: m_p, m_s

  integer          :: i, j, k
  double precision :: r(3), r_p, r_s, grav_force_P, grav_force_S
  double precision :: dm

  ! If the stars have merged, then the Roche radius of the secondary will be zero (or effectively close).
  ! In such a situation, we want to terminate this calculation gracefully. 
  ! The simplest fix is simply by just returning zero for everything;

  ! Now, add to the COM locations and velocities depending on whether we're closer
  ! to the primary or the secondary. Note that in this routine we actually are 
  ! summing mass-weighted quantities for the COM and the velocity; 
  ! we will account for this at the end of the calculation in 
  ! sum_integrated_quantities() by dividing by the mass.

  ! We determine whether a zone is in the region corresponding to one star 
  ! or the other depending on which star is exerting a stronger gravitational 
  ! force on that zone.

  do k = lo(3), hi(3)
     r(3) = problo(3) + (dble(k) + HALF) * dx(3)
     do j = lo(2), hi(2)
        r(2) = problo(2) + (dble(j) + HALF) * dx(2)
        do i = lo(1), hi(1)
           r(1) = problo(1) + (dble(i) + HALF) * dx(1)

           r_P = sqrt(sum((r - com_P)**2))
           r_S = sqrt(sum((r - com_S)**2))

           dm = rho(i,j,k) * vol(i,j,k)

           grav_force_P = mass_P / r_P

           if (single_star .or. mass_S < 1d-12 * mass_P) then
              grav_force_S = ZERO
           else
              grav_force_S = mass_S / r_S
           endif

           if (grav_force_P > grav_force_S) then

              m_p = m_p + dm

              com_p_x = com_p_x + dm * r(1)
              com_p_y = com_p_y + dm * r(2)
              com_p_z = com_p_z + dm * r(3)

              vel_p_x = vel_p_x + xmom(i,j,k) * vol(i,j,k)
              vel_p_y = vel_p_y + ymom(i,j,k) * vol(i,j,k)
              vel_p_z = vel_p_z + zmom(i,j,k) * vol(i,j,k)

           else

              m_s = m_s + dm

              com_s_x = com_s_x + dm * r(1)
              com_s_y = com_s_y + dm * r(2)
              com_s_z = com_s_z + dm * r(3)

              vel_s_x = vel_s_x + xmom(i,j,k) * vol(i,j,k)
              vel_s_y = vel_s_y + ymom(i,j,k) * vol(i,j,k)
              vel_s_z = vel_s_z + zmom(i,j,k) * vol(i,j,k)

           endif

        enddo
     enddo
  enddo

end subroutine wdcom



! This function uses the known center of mass of the two white dwarfs,
! and given a density cutoff, computes the total volume of all zones
! whose density is greater or equal to that density cutoff.
! We also impose a distance requirement so that we only look 
! at zones within the Roche lobe of the white dwarf.

subroutine ca_volumeindensityboundary(rho,r_lo,r_hi,vol,v_lo,v_hi,lo,hi,dx,volp,vols,rho_cutoff) bind(C)

  use bl_constants_module
  use probdata_module, only: mass_P, mass_S, com_P, com_S, single_star
  use prob_params_module, only: problo

  implicit none
  
  integer          :: r_lo(3), r_hi(3)
  integer          :: v_lo(3), v_hi(3)
  integer          :: lo(3), hi(3)
  double precision :: volp, vols, rho_cutoff, dx(3)
  double precision :: rho(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3))
  double precision :: vol(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))

  integer          :: i, j, k
  double precision :: r(3), r_P, r_S, grav_force_P, grav_force_S

  volp = ZERO
  vols = ZERO

  do k = lo(3), hi(3)
     r(3) = problo(3) + (dble(k) + HALF) * dx(3)
     do j = lo(2), hi(2)
        r(2) = problo(2) + (dble(j) + HALF) * dx(2)
        do i = lo(1), hi(1)
           r(1) = problo(1) + (dble(i) + HALF) * dx(1)

           if (rho(i,j,k) > rho_cutoff) then

              r_P = sqrt(sum((r - com_p)**2))
              r_S = sqrt(sum((r - com_s)**2))

              grav_force_P = mass_P / r_P
              if (single_star .or. mass_S < 1.d-12 * mass_P) then
                 grav_force_S = ZERO
              else
                 grav_force_S = mass_S / r_S
              endif
              
              if (grav_force_P > grav_force_S) then
                 volp = volp + vol(i,j,k)
              else
                 vols = vols + vol(i,j,k)
              endif
             
           endif

        enddo
     enddo
  enddo

end subroutine ca_volumeindensityboundary



! Calculate the second time derivative of the quadrupole moment tensor,
! according to the formula in Equation 6.5 of Blanchet, Damour and Schafer 1990.
! It involves integrating the mass distribution and then taking the symmetric 
! trace-free part of the tensor. We can do the latter operation here since the 
! integral is a linear operator and each part of the domain contributes independently.

subroutine quadrupole_tensor_double_dot(rho, r_lo, r_hi, &
                                        xmom, px_lo, px_hi, ymom, py_lo, py_hi, zmom, pz_lo, pz_hi, &
                                        gravx, gx_lo, gx_hi, gravy, gy_lo, gy_hi, gravz, gz_lo, gz_hi, &
                                        vol, vo_lo, vo_hi, &
                                        lo, hi, dx, time, Qtt) bind(C)

  use bl_constants_module, only: ZERO, THIRD, HALF, ONE, TWO
  use prob_params_module, only: center, problo
  use probdata_module, only: inertial_rotation, inertial_velocity
  
  implicit none

  integer         , intent(in   ) :: r_lo(3), r_hi(3)
  integer         , intent(in   ) :: px_lo(3), px_hi(3)
  integer         , intent(in   ) :: py_lo(3), py_hi(3)
  integer         , intent(in   ) :: pz_lo(3), pz_hi(3)
  integer         , intent(in   ) :: gx_lo(3), gx_hi(3)
  integer         , intent(in   ) :: gy_lo(3), gy_hi(3)
  integer         , intent(in   ) :: gz_lo(3), gz_hi(3)
  integer         , intent(in   ) :: vo_lo(3), vo_hi(3)

  double precision, intent(in   ) :: rho(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3))
  double precision, intent(in   ) :: xmom(px_lo(1):px_hi(1),px_lo(2):px_hi(2),px_lo(3):px_hi(3))
  double precision, intent(in   ) :: ymom(py_lo(1):py_hi(1),py_lo(2):py_hi(2),py_lo(3):py_hi(3))
  double precision, intent(in   ) :: zmom(pz_lo(1):pz_hi(1),pz_lo(2):pz_hi(2),pz_lo(3):pz_hi(3))
  double precision, intent(in   ) :: gravx(gx_lo(1):gx_hi(1),gx_lo(2):gx_hi(2),gx_lo(3):gx_hi(3))
  double precision, intent(in   ) :: gravy(gy_lo(1):gy_hi(1),gy_lo(2):gy_hi(2),gy_lo(3):gy_hi(3))
  double precision, intent(in   ) :: gravz(gz_lo(1):gz_hi(1),gz_lo(2):gz_hi(2),gz_lo(3):gz_hi(3))
  double precision, intent(in   ) :: vol(vo_lo(1):vo_hi(1),vo_lo(2):vo_hi(2),vo_lo(3):vo_hi(3))
  integer         , intent(in   ) :: lo(3), hi(3)
  double precision, intent(in   ) :: dx(3), time
  double precision, intent(inout) :: Qtt(3,3)

  integer          :: i, j, k, l, m
  double precision :: r(3), pos(3), vel(3), g(3), rhoInv
  double precision :: dQtt(3,3)

  dQtt(:,:) = ZERO
  
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

           ! Account for rotation, if there is any. These will leave 
           ! r and vel and changed, if not.

           pos = inertial_rotation(r, time)

           ! For constructing the velocity in the inertial frame, we need to 
           ! account for the fact that we have rotated the system already, so that 
           ! the r in omega x r is actually the position in the inertial frame, and 
           ! not the usual position in the rotating frame. It has to be on physical 
           ! grounds, because for binary orbits where the stars aren't moving, that 
           ! r never changes, and so the contribution from rotation would never change.
           ! But it must, since the motion vector of the stars changes in the inertial 
           ! frame depending on where we are in the orbit.

           vel(1) = xmom(i,j,k) * rhoInv
           vel(2) = ymom(i,j,k) * rhoInv
           vel(3) = zmom(i,j,k) * rhoInv           
           
           vel = inertial_velocity(pos, vel, time)

           g(1) = gravx(i,j,k)
           g(2) = gravy(i,j,k)
           g(3) = gravz(i,j,k)

           ! We need to rotate the gravitational field to be consistent with the rotated position.

           g = inertial_rotation(g, time)

           do m = 1, 3
              do l = 1, 3
                 dQtt(l,m) = dQtt(l,m) + TWO * rho(i,j,k) * vol(i,j,k) * (vel(l) * vel(m) + pos(l) * g(m))
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

subroutine gw_strain_tensor(h_plus_1, h_cross_1, h_plus_2, h_cross_2, h_plus_3, h_cross_3, Qtt, time) bind(C)

  use bl_constants_module, only: ZERO, HALF, ONE, TWO
  use fundamental_constants_module, only: Gconst, c_light, parsec
  use probdata_module, only: gw_dist, axis_1, axis_2, axis_3
  use meth_params_module, only: rot_axis

  implicit none

  double precision, intent(inout) :: h_plus_1, h_cross_1, h_plus_2, h_cross_2, h_plus_3, h_cross_3
  double precision, intent(in   ) :: Qtt(3,3)
  double precision, intent(in   ) :: time
  
  integer          :: i, j, k, l, dir
  double precision :: h(3,3), proj(3,3,3,3), delta(3,3), n(3), r
  double precision :: dist(3)
  
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

     if (dir .eq. axis_1) then

        h_plus_1  = h_plus_1  + h(axis_2,axis_2)
        h_cross_1 = h_cross_1 + h(axis_2,axis_3)

     else if (dir .eq. axis_2) then

        h_plus_2  = h_plus_2  + h(axis_3,axis_3)
        h_cross_2 = h_cross_2 + h(axis_3,axis_1)

     else if (dir .eq. axis_3) then

        h_plus_3  = h_plus_3  + h(axis_1,axis_1)
        h_cross_3 = h_cross_3 + h(axis_1,axis_2)

     endif

  enddo

end subroutine gw_strain_tensor



subroutine update_center(time) bind(C)

  use bl_constants_module, only: ZERO  
  use probdata_module, only: bulk_velx, bulk_vely, bulk_velz, &
                             center_fracx, center_fracy, center_fracz
  use prob_params_module, only: center, problo, probhi, dim

  implicit none

  double precision, intent(in) :: time

  ! Determine the original location of the center.

  if (dim .eq. 3) then
  
     center(1) = problo(1) + center_fracx * (probhi(1) - problo(1))
     center(2) = problo(2) + center_fracy * (probhi(2) - problo(2))
     center(3) = problo(3) + center_fracz * (probhi(3) - problo(3))

  else

     center(1) = problo(1)
     center(2) = problo(2) + center_fracz * (probhi(2) - problo(2))
     center(3) = ZERO

  endif

  ! Now update using the time passed since the beginning of the simulation.

  center(1) = center(1) + bulk_velx * time
  center(2) = center(2) + bulk_vely * time
  center(3) = center(3) + bulk_velz * time

end subroutine update_center



! Updates the CASTRO rotational period.

subroutine set_period(period) bind(C)

  use meth_params_module, only: rot_period

  implicit none

  double precision :: period

  rot_period = period

end subroutine set_period



! Returns the CASTRO rotation frequency vector.

subroutine get_omega_vec(omega_in, time) bind(C)

  use rotation_module, only: get_omega

  implicit none

  double precision :: omega_in(3), time

  omega_in = get_omega(time)

end subroutine get_omega_vec
