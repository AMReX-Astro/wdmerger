! Derive momentum, given input vector of the grid momenta.

subroutine ca_derinertialmomentumx(p,p_lo,p_hi,ncomp_p, &
                                   u,u_lo,u_hi,ncomp_u, &
                                   lo,hi,domlo,domhi, &
                                   dx,xlo,time,dt,bc,level,grid_no) &
                                   bind(C,name='ca_derinertialmomentumx')

  use bl_constants_module, only: HALF
  use probdata_module, only: inertial_velocity

  implicit none

  integer          :: p_lo(3), p_hi(3), ncomp_p ! == 1
  integer          :: u_lo(3), u_hi(3), ncomp_u ! == 4
  integer          :: lo(3), hi(3), domlo(3), domhi(3)
  double precision :: p(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3),ncomp_p)
  double precision :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
  double precision :: dx(3), xlo(3), time, dt
  integer          :: bc(3,2,ncomp_u), level, grid_no

  integer          :: i, j, k
  double precision :: loc(3), vel(3), mom(3), rho

  do k = lo(3), hi(3)
     loc(3) = xlo(3) + (dble(k - lo(3)) + HALF) * dx(3)
     do j = lo(2), hi(2)
        loc(2) = xlo(2) + (dble(j - lo(2)) + HALF) * dx(2)
        do i = lo(1), hi(1)
           loc(1) = xlo(1) + (dble(i - lo(1)) + HALF) * dx(1)

           rho = u(i,j,k,1)
           vel = u(i,j,k,2:4) / rho
           mom = rho * inertial_velocity(loc, vel, time)
           p(i,j,k,1) = mom(1)

        enddo
     enddo
  enddo

end subroutine ca_derinertialmomentumx



! Derive momentum, given input vector of the grid momenta.

subroutine ca_derinertialmomentumy(p,p_lo,p_hi,ncomp_p, &
                                   u,u_lo,u_hi,ncomp_u, &
                                   lo,hi,domlo,domhi, &
                                   dx,xlo,time,dt,bc,level,grid_no) &
                                   bind(C,name='ca_derinertialmomentumy')

  use bl_constants_module, only: HALF
  use probdata_module, only: inertial_velocity

  implicit none

  integer          :: p_lo(3), p_hi(3), ncomp_p ! == 1
  integer          :: u_lo(3), u_hi(3), ncomp_u ! == 4
  integer          :: lo(3), hi(3), domlo(3), domhi(3)
  double precision :: p(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3),ncomp_p)
  double precision :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
  double precision :: dx(3), xlo(3), time, dt
  integer          :: bc(3,2,ncomp_u), level, grid_no

  integer          :: i, j, k
  double precision :: loc(3), vel(3), mom(3), rho

  do k = lo(3), hi(3)
     loc(3) = xlo(3) + (dble(k - lo(3)) + HALF) * dx(3)
     do j = lo(2), hi(2)
        loc(2) = xlo(2) + (dble(j - lo(2)) + HALF) * dx(2)
        do i = lo(1), hi(1)
           loc(1) = xlo(1) + (dble(i - lo(1)) + HALF) * dx(1)

           rho = u(i,j,k,1)
           vel = u(i,j,k,2:4) / rho
           mom = rho * inertial_velocity(loc, vel, time)
           p(i,j,k,1) = mom(2)

        enddo
     enddo
  enddo

end subroutine ca_derinertialmomentumy



! Derive momentum, given input vector of the grid momenta.

subroutine ca_derinertialmomentumz(p,p_lo,p_hi,ncomp_p, &
                                   u,u_lo,u_hi,ncomp_u, &
                                   lo,hi,domlo,domhi, &
                                   dx,xlo,time,dt,bc,level,grid_no) &
                                   bind(C,name='ca_derinertialmomentumz')

  use bl_constants_module, only: HALF
  use probdata_module, only: inertial_velocity

  implicit none

  integer          :: p_lo(3), p_hi(3), ncomp_p ! == 1
  integer          :: u_lo(3), u_hi(3), ncomp_u ! == 4
  integer          :: lo(3), hi(3), domlo(3), domhi(3)
  double precision :: p(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3),ncomp_p)
  double precision :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
  double precision :: dx(3), xlo(3), time, dt
  integer          :: bc(3,2,ncomp_u), level, grid_no

  integer          :: i, j, k
  double precision :: loc(3), vel(3), mom(3), rho

  do k = lo(3), hi(3)
     loc(3) = xlo(3) + (dble(k - lo(3)) + HALF) * dx(3)
     do j = lo(2), hi(2)
        loc(2) = xlo(2) + (dble(j - lo(2)) + HALF) * dx(2)
        do i = lo(1), hi(1)
           loc(1) = xlo(1) + (dble(i - lo(1)) + HALF) * dx(1)

           rho = u(i,j,k,1)
           vel = u(i,j,k,2:4) / rho
           mom = rho * inertial_velocity(loc, vel, time)
           p(i,j,k,1) = mom(3)

        enddo
     enddo
  enddo

end subroutine ca_derinertialmomentumz



! Derive angular momentum, given input vector of the grid momenta.

subroutine ca_derinertialangmomx(L,L_lo,L_hi,ncomp_L, &
                                 u,u_lo,u_hi,ncomp_u, &
                                 lo,hi,domlo,domhi, &
                                 dx,xlo,time,dt,bc,level,grid_no) &
                                 bind(C,name='ca_derinertialangmomx')

  use bl_constants_module, only: HALF
  use math_module, only: cross_product
  use probdata_module, only: inertial_velocity

  implicit none

  integer          :: L_lo(3), L_hi(3), ncomp_L ! == 1
  integer          :: u_lo(3), u_hi(3), ncomp_u ! == 4
  integer          :: lo(3), hi(3), domlo(3), domhi(3)
  double precision :: L(L_lo(1):L_hi(1),L_lo(2):L_hi(2),L_lo(3):L_hi(3),ncomp_L)
  double precision :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
  double precision :: dx(3), xlo(3), time, dt
  integer          :: bc(3,2,ncomp_u), level, grid_no

  integer          :: i, j, k
  double precision :: loc(3), vel(3), ang_mom(3), rho

  do k = lo(3), hi(3)
     loc(3) = xlo(3) + (dble(k - lo(3)) + HALF) * dx(3)
     do j = lo(2), hi(2)
        loc(2) = xlo(2) + (dble(j - lo(2)) + HALF) * dx(2)
        do i = lo(1), hi(1)
           loc(1) = xlo(1) + (dble(i - lo(1)) + HALF) * dx(1)

           rho = u(i,j,k,1)
           vel = u(i,j,k,2:4) / rho
           ang_mom = cross_product(loc, rho * inertial_velocity(loc, vel, time))

           L(i,j,k,1) = ang_mom(1)

        enddo
     enddo
  enddo

end subroutine ca_derinertialangmomx



subroutine ca_derinertialangmomy(L,L_lo,L_hi,ncomp_L, &
                                 u,u_lo,u_hi,ncomp_u, &
                                 lo,hi,domlo,domhi, &
                                 dx,xlo,time,dt,bc,level,grid_no) &
                                 bind(C,name='ca_derinertialangmomy')

  use bl_constants_module, only: HALF
  use math_module, only: cross_product
  use probdata_module, only: inertial_velocity

  implicit none

  integer          :: L_lo(3), L_hi(3), ncomp_L ! == 1
  integer          :: u_lo(3), u_hi(3), ncomp_u ! == 4
  integer          :: lo(3), hi(3), domlo(3), domhi(3)
  double precision :: L(L_lo(1):L_hi(1),L_lo(2):L_hi(2),L_lo(3):L_hi(3),ncomp_L)
  double precision :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
  double precision :: dx(3), xlo(3), time, dt
  integer          :: bc(3,2,ncomp_u), level, grid_no

  integer          :: i, j, k
  double precision :: loc(3), vel(3), ang_mom(3), rho

  do k = lo(3), hi(3)
     loc(3) = xlo(3) + (dble(k - lo(3)) + HALF) * dx(3)
     do j = lo(2), hi(2)
        loc(2) = xlo(2) + (dble(j - lo(2)) + HALF) * dx(2)
        do i = lo(1), hi(1)
           loc(1) = xlo(1) + (dble(i - lo(1)) + HALF) * dx(1)

           rho = u(i,j,k,1)
           vel = u(i,j,k,2:4) / rho
           ang_mom = cross_product(loc, rho * inertial_velocity(loc, vel, time))          

           L(i,j,k,1) = ang_mom(2)

        enddo
     enddo
  enddo

end subroutine ca_derinertialangmomy



subroutine ca_derinertialangmomz(L,L_lo,L_hi,ncomp_L, &
                                 u,u_lo,u_hi,ncomp_u, &
                                 lo,hi,domlo,domhi, &
                                 dx,xlo,time,dt,bc,level,grid_no) &
                                 bind(C,name='ca_derinertialangmomz')

  use bl_constants_module, only: HALF
  use math_module, only: cross_product
  use probdata_module, only: inertial_velocity

  implicit none

  integer          :: L_lo(3), L_hi(3), ncomp_L ! == 1
  integer          :: u_lo(3), u_hi(3), ncomp_u ! == 4
  integer          :: lo(3), hi(3), domlo(3), domhi(3)
  double precision :: L(L_lo(1):L_hi(1),L_lo(2):L_hi(2),L_lo(3):L_hi(3),ncomp_L)
  double precision :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
  double precision :: dx(3), xlo(3), time, dt
  integer          :: bc(3,2,ncomp_u), level, grid_no

  integer          :: i, j, k
  double precision :: loc(3), vel(3), ang_mom(3), rho

  do k = lo(3), hi(3)
     loc(3) = xlo(3) + (dble(k - lo(3)) + HALF) * dx(3)
     do j = lo(2), hi(2)
        loc(2) = xlo(2) + (dble(j - lo(2)) + HALF) * dx(2)
        do i = lo(1), hi(1)
           loc(1) = xlo(1) + (dble(i - lo(1)) + HALF) * dx(1)

           rho = u(i,j,k,1)
           vel = u(i,j,k,2:4) / rho
           ang_mom = cross_product(loc, rho * inertial_velocity(loc, vel, time))

           L(i,j,k,1) = ang_mom(3)

        enddo
     enddo
  enddo

end subroutine ca_derinertialangmomz



! Derive the effective potential phiEff = phiGrav + phiRot

subroutine ca_derphieff(phi,phi_lo,phi_hi,ncomp_phi, &
                        u,u_lo,u_hi,ncomp_u, &
                        lo,hi,domlo,domhi, &
                        dx,xlo,time,dt,bc,level,grid_no) bind(C,name='ca_derphieff')

  implicit none

  integer          :: phi_lo(3), phi_hi(3), ncomp_phi ! == 1
  integer          :: u_lo(3), u_hi(3), ncomp_u ! == 2
  integer          :: lo(3), hi(3), domlo(3), domhi(3)
  double precision :: phi(phi_lo(1):phi_hi(1),phi_lo(2):phi_hi(2),phi_lo(3):phi_hi(3),ncomp_phi)
  double precision :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
  double precision :: dx(3), xlo(3), time, dt
  integer          :: bc(3,2,ncomp_u), level, grid_no

  integer          :: i, j, k

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)

           ! Note that we need to account for the fact that phiGrav
           ! is positive in Castro.

           phi(i,j,k,1) = -u(i,j,k,1) + u(i,j,k,2)

        enddo
     enddo
  enddo

end subroutine ca_derphieff



! Derive an approximation to the effective potential of the primary only,
! by treating it as a point-mass at its center of mass.
! The u array contains the rotational potential, so we only need to calculate
! the gravitational potential from the point-mass.

subroutine ca_derphieffpm_p(phi,phi_lo,phi_hi,ncomp_phi, &
                            u,u_lo,u_hi,ncomp_u, &
                            lo,hi,domlo,domhi, &
                            dx,xlo,time,dt,bc,level,grid_no) bind(C,name='ca_derphieffpm_p')

  use bl_constants_module, only: ZERO, HALF
  use probdata_module, only: mass_P, com_P
  use fundamental_constants_module, only: Gconst

  implicit none

  integer          :: phi_lo(3), phi_hi(3), ncomp_phi ! == 1
  integer          :: u_lo(3), u_hi(3), ncomp_u ! == 1
  integer          :: lo(3), hi(3), domlo(3), domhi(3)
  double precision :: phi(phi_lo(1):phi_hi(1),phi_lo(2):phi_hi(2),phi_lo(3):phi_hi(3),ncomp_phi)
  double precision :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
  double precision :: dx(3), xlo(3), time, dt
  integer          :: bc(3,2,ncomp_u), level, grid_no

  integer          :: i, j, k
  double precision :: loc(3), r

  ! Don't do anything here if the star no longer exists

  phi(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) = ZERO

  if (mass_P == ZERO) return

  do k = lo(3), hi(3)
     loc(3) = xlo(3) + (dble(k - lo(3)) + HALF) * dx(3)
     do j = lo(2), hi(2)
        loc(2) = xlo(2) + (dble(j - lo(2)) + HALF) * dx(2)
        do i = lo(1), hi(1)
           loc(1) = xlo(1) + (dble(i - lo(1)) + HALF) * dx(1)

           r = sqrt( sum( (loc - com_P)**2 ) )

           phi(i,j,k,1) = -Gconst * mass_P / r + u(i,j,k,1)

        enddo
     enddo
  enddo

end subroutine ca_derphieffpm_p



! Same as above, but for the secondary.

subroutine ca_derphieffpm_s(phi,phi_lo,phi_hi,ncomp_phi, &
                            u,u_lo,u_hi,ncomp_u, &
                            lo,hi,domlo,domhi, &
                            dx,xlo,time,dt,bc,level,grid_no) bind(C,name='ca_derphieffpm_s')

  use bl_constants_module, only: ZERO, HALF
  use probdata_module, only: mass_S, com_S
  use fundamental_constants_module, only: Gconst

  implicit none

  integer          :: phi_lo(3), phi_hi(3), ncomp_phi ! == 1
  integer          :: u_lo(3), u_hi(3), ncomp_u ! == 1
  integer          :: lo(3), hi(3), domlo(3), domhi(3)
  double precision :: phi(phi_lo(1):phi_hi(1),phi_lo(2):phi_hi(2),phi_lo(3):phi_hi(3),ncomp_phi)
  double precision :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
  double precision :: dx(3), xlo(3), time, dt
  integer          :: bc(3,2,ncomp_u), level, grid_no

  integer          :: i, j, k
  double precision :: loc(3), r

  ! Don't do anything here if the star no longer exists

  phi(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) = ZERO

  if (mass_S == ZERO) return

  do k = lo(3), hi(3)
     loc(3) = xlo(3) + (dble(k - lo(3)) + HALF) * dx(3)
     do j = lo(2), hi(2)
        loc(2) = xlo(2) + (dble(j - lo(2)) + HALF) * dx(2)
        do i = lo(1), hi(1)
           loc(1) = xlo(1) + (dble(i - lo(1)) + HALF) * dx(1)

           r = sqrt( sum( (loc - com_S)**2 ) )

           phi(i,j,k,1) = -Gconst * mass_S / r + u(i,j,k,1)

        enddo
     enddo
  enddo

end subroutine ca_derphieffpm_s
