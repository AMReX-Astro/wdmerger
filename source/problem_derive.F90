module problem_derive_module

  implicit none

  public

contains

  ! Derive momentum, given input vector of the grid momenta.
  
  subroutine ca_derinertialmomentumx(p,p_lo,p_hi,ncomp_p, &
                                     u,u_lo,u_hi,ncomp_u, &
                                     lo,hi,domlo,domhi, &
                                     dx,xlo,time,dt,bc,level,grid_no) bind(C)

    use bl_constants_module, only: HALF
    use meth_params_module, only : URHO, UMX, UMY, UMZ, UEINT
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
                                     dx,xlo,time,dt,bc,level,grid_no) bind(C)

    use bl_constants_module, only: HALF
    use meth_params_module, only : URHO, UMX, UMY, UMZ, UEINT
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
                                     dx,xlo,time,dt,bc,level,grid_no) bind(C)

    use bl_constants_module, only: HALF
    use meth_params_module, only : URHO, UMX, UMY, UMZ, UEINT
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
                                   dx,xlo,time,dt,bc,level,grid_no) bind(C)

    use bl_constants_module, only: HALF
    use meth_params_module, only : URHO, UMX, UMY, UMZ, UEINT
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
                                   dx,xlo,time,dt,bc,level,grid_no) bind(C)

    use bl_constants_module, only: HALF
    use meth_params_module, only : URHO, UMX, UMY, UMZ, UEINT
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
                                   dx,xlo,time,dt,bc,level,grid_no) bind(C)

    use bl_constants_module, only: HALF
    use meth_params_module, only : URHO, UMX, UMY, UMZ, UEINT
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

end module problem_derive_module
