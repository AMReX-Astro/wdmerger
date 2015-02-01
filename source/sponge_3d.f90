module sponge_module

  implicit none

contains

  subroutine sponge(uout,uout_l1,uout_l2,uout_l3,&
                    uout_h1,uout_h2,uout_h3,lo,hi,time,dt, &
                    dx,dy,dz,domlo,domhi)

    use bl_constants_module, only: M_PI, HALF, ZERO, ONE
    use meth_params_module , only: NVAR, URHO, UMX, UMY, UMZ, UEDEN
    use prob_params_module, only: problo, probhi, center

    implicit none
    integer          :: lo(3),hi(3),domlo(3),domhi(3)
    integer          :: uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3
    double precision :: uout(uout_l1:uout_h1, &
                             uout_l2:uout_h2,&
                             uout_l3:uout_h3,NVAR)
    double precision :: time, dt
    double precision :: dx, dy, dz

    double precision :: xx, yy, zz
    double precision :: xctr, yctr, zctr

    double precision :: ke_old, ke_new, sponge_mult, delta_sponge, &
                        radius, r_sponge, fac

    integer :: i,j,k


    r_sponge = 0.25*(probhi(1)-problo(1))
    delta_sponge = 0.025*(probhi(1)-problo(1))

    do k = lo(3),hi(3)
       zz = problo(3) + dble(k + HALF)*dz

       do j = lo(2),hi(2)
          yy = problo(2) + dble(j + HALF)*dy

          do i = lo(1),hi(1)
             xx = problo(1) + dble(i + HALF)*dx

             radius = sqrt((xx-center(1))**2 + (yy-center(2))**2 + (zz-center(3))**2)

             sponge_mult = ONE

             if (radius > r_sponge) then
                if (radius < r_sponge+delta_sponge) then
                   fac = HALF*(ONE - cos(M_PI*(radius-r_sponge)/delta_sponge))
                else
                   fac = ONE
                endif
                sponge_mult = ONE/(ONE+ dt*fac*100.0)
             endif

             ke_old = HALF * (uout(i,j,k,UMX)**2 + &
                              uout(i,j,k,UMY)**2 + &
                              uout(i,j,k,UMZ)**2 )/ uout(i,j,k,URHO)


             uout(i,j,k,UMX  ) = uout(i,j,k,UMX  ) * sponge_mult
             uout(i,j,k,UMY  ) = uout(i,j,k,UMY  ) * sponge_mult
             uout(i,j,k,UMZ  ) = uout(i,j,k,UMZ  ) * sponge_mult

             ke_new = HALF * (uout(i,j,k,UMX)**2 + &
                              uout(i,j,k,UMY)**2 + &
                              uout(i,j,k,UMZ)**2 )/ uout(i,j,k,URHO)

             uout(i,j,k,UEDEN) = uout(i,j,k,UEDEN) + (ke_new-ke_old)

          enddo
       enddo
    enddo


  end subroutine sponge

end module sponge_module
