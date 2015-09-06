module sponge_module

  implicit none

contains

  subroutine sponge(uout,uout_l1,uout_l2,uout_l3,&
                    uout_h1,uout_h2,uout_h3,lo,hi,time,dt, &
                    dx,dy,dz,domlo,domhi,&
                    E_added,xmom_added,ymom_added,zmom_added)

    use bl_constants_module, only: M_PI, HALF, ZERO, ONE
    use meth_params_module , only: NVAR, URHO, UMX, UMY, UMZ, UEDEN
    use prob_params_module, only: problo, probhi, center
    use probdata_module, only: sponge_dist, sponge_width, sponge_timescale

    implicit none
    integer          :: lo(3),hi(3),domlo(3),domhi(3)
    integer          :: uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3
    double precision :: uout(uout_l1:uout_h1, &
                             uout_l2:uout_h2,&
                             uout_l3:uout_h3,NVAR)
    double precision :: time, dt
    double precision :: dx, dy, dz
    double precision :: E_added, xmom_added, ymom_added, zmom_added

    ! Local variables

    double precision :: r(3), radius
    double precision :: ke_old, E_old, xmom_old, ymom_old, zmom_old
    double precision :: sponge_radius, sponge_delta, fac, alpha
    integer :: i,j,k

    sponge_radius = sponge_dist * max( maxval(abs(problo-center)), maxval(abs(probhi-center)) )
    sponge_delta  = sponge_width * max( maxval(abs(problo-center)), maxval(abs(probhi-center)) )
    alpha = dt / sponge_timescale

    do k = lo(3), hi(3)
       r(3) = problo(3) + dble(k + HALF) * dz - center(3)

       do j = lo(2), hi(2)
          r(2) = problo(2) + dble(j + HALF) * dy - center(2)

          do i = lo(1), hi(1)
             r(1) = problo(1) + dble(i + HALF) * dx - center(1)

             ! Starting diagnostic quantities

             E_old    = uout(i,j,k,UEDEN)
             ke_old   = HALF * sum(uout(i,j,k,UMX:UMZ)**2) / uout(i,j,k,URHO)
             xmom_old = uout(i,j,k,UMX)
             ymom_old = uout(i,j,k,UMY)
             zmom_old = uout(i,j,k,UMZ)

             ! Apply sponge

             radius = sqrt(sum(r**2))

             if (radius >= sponge_radius .and. radius < sponge_radius + sponge_delta) then
                fac = HALF * (ONE - cos(M_PI * (radius - sponge_radius) / sponge_delta))
             else if (radius >= sponge_radius + sponge_delta) then
                fac = ONE
             else
                fac = ZERO
             endif

             uout(i,j,k,UMX:UMZ) = uout(i,j,k,UMX:UMZ) / (ONE + alpha * fac)

             uout(i,j,k,UEDEN) = uout(i,j,k,UEDEN) + HALF * sum(uout(i,j,k,UMX:UMZ)**2) / uout(i,j,k,URHO) - ke_old

             ! Ending diagnostic quantities

             E_added    = E_added    + uout(i,j,k,UEDEN) - E_old
             xmom_added = xmom_added + uout(i,j,k,UMX)   - xmom_old
             ymom_added = ymom_added + uout(i,j,k,UMY)   - ymom_old
             zmom_added = zmom_added + uout(i,j,k,UMZ)   - zmom_old

          enddo
       enddo
    enddo


  end subroutine sponge

end module sponge_module
