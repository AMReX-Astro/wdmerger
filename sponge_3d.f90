subroutine sponge(uout,uout_l1,uout_l2,uout_l3,&
                  uout_h1,uout_h2,uout_h3,lo,hi,dt, &
                  dx,dy,dz,domlo,domhi)

  use bl_constants_module, ONLY: M_PI, HALF, ZERO, ONE
  use meth_params_module , only : NVAR, URHO, UMX, UMY, UMZ, UEDEN
      

  implicit none
  integer          :: lo(3),hi(3)
  integer          :: uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3
  double precision :: uout(uout_l1:uout_h1,uout_l2:uout_h2,&
       uout_l3:uout_h3,NVAR)
  double precision :: dt
  double precision :: dx, dy, dz, domlo(3), domhi(3)
  
  double precision :: xx, yy, zz
  double precision :: xctr, yctr, zctr
  
  double precision :: ke_old, ke_new, sponge_mult, delta_sponge, radius, r_sponge, fac

  integer :: i,j,k

  xctr = HALF*(domlo(1) + domhi(1))
  yctr = HALF*(domlo(2) + domhi(2))
  zctr = HALF*(domlo(3) + domhi(3))
  
  r_sponge = 0.25*(domhi(1)-domlo(1))
  delta_sponge = 0.025*(domhi(1)-domlo(1))
  
  do k = lo(3),hi(3)
     zz = domlo(3) + dble(k + HALF)*dz 
     
     do j = lo(2),hi(2)
        yy = domlo(2) + dble(j + HALF)*dy
        
        do i = lo(1),hi(1)
           xx = domlo(1) + dble(i + HALF)*dx
           
           radius = sqrt((xx-xctr)**2 + (yy-yctr)**2 + (zz-zctr)**2)
           
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

