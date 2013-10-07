module grav_sources_module

  implicit none

  private

  public add_grav_source

contains

! :::
! ::: ------------------------------------------------------------------
! :::

    subroutine add_grav_source(uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
                               uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, &
                               grav, gv_l1, gv_l2, gv_l3, gv_h1, gv_h2, gv_h3, &
                               lo,hi,dt,dx,E_added)

      use bl_constants_module, only: M_PI
      use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, grav_source_type, rot_period
      use prob_params_module, only: coord_type
      use probdata_module, only: center

      implicit none

      integer lo(3), hi(3)
      integer uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3
      integer  uout_l1, uout_l2, uout_l3, uout_h1, uout_h2, uout_h3
      integer  gv_l1, gv_l2, gv_l3, gv_h1, gv_h2, gv_h3

      double precision  uin( uin_l1: uin_h1, uin_l2: uin_h2, uin_l3: uin_h3,NVAR)
      double precision uout(uout_l1:uout_h1,uout_l2:uout_h2,uout_l3:uout_h3,NVAR)
      double precision grav(  gv_l1:  gv_h1,  gv_l2:  gv_h2,  gv_l3:  gv_h3,3)
      double precision dt, dx(3)
      double precision E_added

      double precision :: rho
      double precision :: SrU, SrV, SrW, SrE
      double precision :: rhoInv
      double precision :: old_rhoeint, new_rhoeint, old_ke, new_ke, old_re
      double precision :: x, y, z
      integer          :: i, j, k

      double precision, parameter :: TWO_PI = 2.d0 * M_PI
      double precision :: omega

      ! should come in through a module
      double precision :: problo(3)


      if (coord_type == 0) then
         omega = TWO_PI / rot_period
      else
         call bl_error("Error:: Rotate_3d.f90 :: unknown coord_type")
      endif

      problo(:) = 0.0d0

      ! hack: do both gravity and rotation here, so we can couple them.
      ! we need to make sure that Rotate_3d.f90 does nothing.

      do k = lo(3),hi(3)
         z = problo(3) + dx(3)*(float(k)+0.5d0) - center(3)

         do j = lo(2),hi(2)
            y = problo(2) + dx(2)*(float(j)+0.5d0) - center(2)

            do i = lo(1),hi(1)
               x = problo(1) + dx(1)*(float(i)+0.5d0) - center(1)


               rho    = uin(i,j,k,URHO)
               rhoInv = 1.0d0 / rho

               ! momenta

               ! gravity
               SrU = rho * grav(i,j,k,1)
               SrV = rho * grav(i,j,k,2)
               SrW = rho * grav(i,j,k,3)

               ! rotation
               !SrU = SrU + 2.d0*uin(i,j,k,UMY)*omega + rho*omega**2*x
               !SrV = SrV - 2.d0*uin(i,j,k,UMX)*omega + rho*omega**2*y

               ! update the momenta
               uout(i,j,k,UMX)   = uout(i,j,k,UMX) + SrU * dt
               uout(i,j,k,UMY)   = uout(i,j,k,UMY) + SrV * dt
               uout(i,j,k,UMZ)   = uout(i,j,k,UMZ) + SrW * dt

               
               ! energy

               ! Src = rho u dot g, evaluated with all quantities at t^n
               SrE = uin(i,j,k,UMX) * grav(i,j,k,1) + &
                     uin(i,j,k,UMY) * grav(i,j,k,2) + &
                     uin(i,j,k,UMZ) * grav(i,j,k,3)

               ! rotation
               !SrE = SrE + omega**2*(uin(i,j,k,UMX)*x + uin(i,j,k,UMY)*y)

               uout(i,j,k,UEDEN) = uout(i,j,k,UEDEN) + SrE * dt

            enddo
         enddo
      enddo

      end subroutine add_grav_source

end module grav_sources_module

! :::
! ::: ------------------------------------------------------------------
! :::

      subroutine ca_corrgsrc(lo,hi, &
                             gold,gold_l1,gold_l2,gold_l3,gold_h1,gold_h2,gold_h3, &
                             gnew,gnew_l1,gnew_l2,gnew_l3,gnew_h1,gnew_h2,gnew_h3, &
                             uold,uold_l1,uold_l2,uold_l3,uold_h1,uold_h2,uold_h3, &
                             unew,unew_l1,unew_l2,unew_l3,unew_h1,unew_h2,unew_h3, &
                             dt,dx,E_added)

      use bl_constants_module, only: M_PI
      use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, grav_source_type, rot_period
      use prob_params_module, only: coord_type
      use probdata_module, only: center

      implicit none

      integer lo(3),hi(3)
      integer gold_l1,gold_l2,gold_l3,gold_h1,gold_h2,gold_h3
      integer gnew_l1,gnew_l2,gnew_l3,gnew_h1,gnew_h2,gnew_h3
      integer uold_l1,uold_l2,uold_l3,uold_h1,uold_h2,uold_h3
      integer unew_l1,unew_l2,unew_l3,unew_h1,unew_h2,unew_h3
      double precision   gold(gold_l1:gold_h1,gold_l2:gold_h2,gold_l3:gold_h3,3)
      double precision   gnew(gnew_l1:gnew_h1,gnew_l2:gnew_h2,gnew_l3:gnew_h3,3)
      double precision  uold(uold_l1:uold_h1,uold_l2:uold_h2,uold_l3:uold_h3,NVAR)
      double precision  unew(unew_l1:unew_h1,unew_l2:unew_h2,unew_l3:unew_h3,NVAR)
      double precision  dt,dx(3),E_added
      double precision :: x, y, z
      integer i,j,k

      double precision Sx, Sy, Sz
      double precision SrEcorr
      double precision rhoo
      double precision rhon

      double precision rhooinv, rhoninv
      double precision old_ke, old_rhoeint, old_re
      double precision new_ke, new_rhoeint

      double precision, parameter :: TWO_PI = 2.d0 * M_PI
      double precision :: omega

      ! should come in through a module
      double precision :: problo(3)


      if (coord_type == 0) then
         omega = TWO_PI/rot_period
      else
         call bl_error("Error:: Rotate_3d.f90 :: unknown coord_type")
      endif

      problo(:) = 0.0d0

      !print *, 'in corrgsrc, rot_period = ', rot_period

      do k = lo(3),hi(3)
         z = problo(3) + dx(3)*(float(k)+0.5d0) - center(3)

         do j = lo(2),hi(2)
            y = problo(2) + dx(2)*(float(j)+0.5d0) - center(2)

            do i = lo(1),hi(1)
               x = problo(1) + dx(1)*(float(i)+0.5d0) - center(1)


               rhoo    = uold(i,j,k,URHO)
               rhooinv = 1.0d0 / uold(i,j,k,URHO)

               rhon    = unew(i,j,k,URHO)
               rhoninv = 1.0d0 / unew(i,j,k,URHO)


               ! momenta source corrections
               Sx = unew(i,j,k,UMX) + &
                    0.5d0*dt*(rhon*gnew(i,j,k,1) - rhoo*gold(i,j,k,1)) + &
                    0.5d0*dt*(rhon*omega**2*x - rhoo*omega**2*x) - dt*uold(i,j,k,UMY)*omega

               Sy = unew(i,j,k,UMY) + &
                    0.5d0*dt*(rhon*gnew(i,j,k,2) - rhoo*gold(i,j,k,2)) + &
                    0.5d0*dt*(rhon*omega**2*y - rhoo*omega**2*y) + dt*uold(i,j,k,UMX)*omega

               Sz = unew(i,j,k,UMZ) + &
                    0.5d0*dt*(rhon*gnew(i,j,k,3) - rhoo*gold(i,j,k,3)) 


               ! the (rho u) and (rho v) updates are intermingled.  We do it implicitly
               unew(i,j,k,UMY) = (Sy - dt*omega*Sx)/(1.d0 + dt**2 * omega**2)
               unew(i,j,k,UMX) = Sx + dt*omega*unew(i,j,k,UMY)

               unew(i,j,k,UMZ) = Sz



               ! energy
               SrEcorr =  0.5d0 * ( (unew(i,j,k,UMX)*gnew(i,j,k,1) - &
                                     uold(i,j,k,UMX)*gold(i,j,k,1)) + &
                                    (unew(i,j,k,UMY)*gnew(i,j,k,2) - &
                                     uold(i,j,k,UMY)*gold(i,j,k,2)) + &
                                    (unew(i,j,k,UMZ)*gnew(i,j,k,3) - &
                                     uold(i,j,k,UMZ)*gold(i,j,k,3)) )

               SrEcorr = SrEcorr + 0.5d0*omega**2*( (unew(i,j,k,UMX)*x + unew(i,j,k,UMY)*y) - &
                                                    (uold(i,j,k,UMX)*x + uold(i,j,k,UMY)*y) )

               unew(i,j,k,UEDEN) = unew(i,j,k,UEDEN) + SrEcorr*dt


            enddo
         enddo
      enddo


      end subroutine ca_corrgsrc

! :::
! ::: ------------------------------------------------------------------
! :::

     subroutine ca_syncgsrc(lo,hi, &
                            gphi,gphi_l1,gphi_l2,gphi_l3,gphi_h1,gphi_h2,gphi_h3, &
                            gdphi,gdphi_l1,gdphi_l2,gdphi_l3,gdphi_h1,gdphi_h2,gdphi_h3, &
                            state,state_l1,state_l2,state_l3,state_h1,state_h2,state_h3, &
                            dstate,dstate_l1,dstate_l2,dstate_l3, &
                            dstate_h1,dstate_h2,dstate_h3, &
                            sync_src,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3,dt)

     use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ

     implicit none

     integer lo(3),hi(3)
     integer gphi_l1,gphi_l2,gphi_l3,gphi_h1,gphi_h2,gphi_h3
     integer gdphi_l1,gdphi_l2,gdphi_l3,gdphi_h1,gdphi_h2,gdphi_h3
     integer state_l1,state_l2,state_l3,state_h1,state_h2,state_h3
     integer dstate_l1,dstate_l2,dstate_l3,dstate_h1,dstate_h2,dstate_h3
     integer src_l1,src_l2,src_l3,src_h1,src_h2,src_h3
     double precision   gphi(gphi_l1:gphi_h1,gphi_l2:gphi_h2,gphi_l3:gphi_h3,3)
     double precision  gdphi(gdphi_l1:gdphi_h1,gdphi_l2:gdphi_h2,gdphi_l3:gdphi_h3,3)
     double precision  state(state_l1:state_h1,state_l2:state_h2,state_l3:state_h3,NVAR)
     double precision dstate(dstate_l1:dstate_h1,dstate_l2:dstate_h2,dstate_l3:dstate_h3,3+1)
     double precision sync_src(src_l1:src_h1,src_l2:src_h2,src_l3:src_h3,3+1)
     double precision dt

     !    Note that dstate is drho and drhoU, state is the entire state, and src
     !    is S_rhoU and S_rhoE

     integer          :: i,j,k
     double precision :: rho_pre, rhoU_pre, rhoV_pre, rhoW_pre
     double precision :: gx, gy, gz, dgx, dgy, dgz, SrU, SrV, SrW, SrE

     !$OMP PARALLEL DO PRIVATE(i,j,k,rho_pre,rhoU_pre,rhoV_pre,rhoW_pre,gx,gy,gz,dgx,dgy,dgz,SrU,SrV,SrW,SrE)
     do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               rho_pre  = state(i,j,k,URHO) - dstate(i,j,k,1)
               rhoU_pre = state(i,j,k,UMX)  - dstate(i,j,k,2)
               rhoV_pre = state(i,j,k,UMY)  - dstate(i,j,k,3)
               rhoW_pre = state(i,j,k,UMZ)  - dstate(i,j,k,4)

               gx  = gphi(i,j,k,1)
               gy  = gphi(i,j,k,2)
               gz  = gphi(i,j,k,3)

               dgx = gdphi(i,j,k,1)
               dgy = gdphi(i,j,k,2)
               dgz = gdphi(i,j,k,3)

               SrU = dstate(i,j,k,1)*gx + rho_pre*dgx
               SrV = dstate(i,j,k,1)*gy + rho_pre*dgy
               SrW = dstate(i,j,k,1)*gz + rho_pre*dgz

               SrE = ( SrU * (rhoU_pre + (0.5d0*dt)*SrU) + &
                       SrV * (rhoV_pre + (0.5d0*dt)*SrV) + &
                       SrW * (rhoW_pre + (0.5d0*dt)*SrW) ) / rho_pre

               sync_src(i,j,k,1) = SrU
               sync_src(i,j,k,2) = SrV
               sync_src(i,j,k,3) = SrW
               sync_src(i,j,k,4) = SrE

            enddo
         enddo
     enddo
     !$OMP END PARALLEL DO

     end subroutine ca_syncgsrc
