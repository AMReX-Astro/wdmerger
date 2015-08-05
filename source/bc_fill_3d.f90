   ! ::: -----------------------------------------------------------

   subroutine ca_hypfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
                         domlo,domhi,delta,xlo,time,bc)

     use meth_params_module, only: NVAR, URHO, UMX, UMY, UMZ, UTEMP, &
          UEDEN, UEINT, UFS, rot_period
     use prob_params_module, only: center
     use bl_constants_module
     use probdata_module
     use eos_module

     implicit none
     include 'bc_types.fi'
     integer adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3
     integer bc(3,2,*)
     integer domlo(3), domhi(3)
     double precision delta(3), xlo(3), time
     double precision adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3,NVAR)

     integer i, j, k, n
     double precision :: vx, vy, vz
     double precision :: xx, yy

     type (eos_t) :: ambient_state

     call get_ambient(ambient_state)

     do n = 1,NVAR
        call filcc(adv(adv_l1,adv_l2,adv_l3,n), &
                   adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
                   domlo,domhi,delta,xlo,bc(1,1,n))
     enddo

     ! Override the generic routine at the physical boundaries by
     ! setting the material to the ambient state

     if (fill_ambient_bc) then

        ! -x
        if (adv_l1 < domlo(1) .and. bc(1,1,1) .ne. 0) then
           do k = adv_l3, adv_h3
              do j = adv_l2, adv_h2
                 yy = xlo(2) + dble(j - domlo(2) + HALF)*delta(2) - center(2)
                 do i = adv_l1, domlo(1)-1
                    xx = xlo(1) + dble(i - domlo(1) + HALF)*delta(1) - center(1)

                    adv(i,j,k,URHO) = ambient_state % rho
                    adv(i,j,k,UTEMP) = ambient_state % T
                    adv(i,j,k,UFS:UFS-1+nspec) = ambient_state % rho * ambient_state % xn(:)

                    if ( orbital_kick ) then

                      vx = (-2.0d0 * M_PI / rot_period) * yy
                      vy = ( 2.0d0 * M_PI / rot_period) * xx
                      vz = ZERO

                    else

                      vx = adv(domlo(1),j,k,UMX)/adv(domlo(1),j,k,URHO)
                      vy = adv(domlo(1),j,k,UMY)/adv(domlo(1),j,k,URHO)
                      vz = adv(domlo(1),j,k,UMZ)/adv(domlo(1),j,k,URHO)

                    endif

                    adv(i,j,k,UMX) = ambient_state % rho*vx
                    adv(i,j,k,UMY) = ambient_state % rho*vy
                    adv(i,j,k,UMZ) = ambient_state % rho*vz

                    adv(i,j,k,UEINT) = ambient_state % rho*ambient_state % e
                    adv(i,j,k,UEDEN) = ambient_state % rho*ambient_state % e + &
                         HALF*(adv(i,j,k,UMX)**2 + &
                               adv(i,j,k,UMY)**2 + &
                               adv(i,j,k,UMZ)**2)/ambient_state % rho

                 enddo
              enddo
           enddo
        endif

        ! +x
        if (adv_h1 > domhi(1) .and. bc(1,2,1) .ne. 0) then
           do k = adv_l3, adv_h3
              do j = adv_l2, adv_h2
                 yy = xlo(2) + dble(j - domlo(2) + HALF)*delta(2) - center(2)
                 do i = domhi(1)+1, adv_h1
                    xx = xlo(1) + dble(i - domlo(1) + HALF)*delta(1) - center(1)

                    adv(i,j,k,URHO) = ambient_state % rho
                    adv(i,j,k,UTEMP) = ambient_state % T
                    adv(i,j,k,UFS:UFS-1+nspec) = ambient_state % rho*ambient_state % xn(:)

                    if ( orbital_kick ) then

                      vx = (-2.0d0 * M_PI / rot_period) * yy
                      vy = ( 2.0d0 * M_PI / rot_period) * xx
                      vz = ZERO

                    else

                      vx = adv(domhi(1),j,k,UMX)/adv(domhi(1),j,k,URHO)
                      vy = adv(domhi(1),j,k,UMY)/adv(domhi(1),j,k,URHO)
                      vz = adv(domhi(1),j,k,UMZ)/adv(domhi(1),j,k,URHO)

                    endif

                    adv(i,j,k,UMX) = ambient_state % rho*vx
                    adv(i,j,k,UMY) = ambient_state % rho*vy
                    adv(i,j,k,UMZ) = ambient_state % rho*vz

                    adv(i,j,k,UEINT) = ambient_state % rho*ambient_state % e
                    adv(i,j,k,UEDEN) = ambient_state % rho*ambient_state % e + &
                         HALF*(adv(i,j,k,UMX)**2 + &
                               adv(i,j,k,UMY)**2 + &
                               adv(i,j,k,UMZ)**2)/ambient_state % rho

                 enddo
              enddo
           enddo
        endif

        ! -y
        if (adv_l2 < domlo(2) .and. bc(2,1,1) .ne. 0) then
           do k = adv_l3, adv_h3
              do j = adv_l2, domlo(2)-1
                 yy = xlo(2) + dble(j - domlo(2) + HALF)*delta(2) - center(2)
                 do i = adv_l1, adv_h1
                    xx = xlo(1) + dble(i - domlo(1) + HALF)*delta(1) - center(1)

                    adv(i,j,k,URHO) = ambient_state % rho
                    adv(i,j,k,UTEMP) = ambient_state % T
                    adv(i,j,k,UFS:UFS-1+nspec) = ambient_state % rho*ambient_state % xn(:)

                    if ( orbital_kick ) then

                      vx = (-2.0d0 * M_PI / rot_period) * yy
                      vy = ( 2.0d0 * M_PI / rot_period) * xx
                      vz = ZERO

                    else

                      vx = adv(i,domlo(2),k,UMX)/adv(i,domlo(2),k,URHO)
                      vy = adv(i,domlo(2),k,UMY)/adv(i,domlo(2),k,URHO)
                      vz = adv(i,domlo(2),k,UMZ)/adv(i,domlo(2),k,URHO)

                    endif

                    adv(i,j,k,UMX) = ambient_state % rho*vx
                    adv(i,j,k,UMY) = ambient_state % rho*vy
                    adv(i,j,k,UMZ) = ambient_state % rho*vz

                    adv(i,j,k,UEINT) = ambient_state % rho*ambient_state % e
                    adv(i,j,k,UEDEN) = ambient_state % rho*ambient_state % e + &
                         HALF*(adv(i,j,k,UMX)**2 + &
                               adv(i,j,k,UMY)**2 + &
                               adv(i,j,k,UMZ)**2)/ambient_state % rho

                 enddo
              enddo
           enddo
        endif

        ! +y
        if (adv_h2 > domhi(2) .and. bc(2,2,1) .ne. 0) then
           do k = adv_l3, adv_h3
              do j = domhi(2)+1, adv_h2
                 yy = xlo(2) + dble(j - domlo(2) + HALF)*delta(2) - center(2)
                 do i = adv_l1, adv_h1
                    xx = xlo(1) + dble(i - domlo(1) + HALF)*delta(1) - center(1)

                    adv(i,j,k,URHO) = ambient_state % rho
                    adv(i,j,k,UTEMP) = ambient_state % T
                    adv(i,j,k,UFS:UFS-1+nspec) = ambient_state % rho*ambient_state % xn(:)

                    if ( orbital_kick ) then

                      vx = (-2.0d0 * M_PI / rot_period) * yy
                      vy = ( 2.0d0 * M_PI / rot_period) * xx
                      vz = ZERO

                    else

                      vx = adv(i,domhi(2),k,UMX)/adv(i,domhi(2),k,URHO)
                      vy = adv(i,domhi(2),k,UMY)/adv(i,domhi(2),k,URHO)
                      vz = adv(i,domhi(2),k,UMZ)/adv(i,domhi(2),k,URHO)

                    endif

                    adv(i,j,k,UMX) = ambient_state % rho*vx
                    adv(i,j,k,UMY) = ambient_state % rho*vy
                    adv(i,j,k,UMZ) = ambient_state % rho*vz

                    adv(i,j,k,UEINT) = ambient_state % rho*ambient_state % e
                    adv(i,j,k,UEDEN) = ambient_state % rho*ambient_state % e + &
                         HALF*(adv(i,j,k,UMX)**2 + &
                               adv(i,j,k,UMY)**2 + &
                               adv(i,j,k,UMZ)**2)/ambient_state % rho

                 enddo
              enddo
           enddo
        endif

        ! -z
        if (adv_l3 < domlo(3) .and. bc(3,1,1) .ne. 0) then
           do k = adv_l3, domlo(3)-1
              do j = adv_l2, adv_h2
                 yy = xlo(2) + dble(j - domlo(2) + HALF)*delta(2) - center(2)
                 do i = adv_l1, adv_h1
                    xx = xlo(1) + dble(i - domlo(1) + HALF)*delta(1) - center(1)

                    adv(i,j,k,URHO) = ambient_state % rho
                    adv(i,j,k,UTEMP) = ambient_state % T
                    adv(i,j,k,UFS:UFS-1+nspec) = ambient_state % rho*ambient_state % xn(:)

                    if ( orbital_kick ) then

                      vx = (-2.0d0 * M_PI / rot_period) * yy
                      vy = ( 2.0d0 * M_PI / rot_period) * xx
                      vz = ZERO

                    else

                      vx = adv(i,j,domlo(3),UMX)/adv(i,j,domlo(3),URHO)
                      vy = adv(i,j,domlo(3),UMY)/adv(i,j,domlo(3),URHO)
                      vz = adv(i,j,domlo(3),UMZ)/adv(i,j,domlo(3),URHO)

                    endif

                    adv(i,j,k,UMX) = ambient_state % rho*vx
                    adv(i,j,k,UMY) = ambient_state % rho*vy
                    adv(i,j,k,UMZ) = ambient_state % rho*vz

                    adv(i,j,k,UEINT) = ambient_state % rho*ambient_state % e
                    adv(i,j,k,UEDEN) = ambient_state % rho*ambient_state % e + &
                         HALF*(adv(i,j,k,UMX)**2 + &
                               adv(i,j,k,UMY)**2 + &
                               adv(i,j,k,UMZ)**2)/ambient_state % rho

                 enddo
              enddo
           enddo
        endif

        ! +z
        if (adv_h3 > domhi(3) .and. bc(3,2,1) .ne. 0) then
           do k = domhi(3)+1, adv_h3
              do j = adv_l2, adv_h2
                 yy = xlo(2) + dble(j - domlo(2) + HALF)*delta(2) - center(2)
                 do i = adv_l1, adv_h1
                    xx = xlo(1) + dble(i - domlo(1) + HALF)*delta(1) - center(1)

                    adv(i,j,k,URHO) = ambient_state % rho
                    adv(i,j,k,UTEMP) = ambient_state % T
                    adv(i,j,k,UFS:UFS-1+nspec) = ambient_state % rho*ambient_state % xn(:)

                    if ( orbital_kick ) then

                      vx = (-2.0d0 * M_PI / rot_period) * yy
                      vy = ( 2.0d0 * M_PI / rot_period) * xx
                      vz = ZERO

                    else

                      vx = adv(i,j,domhi(3),UMX)/adv(i,j,domhi(3),URHO)
                      vy = adv(i,j,domhi(3),UMY)/adv(i,j,domhi(3),URHO)
                      vz = adv(i,j,domhi(3),UMZ)/adv(i,j,domhi(3),URHO)

                    endif

                    adv(i,j,k,UMX) = ambient_state % rho*vx
                    adv(i,j,k,UMY) = ambient_state % rho*vy
                    adv(i,j,k,UMZ) = ambient_state % rho*vz

                    adv(i,j,k,UEINT) = ambient_state % rho*ambient_state % e
                    adv(i,j,k,UEDEN) = ambient_state % rho*ambient_state % e + &
                         HALF*(adv(i,j,k,UMX)**2 + &
                               adv(i,j,k,UMY)**2 + &
                               adv(i,j,k,UMZ)**2)/ambient_state % rho

                 enddo
              enddo
           enddo
        endif

     endif
        
   end subroutine ca_hypfill

   ! ::: -----------------------------------------------------------

   subroutine ca_denfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2, &
                         adv_h3,domlo,domhi,delta,xlo,time,bc)

     implicit none
     include 'bc_types.fi'

     integer adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3
     integer bc(3,2,*)
     integer domlo(3), domhi(3)
     double precision delta(3), xlo(3), time
     double precision adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3)
     logical rho_only
     integer i,j,k

     call filcc(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3,domlo,domhi,delta,xlo,bc)

   end subroutine ca_denfill

   ! ::: -----------------------------------------------------------

   subroutine ca_gravxfill(grav,grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3, &
                           domlo,domhi,delta,xlo,time,bc)

     implicit none
     include 'bc_types.fi'

     integer :: grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3
     integer :: bc(3,2,*)
     integer :: domlo(3), domhi(3)
     double precision delta(3), xlo(3), time
     double precision grav(grav_l1:grav_h1,grav_l2:grav_h2,grav_l3:grav_h3)

     call filcc(grav,grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3, &
                domlo,domhi,delta,xlo,bc)

   end subroutine ca_gravxfill

   ! ::: -----------------------------------------------------------

   subroutine ca_gravyfill(grav,grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3, &
                           domlo,domhi,delta,xlo,time,bc)

     implicit none
     include 'bc_types.fi'

     integer :: grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3
     integer :: bc(3,2,*)
     integer :: domlo(3), domhi(3)
     double precision delta(3), xlo(3), time
     double precision grav(grav_l1:grav_h1,grav_l2:grav_h2,grav_l3:grav_h3)

     call filcc(grav,grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3, &
                domlo,domhi,delta,xlo,bc)

   end subroutine ca_gravyfill

   ! ::: -----------------------------------------------------------

   subroutine ca_gravzfill(grav,grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3, &
                           domlo,domhi,delta,xlo,time,bc)

     implicit none
     include 'bc_types.fi'

     integer :: grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3
     integer :: bc(3,2,*)
     integer :: domlo(3), domhi(3)
     double precision delta(3), xlo(3), time
     double precision grav(grav_l1:grav_h1,grav_l2:grav_h2,grav_l3:grav_h3)

     call filcc(grav,grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3, &
                domlo,domhi,delta,xlo,bc)

   end subroutine ca_gravzfill

   ! ::: -----------------------------------------------------------

   subroutine ca_phigravfill(phi,phi_l1,phi_l2,phi_l3, &
                             phi_h1,phi_h2,phi_h3,domlo,domhi,delta,xlo,time,bc)

     implicit none
     include 'bc_types.fi'

     integer :: phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3
     integer :: bc(3,2,*)
     integer :: domlo(3), domhi(3)
     double precision delta(3), xlo(3), time
     double precision phi(phi_l1:phi_h1,phi_l2:phi_h2,phi_l3:phi_h3)

     call filcc(phi,phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3, &
                domlo,domhi,delta,xlo,bc)

   end subroutine ca_phigravfill

   ! ::: -----------------------------------------------------------

   subroutine ca_reactfill(react,react_l1,react_l2,react_l3, &
                           react_h1,react_h2,react_h3,domlo,domhi,delta,xlo,time,bc)

     implicit none
     include 'bc_types.fi'

     integer :: react_l1,react_l2,react_l3,react_h1,react_h2,react_h3
     integer :: bc(3,2,*)
     integer :: domlo(3), domhi(3)
     double precision delta(3), xlo(3), time
     double precision react(react_l1:react_h1,react_l2:react_h2,react_l3:react_h3)

     call filcc(react,react_l1,react_l2,react_l3,react_h1,react_h2,react_h3, &
                domlo,domhi,delta,xlo,bc)

   end subroutine ca_reactfill
