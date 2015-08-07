
   ! ::: -----------------------------------------------------------

   subroutine ca_hypfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
                         domlo,domhi,delta,xlo,time,bc)

     use meth_params_module, only: NVAR, URHO, UMX, UMY, UMZ, UTEMP, &
          UEDEN, UEINT, UFS, rot_period, do_rotation
     use prob_params_module, only: center
     use bl_constants_module
     use probdata_module
     use eos_module
     use rot_sources_module, only: cross_product, get_omega

     implicit none

     include 'bc_types.fi'

     integer          :: adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3
     integer          :: bc(3,2,*)
     integer          :: domlo(3), domhi(3)
     double precision :: delta(3), xlo(3), time
     double precision :: adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3,NVAR)

     integer          :: i, j, k, n
     double precision :: loc(3), omega(3)

     type (eos_t)     :: ambient_state

     ! Override the generic routine at the physical boundaries by
     ! setting the material to the ambient state.

     if (fill_ambient_bc) then

        call get_ambient(ambient_state)

        omega = get_omega()

        do k = adv_l3, adv_h3
           loc(3) = xlo(3) + dble(j - adv_l3 + HALF)*delta(3) - center(3)
           do j = adv_l2, adv_h2
              loc(2) = xlo(2) + dble(j - adv_l2 + HALF)*delta(2) - center(2)
              do i = adv_l1, adv_h1
                 loc(1) = xlo(1) + dble(i - adv_l1 + HALF)*delta(1) - center(1)

                 if (k .ge. domlo(3) .and. k .le. domhi(3) .and. &
                     j .ge. domlo(2) .and. j .le. domhi(2) .and. &
                     i .ge. domlo(1) .and. i .le. domhi(1)) then
                    cycle
                 endif

                 adv(i,j,k,URHO) = ambient_state % rho
                 adv(i,j,k,UTEMP) = ambient_state % T
                 adv(i,j,k,UFS:UFS-1+nspec) = ambient_state % rho * ambient_state % xn(:)

                 if ( (do_rotation .ne. 1) .and. (.not. no_orbital_kick) ) then

                    adv(i,j,k,UMX:UMZ) = adv(i,j,k,URHO) * cross_product(omega, loc)

                 else

                    if (adv_l1 < domlo(1) .and. bc(1,1,1) .ne. 0) then

                       adv(i,j,k,UMX:UMZ) = adv(domlo(1),j,k,UMX:UMZ)

                    else if (adv_h1 > domhi(1) .and. bc(1,2,1) .ne. 0) then

                       adv(i,j,k,UMX:UMZ) = adv(domhi(1),j,k,UMX:UMZ)

                    else if (adv_l2 < domlo(2) .and. bc(2,1,1) .ne. 0) then

                       adv(i,j,k,UMX:UMZ) = adv(i,domlo(2),k,UMX:UMZ)

                    else if (adv_h2 > domhi(2) .and. bc(2,2,1) .ne. 0) then

                       adv(i,j,k,UMX:UMZ) = adv(i,domhi(2),k,UMX:UMZ)

                    else if (adv_l3 < domlo(3) .and. bc(3,1,1) .ne. 0) then

                       adv(i,j,k,UMX:UMZ) = adv(i,j,domlo(3),UMX:UMZ)

                    else if (adv_h3 > domhi(3) .and. bc(3,2,1) .ne. 0) then

                       adv(i,j,k,UMX:UMZ) = adv(i,j,domhi(3),UMX:UMZ)

                    endif


                 endif

                 adv(i,j,k,UEINT) = ambient_state % rho * ambient_state % e
                 adv(i,j,k,UEDEN) = adv(i,j,k,UEINT) + &
                      HALF*(adv(i,j,k,UMX)**2 + &
                            adv(i,j,k,UMY)**2 + &
                            adv(i,j,k,UMZ)**2) / adv(i,j,k,URHO)

              enddo
           enddo
        enddo

     else

        do n = 1,NVAR
           call filcc(adv(adv_l1,adv_l2,adv_l3,n), &
                      adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
                      domlo,domhi,delta,xlo,bc(1,1,n))
        enddo

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
