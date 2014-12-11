   subroutine PROBINIT (init,name,namlen,problo,probhi)
     
     use probdata_module
     use bl_constants_module
     use fundamental_constants_module
     use meth_params_module, only: small_temp, small_pres, small_dens, rot_period
     use eos_module

     implicit none

     integer :: init, namlen
     integer :: name(namlen)
     double precision :: problo(3), probhi(3)

     integer :: untin
     integer :: i

     namelist /fortin/ &
          denerr,     dengrad,   max_denerr_lev,   max_dengrad_lev, &
          velerr,     velgrad,   max_velerr_lev,   max_velgrad_lev, &
          presserr, pressgrad, max_presserr_lev, max_pressgrad_lev, &
          temperr,   tempgrad,  max_temperr_lev,  max_tempgrad_lev, &
          density, diameter, ambient_dens, ambient_temp, problem

     integer, parameter :: maxlen=127
     character :: probin*(maxlen)
     character :: model*(maxlen)
     integer :: ipp, ierr, ipp1

     ! Temporary storage variables in case we need to switch the primary and secondary.

     integer :: ioproc

     ! For outputting -- determine if we are the IO processor
     call bl_pd_is_ioproc(ioproc)

     ! Build "probin" filename -- the name of file containing fortin namelist.
     if (namlen .gt. maxlen) then
        call bl_error("ERROR: probin file name too long")
     end if

     do i = 1, namlen
        probin(i:i) = char(name(i))
     end do


     ! Set namelist defaults
     denerr = 1.d20
     dengrad = 1.d20
     max_denerr_lev = 10
     max_dengrad_lev = 10

     presserr = 1.d20
     pressgrad = 1.d20
     max_presserr_lev = -1
     max_pressgrad_lev = -1

     velerr  = 1.d0
     velgrad = 1.d20
     max_velerr_lev = -1
     max_velgrad_lev = -1

     temperr  = 1.d0
     tempgrad = 1.d20
     max_temperr_lev = -1
     max_tempgrad_lev = -1

     density  = 1.0d0
     diameter = 1.0d0

     ambient_dens = 1.0d-8
     ambient_temp = 1.0d6

     problem = 1

     ! Read namelists -- override the defaults
     untin = 9 
     open(untin,file=probin(1:namlen),form='formatted',status='old')
     read(untin,fortin)
     close(unit=untin)

     ! Grid geometry
     center(1) = HALF*(problo(1)+probhi(1))
     center(2) = HALF*(problo(2)+probhi(2))
     center(3) = HALF*(problo(3)+probhi(3))

   end subroutine PROBINIT


   ! ::: -----------------------------------------------------------
   ! ::: This routine is called at problem setup time and is used
   ! ::: to initialize data on each grid.  
   ! ::: 
   ! ::: NOTE:  all arrays have one cell of ghost zones surrounding
   ! :::        the grid interior.  Values in these cells need not
   ! :::        be set here.
   ! ::: 
   ! ::: INPUTS/OUTPUTS:
   ! ::: 
   ! ::: level     => amr level of grid
   ! ::: time      => time at which to init data             
   ! ::: lo,hi     => index limits of grid interior (cell centered)
   ! ::: nstate    => number of state components.  You should know
   ! :::		   this already!
   ! ::: state     <=  Scalar array
   ! ::: delta     => cell size
   ! ::: xlo,xhi   => physical locations of lower left and upper
   ! :::              right hand corner of grid.  (does not include
   ! :::		   ghost region).
   ! ::: -----------------------------------------------------------
   subroutine ca_initdata(level,time,lo,hi,nscal, &
                          state,state_l1,state_l2,state_l3,state_h1,state_h2,state_h3, &
                          delta,xlo,xhi)

     use probdata_module
     use interpolate_module
     use eos_module
     use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UTEMP, &
          UEDEN, UEINT, UFS, UFA, rot_period
     use network, only : nspec
     use bl_constants_module
     use fundamental_constants_module
     use multifab_module

     implicit none

     integer :: level, nscal
     integer :: lo(3), hi(3)
     integer :: state_l1,state_l2,state_l3,state_h1,state_h2,state_h3
     double precision :: xlo(3), xhi(3), time, delta(3)
     double precision :: state(state_l1:state_h1,state_l2:state_h2,state_l3:state_h3,NVAR)

     double precision :: xx,yy,zz
     double precision :: c(0:1,0:2), phi, num1, num2, den1, den2
     integer          :: ii, jj, kk, ll

     type (eos_t) :: eos_state

     integer :: i,j,k,n

     !$OMP PARALLEL DO PRIVATE(i, j, k, xx, yy, zz, eos_state, ii, jj, kk, ll, num1, num2, den1, den2, phi, c)
     do k = lo(3), hi(3)   
        zz = xlo(3) + delta(3)*dble(k-lo(3)+HALF) 

        do j = lo(2), hi(2)     
           yy = xlo(2) + delta(2)*dble(j-lo(2)+HALF)

           do i = lo(1), hi(1)   
              xx = xlo(1) + delta(1)*dble(i-lo(1)+HALF)

              ! Establish the cube or sphere

              if (problem .eq. 1) then

                 if (abs(xx) < diameter/2 .and. abs(yy) < diameter/2 .and. abs(zz) < diameter/2) then
                    state(i,j,k,URHO) = density
                 else
                    state(i,j,k,URHO) = ambient_dens
                 endif

              else if (problem .eq. 2) then

                 if ((xx**2 + yy**2 + zz**2)**0.5 < diameter / 2) then
                    state(i,j,k,URHO) = density
                 else
                    state(i,j,k,URHO) = ambient_dens
                 endif

              else

                 call bl_error("Problem not defined.")
                
              endif

              ! Establish the thermodynamic quantities

              state(i,j,k,UTEMP) = ambient_temp

              state(i,j,k,UFS:UFS-1+nspec) = 1.0d0 / nspec

              eos_state % xn  = state(i,j,k,UFS:UFS-1+nspec)
              eos_state % rho = state(i,j,k,URHO)
              eos_state % T   = state(i,j,k,UTEMP)

              call eos(eos_input_rt, eos_state)

              state(i,j,k,UEINT) = eos_state % e

              state(i,j,k,UEDEN) = state(i,j,k,URHO) * state(i,j,k,UEINT)
              state(i,j,k,UEINT) = state(i,j,k,URHO) * state(i,j,k,UEINT)

              do n = 1,nspec
                 state(i,j,k,UFS+n-1) = state(i,j,k,URHO) * state(i,j,k,UFS+n-1)
              end do

              ! Fill in the true phi state

!              c(0,2) = -diameter/2 - zz
!              c(1,2) =  diameter/2 - zz
!              c(0,1) = -diameter/2 - yy
!              c(1,1) =  diameter/2 - yy
!              c(0,0) = -diameter/2 - xx
!              c(1,0) =  diameter/2 - xx

!              phi = 0.0

!              do ii = 0, 1
!                 do jj = 0, 1
!                    do ll = 0, 2

!                       num1 = ( (c(ii,ll)**2 + c(jj,mod(ll+1,3))**2 + c(1,mod(ll+2,3))**2)**(0.5) + c(1,mod(ll+2,3)) )**3
!                       num2 = ( (c(ii,ll)**2 + c(jj,mod(ll+1,3))**2 + c(0,mod(ll+2,3))**2)**(0.5) - c(0,mod(ll+2,3)) )
!                       den1 = ( (c(ii,ll)**2 + c(jj,mod(ll+1,3))**2 + c(1,mod(ll+2,3))**2)**(0.5) - c(1,mod(ll+2,3)) )
!                       den2 = ( (c(ii,ll)**2 + c(jj,mod(ll+1,3))**2 + c(0,mod(ll+2,3))**2)**(0.5) + c(0,mod(ll+2,3)) )**3

!                       phi = phi + 0.5 * (-1)**(ii+jj) * ( c(ii,ll) * c(jj,mod(ll+1,3)) * &
!                                         log( num1 * num2 / (den1 * den2) ) )

!                    enddo
!                 enddo
!              enddo

!              do ii = 0, 1
!                 do jj = 0, 1
!                    do kk = 0, 1
!                       do ll = 0, 2
!                          phi = phi + (-1)**(ii+jj+kk+1) * c(ii,ll)**2 * &
!                                      atan2( c(ii,ll) * c(kk,mod(ll+2,3)), c(ii,ll)**2 + c(jj,mod(ll+1,3))**2 + &
!                                             c(jj,mod(ll+1,3))*(c(ii,ll)**2 + c(jj,mod(ll+1,3))**2 + c(kk,mod(ll+2,3))**2)**(0.5) )
!                       enddo
!                    enddo
!                 enddo
!              enddo

!              state(i,j,k,UFA) = phi * 0.5 * Gconst * rho

           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO

!     true_mass = 4.0d0 / 3.0d0 * M_PI * (diameter / 2.0d0)**3 * rho

!     call parallel_reduce_d(actual_mass, my_mass, MPI_SUM)

!     do k = lo(3), hi(3)
!        do j = lo(2), hi(2)
!           do i = lo(1), hi(1)
!              if (state(i,j,k,URHO) .gt. 2.d-8) then
!                 state(i,j,k,URHO) = state(i,j,k,URHO) * true_mass / actual_mass
!                 state(i,j,k,UFS:UFS+nspec-1) = state(i,j,k,UFS:UFS+nspec-1) * true_mass / actual_mass
!              endif
!           enddo
!        enddo
!     enddo

   end subroutine ca_initdata


   ! ::: -----------------------------------------------------------

   subroutine ca_hypfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
                         domlo,domhi,delta,xlo,time,bc)

     use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UTEMP, &
          UEDEN, UEINT, UFS, rot_period
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

     do n = 1,NVAR
        call filcc(adv(adv_l1,adv_l2,adv_l3,n), &
                   adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
                   domlo,domhi,delta,xlo,bc(1,1,n))
     enddo

     ! override the generic routine at the physical boundaries by
     ! setting the material to the ambient state

     ! -x
     if (adv_l1 < domlo(1)) then
        do k = adv_l3, adv_h3
           do j = adv_l2, adv_h2
              yy = xlo(2) + dble(j - domlo(2) + HALF)*delta(2) - center(2)
              do i = adv_l1, domlo(1)-1
                 xx = xlo(1) + dble(i - domlo(1) + HALF)*delta(1) - center(1)

                 adv(i,j,k,URHO) = dens_ambient
                 adv(i,j,k,UTEMP) = temp_ambient
                 adv(i,j,k,UFS:UFS-1+nspec) = dens_ambient*xn_ambient(:)

                 vx = adv(domlo(1),j,k,UMX)/adv(domlo(1),j,k,URHO)
                 vy = adv(domlo(1),j,k,UMY)/adv(domlo(1),j,k,URHO)
                 vz = adv(domlo(1),j,k,UMZ)/adv(domlo(1),j,k,URHO)

                 adv(i,j,k,UMX) = dens_ambient*vx
                 adv(i,j,k,UMY) = dens_ambient*vy
                 adv(i,j,k,UMZ) = dens_ambient*vz

                 adv(i,j,k,UEINT) = dens_ambient*eint_ambient
                 adv(i,j,k,UEDEN) = dens_ambient*eint_ambient + &
                      HALF*(adv(i,j,k,UMX)**2 + &
                            adv(i,j,k,UMY)**2 + &
                            adv(i,j,k,UMZ)**2)/dens_ambient

              enddo
           enddo
        enddo
     endif

     ! +x
     if (adv_h1 > domhi(1)) then
        do k = adv_l3, adv_h3
           do j = adv_l2, adv_h2
              yy = xlo(2) + dble(j - domlo(2) + HALF)*delta(2) - center(2)
              do i = domhi(1)+1, adv_h1
                 xx = xlo(1) + dble(i - domlo(1) + HALF)*delta(1) - center(1)

                 adv(i,j,k,URHO) = dens_ambient
                 adv(i,j,k,UTEMP) = temp_ambient
                 adv(i,j,k,UFS:UFS-1+nspec) = dens_ambient*xn_ambient(:)

                 vx = adv(domhi(1),j,k,UMX)/adv(domhi(1),j,k,URHO)
                 vy = adv(domhi(1),j,k,UMY)/adv(domhi(1),j,k,URHO)
                 vz = adv(domhi(1),j,k,UMZ)/adv(domhi(1),j,k,URHO)

                 adv(i,j,k,UMX) = dens_ambient*vx
                 adv(i,j,k,UMY) = dens_ambient*vy
                 adv(i,j,k,UMZ) = dens_ambient*vz

                 adv(i,j,k,UEINT) = dens_ambient*eint_ambient
                 adv(i,j,k,UEDEN) = dens_ambient*eint_ambient + &
                      HALF*(adv(i,j,k,UMX)**2 + &
                            adv(i,j,k,UMY)**2 + &
                            adv(i,j,k,UMZ)**2)/dens_ambient

              enddo
           enddo
        enddo
     endif

     ! -y
     if (adv_l2 < domlo(2)) then
        do k = adv_l3, adv_h3
           do j = adv_l2, domlo(2)-1
              yy = xlo(2) + dble(j - domlo(2) + HALF)*delta(2) - center(2)
              do i = adv_l1, adv_h1
                 xx = xlo(1) + dble(i - domlo(1) + HALF)*delta(1) - center(1)

                 adv(i,j,k,URHO) = dens_ambient
                 adv(i,j,k,UTEMP) = temp_ambient
                 adv(i,j,k,UFS:UFS-1+nspec) = dens_ambient*xn_ambient(:)

                 vx = adv(i,domlo(2),k,UMX)/adv(i,domlo(2),k,URHO)
                 vy = adv(i,domlo(2),k,UMY)/adv(i,domlo(2),k,URHO)
                 vz = adv(i,domlo(2),k,UMZ)/adv(i,domlo(2),k,URHO)

                 adv(i,j,k,UMX) = dens_ambient*vx
                 adv(i,j,k,UMY) = dens_ambient*vy
                 adv(i,j,k,UMZ) = dens_ambient*vz

                 adv(i,j,k,UEINT) = dens_ambient*eint_ambient
                 adv(i,j,k,UEDEN) = dens_ambient*eint_ambient + &
                      HALF*(adv(i,j,k,UMX)**2 + &
                            adv(i,j,k,UMY)**2 + &
                            adv(i,j,k,UMZ)**2)/dens_ambient

              enddo
           enddo
        enddo
     endif

     ! +y
     if (adv_h2 > domhi(2)) then
        do k = adv_l3, adv_h3
           do j = domhi(2)+1, adv_h2
              yy = xlo(2) + dble(j - domlo(2) + HALF)*delta(2) - center(2)
              do i = adv_l1, adv_h1
                 xx = xlo(1) + dble(i - domlo(1) + HALF)*delta(1) - center(1)

                 adv(i,j,k,URHO) = dens_ambient
                 adv(i,j,k,UTEMP) = temp_ambient
                 adv(i,j,k,UFS:UFS-1+nspec) = dens_ambient*xn_ambient(:)

                 vx = adv(i,domhi(2),k,UMX)/adv(i,domhi(2),k,URHO)
                 vy = adv(i,domhi(2),k,UMY)/adv(i,domhi(2),k,URHO)
                 vz = adv(i,domhi(2),k,UMZ)/adv(i,domhi(2),k,URHO)

                 adv(i,j,k,UMX) = dens_ambient*vx
                 adv(i,j,k,UMY) = dens_ambient*vy
                 adv(i,j,k,UMZ) = dens_ambient*vz

                 adv(i,j,k,UEINT) = dens_ambient*eint_ambient
                 adv(i,j,k,UEDEN) = dens_ambient*eint_ambient + &
                      HALF*(adv(i,j,k,UMX)**2 + &
                            adv(i,j,k,UMY)**2 + &
                            adv(i,j,k,UMZ)**2)/dens_ambient

              enddo
           enddo
        enddo
     endif

     ! -z
     if (adv_l3 < domlo(3)) then
        do k = adv_l3, domlo(3)-1
           do j = adv_l2, adv_h2
              yy = xlo(2) + dble(j - domlo(2) + HALF)*delta(2) - center(2)
              do i = adv_l1, adv_h1
                 xx = xlo(1) + dble(i - domlo(1) + HALF)*delta(1) - center(1)

                 adv(i,j,k,URHO) = dens_ambient
                 adv(i,j,k,UTEMP) = temp_ambient
                 adv(i,j,k,UFS:UFS-1+nspec) = dens_ambient*xn_ambient(:)

                 vx = adv(i,j,domlo(3),UMX)/adv(i,j,domlo(3),URHO)
                 vy = adv(i,j,domlo(3),UMY)/adv(i,j,domlo(3),URHO)
                 vz = adv(i,j,domlo(3),UMZ)/adv(i,j,domlo(3),URHO)

                 adv(i,j,k,UMX) = dens_ambient*vx
                 adv(i,j,k,UMY) = dens_ambient*vy
                 adv(i,j,k,UMZ) = dens_ambient*vz

                 adv(i,j,k,UEINT) = dens_ambient*eint_ambient
                 adv(i,j,k,UEDEN) = dens_ambient*eint_ambient + &
                      HALF*(adv(i,j,k,UMX)**2 + &
                            adv(i,j,k,UMY)**2 + &
                            adv(i,j,k,UMZ)**2)/dens_ambient

              enddo
           enddo
        enddo
     endif

     ! +z
     if (adv_h3 > domhi(3)) then
        do k = domhi(3)+1, adv_h3
           do j = adv_l2, adv_h2
              yy = xlo(2) + dble(j - domlo(2) + HALF)*delta(2) - center(2)
              do i = adv_l1, adv_h1
                 xx = xlo(1) + dble(i - domlo(1) + HALF)*delta(1) - center(1)

                 adv(i,j,k,URHO) = dens_ambient
                 adv(i,j,k,UTEMP) = temp_ambient
                 adv(i,j,k,UFS:UFS-1+nspec) = dens_ambient*xn_ambient(:)

                 vx = adv(i,j,domhi(3),UMX)/adv(i,j,domhi(3),URHO)
                 vy = adv(i,j,domhi(3),UMY)/adv(i,j,domhi(3),URHO)
                 vz = adv(i,j,domhi(3),UMZ)/adv(i,j,domhi(3),URHO)

                 adv(i,j,k,UMX) = dens_ambient*vx
                 adv(i,j,k,UMY) = dens_ambient*vy
                 adv(i,j,k,UMZ) = dens_ambient*vz

                 adv(i,j,k,UEINT) = dens_ambient*eint_ambient
                 adv(i,j,k,UEDEN) = dens_ambient*eint_ambient + &
                      HALF*(adv(i,j,k,UMX)**2 + &
                            adv(i,j,k,UMY)**2 + &
                            adv(i,j,k,UMZ)**2)/dens_ambient

              enddo
           enddo
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

     use probdata_module
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

     use probdata_module
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

     use probdata_module
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

   subroutine ca_reactfill(react,react_l1,react_l2,react_l3, &
                           react_h1,react_h2,react_h3,domlo,domhi,delta,xlo,time,bc)

     use probdata_module
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
