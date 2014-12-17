   subroutine PROBINIT (init,name,namlen,problo,probhi)
     
     use probdata_module
     use model_parser_module, only: idens_model, itemp_model, ipres_model, ispec_model
     use bl_constants_module
     use fundamental_constants_module
     use meth_params_module, only: small_temp, small_pres, small_dens, rot_period
     use network
     use eos_module
     use initial_model_module

     implicit none

     integer :: init, namlen
     integer :: name(namlen)
     double precision :: problo(3), probhi(3)

     integer :: untin
     integer :: i

     namelist /fortin/ &
          mass_p, mass_s, &
          nsub, &
          inertial, &
          interp_temp, &
          damping, damping_alpha, &
          do_relax, relax_tau, &
          denerr,     dengrad,   max_denerr_lev,   max_dengrad_lev, &
          velerr,     velgrad,   max_velerr_lev,   max_velgrad_lev, &
          presserr, pressgrad, max_presserr_lev, max_pressgrad_lev, &
          temperr,   tempgrad,  max_temperr_lev,  max_tempgrad_lev, &
          starBuffer, boundaryBuffer, star_axis

     integer, parameter :: maxlen=127
     character :: probin*(maxlen)
     character :: model*(maxlen)
     integer :: ipp, ierr, ipp1

     double precision :: r_l, r_r
     double precision :: a
     double precision :: length

     type (eos_t) :: eos_state

     ! Temporary storage variables in case we need to switch the primary and secondary.

     double precision, allocatable :: temp_model_state(:,:), temp_model_r(:)
     double precision   :: temp_mass
     integer            :: temp_npts
     character (len=80) :: temp_name

     double precision :: temp_core, xn_core(nspec)
     double precision :: smallx, dx

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

     starBuffer = 1.5d0
     boundaryBuffer = 0.6d0

     nsub = 1

     mass_p = 1.0
     mass_s = 1.0

     inertial = .false.
     interp_temp = .false.
     damping  = .false.
     do_relax = .false.
     star_axis = 1

     ! Read namelists -- override the defaults
     untin = 9 
     open(untin,file=probin(1:namlen),form='formatted',status='old')
     read(untin,fortin)
     close(unit=untin)


     ! Grid geometry
     center(1) = HALF*(problo(1)+probhi(1))
     center(2) = HALF*(problo(2)+probhi(2))
     center(3) = HALF*(problo(3)+probhi(3))

     npts_model = 1024

     dx = 1.0d6

     smallx = 1.0d-10

     ! Define stellar interior quantities

     temp_core = 1.0d7

     xn_core(:) = smallx
     xn_core(2) = HALF - TWO * smallx
     xn_core(3) = HALF - TWO * smallx

     ! Define ambient state and call EOS to get eint and pressure

     ambient_state % rho = 1.0d-4
     ambient_state % T   = 1.0d7
     ambient_state % xn  = xn_core

     call eos(eos_input_rt, ambient_state, .false.)

     ! Allocate arrays to hold the stellar models

     allocate(model_P_state(npts_model,3+nspec))
     allocate(model_S_state(npts_model,3+nspec))

     allocate(model_P_r(npts_model))
     allocate(model_S_r(npts_model))

     ! We want the primary WD to be more massive. If what we're calling
     ! the primary is less massive, switch the stars.

     if ( mass_p < mass_s ) then

       if (ioproc == 1) then
         print *, "Primary mass is less than secondary mass; switching the stars."
       endif

       temp_mass = mass_p
       mass_p = mass_s
       mass_s = temp_mass

     endif

     ! Generate primary and secondary WD

     call init_1d(model_P_r, model_P_state, npts_model, dx, mass_P, radius_P_initial, &
                  temp_core, xn_core, ambient_state)

     call init_1d(model_S_r, model_S_state, npts_model, dx, mass_S, radius_S_initial, &
                  temp_core, xn_core, ambient_state)

     ! Given the inputs of small_dens and small_temp, figure out small_pres.
     ! Use the ambient gas composition (we don't need to worry about saving
     ! our result in eint_ambient since it will be reset in the next call).
 
     eos_state % xn  = ambient_state % xn
     eos_state % rho = small_dens
     eos_state % T   = small_temp
 
     call eos(eos_input_rt, eos_state, .false.)

     small_pres = eos_state % p

     ! Get the orbit from Kepler's third law

     call kepler_third_law(radius_P_initial, mass_P, radius_S_initial, mass_S, &
                           rot_period, a, a_P_initial, a_S_initial, problo, probhi)

     ! Star center positions -- we'll put them in the midplane on the
     ! axis specified by star_axis, with the center of mass at the center of the domain.

     center_P_initial = center
     center_S_initial = center

     center_P_initial(star_axis) = center_P_initial(star_axis) - a_P_initial
     center_S_initial(star_axis) = center_S_initial(star_axis) + a_S_initial

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
          UEDEN, UEINT, UFS, rot_period
     use network, only : nspec
     use bl_constants_module
     use model_parser_module, only: idens_model, itemp_model, ipres_model, ispec_model

     implicit none

     integer :: level, nscal
     integer :: lo(3), hi(3)
     integer :: state_l1,state_l2,state_l3,state_h1,state_h2,state_h3
     double precision :: xlo(3), xhi(3), time, delta(3)
     double precision :: state(state_l1:state_h1,state_l2:state_h2,state_l3:state_h3,NVAR)

     double precision :: xl,yl,zl,xx,yy,zz
     double precision :: pres_zone, temp_zone
     double precision :: dist_P, dist_S

     integer :: pt_index(3)

     type (eos_t) :: eos_state

     integer :: i,j,k,ii,jj,kk,n

     !$OMP PARALLEL DO PRIVATE(i, j, k, xl, yl, zl, xx, yy, zz) &
     !$OMP PRIVATE(ii, jj, kk, n) &
     !$OMP PRIVATE(pres_zone, temp_zone, dist_P, dist_S, eos_state, pt_index)
     do k = lo(3), hi(3)   
        zl = xlo(3) + delta(3)*dble(k-lo(3)) 

        do j = lo(2), hi(2)     
           yl = xlo(2) + delta(2)*dble(j-lo(2))

           do i = lo(1), hi(1)   
              xl = xlo(1) + delta(1)*dble(i-lo(1))

              state(i,j,k,URHO) = 0.0d0
              pres_zone = 0.0d0
              temp_zone = 0.0d0
              state(i,j,k,UFS:UFS-1+nspec) = 0.0d0

              do kk = 0, nsub-1
                 zz = zl + dble(kk + HALF)*delta(3)/nsub

                 do jj = 0, nsub-1
                    yy = yl + dble(jj + HALF)*delta(2)/nsub

                    do ii = 0, nsub-1
                       xx = xl + dble(ii + HALF)*delta(1)/nsub

                       dist_P = sqrt(( xx - center_P_initial(1) )**2 + &
                                     ( yy - center_P_initial(2) )**2 + &
                                     ( zz - center_P_initial(3) )**2)

                       dist_S = sqrt(( xx - center_S_initial(1) )**2 + &
                                     ( yy - center_S_initial(2) )**2 + &
                                     ( zz - center_S_initial(3) )**2)

                       ! Are we inside the primary WD?
                       if (dist_P < radius_P_initial) then

                          state(i,j,k,URHO) = state(i,j,k,URHO) + &
                               interpolate(dist_P,npts_model, &
                                           model_P_r,model_P_state(:,idens_model))

                          if (interp_temp) then
                             temp_zone = temp_zone + &
                                  interpolate(dist_P,npts_model, &
                                              model_P_r,model_P_state(:,itemp_model))

                          else
                             pres_zone = pres_zone + &
                                  interpolate(dist_P,npts_model, &
                                              model_P_r,model_P_state(:,ipres_model))

                          endif

                          do n = 1, nspec
                             state(i,j,k,UFS-1+n) = state(i,j,k,UFS-1+n) + &
                                  interpolate(dist_P,npts_model, &
                                              model_P_r,model_P_state(:,ispec_model-1+n))
                          enddo


                       ! Are we inside the secondary WD?
                       else if (dist_S < radius_S_initial) then

                          state(i,j,k,URHO) = state(i,j,k,URHO) + &
                               interpolate(dist_S,npts_model, &
                                           model_S_r,model_S_state(:,idens_model))

                          if (interp_temp) then
                             temp_zone = temp_zone + &
                                  interpolate(dist_S,npts_model, &
                                              model_S_r,model_S_state(:,itemp_model))

                          else
                             pres_zone = pres_zone + &
                                  interpolate(dist_S,npts_model, &
                                              model_S_r,model_S_state(:,ipres_model))

                          endif

                          do n = 1, nspec
                             state(i,j,k,UFS-1+n) = state(i,j,k,UFS-1+n) + &
                                  interpolate(dist_S,npts_model, &
                                              model_S_r,model_S_state(:,ispec_model-1+n))
                          enddo


                       ! ambient medium
                       else
                          state(i,j,k,URHO) = state(i,j,k,URHO) + ambient_state % rho
                          if (interp_temp) then
                             temp_zone = temp_zone + ambient_state % T
                          else
                             pres_zone = pres_zone + ambient_state % p
                          endif
                          state(i,j,k,UFS:UFS-1+nspec) = state(i,j,k,UFS:UFS-1+nspec) + ambient_state % xn(:)

                       endif

                    enddo
                 enddo
              enddo


              ! normalize
              state(i,j,k,URHO) = state(i,j,k,URHO)/(nsub*nsub*nsub)
              pres_zone = pres_zone/(nsub*nsub*nsub)
              temp_zone = temp_zone/(nsub*nsub*nsub)
              state(i,j,k,UFS:UFS-1+nspec) = state(i,j,k,UFS:UFS-1+nspec)/(nsub*nsub*nsub)

              ! thermodynamics              

              eos_state % xn  = state(i,j,k,UFS:UFS-1+nspec)
              eos_state % p   = pres_zone
              eos_state % rho = state(i,j,k,URHO)
              eos_state % T   = temp_zone

              pt_index(1) = i
              pt_index(2) = j
              pt_index(3) = k

              if (interp_temp) then
                 state(i,j,k,UTEMP) = eos_state % T
                 call eos(eos_input_rt, eos_state, .false., pt_index = pt_index)
              else
                 call eos(eos_input_rp, eos_state, .false., pt_index = pt_index)
                 state(i,j,k,UTEMP) = eos_state % T
              endif

              state(i,j,k,UEINT) = eos_state % e

              state(i,j,k,UEDEN) = state(i,j,k,URHO) * state(i,j,k,UEINT)
              state(i,j,k,UEINT) = state(i,j,k,URHO) * state(i,j,k,UEINT)

              do n = 1,nspec
                 state(i,j,k,UFS+n-1) = state(i,j,k,URHO) * state(i,j,k,UFS+n-1)
              end do

           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO

     ! Initial velocities = 0

     state(:,:,:,UMX:UMZ) = ZERO

     ! If we're in the inertial reference frame, 
     ! set counter-clockwise rigid body rotation

     if ( inertial ) then
 
       !$OMP PARALLEL DO PRIVATE(i, j, k, xx, yy)
       do k = lo(3), hi(3)

         do j = lo(2), hi(2)
           yy = xlo(2) + dble(j - lo(2) + HALF)*delta(2) - center(2)
             
           do i = lo(1), hi(1)
             xx = xlo(1) + dble(i - lo(1) + HALF)*delta(1) - center(1)

             ! x velocity is -omega * r * sin(theta) == -omega * y

             state(i,j,k,UMX) = state(i,j,k,URHO) * (-2.0d0 * M_PI / rot_period) * yy
           
             ! y velocity is +omega * r * cos(theta) == +omega * x

             state(i,j,k,UMY) = state(i,j,k,URHO) * ( 2.0d0 * M_PI / rot_period) * xx

             ! Add corresponding kinetic energy

             state(i,j,k,UEDEN) = state(i,j,k,UEDEN) + &
               ( state(i,j,k,UMX)**2 + state(i,j,k,UMY)**2 ) / &
               ( 2.0 * state(i,j,k,URHO) )

           enddo

         enddo
       
       enddo
       !$OMP END PARALLEL DO

     endif

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

                 adv(i,j,k,URHO) = ambient_state % rho
                 adv(i,j,k,UTEMP) = ambient_state % T
                 adv(i,j,k,UFS:UFS-1+nspec) = ambient_state % rho * ambient_state % xn(:)

                 if ( inertial ) then

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
     if (adv_h1 > domhi(1)) then
        do k = adv_l3, adv_h3
           do j = adv_l2, adv_h2
              yy = xlo(2) + dble(j - domlo(2) + HALF)*delta(2) - center(2)
              do i = domhi(1)+1, adv_h1
                 xx = xlo(1) + dble(i - domlo(1) + HALF)*delta(1) - center(1)

                 adv(i,j,k,URHO) = ambient_state % rho
                 adv(i,j,k,UTEMP) = ambient_state % T
                 adv(i,j,k,UFS:UFS-1+nspec) = ambient_state % rho*ambient_state % xn(:)

                 if ( inertial ) then

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
     if (adv_l2 < domlo(2)) then
        do k = adv_l3, adv_h3
           do j = adv_l2, domlo(2)-1
              yy = xlo(2) + dble(j - domlo(2) + HALF)*delta(2) - center(2)
              do i = adv_l1, adv_h1
                 xx = xlo(1) + dble(i - domlo(1) + HALF)*delta(1) - center(1)

                 adv(i,j,k,URHO) = ambient_state % rho
                 adv(i,j,k,UTEMP) = ambient_state % T
                 adv(i,j,k,UFS:UFS-1+nspec) = ambient_state % rho*ambient_state % xn(:)

                 if ( inertial ) then

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
     if (adv_h2 > domhi(2)) then
        do k = adv_l3, adv_h3
           do j = domhi(2)+1, adv_h2
              yy = xlo(2) + dble(j - domlo(2) + HALF)*delta(2) - center(2)
              do i = adv_l1, adv_h1
                 xx = xlo(1) + dble(i - domlo(1) + HALF)*delta(1) - center(1)

                 adv(i,j,k,URHO) = ambient_state % rho
                 adv(i,j,k,UTEMP) = ambient_state % T
                 adv(i,j,k,UFS:UFS-1+nspec) = ambient_state % rho*ambient_state % xn(:)

                 if ( inertial ) then

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
     if (adv_l3 < domlo(3)) then
        do k = adv_l3, domlo(3)-1
           do j = adv_l2, adv_h2
              yy = xlo(2) + dble(j - domlo(2) + HALF)*delta(2) - center(2)
              do i = adv_l1, adv_h1
                 xx = xlo(1) + dble(i - domlo(1) + HALF)*delta(1) - center(1)

                 adv(i,j,k,URHO) = ambient_state % rho
                 adv(i,j,k,UTEMP) = ambient_state % T
                 adv(i,j,k,UFS:UFS-1+nspec) = ambient_state % rho*ambient_state % xn(:)

                 if ( inertial ) then

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
     if (adv_h3 > domhi(3)) then
        do k = domhi(3)+1, adv_h3
           do j = adv_l2, adv_h2
              yy = xlo(2) + dble(j - domlo(2) + HALF)*delta(2) - center(2)
              do i = adv_l1, adv_h1
                 xx = xlo(1) + dble(i - domlo(1) + HALF)*delta(1) - center(1)

                 adv(i,j,k,URHO) = ambient_state % rho
                 adv(i,j,k,UTEMP) = ambient_state % T
                 adv(i,j,k,UFS:UFS-1+nspec) = ambient_state % rho*ambient_state % xn(:)

                 if ( inertial ) then

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
