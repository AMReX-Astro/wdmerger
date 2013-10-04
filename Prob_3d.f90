   subroutine PROBINIT (init,name,namlen,problo,probhi)
     
     use probdata_module
     use model_parser_module
     use bl_constants_module
     use fundamental_constants_module
     use meth_params_module, only: small_temp, small_pres, small_dens, rot_period
     use eos_module
     use com, only: mass_p, mass_s, com_loc_p, com_loc_s

     implicit none

     integer :: init, namlen
     integer :: name(namlen)
     double precision :: problo(3), probhi(3)

     integer :: untin
     integer :: i

     namelist /fortin/ &
          model_P_name, model_S_name, &
          nsub, &
          inertial, &
          interp_temp, &
          damping, damping_alpha, &
          denerr,     dengrad,   max_denerr_lev,   max_dengrad_lev, &
          velerr,     velgrad,   max_velerr_lev,   max_velgrad_lev, &
          presserr, pressgrad, max_presserr_lev, max_pressgrad_lev, &
          temperr,   tempgrad,  max_temperr_lev,  max_tempgrad_lev

     integer, parameter :: maxlen=127
     character :: probin*(maxlen)
     character :: model*(maxlen)
     integer :: ipp, ierr, ipp1

     double precision :: r_l, r_r
     double precision :: a
     double precision :: length

     ! Temporary storage variables in case we need to switch the primary and secondary.

     double precision, allocatable :: temp_model_state(:,:), temp_model_r(:)
     double precision   :: temp_mass
     integer            :: temp_npts
     character (len=80) :: temp_name

     integer :: ioproc

     if (init == 0) return

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

     model_P_name = "P"
     model_S_name = "S"

     nsub = 1

     inertial = .false.

     interp_temp = .false.

     ! Read namelists -- override the defaults
     untin = 9 
     open(untin,file=probin(1:namlen),form='formatted',status='old')
     read(untin,fortin)
     close(unit=untin)


     ! Grid geometry
     center(1) = HALF*(problo(1)+probhi(1))
     center(2) = HALF*(problo(2)+probhi(2))
     center(3) = HALF*(problo(3)+probhi(3))

     ! Read in model for the primary WD
     call read_model_file(model_P_name)
     
     ! Copy the data over to P's arrays
     allocate(model_P_r(npts_model))
     allocate(model_P_state(npts_model, nvars_model))

     model_P_r(:) = model_r(:)
     model_P_state(:,:) = model_state(:,:)
     npts_model_P = npts_model

     call close_model_file()

     ! Compute the mass of the primary
     mass_p_initial = ZERO
     do i = 1, npts_model_P-1
        if (i == 1) then
           r_l = ZERO
        else
           r_l = HALF*(model_P_r(i) + model_P_r(i-1))
        endif

        r_r = HALF*(model_P_r(i) + model_P_r(i+1))

        mass_p_initial = mass_p_initial + FOUR3RD*M_PI*(r_r - r_l)* &
             (r_r**2 + r_l*r_r + r_l**2)*model_P_state(i,idens_model)

     enddo


     ! Read in model for the secondary WD
     call read_model_file(model_S_name)
     
     ! Copy the data over to B's arrays
     allocate(model_S_r(npts_model))
     allocate(model_S_state(npts_model, nvars_model))

     model_S_r(:) = model_r(:)
     model_S_state(:,:) = model_state(:,:)
     npts_model_S = npts_model

     call close_model_file()

     ! Compute the mass of the secondary WD.

     mass_s_initial = ZERO
     do i = 1, npts_model_S-1
        if (i == 1) then
           r_l = ZERO
        else
           r_l = HALF*(model_S_r(i) + model_S_r(i-1))
        endif

        r_r = HALF*(model_S_r(i) + model_S_r(i+1))

        mass_s_initial = mass_s_initial + FOUR3RD*M_PI*(r_r - r_l)* &
             (r_r**2 + r_l*r_r + r_l**2)*model_S_state(i,idens_model)

     enddo

     ! We want the primary WD to be more massive. If what we're calling
     ! the primary ends up being less massive, switch the stars and all
     ! the relevant data.

     if ( mass_p < mass_s ) then

       if (ioproc == 1) then
         print *, "Primary mass is less than secondary mass; switching the stars."
       endif

       allocate(temp_model_r(npts_model_p))
       allocate(temp_model_state(npts_model_p,nvars_model))

       temp_model_state(:,:) = model_P_state(:,:)
       temp_model_r(:)       = model_P_r(:)
 
       model_P_state(:,:) = model_S_state(:,:)
       model_P_r(:)       = model_S_r(:)

       model_S_state(:,:) = temp_model_state(:,:)
       model_S_r(:)       = temp_model_r(:)

       deallocate(temp_model_r)
       deallocate(temp_model_state)

       temp_mass = mass_p
       mass_p_initial = mass_s_initial
       mass_s_initial = temp_mass

       temp_npts = npts_model_p
       npts_model_p = npts_model_s
       npts_model_s = temp_npts

       temp_name = model_p_name
       model_p_name = model_s_name
       model_s_name = temp_name

     endif
     
     if (ioproc == 1) then
        print *, "Mass of the primary = ",   mass_p_initial
        print *, "Mass of the secondary = ", mass_s_initial
     endif


     ! Check that the cutoff densities match
     if (model_S_state(npts_model_S,idens_model) /= &
         model_P_state(npts_model_P,idens_model)) then
        call bl_error("ERROR: Cutoff densities for primary and secondary models do not agree.")
     endif

     dens_ambient = model_S_state(npts_model_S,idens_model)

     ! Check that the cutoff temperature match
     if (model_S_state(npts_model_S,itemp_model) /= &
         model_P_state(npts_model_P,itemp_model)) then
        call bl_error("ERROR: Cutoff temperatures for primary and secondary models do not agree.")
     endif

     temp_ambient = model_S_state(npts_model_S,itemp_model)


     ! Ambient X
     xn_ambient(:) = model_S_state(npts_model_S,ispec_model:ispec_model-1+nspec)


     ! Given the inputs of small_dens and small_temp, figure out small_pres.
     ! Use the ambient gas composition (we don't need to worry about saving
     ! our result in eint_ambient since it will be reset in the next call).

     call eos_given_RTX(eint_ambient,small_pres,small_dens, &
                        small_temp,xn_ambient)

     ! Get the rest of the ambient thermodynamic state
     call eos_given_RTX(eint_ambient,pres_ambient,dens_ambient, &
                        temp_ambient,xn_ambient)

     ! Compute the radius of the primary
     radius_P_initial = -1.0d0
     do i = 1, npts_model_P
        if (model_P_state(i,idens_model) <= 1.1*dens_ambient) then
           radius_P_initial = model_P_r(i)
           exit
        endif
     enddo

     if (radius_P_initial < 0.0d0) then
        call bl_error("ERROR: unable to compute radius of the primary")
     endif

     ! compute the radius of the secondary
     radius_S_initial = -1.0d0
     do i = 1, npts_model_S
        if (model_S_state(i,idens_model) <= 1.1*dens_ambient) then
           radius_S_initial = model_S_r(i)
           exit
        endif
     enddo

     if (radius_S_initial < ZERO) then
        call bl_error("ERROR: unable to compute radius of the secondary")
     endif

     
     ! Orbital properties

     ! Kepler's third law tells us the separation
     a = (Gconst*(mass_P_initial + mass_s_initial)*rot_period**2/(FOUR*M_PI**2))**THIRD

     ! Semi-major axes
     a_S_initial = a/(ONE + mass_S_initial/mass_P_initial)
     a_P_initial = (mass_S_initial/mass_P_initial)*a_S_initial

     if (ioproc == 1) then
        print *, "a = ", a
        print *, "(a_P, a_S) = ", a_P_initial, a_S_initial
     endif

     ! Make sure the stars are not touching.
     if (radius_P_initial + radius_S_initial > a) then
        call bl_error("ERROR: Stars are touching!")
     endif

     ! Make sure the domain is big enough

     ! For simplicity, make sure the x and y sizes are 2x the greatest a + R
     length = max(a_P_initial + radius_P_initial, a_S_initial + radius_S_initial)
     if (length > HALF*(probhi(1) - problo(1)) .or. &
         length > HALF*(probhi(2) - problo(2))) then
        call bl_error("ERROR: The domain width is too small to include the stars (along their axis).")
     endif

     ! For the height, let's take 2x the radius
     if (TWO*max(radius_P_initial, radius_S_initial) > HALF*(probhi(3) - problo(3))) then
        call bl_error("ERROR: The domain height is too small to include the stars (perpendicular to their axis).")
     endif

     
     ! Star center positions -- we'll put them in the midplane on the
     ! x-axis, with the CM at the center of the domain
     center_P_initial(1) = center(1) - a_P_initial
     center_P_initial(2) = center(2)
     center_P_initial(3) = center(3)

     center_S_initial(1) = center(1) + a_S_initial
     center_S_initial(2) = center(2)
     center_S_initial(3) = center(3)

     ! Initialize COM quantities

     mass_p = mass_P_initial
     mass_s = mass_S_initial
     com_loc_p = center_P_initial
     com_loc_s = center_S_initial
     
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

     integer :: i,j,k,ii,jj,kk,n

     !$OMP PARALLEL DO PRIVATE(i, j, k, xl, yl, zl, xx, yy, zz) &
     !$OMP PRIVATE(pres_zone, temp_zone, dist_P, dist_S)
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
                               interpolate(dist_P,npts_model_P, &
                                           model_P_r,model_P_state(:,idens_model))

                          if (interp_temp) then
                             temp_zone = temp_zone + &
                                  interpolate(dist_P,npts_model_P, &
                                              model_P_r,model_P_state(:,itemp_model))

                          else
                             pres_zone = pres_zone + &
                                  interpolate(dist_P,npts_model_P, &
                                              model_P_r,model_P_state(:,ipres_model))

                          endif

                          do n = 1, nspec
                             state(i,j,k,UFS-1+n) = state(i,j,k,UFS-1+n) + &
                                  interpolate(dist_P,npts_model_P, &
                                              model_P_r,model_P_state(:,ispec_model-1+n))
                          enddo


                       ! Are we inside the secondary WD?
                       else if (dist_S < radius_S_initial) then

                          state(i,j,k,URHO) = state(i,j,k,URHO) + &
                               interpolate(dist_S,npts_model_S, &
                                           model_S_r,model_S_state(:,idens_model))

                          if (interp_temp) then
                             temp_zone = temp_zone + &
                                  interpolate(dist_S,npts_model_S, &
                                              model_S_r,model_S_state(:,itemp_model))

                          else
                             pres_zone = pres_zone + &
                                  interpolate(dist_S,npts_model_S, &
                                              model_S_r,model_S_state(:,ipres_model))

                          endif

                          do n = 1, nspec
                             state(i,j,k,UFS-1+n) = state(i,j,k,UFS-1+n) + &
                                  interpolate(dist_S,npts_model_S, &
                                              model_S_r,model_S_state(:,ispec_model-1+n))
                          enddo


                       ! ambient medium
                       else
                          state(i,j,k,URHO) = state(i,j,k,URHO) + dens_ambient
                          if (interp_temp) then
                             temp_zone = temp_zone + temp_ambient
                          else
                             pres_zone = pres_zone + pres_ambient
                          endif
                          state(i,j,k,UFS:UFS-1+nspec) = state(i,j,k,UFS:UFS-1+nspec) + xn_ambient(:)

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
              if (interp_temp) then
                 state(i,j,k,UTEMP) = temp_zone
                 call eos_given_RTX(state(i,j,k,UEINT),pres_zone,state(i,j,k,URHO),state(i,j,k,UTEMP), &
                                      state(i,j,k,UFS:UFS-1+nspec))
              else

                 call eos_e_given_RPX(state(i,j,k,UEINT),state(i,j,k,UTEMP),state(i,j,k,URHO), &
                                      pres_zone,state(i,j,k,UFS:UFS-1+nspec))

              endif

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

                 adv(i,j,k,URHO) = dens_ambient
                 adv(i,j,k,UTEMP) = temp_ambient
                 adv(i,j,k,UFS:UFS-1+nspec) = dens_ambient*xn_ambient(:)

                 if ( inertial ) then

                   vx = (-2.0d0 * M_PI / rot_period) * yy
                   vy = ( 2.0d0 * M_PI / rot_period) * xx
                   vz = ZERO

                 else

                   vx = adv(domlo(1),j,k,UMX)/adv(domlo(1),j,k,URHO)
                   vy = adv(domlo(1),j,k,UMY)/adv(domlo(1),j,k,URHO)
                   vz = adv(domlo(1),j,k,UMZ)/adv(domlo(1),j,k,URHO)
 
                 endif

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

                 if ( inertial ) then

                   vx = (-2.0d0 * M_PI / rot_period) * yy
                   vy = ( 2.0d0 * M_PI / rot_period) * xx
                   vz = ZERO

                 else

                   vx = adv(domhi(1),j,k,UMX)/adv(domhi(1),j,k,URHO)
                   vy = adv(domhi(1),j,k,UMY)/adv(domhi(1),j,k,URHO)
                   vz = adv(domhi(1),j,k,UMZ)/adv(domhi(1),j,k,URHO)
 
                 endif

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

                 if ( inertial ) then

                   vx = (-2.0d0 * M_PI / rot_period) * yy
                   vy = ( 2.0d0 * M_PI / rot_period) * xx
                   vz = ZERO

                 else

                   vx = adv(i,domlo(2),k,UMX)/adv(i,domlo(2),k,URHO)
                   vy = adv(i,domlo(2),k,UMY)/adv(i,domlo(2),k,URHO)
                   vz = adv(i,domlo(2),k,UMZ)/adv(i,domlo(2),k,URHO)
 
                 endif

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

                 if ( inertial ) then

                   vx = (-2.0d0 * M_PI / rot_period) * yy
                   vy = ( 2.0d0 * M_PI / rot_period) * xx
                   vz = ZERO

                 else

                   vx = adv(i,domhi(2),k,UMX)/adv(i,domhi(2),k,URHO)
                   vy = adv(i,domhi(2),k,UMY)/adv(i,domhi(2),k,URHO)
                   vz = adv(i,domhi(2),k,UMZ)/adv(i,domhi(2),k,URHO)
 
                 endif

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

                 if ( inertial ) then

                   vx = (-2.0d0 * M_PI / rot_period) * yy
                   vy = ( 2.0d0 * M_PI / rot_period) * xx
                   vz = ZERO

                 else

                   vx = adv(i,j,domlo(3),UMX)/adv(i,j,domlo(3),URHO)
                   vy = adv(i,j,domlo(3),UMY)/adv(i,j,domlo(3),URHO)
                   vz = adv(i,j,domlo(3),UMZ)/adv(i,j,domlo(3),URHO)
 
                 endif

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

                 if ( inertial ) then

                   vx = (-2.0d0 * M_PI / rot_period) * yy
                   vy = ( 2.0d0 * M_PI / rot_period) * xx
                   vz = ZERO

                 else

                   vx = adv(i,j,domhi(3),UMX)/adv(i,j,domhi(3),URHO)
                   vy = adv(i,j,domhi(3),UMY)/adv(i,j,domhi(3),URHO)
                   vz = adv(i,j,domhi(3),UMZ)/adv(i,j,domhi(3),URHO)
 
                 endif

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
