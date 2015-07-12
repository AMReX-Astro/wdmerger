   subroutine PROBINIT (init,name,namlen,problo,probhi)
     
     use probdata_module, only: initialize

     implicit none

     integer :: init, namlen
     integer :: name(namlen)
     double precision :: problo(3), probhi(3)

     call initialize(name, namlen, init)

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
     use prob_params_module, only: center
     use eos_module
     use meth_params_module, only: NVAR, URHO, UMX, UMY, UMZ, UTEMP, &
          UEDEN, UEINT, UFS, rot_period
     use network, only: nspec
     use bl_constants_module
     use model_parser_module, only: idens_model, itemp_model, ipres_model, ispec_model
     use initial_model_module, only: interpolate_3d_from_1d

     implicit none

     integer :: level, nscal
     integer :: lo(3), hi(3)
     integer :: state_l1,state_l2,state_l3,state_h1,state_h2,state_h3
     double precision :: xlo(3), xhi(3), time, delta(3)
     double precision :: state(state_l1:state_h1,state_l2:state_h2,state_l3:state_h3,NVAR)

     double precision :: loc(3)
     double precision :: dist_P(3), dist_S(3)

     type (eos_t) :: zone_state, ambient_state

     integer :: i,j,k,ii,jj,kk,n

     ! Loop through the zones and set the zone state depending on whether we are
     ! inside the primary or secondary (in which case interpolate from the respective model)
     ! or if we are in an ambient zone.

     call get_ambient(ambient_state)

     !$OMP PARALLEL DO PRIVATE(i, j, k, loc) &
     !$OMP PRIVATE(dist_P, dist_S, zone_state)
     do k = lo(3), hi(3)
        loc(3) = xlo(3) + delta(3)*dble(k+HALF-lo(3)) 

        do j = lo(2), hi(2)
           loc(2) = xlo(2) + delta(2)*dble(j+HALF-lo(2))

           do i = lo(1), hi(1)
              loc(1) = xlo(1) + delta(1)*dble(i+HALF-lo(1))

              dist_P = loc - center_P_initial
              dist_S = loc - center_S_initial

              if (sum(dist_P**2) < radius_P_initial**2) then
                 call interpolate_3d_from_1d(model_P_r, model_P_state, initial_model_npts, dist_P, delta, zone_state, nsub)
              else if (sum(dist_S**2) < radius_S_initial**2) then
                 call interpolate_3d_from_1d(model_S_r, model_S_state, initial_model_npts, dist_S, delta, zone_state, nsub)
              else
                 zone_state = ambient_state
              endif

              state(i,j,k,URHO)  = zone_state % rho
              state(i,j,k,UTEMP) = zone_state % T
              state(i,j,k,UEINT) = zone_state % e * zone_state % rho
              state(i,j,k,UEDEN) = zone_state % e * zone_state % rho
              state(i,j,k,UFS:UFS+nspec-1) = zone_state % rho * zone_state % xn

           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO

     ! Set the velocities in each direction equal to the bulk
     ! velocity of the system. By default this is zero so that
     ! the system is at rest in our reference frame.

     do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              state(i,j,k,UMX) = state(i,j,k,URHO) * (bulk_velx + smallu)
              state(i,j,k,UMY) = state(i,j,k,URHO) * (bulk_vely + smallu)
              state(i,j,k,UMZ) = state(i,j,k,URHO) * (bulk_velz + smallu)
           enddo
        enddo
     enddo

    !$OMP PARALLEL DO PRIVATE(i, j, k, loc, dist_P, dist_S)
    do k = lo(3), hi(3)
       loc(3) = xlo(3) + dble(k - lo(3) + HALF)*delta(3) - center(3)
       do j = lo(2), hi(2)
          loc(2) = xlo(2) + dble(j - lo(2) + HALF)*delta(2) - center(2)
          do i = lo(1), hi(1)
             loc(1) = xlo(1) + dble(i - lo(1) + HALF)*delta(1) - center(1)

             ! If we're in the inertial reference frame, and we want to provide an
             ! initial orbital kick, use counter-clockwise rigid body rotation.

             if (orbital_kick .and. (.not. single_star)) then

                ! x velocity is -omega * r * sin(theta) == -omega * y
                state(i,j,k,UMX+star_axis-1) = state(i,j,k,UMX+star_axis-1) &
                                             + state(i,j,k,URHO) * (-TWO * M_PI / rot_period) * loc(initial_motion_dir)

                ! y velocity is +omega * r * cos(theta) == +omega * x
                state(i,j,k,UMX+initial_motion_dir-1) = state(i,j,k,UMX+initial_motion_dir-1) &
                                                      + state(i,j,k,URHO) * ( TWO * M_PI / rot_period) * loc(star_axis)

             ! If we want a collision calculation, set the stars in 
             ! motion with the respective free-fall velocities.

             else if (collision) then

                dist_P = loc - center_P_initial
                dist_S = loc - center_S_initial

                if (sum(dist_P**2) < radius_P_initial**2) then
                   state(i,j,k,UMX:UMZ) = state(i,j,k,UMX:UMZ) + vel_P(:) * state(i,j,k,URHO)
                else if (sum(dist_S**2) < radius_S_initial**2) then
                   state(i,j,k,UMX:UMZ) = state(i,j,k,UMX:UMZ) + vel_S(:) * state(i,j,k,URHO)
                endif

             endif

        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

     ! Add corresponding kinetic energy from the system motion

     do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)

              state(i,j,k,UEDEN) = state(i,j,k,UEDEN) + &
                ( state(i,j,k,UMX)**2 + state(i,j,k,UMY)**2 + state(i,j,k,UMZ)**2 ) / &
                ( 2.0 * state(i,j,k,URHO) )

           enddo
        enddo
     enddo

   end subroutine ca_initdata


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
