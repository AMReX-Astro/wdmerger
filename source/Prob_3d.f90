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
          UEDEN, UEINT, UFS, rot_period, do_rotation
     use network, only: nspec
     use bl_constants_module
     use model_parser_module, only: idens_model, itemp_model, ipres_model, ispec_model
     use initial_model_module, only: interpolate_3d_from_1d
     use rot_sources_module, only: cross_product, get_omega

     implicit none

     integer :: level, nscal
     integer :: lo(3), hi(3)
     integer :: state_l1,state_l2,state_l3,state_h1,state_h2,state_h3
     double precision :: xlo(3), xhi(3), time, delta(3)
     double precision :: state(state_l1:state_h1,state_l2:state_h2,state_l3:state_h3,NVAR)

     double precision :: loc(3), omega(3)
     double precision :: dist_P(3), dist_S(3)

     type (eos_t) :: zone_state, ambient_state

     integer :: i,j,k,ii,jj,kk,n

     ! Loop through the zones and set the zone state depending on whether we are
     ! inside the primary or secondary (in which case interpolate from the respective model)
     ! or if we are in an ambient zone.

     call get_ambient(ambient_state)

     omega = get_omega()

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

             ! If we want a collision calculation, set the stars in 
             ! motion with the respective free-fall velocities.

             if (collision) then

                dist_P = loc - center_P_initial
                dist_S = loc - center_S_initial

                if (sum(dist_P**2) < radius_P_initial**2) then
                   state(i,j,k,UMX:UMZ) = state(i,j,k,UMX:UMZ) + vel_P(:) * state(i,j,k,URHO)
                else if (sum(dist_S**2) < radius_S_initial**2) then
                   state(i,j,k,UMX:UMZ) = state(i,j,k,UMX:UMZ) + vel_S(:) * state(i,j,k,URHO)
                endif

             ! If we're in the inertial reference frame, and we want to provide an
             ! initial orbital kick, use rigid body rotation with velocity omega x r.

             else if ((do_rotation .ne. 1) .and. (.not. no_orbital_kick) .and. (.not. single_star)) then

                state(i,j,k,UMX:UMZ) = state(i,j,k,UMX:UMZ) + state(i,j,k,URHO) * cross_product(omega, loc)

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
