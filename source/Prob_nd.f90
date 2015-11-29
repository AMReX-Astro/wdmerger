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
                          state,state_lo,state_hi, &
                          dx,xlo,xhi)

     use probdata_module
     use prob_params_module, only: center, dim
     use eos_module
     use meth_params_module, only: NVAR, URHO, UMX, UMY, UMZ, UTEMP, &
          UEDEN, UEINT, UFS, do_rotation
     use network, only: nspec
     use bl_constants_module
     use model_parser_module, only: idens_model, itemp_model, ipres_model, ispec_model
     use initial_model_module, only: interpolate_3d_from_1d
     use rotation_module, only: cross_product, get_omega

     implicit none

     integer          :: level, nscal
     integer          :: lo(3), hi(3)
     integer          :: state_lo(3), state_hi(3)
     double precision :: xlo(3), xhi(3), time, dx(3)
     double precision :: state(state_lo(1):state_hi(1),state_lo(2):state_hi(2),state_lo(3):state_hi(3),NVAR)

     double precision :: loc(3), omega(3)
     double precision :: dist_P, dist_S

     type (eos_t) :: zone_state, ambient_state

     integer :: i, j, k, n

     double precision :: rho_P(model_P % npts), rho_S(model_S % npts)
     double precision :: T_P(model_P % npts), T_S(model_S % npts)
     double precision :: xn_P(model_P % npts, nspec), xn_S(model_S % npts, nspec)
     double precision :: r_P(model_P % npts), r_S(model_S % npts)

     ! Loop through the zones and set the zone state depending on whether we are
     ! inside the primary or secondary (in which case interpolate from the respective model)
     ! or if we are in an ambient zone.

     call get_ambient(ambient_state)

     omega = get_omega(time)

     rho_P = model_P % state(:) % rho
     rho_S = model_S % state(:) % rho

     T_P = model_P % state(:) % T
     T_S = model_S % state(:) % T

     do n = 1, nspec
        xn_P(:,n) = model_P % state(:) % xn(n)
        xn_S(:,n) = model_S % state(:) % xn(n)
     enddo

     r_P = model_P % r
     r_S = model_S % r
     
     !$OMP PARALLEL DO PRIVATE(i, j, k, loc) &
     !$OMP PRIVATE(dist_P, dist_S, zone_state)
     do k = lo(3), hi(3)
        loc(3) = xlo(3) + dx(3) * dble(k + HALF - lo(3)) 

        do j = lo(2), hi(2)
           loc(2) = xlo(2) + dx(2) * dble(j + HALF - lo(2))

           do i = lo(1), hi(1)
              loc(1) = xlo(1) + dx(1) * dble(i + HALF - lo(1))

              dist_P = sum((loc - center_P_initial)**2)**HALF
              dist_S = sum((loc - center_S_initial)**2)**HALF

              if (dist_P < model_P % radius) then
                 call interpolate_3d_from_1d(rho_P, T_P, xn_P, r_P, model_P % npts, loc - center_P_initial, dx, zone_state, nsub)
              else if (dist_S < model_S % radius) then
                 call interpolate_3d_from_1d(rho_S, T_S, xn_S, r_S, model_S % npts, loc - center_S_initial, dx, zone_state, nsub)
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
        loc(3) = xlo(3) + dble(k - lo(3) + HALF)*dx(3) - center(3)
        do j = lo(2), hi(2)
           loc(2) = xlo(2) + dble(j - lo(2) + HALF)*dx(2) - center(2)
           do i = lo(1), hi(1)
              loc(1) = xlo(1) + dble(i - lo(1) + HALF)*dx(1) - center(1)

              ! Add any additional velocity imparted to the stars, usually
              ! from an eccentric orbit or from a collision calculation.
              
              dist_P = sum((loc - center_P_initial)**2)**HALF
              dist_S = sum((loc - center_S_initial)**2)**HALF              
              
              if (dist_P < model_P % radius) then
                 state(i,j,k,UMX:UMZ) = state(i,j,k,UMX:UMZ) + vel_P(:) * state(i,j,k,URHO)
              else if (dist_S < model_S % radius) then
                 state(i,j,k,UMX:UMZ) = state(i,j,k,UMX:UMZ) + vel_S(:) * state(i,j,k,URHO)
              endif

              ! If we're in the inertial reference frame, use rigid body rotation with velocity omega x r.
              ! In 2D we have to be careful, though: the third coordinate is an angular
              ! coordinate, whose unit vector is tangent to the unit circle, so we should
              ! have the same velocity everywhere along that coordinate to begin with.

              if ( (do_rotation .ne. 1) .and. (.not. no_orbital_kick) .and. (.not. collision) ) then

                 state(i,j,k,UMX:UMZ) = state(i,j,k,UMX:UMZ) + state(i,j,k,URHO) * cross_product(omega, loc)
                 if (dim .eq. 2) state(i,j,k,UMZ) = abs(state(i,j,k,UMZ))

              endif

           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO

     ! Add corresponding kinetic energy from the velocity on the grid.

     do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)

              state(i,j,k,UEDEN) = state(i,j,k,UEDEN) + HALF * sum(state(i,j,k,UMX:UMZ)**2) / state(i,j,k,URHO)

           enddo
        enddo
     enddo

   end subroutine ca_initdata