   subroutine PROBINIT (init,name,namlen,problo,probhi)
     
     use probdata_module
     use model_parser_module
     use bl_constants_module
     use fundamental_constants_module
     use eos_module
     use com

     implicit none

     integer :: init, namlen
     integer :: name(namlen)
     double precision :: problo(3), probhi(3)

     integer :: untin
     integer :: i

     namelist /fortin/ &
          model_A_name, model_B_name, &
          period, &
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

     integer :: ioproc

     ! for outputting -- determine if we are the IO processor
     call bl_pd_is_ioproc(ioproc)

     ! build "probin" filename -- the name of file containing fortin namelist.
     if (namlen .gt. maxlen) then
        call bl_error("ERROR: probin file name too long")
     end if

     do i = 1, namlen
        probin(i:i) = char(name(i))
     end do


     ! set namelist defaults
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

     period = 1.0
     model_A_name = "A"
     model_B_name = "B"


     ! read namelists -- override the defaults
     untin = 9 
     open(untin,file=probin(1:namlen),form='formatted',status='old')
     read(untin,fortin)
     close(unit=untin)


     ! grid geometry
     center(1) = HALF*(problo(1)+probhi(1))
     center(2) = HALF*(problo(2)+probhi(2))
     center(3) = HALF*(problo(3)+probhi(3))

     
     ! read in model A
     call read_model_file(model_A_name)
     
     ! copy the data over to A's arrays
     allocate(model_A_r(npts_model))
     allocate(model_A_state(npts_model, nvars_model))

     model_A_r(:) = model_r(:)
     model_A_state(:,:) = model_state(:,:)
     npts_model_A = npts_model

     call close_model_file()

     ! compute the mass of model A
     mass_A = ZERO
     do i = 1, npts_model_A-1
        if (i == 1) then
           r_l = ZERO
        else
           r_l = HALF*(model_A_r(i) + model_A_r(i-1))
        endif

        r_r = HALF*(model_A_r(i) + model_A_r(i+1))

        mass_A = mass_A + FOUR3RD*M_PI*(r_r - r_l)* &
             (r_r**2 + r_l*r_r + r_l**2)*model_A_state(i,idens_model)

     enddo


     ! read in model B
     call read_model_file(model_B_name)
     
     ! copy the data over to B's arrays
     allocate(model_B_r(npts_model))
     allocate(model_B_state(npts_model, nvars_model))

     model_B_r(:) = model_r(:)
     model_B_state(:,:) = model_state(:,:)
     npts_model_B = npts_model

     call close_model_file()

     ! compute the mass of model B
     mass_B = ZERO
     do i = 1, npts_model_B-1
        if (i == 1) then
           r_l = ZERO
        else
           r_l = HALF*(model_B_r(i) + model_B_r(i-1))
        endif

        r_r = HALF*(model_B_r(i) + model_B_r(i+1))

        mass_B = mass_B + FOUR3RD*M_PI*(r_r - r_l)* &
             (r_r**2 + r_l*r_r + r_l**2)*model_B_state(i,idens_model)

     enddo

     if (ioproc == 1) then
        print *, "mass of model A = ", mass_A
        print *, "mass of model B = ", mass_B
     endif


     ! check that the cutoff densities match
     if (model_B_state(npts_model_B,idens_model) /= &
         model_A_state(npts_model_A,idens_model)) then
        call bl_error("ERROR: cutoff densities for models A and B don't agree")
     endif

     dens_ambient = model_B_state(npts_model_B,idens_model)

     ! check that the cutoff temperature match
     if (model_B_state(npts_model_B,itemp_model) /= &
         model_A_state(npts_model_A,itemp_model)) then
        call bl_error("ERROR: cutoff temps for models A and B don't agree")
     endif

     temp_ambient = model_B_state(npts_model_B,itemp_model)


     ! ambient X
     xn_ambient(:) = model_B_state(npts_model_B,ispec_model:ispec_model-1+nspec)


     ! get the rest of the ambient thermodynamic state
     call eos_given_RTX(eint_ambient,pres_ambient,dens_ambient, &
                        temp_ambient,xn_ambient)


     ! compute the radius of model A
     radius_A = -1.0d0
     do i = 1, npts_model_A
        if (model_A_state(i,idens_model) <= 1.1*dens_ambient) then
           radius_A = model_A_r(i)
           exit
        endif
     enddo

     if (radius_A < 0.0d0) then
        call bl_error("ERROR: unable to compute radius of model A")
     endif

     ! compute the radius of model B
     radius_B= -1.0d0
     do i = 1, npts_model_B
        if (model_B_state(i,idens_model) <= 1.1*dens_ambient) then
           radius_B = model_B_r(i)
           exit
        endif
     enddo

     if (radius_B < ZERO) then
        call bl_error("ERROR: unable to compute radius of model B")
     endif

     
     ! orbital properties
     ! Kepler tells us the separation
     a = (Gconst*(mass_A + mass_B)*period**2/(FOUR*M_PI**2))**THIRD

     ! semi-major axes
     a_B = a/(ONE + mass_B/mass_A)
     a_A = (mass_B/mass_A)*a_B

     if (ioproc == 1) then
        print *, "a = ", a
        print *, "(a_A, a_B) = ", a_A, a_B
     endif

     ! make sure the stars are not touching
     if (radius_A + radius_B > a) then
        call bl_error("ERROR: stars are touching!")
     endif


     ! make sure the domain is big enough

     ! for simplicity, make sure the x and y sizes are 2x the greatest a + R
     length = max(a_A + radius_A, a_B + radius_B)
     if (length > HALF*(probhi(1) - problo(1)) .or. &
         length > HALF*(probhi(2) - problo(2))) then
        call bl_error("ERROR: domain width too small", length)
     endif

     ! for the height, let's take 2x the radius
     if (TWO*max(radius_A,radius_B) > HALF*(probhi(3) - problo(3))) then
        call bl_error("ERROR: domain height too small")
     endif

     
     ! star center positions -- we'll put them in the midplane on the
     ! x-axis, with the CM at the center of the domain
     x_cen_A = center(1) - a_A
     y_cen_A = center(2)
     z_cen_A = center(3)

     x_cen_B = center(1) + a_B
     y_cen_B = center(2)
     z_cen_B = center(3)

     ! Initialize COM quantities

     call com_save(mass_A,mass_B,x_cen_A,x_cen_B,y_cen_A,y_cen_B,z_cen_A,z_cen_B)
     

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
          UEDEN, UEINT, UFS
     use network, only : nspec
     use bl_constants_module
     use model_parser_module, only: idens_model, itemp_model, ispec_model

     implicit none

     integer :: level, nscal
     integer :: lo(3), hi(3)
     integer :: state_l1,state_l2,state_l3,state_h1,state_h2,state_h3
     double precision :: xlo(3), xhi(3), time, delta(3)
     double precision :: state(state_l1:state_h1,state_l2:state_h2,state_l3:state_h3,NVAR)

     double precision :: x,y,z
     double precision :: pres
     double precision :: dist_A, dist_B

     integer :: i,j,k,n


     do k = lo(3), hi(3)   
        z = xlo(3) + delta(3)*(dble(k-lo(3)) + HALF) 

        do j = lo(2), hi(2)     
           y = xlo(2) + delta(2)*(dble(j-lo(2)) + HALF) 

           do i = lo(1), hi(1)   
              x = xlo(1) + delta(1)*(dble(i-lo(1)) + HALF) 

              dist_A = sqrt((x - x_cen_A)**2 + &
                            (y - y_cen_A)**2 + &
                            (z - z_cen_A)**2)

              dist_B = sqrt((x - x_cen_B)**2 + &
                            (y - y_cen_B)**2 + &
                            (z - z_cen_B)**2)

              ! are we "inside" star A?
              if (dist_A < radius_A) then

                 state(i,j,k,URHO) = &
                      interpolate(dist_A,npts_model_A, &
                                  model_A_r,model_A_state(:,idens_model))

                 state(i,j,k,UTEMP) = &
                      interpolate(dist_A,npts_model_A, &
                                  model_A_r,model_A_state(:,itemp_model))

                 do n = 1, nspec
                    state(i,j,k,UFS-1+n) = &
                         interpolate(dist_A,npts_model_A, &
                                     model_A_r,model_A_state(:,ispec_model-1+n))
                 enddo


              ! inside B?
              else if (dist_B < radius_B) then

                 state(i,j,k,URHO) = &
                      interpolate(dist_B,npts_model_B, &
                                  model_B_r,model_B_state(:,idens_model))

                 state(i,j,k,UTEMP) = &
                      interpolate(dist_B,npts_model_B, &
                                  model_B_r,model_B_state(:,itemp_model))

                 do n = 1, nspec
                    state(i,j,k,UFS-1+n) = &
                         interpolate(dist_B,npts_model_B, &
                                     model_B_r,model_B_state(:,ispec_model-1+n))
                 enddo

              ! ambient medium
              else
                 state(i,j,k,URHO)  = dens_ambient
                 state(i,j,k,UTEMP) = temp_ambient
                 state(i,j,k,UFS:UFS-1+nspec) = xn_ambient(:)

              endif

           enddo
        enddo
     enddo

     ! thermodynamics
     do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              call eos_given_RTX(state(i,j,k,UEINT),pres,state(i,j,k,URHO), &
                                 state(i,j,k,UTEMP),state(i,j,k,UFS:UFS-1+nspec))

              state(i,j,k,UEDEN) = state(i,j,k,URHO) * state(i,j,k,UEINT)
              state(i,j,k,UEINT) = state(i,j,k,URHO) * state(i,j,k,UEINT)

              do n = 1,nspec
                 state(i,j,k,UFS+n-1) = state(i,j,k,URHO) * state(i,j,k,UFS+n-1)
              end do

           enddo
        enddo
     enddo

     ! Initial velocities = 0
     state(:,:,:,UMX:UMZ) = ZERO

   end subroutine ca_initdata


   ! ::: -----------------------------------------------------------

   subroutine ca_hypfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
                         domlo,domhi,delta,xlo,time,bc)

     use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UTEMP, &
          UEDEN, UEINT, UFS
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
              do i = adv_l1, domlo(1)-1
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
              do i = domhi(1)+1, adv_h1
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
              do i = adv_l1, adv_h1
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
              do i = adv_l1, adv_h1
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
              do i = adv_l1, adv_h1
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
              do i = adv_l1, adv_h1
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
