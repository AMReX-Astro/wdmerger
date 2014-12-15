   subroutine PROBINIT (init,name,namlen,problo,probhi)
     
     use probdata_module
     use model_parser_module
     use bl_constants_module
     use fundamental_constants_module
     use eos_module
     use eos_type_module
     use meth_params_module, only : small_dens, small_temp, small_pres
     implicit none

     integer :: init, namlen
     integer :: name(namlen)
     double precision :: problo(3), probhi(3)

     integer :: untin
     integer :: i

     namelist /fortin/ &
          model_name, &
          denerr,     dengrad,   max_denerr_lev,   max_dengrad_lev, &
          velerr,     velgrad,   max_velerr_lev,   max_velgrad_lev, &
          presserr, pressgrad, max_presserr_lev, max_pressgrad_lev, &
          temperr,   tempgrad,  max_temperr_lev,  max_tempgrad_lev

     integer, parameter :: maxlen=127
     character :: probin*(maxlen)
     character :: model*(maxlen)
     integer :: ipp, ierr, ipp1

     double precision :: r_l, r_r

     type (eos_t) :: eos_state

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

     model_name = "A"


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
     call read_model_file(model_name)
     
     ! copy the data over to A's arrays
     allocate(model_star_r(npts_model))
     allocate(model_star_state(npts_model, nvars_model))

     model_star_r(:) = model_r(:)
     model_star_state(:,:) = model_state(:,:)
     npts_model_star = npts_model

     call close_model_file()

     ! compute the mass of the model
     total_mass = ZERO
     do i = 1, npts_model_star-1
        if (i == 1) then
           r_l = ZERO
        else
           r_l = HALF*(model_star_r(i) + model_star_r(i-1))
        endif

        r_r = HALF*(model_star_r(i) + model_star_r(i+1))

        total_mass = total_mass + FOUR3RD*M_PI*(r_r - r_l)* &
             (r_r**2 + r_l*r_r + r_l**2)*model_star_state(i,idens_model)

     enddo


     if (ioproc == 1) then
        print *, "mass of model = ", total_mass
     endif

     dens_ambient = model_star_state(npts_model_star,idens_model)
     temp_ambient = model_star_state(npts_model_star,itemp_model)
     xn_ambient(:) = model_star_state(npts_model_star,ispec_model:ispec_model-1+nspec)


     ! get the rest of the ambient thermodynamic state
     eos_state%rho = dens_ambient
     eos_state%T   = temp_ambient
     eos_state%xn(:) = xn_ambient(:)

     call eos(eos_input_rt, eos_state)

     pres_ambient = eos_state%p
     eint_ambient = eos_state%e


     ! set small_press
     eos_state%rho = small_dens
     eos_state%T = small_temp
     eos_state%xn(:) = xn_ambient(:)

     call eos(eos_input_rt, eos_state)
     
     small_pres = eos_state%p

     if (ioproc == 1) then
        print *, 'small_pres set to: ', small_pres
     endif

     ! compute the radius of the model
     radius = -1.0d0
     do i = 1, npts_model_star
        if (model_star_state(i,idens_model) <= 1.1*dens_ambient) then
           radius = model_star_r(i)
           exit
        endif
     enddo

     if (radius < 0.0d0) then
        call bl_error("ERROR: unable to compute radius of the model star")
     endif


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
     use eos_type_module
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
     double precision :: dist

     type (eos_t) :: eos_state

     integer :: i,j,k,n


     do k = lo(3), hi(3)   
        z = xlo(3) + delta(3)*(dble(k-lo(3)) + HALF) 

        do j = lo(2), hi(2)     
           y = xlo(2) + delta(2)*(dble(j-lo(2)) + HALF) 

           do i = lo(1), hi(1)   
              x = xlo(1) + delta(1)*(dble(i-lo(1)) + HALF) 

              dist = sqrt((x - center(1))**2 + &
                          (y - center(2))**2 + &
                          (z - center(3))**2)


              state(i,j,k,URHO) = &
                   interpolate(dist,npts_model_star, &
                               model_star_r,model_star_state(:,idens_model))

              state(i,j,k,UTEMP) = &
                   interpolate(dist,npts_model_star, &
                               model_star_r,model_star_state(:,itemp_model))

              do n = 1, nspec
                 state(i,j,k,UFS-1+n) = &
                      interpolate(dist,npts_model_star, &
                                  model_star_r,model_star_state(:,ispec_model-1+n))
              enddo

           enddo
        enddo
     enddo

     ! thermodynamics
     do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)

              eos_state%rho = state(i,j,k,URHO)
              eos_state%T = state(i,j,k,UTEMP)
              eos_state%xn(:) = state(i,j,k,UFS:UFS-1+nspec)

              call eos(eos_input_rt, eos_state)

              state(i,j,k,UEINT) = eos_state%e
              pres = eos_state%p

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
