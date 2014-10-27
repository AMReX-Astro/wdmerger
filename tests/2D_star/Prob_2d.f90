subroutine PROBINIT (init,name,namlen,problo,probhi)

  use eos_module
  use eos_type_module
  use bl_error_module
  use network
  use probdata_module

  implicit none

  integer init, namlen
  integer name(namlen)
  double precision problo(2), probhi(2)
  double precision xn(nspec)
  
  integer untin,i

  type (eos_t) :: eos_state

  namelist /fortin/ p_s, u_s, rho_s, p_a, u_a, rho_a, T_s, T_a, width, idir, &
       use_Tinit, &
       denerr,  dengrad,  max_denerr_lev,  max_dengrad_lev, &
       velgrad,  max_velgrad_lev, &
       presserr,pressgrad,max_presserr_lev,max_pressgrad_lev

  !
  !     Build "probin" filename -- the name of file containing fortin namelist.
  !     
  integer maxlen
  parameter (maxlen=256)
  character probin*(maxlen)

  if (namlen .gt. maxlen) then
     call bl_error("probin file name too long")
  end if

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do
         
  ! set namelist defaults

  p_s = 1.0               ! left pressure (erg/cc)
  u_s = 0.0               ! left velocity (cm/s)
  rho_s = 1.0             ! left density (g/cc)
  T_s = 1.0
  
  p_a = 0.1               ! right pressure (erg/cc)
  u_a = 0.0               ! right velocity (cm/s)
  rho_a = 0.125           ! right density (g/cc)
  T_a = 1.0

  width = 0.25            ! fraction of the domain for the interface

  use_Tinit = .false.     ! optionally use T_l/r instead of p_l/r for initialization

  denerr = 1.d20
  dengrad = 1.d20
  max_denerr_lev = -1
  max_dengrad_lev = -1

  presserr = 1.d20
  pressgrad = 1.d20
  max_presserr_lev = -1
  max_pressgrad_lev = -1

  velgrad = 1.d20
  max_velgrad_lev = -1

  !     Read namelists
  untin = 9
  open(untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin,fortin)
  close(unit=untin)

  print *, u_s, u_a
  
  center(1) = 0.5*(problo(1)+probhi(1))
  center(2) = 0.5*(problo(2)+probhi(2))

  !     compute the internal energy (erg/cc) for the left and right state
  xn(:) = 0.0d0
  xn(1) = 1.0d0

  if (use_Tinit) then

     eos_state%rho = rho_s
     eos_state%T = T_s
     eos_state%xn(:) = xn(:)

     call eos(eos_input_rt, eos_state, .false.)

     rhoe_s = rho_s*eos_state%e
     p_s = eos_state%p

     eos_state%rho = rho_a
     eos_state%T = T_a
     eos_state%xn(:) = xn(:)

     call eos(eos_input_rt, eos_state, .false.)

     rhoe_a = rho_a*eos_state%e
     p_a = eos_state%p

  else

     eos_state%rho = rho_s
     eos_state%p = p_s
     eos_state%T = 100000.d0   ! initial guess
     eos_state%xn(:) = xn(:)

     call eos(eos_input_rp, eos_state, .false.)

     rhoe_s = rho_s*eos_state%e
     T_s = eos_state%T

     eos_state%rho = rho_a
     eos_state%p = p_a
     eos_state%T = 100000.d0   ! initial guess
     eos_state%xn(:) = xn(:)

     call eos(eos_input_rp, eos_state, .false.)

     rhoe_a = rho_a*eos_state%e
     T_a = eos_state%T

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
                      state,state_l1,state_l2,state_h1,state_h2,delta,xlo,xhi)

  use network, only: nspec
  use probdata_module
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UEDEN, UEINT, UTEMP, UFS

  implicit none
  integer level, nscal
  integer lo(2), hi(2)
  integer state_l1,state_h1,state_l2,state_h2
  double precision state(state_l1:state_h1,state_l2:state_h2,NVAR)
  double precision time, delta(2)
  double precision xlo(2), xhi(2)
  
  double precision xcen, ycen, r
  integer i, j

  do j = lo(2), hi(2)
     do i = lo(1), hi(1)
        xcen = xlo(1) + delta(1)*(float(i-lo(1)) + 0.5d0) - center(1)
        ycen = xlo(2) + delta(2)*(float(j-lo(2)) + 0.5d0) - center(2)

        r = (xcen**2 + ycen**2)**0.5d0

        if (r .le. (width * center(1))) then
           state(i,j,URHO ) = rho_s
           state(i,j,UMX  ) = rho_s*u_s
           state(i,j,UMY  ) = 0.0d0
           state(i,j,UEDEN) = rhoe_s + 0.5*rho_s*u_s*u_s
           state(i,j,UEINT) = rhoe_s
           state(i,j,UTEMP) = T_s
        else
           state(i,j,URHO ) = rho_a
           state(i,j,UMX  ) = rho_a*u_a
           state(i,j,UMY  ) = 0.0d0
           state(i,j,UEDEN) = rhoe_a + 0.5*rho_a*u_a*u_a
           state(i,j,UEINT) = rhoe_a
           state(i,j,UTEMP) = T_a
        endif

        state(i,j,UFS:UFS-1+nspec) = 0.0d0
        state(i,j,UFS  ) = state(i,j,URHO)

     enddo
  enddo

end subroutine ca_initdata


! ::: -----------------------------------------------------------

subroutine ca_hypfill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
                      domlo,domhi,delta,xlo,time,bc)

  use bl_error_module
  use meth_params_module, only : NVAR

  implicit none

  include 'bc_types.fi'
  integer adv_l1,adv_l2,adv_h1,adv_h2
  integer bc(2,2,*)
  integer domlo(2), domhi(2)
  double precision delta(2), xlo(2), time
  double precision adv(adv_l1:adv_h1,adv_l2:adv_h2,NVAR)
  
  integer n

  do n = 1,NVAR
     call filcc(adv(adv_l1,adv_l2,n), &
                adv_l1,adv_l2,adv_h1,adv_h2, &
                domlo,domhi,delta,xlo,bc(1,1,n))
  enddo

  do n = 1,NVAR

     !        XLO
     if ( bc(1,1,n).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then
        call bl_error('SHOULD NEVER GET HERE bc(1,1,n) .eq. EXT_DIR) ')
     end if

     !        XHI
     if ( bc(1,2,n).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then
        call bl_error('SHOULD NEVER GET HERE bc(1,2,n) .eq. EXT_DIR) ')
     end if

     !        YLO
     if ( bc(2,1,n).eq.EXT_DIR .and. adv_l2.lt.domlo(2)) then
        call bl_error('SHOULD NEVER GET HERE bc(2,1,n) .eq. EXT_DIR) ')
     end if

     !        YHI
     if ( bc(2,2,n).eq.EXT_DIR .and. adv_l2.lt.domlo(2)) then
        call bl_error('SHOULD NEVER GET HERE bc(2,2,n) .eq. EXT_DIR) ')
     end if

  end do

end subroutine ca_hypfill


! ::: -----------------------------------------------------------

subroutine ca_denfill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
                      domlo,domhi,delta,xlo,time,bc)

  use bl_error_module

  implicit none

  include 'bc_types.fi'
  integer adv_l1,adv_l2,adv_h1,adv_h2
  integer bc(2,2,*)
  integer domlo(2), domhi(2)
  double precision delta(2), xlo(2), time
  double precision adv(adv_l1:adv_h1,adv_l2:adv_h2)

  !     Note: this function should not be needed, technically, but is provided
  !     to filpatch because there are many times in the algorithm when just
  !     the density is needed.  We try to rig up the filling so that the same
  !     function is called here and in hypfill where all the states are filled.
  
  call filcc(adv,adv_l1,adv_l2,adv_h1,adv_h2,domlo,domhi,delta,xlo,bc)

  !     XLO
  if ( bc(1,1,1).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then
     call bl_error('SHOULD NEVER GET HERE bc(1,1,1) .eq. EXT_DIR) ')
  end if

  !     XHI
  if ( bc(1,2,1).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then
     call bl_error('SHOULD NEVER GET HERE bc(1,2,1) .eq. EXT_DIR) ')
  end if

  !     YLO
  if ( bc(2,1,1).eq.EXT_DIR .and. adv_l2.lt.domlo(2)) then
     call bl_error('SHOULD NEVER GET HERE bc(2,1,1) .eq. EXT_DIR) ')
  end if

  !     YHI
  if ( bc(2,2,1).eq.EXT_DIR .and. adv_l2.lt.domlo(2)) then
     call bl_error('SHOULD NEVER GET HERE bc(2,2,1) .eq. EXT_DIR) ')
  end if

end subroutine ca_denfill
