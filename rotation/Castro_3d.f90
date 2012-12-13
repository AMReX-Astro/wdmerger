! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine ca_umdrv(is_finest_level,time,lo,hi,domlo,domhi, &
           uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
           uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, &
           ugdnvx_out,ugdnvx_l1,ugdnvx_l2,ugdnvx_l3,ugdnvx_h1,ugdnvx_h2,ugdnvx_h3, &
           ugdnvy_out,ugdnvy_l1,ugdnvy_l2,ugdnvy_l3,ugdnvy_h1,ugdnvy_h2,ugdnvy_h3, &
           ugdnvz_out,ugdnvz_l1,ugdnvz_l2,ugdnvz_l3,ugdnvz_h1,ugdnvz_h2,ugdnvz_h3, &
           src ,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3, &
           grav,gv_l1,gv_l2,gv_l3,gv_h1,gv_h2,gv_h3, &
           delta,dt, &
           flux1,flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3, &
           flux2,flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3, &
           flux3,flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3, &
           area1,area1_l1,area1_l2,area1_l3,area1_h1,area1_h2,area1_h3, &
           area2,area2_l1,area2_l2,area2_l3,area2_h1,area2_h2,area2_h3, &
           area3,area3_l1,area3_l2,area3_l3,area3_h1,area3_h2,area3_h3, &
           vol,vol_l1,vol_l2,vol_l3,vol_h1,vol_h2,vol_h3, &
           courno,verbose,mass_added,eint_added,eden_added,&
           E_added_flux,E_added_grav)

      use meth_params_module, only : QVAR, NVAR, NHYP, do_sponge, &
                                     normalize_species
      use advection_module, only : umeth3d, ctoprim, divu, consup, enforce_minimum_density, &
           normalize_new_species
      use sponge_module, only : sponge
      use grav_sources_module, only : add_grav_source

      ! This is used for IsoTurb only
      ! use probdata_module   , only : radiative_cooling_type

      implicit none

      integer is_finest_level
      integer lo(3),hi(3),verbose
      integer domlo(3),domhi(3)
      integer uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3
      integer uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3
      integer ugdnvx_l1,ugdnvx_l2,ugdnvx_l3,ugdnvx_h1,ugdnvx_h2,ugdnvx_h3
      integer ugdnvy_l1,ugdnvy_l2,ugdnvy_l3,ugdnvy_h1,ugdnvy_h2,ugdnvy_h3
      integer ugdnvz_l1,ugdnvz_l2,ugdnvz_l3,ugdnvz_h1,ugdnvz_h2,ugdnvz_h3
      integer flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3
      integer flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3
      integer flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3
      integer area1_l1,area1_l2,area1_l3,area1_h1,area1_h2,area1_h3
      integer area2_l1,area2_l2,area2_l3,area2_h1,area2_h2,area2_h3
      integer area3_l1,area3_l2,area3_l3,area3_h1,area3_h2,area3_h3
      integer vol_l1,vol_l2,vol_l3,vol_h1,vol_h2,vol_h3
      integer src_l1,src_l2,src_l3,src_h1,src_h2,src_h3
      integer gv_l1,gv_l2,gv_l3,gv_h1,gv_h2,gv_h3
      double precision   uin(  uin_l1:uin_h1,    uin_l2:uin_h2,     uin_l3:uin_h3,  NVAR)
      double precision  uout( uout_l1:uout_h1,  uout_l2:uout_h2,   uout_l3:uout_h3, NVAR)
      double precision ugdnvx_out(ugdnvx_l1:ugdnvx_h1,ugdnvx_l2:ugdnvx_h2,ugdnvx_l3:ugdnvx_h3)
      double precision ugdnvy_out(ugdnvy_l1:ugdnvy_h1,ugdnvy_l2:ugdnvy_h2,ugdnvy_l3:ugdnvy_h3)
      double precision ugdnvz_out(ugdnvz_l1:ugdnvz_h1,ugdnvz_l2:ugdnvz_h2,ugdnvz_l3:ugdnvz_h3)
      double precision   src(  src_l1:src_h1,    src_l2:src_h2,     src_l3:src_h3,  NVAR)
      double precision  grav( gv_l1:gv_h1,  gv_l2:gv_h2,   gv_l3:gv_h3,    3)
      double precision flux1(flux1_l1:flux1_h1,flux1_l2:flux1_h2, flux1_l3:flux1_h3,NVAR)
      double precision flux2(flux2_l1:flux2_h1,flux2_l2:flux2_h2, flux2_l3:flux2_h3,NVAR)
      double precision flux3(flux3_l1:flux3_h1,flux3_l2:flux3_h2, flux3_l3:flux3_h3,NVAR)
      double precision area1(area1_l1:area1_h1,area1_l2:area1_h2, area1_l3:area1_h3)
      double precision area2(area2_l1:area2_h1,area2_l2:area2_h2, area2_l3:area2_h3)
      double precision area3(area3_l1:area3_h1,area3_l2:area3_h2, area3_l3:area3_h3)
      double precision vol(vol_l1:vol_h1,vol_l2:vol_h2, vol_l3:vol_h3)
      double precision delta(3),dt,time,courno,E_added_flux,E_added_grav
      double precision mass_added,eint_added,eden_added

      ! Automatic arrays for workspace
      double precision, allocatable:: q(:,:,:,:)
      double precision, allocatable:: gamc(:,:,:)
      double precision, allocatable:: flatn(:,:,:)
      double precision, allocatable:: c(:,:,:)
      double precision, allocatable:: csml(:,:,:)
      double precision, allocatable:: div(:,:,:)
      double precision, allocatable:: pdivu(:,:,:)
      double precision, allocatable:: srcQ(:,:,:,:)

      double precision dx,dy,dz
      integer ngq,ngf,iflaten

      allocate(     q(uin_l1:uin_h1,uin_l2:uin_h2,uin_l3:uin_h3,QVAR))
      allocate(  gamc(uin_l1:uin_h1,uin_l2:uin_h2,uin_l3:uin_h3))
      allocate( flatn(uin_l1:uin_h1,uin_l2:uin_h2,uin_l3:uin_h3))
      allocate(     c(uin_l1:uin_h1,uin_l2:uin_h2,uin_l3:uin_h3))
      allocate(  csml(uin_l1:uin_h1,uin_l2:uin_h2,uin_l3:uin_h3))
      allocate(   div(lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1))

      allocate( pdivu(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))

      allocate(  srcQ(src_l1:src_h1,src_l2:src_h2,src_l3:src_h3,QVAR))

      dx = delta(1)
      dy = delta(2)
      dz = delta(3)

      ngq = NHYP
      ngf = 1
      iflaten = 1

      ! 1) Translate conserved variables (u) to primitive variables (q).
      ! 2) Compute sound speeds (c) and gamma (gamc).
      !    Note that (q,c,gamc,csml,flatn) are all dimensioned the same
      !    and set to correspond to coordinates of (lo:hi)
      ! 3) Translate source terms
      call ctoprim(lo,hi,uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
                   q,c,gamc,csml,flatn,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
                   src,srcQ,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3, &
                   courno,dx,dy,dz,dt,ngq,ngf,iflaten)


      ! Compute hyperbolic fluxes using unsplit Godunov
      call umeth3d(q,c,gamc,csml,flatn,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
                   srcQ,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3, &
                   grav,gv_l1,gv_l2,gv_l3,gv_h1,gv_h2,gv_h3, &
                   lo(1),lo(2),lo(3),hi(1),hi(2),hi(3),dx,dy,dz,dt, &
                   flux1,flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3, &
                   flux2,flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3, &
                   flux3,flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3, &
                   ugdnvx_out,ugdnvx_l1,ugdnvx_l2,ugdnvx_l3,ugdnvx_h1,ugdnvx_h2,ugdnvx_h3, &
                   ugdnvy_out,ugdnvy_l1,ugdnvy_l2,ugdnvy_l3,ugdnvy_h1,ugdnvy_h2,ugdnvy_h3, &
                   ugdnvz_out,ugdnvz_l1,ugdnvz_l2,ugdnvz_l3,ugdnvz_h1,ugdnvz_h2,ugdnvz_h3, &
                   pdivu)

      ! Compute divergence of velocity field (on surroundingNodes(lo,hi))
      call divu(lo,hi,q,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
                dx,dy,dz,div,lo(1),lo(2),lo(3),hi(1)+1,hi(2)+1,hi(3)+1)

      ! Conservative update
      call consup(uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
                  uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, &
                  src ,  src_l1,  src_l2,  src_l3,  src_h1,  src_h2,  src_h3, &
                  flux1,flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3, &
                  flux2,flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3, &
                  flux3,flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3, &
                  area1,area1_l1,area1_l2,area1_l3,area1_h1,area1_h2,area1_h3, &
                  area2,area2_l1,area2_l2,area2_l3,area2_h1,area2_h2,area2_h3, &
                  area3,area3_l1,area3_l2,area3_l3,area3_h1,area3_h2,area3_h3, &
                  vol,vol_l1,vol_l2,vol_l3,vol_h1,vol_h2,vol_h3, &
                  div,pdivu,lo,hi,dx,dy,dz,dt,E_added_flux)

      ! Add the radiative cooling -- for SGS only.
      ! if (radiative_cooling_type.eq.2) then
      !    call post_step_radiative_cooling(lo,hi,dt, &
      !         uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3)
      ! endif

      ! Enforce the density >= small_dens.
      call enforce_minimum_density(uin, uin_l1, uin_l2, uin_l3, uin_h1, uin_h2, uin_h3, &
                                   uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, &
                                   lo,hi,mass_added,eint_added,eden_added,verbose)

      ! Enforce species >= 0
      call ca_enforce_nonnegative_species(uout,uout_l1,uout_l2,uout_l3, &
                                          uout_h1,uout_h2,uout_h3,lo,hi)
 
      ! Re-normalize the species
      if (normalize_species .eq. 1) then
         call normalize_new_species(uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, &
                                    lo,hi)
      end if

      call add_grav_source(uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
                           uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, &
                           grav, gv_l1, gv_l2, gv_l3, gv_h1, gv_h2, gv_h3, &
                           lo,hi,dt,delta,E_added_grav)

      ! Impose sponge
      if (do_sponge .eq. 1) then
         call sponge(uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3,lo,hi, &
                     time,dt, &
                     dx,dy,dz,domlo,domhi)
      end if

      deallocate(q,gamc,flatn,c,csml,div,srcQ,pdivu)

      end subroutine ca_umdrv

! ::
! :: ----------------------------------------------------------
! ::

      subroutine ca_check_initial_species(lo,hi,&
                          state,state_l1,state_l2,state_l3,state_h1,state_h2,state_h3)

      use network           , only : nspec
      use meth_params_module, only : NVAR, URHO, UFS

      implicit none

      integer          :: lo(3), hi(3)
      integer          :: state_l1,state_l2,state_l3,state_h1,state_h2,state_h3
      double precision :: state(state_l1:state_h1,state_l2:state_h2,state_l3:state_h3,NVAR)

      ! Local variables
      integer          :: i,j,k,n
      double precision :: sum

      do k = lo(3), hi(3)
      do j = lo(2), hi(2)
      do i = lo(1), hi(1)

         sum = 0.d0
         do n = 1, nspec
             sum = sum + state(i,j,k,UFS+n-1)
         end do
         if (abs(state(i,j,k,URHO)-sum).gt. 1.d-8 * state(i,j,k,URHO)) then
            print *,'Sum of (rho X)_i vs rho at (i,j,k): ',i,j,k,sum,state(i,j,k,URHO)
            call bl_error("Error:: Failed check of initial species summing to 1")
         end if

      enddo
      enddo
      enddo

      end subroutine ca_check_initial_species

! :: ----------------------------------------------------------
! :: Volume-weight average the fine grid data onto the coarse
! :: grid.  Overlap is given in coarse grid coordinates.
! ::
! :: INPUTS / OUTPUTS:
! ::  crse      <=  coarse grid data
! ::  clo,chi    => index limits of crse array interior
! ::  nvar	 => number of components in arrays
! ::  fine       => fine grid data
! ::  flo,fhi    => index limits of fine array interior
! ::  lo,hi      => index limits of overlap (crse grid)
! ::  lrat       => refinement ratio
! ::
! :: NOTE:
! ::  Assumes all data cell centered
! :: ----------------------------------------------------------
! ::
      subroutine ca_avgdown(crse,c_l1,c_l2,c_l3,c_h1,c_h2,c_h3,nvar, &
                            cv,cv_l1,cv_l2,cv_l3,cv_h1,cv_h2,cv_h3, &
                            fine,f_l1,f_l2,f_l3,f_h1,f_h2,f_h3, &
                            fv,fv_l1,fv_l2,fv_l3,fv_h1,fv_h2,fv_h3,lo,hi,lrat)

      implicit none
      integer c_l1,c_l2,c_l3,c_h1,c_h2,c_h3
      integer cv_l1,cv_l2,cv_l3,cv_h1,cv_h2,cv_h3
      integer f_l1,f_l2,f_l3,f_h1,f_h2,f_h3
      integer fv_l1,fv_l2,fv_l3,fv_h1,fv_h2,fv_h3
      integer lo(3), hi(3)
      integer nvar, lrat(3)
      double precision crse(c_l1:c_h1,c_l2:c_h2,c_l3:c_h3,nvar)
      double precision cv(cv_l1:cv_h1,cv_l2:cv_h2,cv_l3:cv_h3)
      double precision fine(f_l1:f_h1,f_l2:f_h2,f_l3:f_h3,nvar)
      double precision fv(fv_l1:fv_h1,fv_l2:fv_h2,fv_l3:fv_h3)

      integer i, j, k, n, ic, jc, kc, ioff, joff, koff
      integer lratx, lraty, lratz
      double precision   volfrac

      lratx   = lrat(1)
      lraty   = lrat(2)
      lratz   = lrat(3)
      volfrac = 1.d0/float(lrat(1)*lrat(2)*lrat(3))

      do n = 1, nvar
         !
         ! Set coarse grid to zero on overlap.
         !
         do kc = lo(3), hi(3)
            do jc = lo(2), hi(2)
               do ic = lo(1), hi(1)
                  crse(ic,jc,kc,n) = 0.d0
               enddo
            enddo
         enddo
         !
         ! Sum fine data.
         !
         do koff = 0, lratz-1
            !$OMP PARALLEL DO PRIVATE(i,j,k,ic,jc,kc,ioff,joff)
            do kc = lo(3),hi(3)
               k = kc*lratz + koff
               do joff = 0, lraty-1
                  do jc = lo(2), hi(2)
                     j = jc*lraty + joff
                     do ioff = 0, lratx-1
                        do ic = lo(1), hi(1)
                           i = ic*lratx + ioff
                           crse(ic,jc,kc,n) = crse(ic,jc,kc,n) + fine(i,j,k,n)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
            !$OMP END PARALLEL DO
         enddo
         !
         ! Divide out by volume weight.
         !
         !$OMP PARALLEL DO PRIVATE(ic,jc,kc)
         do kc = lo(3), hi(3)
            do jc = lo(2), hi(2)
               do ic = lo(1), hi(1)
                  crse(ic,jc,kc,n) = volfrac*crse(ic,jc,kc,n)
               enddo
            enddo
         enddo
         !$OMP END PARALLEL DO

      enddo

      end subroutine ca_avgdown

! ::
! :: ----------------------------------------------------------
! ::

      subroutine ca_compute_avgstate(lo,hi,dx,dr,nc,&
                                     state,s_l1,s_l2,s_l3,s_h1,s_h2,s_h3,radial_state, &
                                     vol,v_l1,v_l2,v_l3,v_h1,v_h2,v_h3,radial_vol, &
                                     problo,numpts_1d)

      use meth_params_module, only : URHO, UMX, UMY, UMZ
      use probdata_module
      implicit none

      integer          :: lo(3),hi(3),nc
      double precision :: dx(3),dr,problo(3)

      integer          :: numpts_1d
      double precision :: radial_state(nc,0:numpts_1d-1)
      double precision :: radial_vol(0:numpts_1d-1)

      integer          :: s_l1,s_l2,s_l3,s_h1,s_h2,s_h3
      double precision :: state(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3,nc)

      integer          :: v_l1,v_l2,v_l3,v_h1,v_h2,v_h3
      double precision :: vol(v_l1:v_h1,v_l2:v_h2,v_l3:v_h3)

      integer          :: i,j,k,n,index
      double precision :: x,y,z,r
      double precision :: x_mom,y_mom,z_mom,radial_mom
      !
      ! Do not OMP this.
      !
      do k = lo(3), hi(3)
         z = problo(3) + (dble(k)+0.50d0) * dx(3) - center(3)
         do j = lo(2), hi(2)
            y = problo(2) + (dble(j)+0.50d0) * dx(2) - center(2)
            do i = lo(1), hi(1)
               x = problo(1) + (dble(i)+0.50d0) * dx(1) - center(1)
               r = sqrt(x**2 + y**2 + z**2)
               index = int(r/dr)
               if (index .gt. numpts_1d-1) then
                   print *,'COMPUTE_AVGSTATE: INDEX TOO BIG ',index,' > ',numpts_1d-1
                   print *,'AT (i,j,k) ',i,j,k
                   print *,'R / DR ',r,dr
                   call bl_error("Error:: Castro_3d.f90 :: ca_compute_avgstate")
               end if
               radial_state(URHO,index) = radial_state(URHO,index) &
                    + vol(i,j,k)*state(i,j,k,URHO)
               !
               ! Store the radial component of the momentum in the 
               ! UMX, UMY and UMZ components for now.
               !
               x_mom = state(i,j,k,UMX)
               y_mom = state(i,j,k,UMY)
               z_mom = state(i,j,k,UMZ)
               radial_mom = x_mom * (x/r) + y_mom * (y/r) + z_mom * (z/r)
               radial_state(UMX,index) = radial_state(UMX,index) + vol(i,j,k)*radial_mom
               radial_state(UMY,index) = radial_state(UMY,index) + vol(i,j,k)*radial_mom
               radial_state(UMZ,index) = radial_state(UMZ,index) + vol(i,j,k)*radial_mom

               do n = UMZ+1,nc
                  radial_state(n,index) = radial_state(n,index) + vol(i,j,k)*state(i,j,k,n)
               end do
               radial_vol(index) = radial_vol(index) + vol(i,j,k)
            enddo
         enddo
      enddo

      end subroutine ca_compute_avgstate

! ::
! :: ----------------------------------------------------------
! ::

      subroutine ca_enforce_nonnegative_species(uout,uout_l1,uout_l2,uout_l3, &
                                                uout_h1,uout_h2,uout_h3,lo,hi)

      use network, only : nspec
      use meth_params_module, only : NVAR, URHO, UFS

      implicit none

      integer          :: lo(3), hi(3)
      integer          :: uout_l1, uout_l2, uout_l3, uout_h1, uout_h2, uout_h3
      double precision :: uout(uout_l1:uout_h1,uout_l2:uout_h2,uout_l3:uout_h3,NVAR)

      ! Local variables
      integer          :: i,j,k,n
      integer          :: int_dom_spec
      logical          :: any_negative
      double precision :: dom_spec,x

      double precision, parameter :: eps = -1.0d-16

      !$OMP PARALLEL DO PRIVATE(i,j,k,n,dom_spec,int_dom_spec,any_negative,x)
      do k = lo(3),hi(3)
      do j = lo(2),hi(2)
      do i = lo(1),hi(1)

         any_negative = .false.
         !
         ! First deal with tiny undershoots by just setting them to zero.
         !
         do n = UFS, UFS+nspec-1
           if (uout(i,j,k,n) .lt. 0.d0) then
              x = uout(i,j,k,n)/uout(i,j,k,URHO)
              if (x .gt. eps) then
                 uout(i,j,k,n) = 0.d0
              else
                 any_negative = .true.
              end if
           end if
         end do
         !
         ! We know there are one or more undershoots needing correction.
         !
         if (any_negative) then
            !
            ! Find the dominant species.
            !
            int_dom_spec = UFS
            dom_spec     = uout(i,j,k,int_dom_spec)

            do n = UFS,UFS+nspec-1
              if (uout(i,j,k,n) .gt. dom_spec) then
                dom_spec     = uout(i,j,k,n)
                int_dom_spec = n
              end if
           end do
           !
           ! Now take care of undershoots greater in magnitude than 1e-16.
           !
           do n = UFS, UFS+nspec-1
 
              if (uout(i,j,k,n) .lt. 0.d0) then

                 x = uout(i,j,k,n)/uout(i,j,k,URHO)
                 !
                 ! Here we only print the bigger negative values.
                 !
                 if (x .lt. -1.d-2) then
                    print *,'Correcting nth negative species ',n-UFS+1
                    print *,'   at cell (i,j,k)              ',i,j,k
                    print *,'Negative (rho*X) is             ',uout(i,j,k,n)
                    print *,'Negative      X  is             ',x
                    print *,'Filling from dominant species   ',int_dom_spec-UFS+1
                    print *,'  which had X =                 ',&
                             uout(i,j,k,int_dom_spec) / uout(i,j,k,URHO)
                 end if
                 !
                 ! Take enough from the dominant species to fill the negative one.
                 !
                 uout(i,j,k,int_dom_spec) = uout(i,j,k,int_dom_spec) + uout(i,j,k,n)
                 !
                 ! Test that we didn't make the dominant species negative.
                 !
                 if (uout(i,j,k,int_dom_spec) .lt. 0.d0) then 
                    print *,' Just made nth dominant species negative ',int_dom_spec-UFS+1,' at ',i,j,k 
                    print *,'We were fixing species ',n-UFS+1,' which had value ',x
                    print *,'Dominant species became ',uout(i,j,k,int_dom_spec) / uout(i,j,k,URHO)
                    call bl_error("Error:: Castro_3d.f90 :: ca_enforce_nonnegative_species")
                 end if
                 !
                 ! Now set the negative species to zero.
                 !
                 uout(i,j,k,n) = 0.d0

              end if

           enddo
         end if
      enddo
      enddo
      enddo
      !$OMP END PARALLEL DO

      end subroutine ca_enforce_nonnegative_species

! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine get_center(center_out)

        use probdata_module, only : center

        implicit none

        double precision, intent(inout) :: center_out(3)

        center_out(1:3) = center(1:3)

      end subroutine get_center

! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine set_center(center_in)

        use probdata_module, only : center

        implicit none

        double precision :: center_in(3)

        center(1:3) = center_in(1:3)

      end subroutine set_center

! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine find_center(data,new_center,icen,dx,problo)

        implicit none

        double precision :: data(-1:1,-1:1,-1:1)
        double precision :: new_center(3)
        double precision :: dx(3),problo(3)
        double precision :: a,b,x,y,z,cen
        integer          :: icen(3)
        integer          :: i,j,k

        ! We do this to take care of precision issues
        cen = data(0,0,0)
        do k = -1,1
        do j = -1,1
        do i = -1,1
           data(i,j,k) = data(i,j,k) - cen 
        end do
        end do
        end do

!       This puts the "center" at the cell center
        new_center(1) = problo(1) +  (icen(1)+0.5d0) * dx(1)
        new_center(2) = problo(2) +  (icen(2)+0.5d0) * dx(2)
        new_center(3) = problo(3) +  (icen(3)+0.5d0) * dx(3)
   
        ! Fit parabola y = a x^2  + b x + c through three points
        ! a = 1/2 ( y_1 + y_-1)
        ! b = 1/2 ( y_1 - y_-1)
        ! x_vertex = -b / 2a

        ! ... in x-direction
        a = 0.5d0 * (data(1,0,0) + data(-1,0,0)) - data(0,0,0)
        b = 0.5d0 * (data(1,0,0) - data(-1,0,0)) - data(0,0,0)
        x = -b / (2.d0*a)
        new_center(1) = new_center(1) +  x*dx(1)

        ! ... in y-direction
        a = 0.5d0 * (data(0,1,0) + data(0,-1,0)) - data(0,0,0)
        b = 0.5d0 * (data(0,1,0) - data(0,-1,0)) - data(0,0,0)
        y = -b / (2.d0*a)
        new_center(2) = new_center(2) +  y*dx(2)

        ! ... in z-direction
        a = 0.5d0 * (data(0,0,1) + data(0,0,-1)) - data(0,0,0)
        b = 0.5d0 * (data(0,0,1) - data(0,0,-1)) - data(0,0,0)
        z = -b / (2.d0*a)
        new_center(3) = new_center(3) +  z*dx(3)

      end subroutine find_center
