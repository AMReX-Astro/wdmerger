
     subroutine ca_ext_src(lo,hi,&
                           old_state,os_l1,os_l2,os_l3,os_h1,os_h2,os_h3,&
                           new_state,ns_l1,ns_l2,ns_l3,ns_h1,ns_h2,ns_h3,&
                           src,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3,problo,dx,time,dt)

       use meth_params_module,  only: NVAR, URHO, UMX, UMZ, UEDEN
       use probdata_module   ,  only: orbital_damping, orbital_damping_alpha
       use prob_params_module,  only: center
       use bl_constants_module, only: ZERO, HALF, ONE, TWO
       use rot_sources_module,  only: get_omega, cross_product

       implicit none

       integer         , intent(in   ) :: lo(3),hi(3)
       integer         , intent(in   ) :: os_l1,os_l2,os_l3,os_h1,os_h2,os_h3
       integer         , intent(in   ) :: ns_l1,ns_l2,ns_l3,ns_h1,ns_h2,ns_h3
       integer         , intent(in   ) :: src_l1,src_l2,src_l3,src_h1,src_h2,src_h3
       double precision, intent(in   ) :: old_state(os_l1:os_h1,os_l2:os_h2,os_l3:os_h3,NVAR)
       double precision, intent(in   ) :: new_state(ns_l1:ns_h1,ns_l2:ns_h2,ns_l3:ns_h3,NVAR)
       double precision, intent(inout) :: src(src_l1:src_h1,src_l2:src_h2,src_l3:src_h3,NVAR)
       double precision, intent(in   ) :: problo(3),dx(3),time,dt

       ! Local variables

       double precision :: rho, rhoInv, omegacrossv(3), omegacrossr(3)
       double precision :: r(3), v(3), Sr(3), omega(3)
       integer          :: i, j, k

       omega = get_omega()

       ! Note that this function exists in a tiling region so we should only 
       ! modify the zones between lo and hi. 

       src(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) = ZERO

       do k = lo(3), hi(3)
          r(3) = problo(3) + (dble(k) + HALF) * dx(3) - center(3)

          do j = lo(2), hi(2)
             r(2) = problo(2) + (dble(j) + HALF) * dx(2) - center(2)

             do i = lo(1), hi(1)
                r(1) = problo(1) + (dble(i) + HALF) * dx(1) - center(1)

                if ( orbital_damping ) then

                   rho = new_state(i,j,k,URHO)
                   rhoInv = ONE / rho

                   v = new_state(i,j,k,UMX:UMZ) * rhoInv
                   omegacrossv = cross_product(omega, v)
                   omegacrossr = cross_product(omega, r)

                   ! For the orbital damping term, our strategy is simply to subtract some fraction
                   ! of the rotational source term from the momentum update so that the system isn't
                   ! fully rotationally supported.

                   Sr = -TWO * rho * omegacrossv - rho * cross_product(omega, omegacrossr)
                   Sr = -orbital_damping_alpha * Sr

                   src(i,j,k,UMX:UMZ) = Sr

                   ! Kinetic energy source term is v . Sr

                   src(i,j,k,UEDEN) = dot_product(v, sr)

                endif

             enddo
          enddo
       enddo
 
     end subroutine ca_ext_src

