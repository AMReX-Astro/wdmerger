
     subroutine ca_ext_src(lo,hi,&
                           old_state,os_l1,os_l2,os_l3,os_h1,os_h2,os_h3,&
                           new_state,ns_l1,ns_l2,ns_l3,ns_h1,ns_h2,ns_h3,&
                           src,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3,problo,dx,time,dt)

       use meth_params_module,  only: NVAR, URHO, UMX, UMZ, UEDEN
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

       double precision :: rho, rhoInv
       double precision :: r(3)
       integer          :: i, j, k

       ! Note that this function exists in a tiling region so we should only 
       ! modify the zones between lo and hi. 

       src(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) = ZERO

       do k = lo(3), hi(3)
          r(3) = problo(3) + (dble(k) + HALF) * dx(3) - center(3)

          do j = lo(2), hi(2)
             r(2) = problo(2) + (dble(j) + HALF) * dx(2) - center(2)

             do i = lo(1), hi(1)
                r(1) = problo(1) + (dble(i) + HALF) * dx(1) - center(1)

             enddo
          enddo
       enddo
 
     end subroutine ca_ext_src

