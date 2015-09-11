
     subroutine ca_ext_src(lo,hi,&
                           old_state,os_lo,os_hi,&
                           new_state,ns_lo,ns_hi,&
                           src,src_lo,src_hi,problo,dx,time,dt)

       use meth_params_module,  only: NVAR, URHO, UMX, UMZ, UEDEN
       use prob_params_module,  only: center
       use bl_constants_module, only: ZERO, HALF, ONE, TWO
       use probdata_module,     only: do_initial_relaxation, relaxation_timescale

       implicit none

       integer          :: lo(3),hi(3)
       integer          :: os_lo(3),os_hi(3)
       integer          :: ns_lo(3),ns_hi(3)
       integer          :: src_lo(3),src_hi(3)
       double precision :: old_state(os_lo(1):os_hi(1),os_lo(2):os_hi(2),os_lo(3):os_hi(3),NVAR)
       double precision :: new_state(ns_lo(1):ns_hi(1),ns_lo(2):ns_hi(2),ns_lo(3):ns_hi(3),NVAR)
       double precision :: src(src_lo(1):src_hi(1),src_lo(2):src_hi(2),src_lo(3):src_hi(3),NVAR)
       double precision :: problo(3),dx(3),time,dt

       ! Local variables

       double precision :: damping_factor
       double precision :: r(3)
       integer          :: i, j, k

       ! Note that this function exists in a tiling region so we should only 
       ! modify the zones between lo and hi. 

       src(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) = ZERO

       if (do_initial_relaxation) then

          damping_factor = ONE / relaxation_timescale

       else

          damping_factor = ZERO
          
       endif

       do k = lo(3), hi(3)
          r(3) = problo(3) + (dble(k) + HALF) * dx(3) - center(3)

          do j = lo(2), hi(2)
             r(2) = problo(2) + (dble(j) + HALF) * dx(2) - center(2)

             do i = lo(1), hi(1)
                r(1) = problo(1) + (dble(i) + HALF) * dx(1) - center(1)

                   src(i,j,k,UMX:UMZ) = -new_state(i,j,k,UMX:UMZ) * damping_factor

                   ! Kinetic energy source term is v . momentum source term

                   src(i,j,k,UEDEN) = dot_product(new_state(i,j,k,UMX:UMZ) / new_state(i,j,k,URHO), src(i,j,k,UMX:UMZ))

             enddo
          enddo
       enddo
 
     end subroutine ca_ext_src

