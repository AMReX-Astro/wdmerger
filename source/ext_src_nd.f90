
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
       integer          :: i, j, k
       double precision :: new_mom(3), old_mom(3), new_ke, old_ke, rhoInv, dtInv

       ! Note that this function exists in a tiling region so we should only 
       ! modify the zones between lo and hi. 

       src(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) = ZERO

       if (do_initial_relaxation) then

          damping_factor = ONE / (ONE + HALF * dt / relaxation_timescale)

       else

          damping_factor = ONE
          
       endif

       dtInv = ONE / dt

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                
                rhoInv = ONE / new_state(i,j,k,URHO)

                ! The form of the source term is d(rho * v) / dt = - (rho * v) / tau.
                ! The solution to this is an exponential decline, (rho * v) = (rho * v)_0 * exp(-t / tau).
                ! However, an explicit discretized update of mom -> -mom * (dt / 2) / tau does not provide
                ! a good approximation to this if tau << t. So we want to implicitly solve
                ! that update: (rho * v) -> (rho * v) / (1 + (dt / 2) / tau). Then, since we will be
                ! multiplying this source term by dt / 2, divide by that here to get the scaling right.

                old_mom = new_state(i,j,k,UMX:UMZ)
                new_mom = new_state(i,j,k,UMX:UMZ) * damping_factor

                src(i,j,k,UMX:UMZ) = (TWO * dtInv) * (new_mom - old_mom)

                ! Do the same thing for the kinetic energy update.

                old_ke = HALF * sum(old_mom**2) * rhoInv
                new_ke = HALF * sum(new_mom**2) * rhoInv

                src(i,j,k,UEDEN) = (TWO * dtInv) * (new_ke - old_ke)
                   
             enddo
          enddo
       enddo
 
     end subroutine ca_ext_src

