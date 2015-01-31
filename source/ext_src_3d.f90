
     subroutine ca_ext_src(lo,hi,&
                           old_state,old_state_l1,old_state_l2,old_state_l3,old_state_h1,old_state_h2,old_state_h3,&
                           new_state,new_state_l1,new_state_l2,new_state_l3,new_state_h1,new_state_h2,new_state_h3,&
                           src,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3,problo,dx,time,dt)

     use meth_params_module, only : NVAR, URHO, UMX, UMY, UEDEN, rot_period
     use probdata_module   , only : center, damping, damping_alpha
     use bl_constants_module

     implicit none

     integer         ,intent(in   ) :: lo(3),hi(3)
     integer         ,intent(in   ) :: old_state_l1,old_state_l2,old_state_l3,old_state_h1,old_state_h2,old_state_h3
     integer         ,intent(in   ) :: new_state_l1,new_state_l2,new_state_l3,new_state_h1,new_state_h2,new_state_h3
     integer         ,intent(in   ) :: src_l1,src_l2,src_l3,src_h1,src_h2,src_h3
     double precision,intent(in   ) :: old_state(old_state_l1:old_state_h1,old_state_l2:old_state_h2, &
                                                 old_state_l3:old_state_h3,NVAR)
     double precision,intent(in   ) :: new_state(new_state_l1:new_state_h1,new_state_l2:new_state_h2, &
                                                 new_state_l3:new_state_h3,NVAR)
     double precision,intent(  out) :: src(src_l1:src_h1,src_l2:src_h2,src_l3:src_h3,NVAR)
     double precision,intent(in   ) :: problo(3),dx(3),time,dt

     double precision :: x, y
     double precision :: omega
     double precision :: rho, px, py
     integer          :: i, j, k

     ! lo and hi specify work region     
     src(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) = ZERO ! Fill work region only
     
     if ( damping ) then

       omega = 2.d0 * M_PI / rot_period

       do i = lo(1), hi(1)
         x = problo(1) + (dble(i) + HALF) * dx(1) - center(1)

         do j = lo(2), hi(2)
           y = problo(2) + (dble(j) + HALF) * dx(2) - center(2)

           do k = lo(3), hi(3)
             rho = new_state(i,j,k,URHO)
             px  = new_state(i,j,k,UMX )
             py  = new_state(i,j,k,UMY )
             src(i,j,k,UMX)   =         damping_alpha / rot_period * rho * omega * y
             src(i,j,k,UMY)   = -1.d0 * damping_alpha / rot_period * rho * omega * x
             src(i,j,k,UEDEN) =         damping_alpha / rot_period * omega * ( px * y - py * x ) / (2.d0 * rho)
           enddo

         enddo

       enddo 
 
     endif

     end subroutine ca_ext_src

