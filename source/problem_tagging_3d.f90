subroutine set_problem_tags(tag,tagl1,tagl2,tagl3,tagh1,tagh2,tagh3, &
                            state,state_l1,state_l2,state_l3,state_h1,state_h2,state_h3,&
                            set,clear,&
                            lo,hi,&
                            dx,problo,time,level)

  use bl_constants_module, only: ZERO, HALF
  use meth_params_module, only: URHO, UMX, UMY, UMZ, UEDEN, NVAR
  use prob_params_module, only: center, probhi
  use probdata_module, only: maxTaggingRadius
 
  implicit none
  
  integer         ,intent(in   ) :: lo(3),hi(3)
  integer         ,intent(in   ) :: state_l1,state_l2,state_l3, &
                                    state_h1,state_h2,state_h3
  integer         ,intent(in   ) :: tagl1,tagl2,tagl3,tagh1,tagh2,tagh3
  double precision,intent(in   ) :: state(state_l1:state_h1, &
                                          state_l2:state_h2, &
                                          state_l3:state_h3,NVAR)
  integer         ,intent(inout) :: tag(tagl1:tagh1,tagl2:tagh2,tagl3:tagh3)
  double precision,intent(in   ) :: problo(3),dx(3),time
  integer         ,intent(in   ) :: level,set,clear

  integer          :: i, j, k
  double precision :: x,y,z,r

  ! Clear all tagging that occurs outside the radius set by maxTaggingRadius.

  do k = lo(3), hi(3)
     z = problo(3) + dble(k + HALF)*dx(3) - center(3)

     do j = lo(2), hi(2)
        y = problo(2) + dble(j + HALF)*dx(2) - center(2)

        do i = lo(1), hi(1)
           x = problo(1) + dble(i + HALF)*dx(1) - center(1)

           r = (x**2 + y**2 + z**2)**HALF

           if (r .gt. maxTaggingRadius * maxval(abs(problo-center)) .or. &
               r .gt. maxTaggingRadius * maxval(abs(probhi-center)) ) then

             tag(i,j,k) = clear

           endif

        enddo
     enddo
  enddo

end subroutine set_problem_tags
