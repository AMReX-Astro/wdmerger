subroutine set_problem_tags(tag,tagl1,tagl2,tagl3,tagh1,tagh2,tagh3, &
                            state,state_l1,state_l2,state_l3,state_h1,state_h2,state_h3,&
                            set,clear,&
                            lo,hi,&
                            dx,problo,time,level)

  use bl_constants_module, only: ZERO, HALF
  use meth_params_module, only: URHO, UMX, UMY, UMZ, UEDEN, NVAR
  use prob_params_module, only: center
  use probdata_module, only: radius_P_initial, radius_S_initial, &
                             a_P_initial, a_S_initial, &
                             starBuffer, boundaryBuffer
 
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
  double precision :: x,y,z,rSq
  double precision :: xMax, yMax, zMax, rP, rS, rMax, rMaxSq, rMin, rMinSq

  ! We define two spherical volumes that define boundaries for refinement.
  ! The first is the volume of the domain outside of which no refinement may happen.
  ! boundaryBuffer = fraction of the distance from the center to the
  ! edges where refinement is permitted. Set this slightly lower than where you
  ! want the refinement to stop because the error buffers will tag zones outside this
  ! boundary for refinement.

  xMax = boundaryBuffer * (center(1) - problo(1))
  yMax = boundaryBuffer * (center(2) - problo(2))
  zMax = boundaryBuffer * (center(3) - problo(3))

  rMax = xMax

  if (yMax < rMax) rMax = yMax
  if (zMax < rMax) rMax = zMax

  ! We also set some inner portion of the domain that should always be refined.
  ! This is a spherical volume that extends past the stars by a fraction starBuffer
  ! of the radius of the stars. For a WD the secondary is larger so we'll use that radius
  ! for both sides of the domain to retain parity along the axis joining the stars.

  rP = starBuffer * radius_P_initial + a_P_initial
  rS = starBuffer * radius_S_initial + a_S_initial

  rMin = rS

  if (rP > rMin) rMin = rP

  rMaxSq = rMax * rMax
  rMinSq = rMin * rMin

  if (level .eq. 0) then
     do k = lo(3), hi(3)
        z = problo(3) + dble(k + HALF)*dx(3) - center(3)

        do j = lo(2), hi(2)
           y = problo(2) + dble(j + HALF)*dx(2) - center(2)

           do i = lo(1), hi(1)
              x = problo(1) + dble(i + HALF)*dx(1) - center(1)

              rSq = x**2 + y**2 + z**2

              if ( rSq .lt. rMinSq ) then

                tag(i,j,k) = set

              elseif (rSq .gt. rMaxSq ) then

                tag(i,j,k) = clear

              endif

           enddo
        enddo
     enddo
  endif

end subroutine set_problem_tags
