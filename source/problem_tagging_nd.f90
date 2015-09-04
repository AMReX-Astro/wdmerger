subroutine set_problem_tags(tag,tag_lo,tag_hi, &
                            state,state_lo,state_hi, &
                            set,clear,&
                            lo,hi,&
                            dx,problo,time,level)

  use bl_constants_module, only: ZERO, HALF, TWO
  use meth_params_module, only: URHO, UMX, UMY, UMZ, UEDEN, NVAR
  use prob_params_module, only: center, probhi
  use probdata_module, only: maxTaggingRadius, com_P, com_S, roche_rad_P, roche_rad_S  
  
  implicit none
  
  integer          :: lo(3),hi(3)
  integer          :: state_lo(3),state_hi(3)
  integer          :: tag_lo(3),tag_hi(3)
  double precision :: state(state_lo(1):state_hi(1), &
                            state_lo(2):state_hi(2), &
                            state_lo(3):state_hi(3),NVAR)
  integer          :: tag(tag_lo(1):tag_hi(1),tag_lo(2):tag_hi(2),tag_lo(3):tag_hi(3))
  double precision :: problo(3),dx(3),time
  integer          :: level,set,clear

  integer          :: i, j, k
  double precision :: x,y,z,r,r_P,r_S

  do k = lo(3), hi(3)
     z = problo(3) + (dble(k) + HALF)*dx(3)

     do j = lo(2), hi(2)
        y = problo(2) + (dble(j) + HALF)*dx(2)

        do i = lo(1), hi(1)
           x = problo(1) + (dble(i) + HALF)*dx(1)

           ! Tag all regions within the Roche radii of each star.
           ! We'll add a buffer around each star to double the Roche
           ! radius to ensure there aren't any sharp gradients.

           r_P = ( (x-com_P(1))**2 + (y-com_P(2))**2 + (z-com_P(3))**2 )**HALF
           r_S = ( (x-com_S(1))**2 + (y-com_S(2))**2 + (z-com_S(3))**2 )**HALF

           if (r_P <= TWO * roche_rad_P) then
              tag(i,j,k) = set
           endif
           
           if (r_S <= TWO * roche_rad_S) then
              tag(i,j,k) = set
           endif

           ! Clear all tagging that occurs outside the radius set by maxTaggingRadius.

           r = ((x-center(1))**2 + (y-center(2))**2 + (z-center(3))**2)**HALF

           if (r .gt. maxTaggingRadius * maxval(abs(problo-center)) .or. &
               r .gt. maxTaggingRadius * maxval(abs(probhi-center)) ) then

              tag(i,j,k) = clear

           endif

           ! We must ensure that the outermost zones are untagged due to the Poisson equation boundary conditions.
           ! This is necessary in addition to the above step because the problem center doesn't always coincide with the geometric center.

           if (x .lt. problo(1) + TWO * dx(1) .or. x .gt. probhi(1) - TWO * dx(1) .or. &
               y .lt. problo(2) + TWO * dx(2) .or. y .gt. probhi(2) - TWO * dx(2) .or. &
               z .lt. problo(3) + TWO * dx(3) .or. z .gt. probhi(3) - TWO * dx(3)) then

              tag(i,j,k) = clear

           endif

        enddo
     enddo
  enddo

end subroutine set_problem_tags
