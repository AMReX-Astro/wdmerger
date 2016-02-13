module problem_tagging_module

  implicit none

  public

contains

  subroutine set_problem_tags(tag,tag_lo,tag_hi, &
                              state,state_lo,state_hi, &
                              set,clear,&
                              lo,hi,&
                              dx,problo,time,level) bind(C,name='set_problem_tags')

    use bl_constants_module, only: ZERO, HALF, TWO
    use meth_params_module, only: NVAR, URHO, UTEMP
    use prob_params_module, only: center, probhi, dim, Symmetry, physbc_lo, physbc_hi
    use probdata_module, only: max_tagging_radius, &
                               max_stellar_tagging_level, &
                               max_temperature_tagging_level, &
                               roche_tagging_factor, &
                               stellar_density_threshold, &
                               temperature_tagging_threshold, &
                               com_P, com_S, roche_rad_P, roche_rad_S

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

    integer          :: i, j, k, n
    double precision :: loc(3), r, r_P, r_S

    logical          :: outer_boundary_test(3)

    do k = lo(3), hi(3)
       loc(3) = problo(3) + (dble(k) + HALF) * dx(3)

       do j = lo(2), hi(2)
          loc(2) = problo(2) + (dble(j) + HALF) * dx(2)

          do i = lo(1), hi(1)
             loc(1) = problo(1) + (dble(i) + HALF) * dx(1)

             if (level < max_stellar_tagging_level) then

                if (level == 0) then

                   ! On the coarse grid, tag all regions within the Roche radii of each star.
                   ! We'll add a buffer around each star to double the Roche
                   ! radius to ensure there aren't any sharp gradients in regions of
                   ! greater than ambient density.

                   r_P = ( sum((loc-com_P)**2) )**HALF
                   r_S = ( sum((loc-com_S)**2) )**HALF

                   if (r_P <= roche_tagging_factor * roche_rad_P) then
                      tag(i,j,k) = set
                   endif

                   if (r_S <= roche_tagging_factor * roche_rad_S) then
                      tag(i,j,k) = set
                   endif

                else if (level >= 1) then

                   ! On more refined levels, tag all regions within the stars themselves (defined as 
                   ! areas where the density is greater than some threshold).

                   if (state(i,j,k,URHO) > stellar_density_threshold) then

                      tag(i,j,k) = set

                   endif

                endif

             endif

             ! Tag all zones at all levels that are hotter than a specified temperature threshold.

             if (level < max_temperature_tagging_level) then

                if (state(i,j,k,UTEMP) > temperature_tagging_threshold) then

                   tag(i,j,k) = set

                endif

             endif

             ! Clear all tagging that occurs outside the radius set by max_tagging_radius.

             r = ( sum((loc-center)**2) )**HALF

             if (r .gt. max_tagging_radius * maxval(abs(problo-center)) .or. &
                 r .gt. max_tagging_radius * maxval(abs(probhi-center)) ) then

                tag(i,j,k) = clear

             endif

             ! We must ensure that the outermost zones are untagged
             ! due to the Poisson equation boundary conditions.
             ! This is necessary in addition to the above step
             ! because the problem center doesn't always coincide
             ! with the geometric center. However we can avoid this
             ! restriction on symmetric boundaries.

             outer_boundary_test = .false.

             do n = 1, dim

                if ((physbc_lo(n) .ne. Symmetry) .and. (loc(n) .lt. problo(n) + TWO * dx(1))) then
                   outer_boundary_test = .true.
                endif

                if ((physbc_hi(n) .ne. Symmetry) .and. (loc(n) .gt. probhi(n) - TWO * dx(1))) then
                   outer_boundary_test = .true.
                endif

             enddo

             if ( any(outer_boundary_test) ) then

                tag(i,j,k) = clear

             endif

          enddo
       enddo
    enddo

  end subroutine set_problem_tags

end module problem_tagging_module
