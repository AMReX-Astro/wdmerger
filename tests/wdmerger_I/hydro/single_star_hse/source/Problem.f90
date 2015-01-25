! problem-specific Fortran stuff goes here

subroutine problem_checkpoint(int_dir_name, len)

  ! called by the IO processor during checkpoint

  implicit none

  integer :: len
  integer :: int_dir_name(len)
  character (len=len) :: dir

  integer :: i

  ! dir will be the string name of the checkpoint directory
  do i = 1, len
     dir(i:i) = char(int_dir_name(i))
  enddo



end subroutine problem_checkpoint


subroutine problem_restart(int_dir_name, len)

  ! called by ALL processors during restart 

  implicit none

  integer :: len
  integer :: int_dir_name(len)
  character (len=len) :: dir

  integer :: i

  ! dir will be the string name of the checkpoint directory
  do i = 1, len
     dir(i:i) = char(int_dir_name(i))
  enddo

end subroutine problem_restart

! ::
! :: ----------------------------------------------------------
! ::  volumeInDensityBoundary
! ::
! :: INPUTS / OUTPUTS:
! ::  rho        => density field
! ::  rlo,rhi    => index limits of rho array
! ::  lo,hi      => index limits of grid interior
! ::  dx         => cell size
! ::  rho_cutoff => cutoff density
! ::  vol        <= total volume with density greater than the cutoff
! :: ----------------------------------------------------------
! ::

subroutine ca_volumeindensityboundary(rho,r_l1,r_l2,r_l3,r_h1,r_h2,r_h3,lo,hi,dx,vol,rho_cutoff)

  use probdata_module, only : center
  use bl_constants_module

  implicit none
  integer          :: r_l1,r_l2,r_l3,r_h1,r_h2,r_h3
  integer          :: lo(3), hi(3)
  double precision :: vol, rho_cutoff, dx(3)
  double precision :: rho(r_l1:r_h1,r_l2:r_h2,r_l3:r_h3)

  integer          :: i, j, k
  double precision :: dV

  dV  = dx(1)*dx(2)*dx(3)
  vol = ZERO

  !$OMP PARALLEL DO PRIVATE(i,j,k) REDUCTION(+:vol)
  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           if (rho(i,j,k) > rho_cutoff) vol = vol + dV
        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO

end subroutine ca_volumeindensityboundary
