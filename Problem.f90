! problem-specific Fortran stuff goes here

subroutine problem_checkpoint(int_dir_name, len)

  ! called by the IO processor during checkpoint

  use bl_IO_module
  use com

  implicit none

  integer :: len
  integer :: int_dir_name(len)
  character (len=len) :: dir

  integer :: i, un

  ! dir will be the string name of the checkpoint directory
  do i = 1, len
     dir(i:i) = char(int_dir_name(i))
  enddo

  un = unit_new()
  open (unit=un, file=trim(dir)//"/COM", status="unknown")

100 format(1x, g30.20, 1x, g30.20)

  write (un,100) mass_left, mass_right
  write (un,100) com_xloc_l, com_xloc_r
  write (un,100) com_yloc_l, com_yloc_r
  write (un,100) com_zloc_l, com_zloc_r

  close (un)

end subroutine problem_checkpoint


subroutine problem_restart(int_dir_name, len)

  ! called by ALL processors during restart 

  use bl_IO_module
  use com

  implicit none

  integer :: len
  integer :: int_dir_name(len)
  character (len=len) :: dir

  integer :: i, un

  ! dir will be the string name of the checkpoint directory
  do i = 1, len
     dir(i:i) = char(int_dir_name(i))
  enddo

  un = unit_new()
  open (unit=un, file=trim(dir)//"/COM", status="old")

100 format(1x, g30.20, 1x, g30.20)

  ! here we read in, and this sets the values in the com module
  read (un,100) mass_left, mass_right
  read (un,100) com_xloc_l, com_xloc_r
  read (un,100) com_yloc_l, com_yloc_r
  read (un,100) com_zloc_l, com_zloc_r

  close (un)

end subroutine problem_restart
