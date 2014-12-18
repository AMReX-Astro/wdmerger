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

  write (un,100) mass_p, mass_s
  write (un,100) com_loc_p(1), com_loc_s(1)
  write (un,100) com_loc_p(2), com_loc_s(2)
  write (un,100) com_loc_p(3), com_loc_s(3)

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
  read (un,100) mass_p, mass_s
  read (un,100) com_loc_p(1), com_loc_s(1)
  read (un,100) com_loc_p(2), com_loc_s(2)
  read (un,100) com_loc_p(3), com_loc_s(3)

  close (un)

end subroutine problem_restart
