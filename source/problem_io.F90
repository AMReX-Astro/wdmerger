module problem_io_module

  implicit none

  ! For determining if we are the I/O processor.
  
  integer, save :: ioproc

  ! Probin file

  character (len=:), allocatable, save :: probin  

contains

  subroutine initialize_io(name, namlen)
    
    implicit none

    integer :: namlen, i
    integer :: name(namlen)
    
    ! Build "probin" filename -- the name of the file containing the fortin namelist.
    
    allocate(character(len=namlen) :: probin)
    do i = 1, namlen
       probin(i:i) = char(name(i))
    enddo

    ! Determine whether we are the I/O procoessor.
    
    call bl_pd_is_ioproc(ioproc)
    
  end subroutine initialize_io
  
end module problem_io_module
