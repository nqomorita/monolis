program main
  use mod_monolis
  implicit none
  integer(4) :: nnode, nelem
  integer(4), allocatable :: elem(:,:)
  character :: filename*128

  call monolis_get_mtx_arg(filename)
  write(*,*) "filename: ", trim(filename)

end program main
