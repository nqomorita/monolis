program main
  use mod_monolis
  implicit none
  integer(4) :: nnode, nelem, i
  integer(4), allocatable :: elem(:,:)
  real(kdouble), allocatable :: coef(:)
  character :: filename*128

  call monolis_std_debug_log_flag(.true.)

  call monolis_get_mtx_arg(filename)
  write(*,*) "filename: ", trim(filename)

  call monolis_input_mtx(filename, nnode, nelem, elem, coef)

  open(20, file = "node.dat", status = "replace")
    write(20,"(i8,i8)")nnode
    do i = 1, nnode
      write(20,*)"0., 0., 0."
    enddo
  close(20)

  open(20, file = "elem.dat", status = "replace")
    write(20,"(i8,i8)")nelem, 2
    do i = 1, nelem
      write(20,"(2i10)")elem(1,i), elem(2,i)
    enddo
  close(20)

  open(20, file = "coef.dat", status = "replace")
    write(20,"(i8,i8)")nelem
    do i = 1, nelem
      write(20,"(1pe22.15)")coef(i)
    enddo
  close(20)
end program main
