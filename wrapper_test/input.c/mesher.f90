program main
  use mod_monolis
  use mod_monolis_utils
  implicit none
  integer(kint) :: n_node, n_elem, i
  integer(kint), allocatable :: elem(:,:)
  real(kdouble), allocatable :: coef(:)
  character(monolis_charlen) :: filename
  logical :: is_get

  call monolis_global_initialize()

  call monolis_std_log_string("monolis matrix to mesh")

  call monolis_get_arg_input_i_tag(filename, is_get)

  call monolis_input_mtx_R(filename, n_node, n_elem, elem, coef)

  open(20, file = "node.dat", status = "replace")
    write(20,"(i8,i8)")n_node
    do i = 1, n_node
      write(20,*)"0., 0., 0."
    enddo
  close(20)

  open(20, file = "elem.dat", status = "replace")
    write(20,"(i8,i8)")n_elem, 2
    do i = 1, n_elem
      write(20,"(2i10)")elem(1,i) - 1, elem(2,i) - 1
    enddo
  close(20)

  open(20, file = "coef.dat", status = "replace")
    write(20,"(i8,i8)")n_elem
    do i = 1, n_elem
      write(20,"(1p2e23.15)")coef(i)
    enddo
  close(20)

  call monolis_global_finalize()

  contains
end program main
