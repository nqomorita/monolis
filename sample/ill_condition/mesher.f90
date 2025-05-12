program main
  use mod_monolis
  use mod_monolis_utils
  implicit none
  integer(kint) :: n_node, n_elem, i, in
  real(kdouble) :: gamma
  integer(kint), allocatable :: elem(:,:)
  real(kdouble), allocatable :: coef(:)
  logical :: is_get

  call monolis_global_initialize()

  call monolis_std_log_string("monolis ill-condition matrix")

  call monolis_get_arg_input_n_tag(n_node, is_get)
  call monolis_get_arg_input_R("-g", gamma, is_get)

  if(.not. is_get) n_node = 5
  if(n_node < 5) stop "n > 5"
  n_elem = n_node*3 - 3

  call monolis_alloc_I_2d(elem, 2, n_elem)
  call monolis_alloc_R_1d(coef, n_elem)

  !> make matrix
  in = 0
  do i = 1, n_node
    in = in + 1
    elem(1,in) = i
    elem(2,in) = i
    coef(in) = 2.0d0
  enddo

  do i = 1, n_node - 1
    in = in + 1
    elem(1,in) = i
    elem(2,in) = i + 1
    coef(in) = gamma
  enddo

  do i = 1, n_node - 2
    in = in + 1
    elem(1,in) = i + 2
    elem(2,in) = i
    coef(in) = 1.0d0
  enddo

  !> output
  open(20, file = "node.f.dat", status = "replace")
    write(20,"(i8,i8)")n_node
    do i = 1, n_node
      write(20,*)"0., 0., 0."
    enddo
  close(20)

  open(20, file = "node.c.dat", status = "replace")
    write(20,"(i8,i8)")n_node
    do i = 1, n_node
      write(20,*)"0., 0., 0."
    enddo
  close(20)

  open(20, file = "elem.f.dat", status = "replace")
    write(20,"(i8,i8)")n_elem, 2
    do i = 1, n_elem
      write(20,"(2i10)")elem(1,i), elem(2,i)
    enddo
  close(20)

  open(20, file = "elem.c.dat", status = "replace")
    write(20,"(i8,i8)")n_elem, 2
    do i = 1, n_elem
      write(20,"(2i10)")elem(1,i) - 1, elem(2,i) - 1
    enddo
  close(20)

  open(20, file = "coef.f.dat", status = "replace")
    write(20,"(i8)")n_elem
    do i = 1, n_elem
      write(20,"(i8,i8,1pe12.5)")elem(1,i), elem(2,i), coef(i)
    enddo
  close(20)

  open(20, file = "coef.c.dat", status = "replace")
    write(20,"(i8)")n_elem
    do i = 1, n_elem
      write(20,"(i8,i8,1pe12.5)")elem(1,i) - 1, elem(2,i) - 1, coef(i)
    enddo
  close(20)

  call monolis_global_finalize()

end program main
