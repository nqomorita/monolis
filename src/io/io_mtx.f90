module mod_monolis_io_mtx
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_util
  use mod_monolis_util_debug
  use mod_monolis_mesh
  use mod_monolis_stdlib
  implicit none

contains

  subroutine monolis_input_mtx(fname, nnode, nelem, elem, coef)
    implicit none
    integer(kint) :: nnode, nelem, i, in, ierr
    integer(kint), allocatable :: elem(:,:)
    real(kdouble), allocatable :: coef(:)
    character :: fname*100, ctemp*100
    logical :: is_first

    !> count lines
    is_first = .true.
    open(20, file = trim(fname), status = "old")
    in = 0

    do
      read(20, "(a)", iostat = ierr) ctemp
      if(0 /= ierr) exit

      if(ctemp(1:1) == "%") cycle

      if(is_first)then
        backspace(20)
        read(20,*) nnode, i, nelem
        call monolis_debug_int("nnode", nnode)
        call monolis_debug_int("nelem", nelem)

        allocate(elem(2,nelem), source = 0)
        allocate(coef(nelem), source = 0.0d0)
        is_first = .false.
      else
        in = in + 1
        backspace(20)
        read(20,*) elem(1,in), elem(2,in), coef(in)
      endif
      if(in == nelem) exit
    enddo
    !write(*,"(a,i8)")"nelem: ", nelem
    close(20)

  end subroutine monolis_input_mtx

  subroutine monolis_get_mtx_arg(filename)
    implicit none
    integer(kint) :: i, count
    character :: argc1*128, argc2*128, filename*128

    call monolis_debug_header("monolis_get_mtx_arg")
    count = iargc()

    if(mod(count,2) /= 0) stop "* monolis monolis_get_mtx_arg"
    do i = 1, count/2
      call getarg(2*i-1, argc1)
      call getarg(2*i  , argc2)
      if(trim(argc1) == "-i")then
        filename = trim(argc2)
      else
        write(*,"(a)")"* monolis input arg error"
        stop
      endif
    enddo
  end subroutine monolis_get_mtx_arg

end module mod_monolis_io_mtx
