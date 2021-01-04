module mod_monolis_io_mtx
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_util
  use mod_monolis_mesh
  use mod_monolis_stdlib
  implicit none

contains

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
