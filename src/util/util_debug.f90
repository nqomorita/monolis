module mod_monolis_util_debug
  use mod_monolis_util
  use mod_monolis_prm
  use mod_monolis_com
  implicit none

  private
  public :: monolis_print_char
  public :: monolis_set_debug
  public :: monolis_debug_header
  public :: monolis_debug_int
  public :: monolis_debug_real
  public :: monolis_debug_char
  public :: monolis_debug_logical
  public :: monolis_warning_header
  public :: monolis_error_stop

  logical, save :: is_debug

contains

  !> debug section
  subroutine monolis_print_char(header, char)
    implicit none
    character(*) :: header, char
    if(monolis_global_myrank() == 0) write(*,"(a,a)")">> "//trim(header)//": ", trim(char)
  end subroutine monolis_print_char

  subroutine monolis_set_debug(flag)
    implicit none
    logical :: flag
    is_debug = flag
  end subroutine monolis_set_debug

  subroutine monolis_debug_header(header)
    implicit none
    character(*) :: header

    if(.not. is_debug) return
    if(monolis_global_myrank() == 0) write(*,"(a)")"** monolis debug: "//trim(header)
  end subroutine monolis_debug_header

  subroutine monolis_debug_int(header, n)
    implicit none
    integer(kint) :: n
    character(*) :: header

    if(.not. is_debug) return
    if(monolis_global_myrank() == 0) write(*,"(a,i12)")"** monolis debug: "//trim(header)//": ", n
  end subroutine monolis_debug_int

  subroutine monolis_debug_real(header, v)
    implicit none
    real(kdouble) :: v
    character(*) :: header

    if(.not. is_debug) return
    if(monolis_global_myrank() == 0) write(*,"(a,1pe12.4)")"** monolis debug: "//trim(header)//": ", v
  end subroutine monolis_debug_real

  subroutine monolis_debug_char(header, char)
    implicit none
    character(*) :: header, char

    if(.not. is_debug) return
    if(monolis_global_myrank() == 0) write(*,"(a,a)")"** monolis debug: "//trim(header)//": ", trim(char)
  end subroutine monolis_debug_char

  subroutine monolis_debug_logical(header, l)
    implicit none
    logical :: l
    character(*) :: header

    if(.not. is_debug) return
    if(monolis_global_myrank() == 0) write(*,"(a,l)")"** monolis debug: "//trim(header)//": ", l
  end subroutine monolis_debug_logical

  subroutine monolis_warning_header(header)
    implicit none
    character(*) :: header

    if(monolis_global_myrank() == 0) write(*,"(a)")"** monolis warning: "//trim(header)
  end subroutine monolis_warning_header

  subroutine monolis_error_stop()
    implicit none
    error stop monolis_fail
  end subroutine monolis_error_stop

end module mod_monolis_util_debug
