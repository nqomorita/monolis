module mod_monolis_util
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat

  implicit none

  private
  public :: monolis_initialize
  public :: monolis_finalize
  public :: monolis_timer_initialize
  public :: monolis_timer_finalize
  public :: monolis_check_diagonal
  public :: monolis_debug_header

contains

  subroutine monolis_initialize(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT

    call monolis_prm_initialize(monoPRM)
    call monolis_com_initialize(monoCOM)
    call monolis_mat_initialize(monoMAT)
  end subroutine monolis_initialize

  subroutine monolis_finalize(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT

    call monolis_prm_finalize(monoPRM)
    call monolis_com_finalize(monoCOM)
    call monolis_mat_finalize(monoMAT)
  end subroutine monolis_finalize

  subroutine monolis_timer_initialize(monoPRM)
    implicit none
    type(monolis_prm) :: monoPRM

    if(monoPRM%is_debug) call monolis_debug_header("monolis_timer_initialize")

    monoPRM%tsol  = 0.0d0
    monoPRM%tspmv = 0.0d0
    monoPRM%tprec = 0.0d0
    monoPRM%tcomm = 0.0d0
  end subroutine monolis_timer_initialize

  subroutine monolis_timer_finalize(monoPRM, monoCOM)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM

    if(monoPRM%is_debug) call monolis_debug_header("monolis_timer_finalize")

    !if(monoCOM%myrank == 0) write(*,"(a,i8,1p4e12.5)")" ** monolis solved:", 0, tsol, tspmv, tprec, tcomm
  end subroutine monolis_timer_finalize

  subroutine monolis_check_diagonal(monoPRM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_mat) :: monoMAT
    integer(kind=kint) :: i, j, k, jS, jE, in, kn, NP, NDOF, NDOF2

    if(monoPRM%is_debug) call monolis_debug_header("monolis_check_diagonal")

    NP =  monoMAT%NP
    NDOF  = monoMAT%NDOF
    NDOF2 = NDOF*NDOF

    do i = 1, NP
      jS = monoMAT%index(i-1) + 1
      jE = monoMAT%index(i)
      do j = jS, jE
        in = monoMAT%item(j)
        if(i == in)then
          do k = 1, NDOF
            kn = NDOF2*(j-1) + (NDOF+1)*(k-1) + 1
            if(monoMAT%A(kn) == 0.0d0) stop " ** monolis error: zero diagonal"
          enddo
        endif
      enddo
    enddo
  end subroutine monolis_check_diagonal

  subroutine monolis_debug_header(header)
    implicit none
    character(*) :: header

    write(*,"(a)")"** monolis debug: "//trim(header)
  end subroutine monolis_debug_header

  subroutine monolis_debug_equal_R(header, a, b)
    implicit none
    character(*) :: header
    real(kind=kdouble) :: a, b

    if(a /= b)then
      write(*,"(a,1pe12.5,a,1pe12.5,a)")"** monolis debug: ", a, "is not equal", b, " at "//trim(header)
    endif
  end subroutine monolis_debug_equal_R

end module mod_monolis_util