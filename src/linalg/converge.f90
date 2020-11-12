module mod_monolis_converge
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_linalg
  use mod_monolis_linalg_util
  use mod_monolis_linalg_com
  use mod_monolis_util
  implicit none

contains

  subroutine monolis_set_converge(monoPRM, monoCOM, monoMAT, B, B2, is_converge, tdotp, tcomm)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    real(kind=kdouble) :: B(:), B2
    real(kind=kdouble) :: tdotp, tcomm
    logical :: is_converge

    if(monoPRM%is_debug) call monolis_debug_header("monolis_set_converge")

    is_converge = .false.
    call monolis_inner_product_R(monoCOM, monoMAT%N, monoMAT%NDOF, B, B, B2, tdotp, tcomm)

    if(B2 == 0.0d0)then
      if(monoCOM%myrank == 0) write (*,"(a,1pe16.6)")" ** monolis warning: bnorm ", B2
      monoMAT%X = 0.0d0
      is_converge = .true.
    endif

  end subroutine monolis_set_converge

  subroutine monolis_check_converge(monoPRM, monoCOM, monoMAT, R, B2, iter, is_converge, tdotp, tcomm)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kind=kint) :: iter
    real(kind=kdouble) :: R(:), R2, B2, resid
    real(kind=kdouble) :: tdotp, tcomm
    logical :: is_converge

    is_converge = .false.
    call monolis_inner_product_R(monoCOM, monoMAT%N, monoMAT%NDOF, R, R, R2, tdotp, tcomm)
    resid = dsqrt(R2/B2)

    monoPRM%curiter = iter
    monoPRM%curresid = resid
    if(monoCOM%myrank == 0 .and. monoPRM%show_iterlog) write (*,"(i7, 1pe16.6)") iter, resid
    if(resid < monoPRM%tol .or. iter == monoPRM%maxiter)then
      is_converge = .true.
      if(monoCOM%myrank == 0 .and. monoPRM%is_debug) write (*,"(i7, 1pe16.6)") iter, resid
    endif

  end subroutine monolis_check_converge

  subroutine monolis_check_converge_2(monoPRM, monoCOM, monoMAT, R2, B2, iter, is_converge)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kind=kint) :: iter
    real(kind=kdouble) :: R2, B2, resid
    logical :: is_converge

    is_converge = .false.
    resid = dsqrt(R2/B2)
    if(monoCOM%myrank == 0 .and. monoPRM%show_iterlog) write (*,"(i7, 1pe16.6)") iter, resid
    if(resid < monoPRM%tol) is_converge = .true.

  end subroutine monolis_check_converge_2
end module mod_monolis_converge