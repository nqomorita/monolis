module mod_monolis_matvec
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_linalg_com
  implicit none

contains

  subroutine monolis_residual(monoCOM, monoMAT, X, B, R, tcomm)
    use mod_monolis_com
    use mod_monolis_mat
    use mod_monolis_linalg_com
    implicit none
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kind=kint) :: ndof
    integer(kind=kint) :: i
    real(kind=kdouble) :: X(:), B(:), R(:)
    real(kind=kdouble) :: t1, t2
    real(kind=kdouble), optional :: tcomm

  end subroutine monolis_residual

  subroutine monolis_matvec(monoCOM, monoMAT, X, Y, tcomm)
    use mod_monolis_com
    use mod_monolis_mat
    use mod_monolis_linalg_com
    implicit none
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kind=kint) :: ndof
    integer(kind=kint) :: i
    real(kind=kdouble) :: X(:), Y(:)
    real(kind=kdouble) :: t1, t2
    real(kind=kdouble), optional :: tcomm

  end subroutine monolis_matvec

end module mod_monolis_matvec