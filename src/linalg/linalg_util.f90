module mod_monolis_linalg_util
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_linalg_com
  implicit none

contains

  subroutine monolis_update_R(monoCOM, monoMAT, ndof, X, tcomm)
    use mod_monolis_com
    use mod_monolis_mat
    use mod_monolis_linalg_com
    implicit none
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kind=kint) :: ndof
    integer(kind=kint) :: i
    real(kind=kdouble) :: X(:)
    real(kind=kdouble) :: t1, t2
    real(kind=kdouble), optional :: tcomm

  end subroutine monolis_update_R

end module mod_monolis_linalg_util