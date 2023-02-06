module mod_monolis_eigen_LOBPCG
  use mod_monolis_prm
  use mod_monolis_mat
  use mod_monolis_solve
  use mod_monolis_converge
  use mod_monolis_eigen_lanczos_util

  implicit none

contains

  subroutine monolis_eigen_inverted_lobpcg(monolis, n_get_eigen, ths, maxiter, vec)
    implicit none
    type(monolis_structure) :: monolis
    integer(kint) :: n_get_eigen, maxiter
    real(kdouble) :: ths, vec(:,:)

    call monolis_eigen_inverted_lobpcg_mat (monolis%PRM, monolis%COM, monolis%MAT, &
      & n_get_eigen, ths, maxiter, vec)

    !call monolis_eigen_inverted_lobpcg_ (monolis%PRM, monolis%COM, monolis%MAT, n_get_eigen, ths)
  end subroutine monolis_eigen_inverted_lobpcg

  subroutine monolis_eigen_inverted_lobpcg_mat(monoPRM, monoCOM, monoMAT, n_get_eigen, ths, maxiter, vec)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kint), parameter :: XX = 0 !> A
    integer(kint), parameter :: RR = 1 !> A
    integer(kint), parameter :: PP = 2 !> A
    integer(kint) :: N, NP, NDOF, NG, total_dof
    integer(kint) :: i, j, iter, maxiter,  n_get_eigen
    real(kdouble) :: ths, mu
    real(kdouble) :: vec(:,:)
    real(kdouble), allocatable :: X(:,:), R(:,:), P(:,:), XAX(:,:), XBX(:,:), evec(:,:), eval(:), T(:)
    real(kdouble), allocatable :: A(:,:), B(:,:), R0(:), R2(:), resid(:)
  end subroutine monolis_eigen_inverted_lobpcg_mat

  subroutine monolis_eigen_inverted_lobpcg_(monoPRM, monoCOM, monoMAT, n_get_eigen, ths)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kint), parameter :: W = 1 !> A
    integer(kint), parameter :: X = 2 !> A
    integer(kint), parameter :: P = 3 !> A
    integer(kint), parameter :: V = 1 !> B
    integer(kint), parameter :: Y = 2 !> B
    integer(kint), parameter :: Q = 3 !> B
    integer(kint) :: N, NP, NDOF, total_dof
    integer(kint) :: i, iter, n_get_eigen
    real(kdouble) :: ths, mu, Sa(3,3), Sb(3,3), lambda, coef(3), norm_x, norm_p, R0, R2, resid
    real(kdouble), allocatable :: A(:,:), B(:,:)
  end subroutine monolis_eigen_inverted_lobpcg_

end module mod_monolis_eigen_LOBPCG
