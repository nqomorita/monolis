module mod_monolis_precond_ROM
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat

  implicit none

  integer(kint), save :: n_get_eigen
  real(kdouble), save, allocatable :: val(:), vec(:,:)

contains

  subroutine monolis_precond_ROM_setup(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kint) :: maxiter
    real(kdouble) :: ths
    !logical, optional :: is_bc(:)

    maxiter = 10
    n_get_eigen = 100
    ths = 1.0d-4
    allocate(val(n_get_eigen), source = 0.0d0)
    allocate(vec(monoMAT%NP*monoMAT%NDOF, n_get_eigen), source = 0.0d0)

    call monolis_eigen_inverted_standard_lanczos_ &
    & (monoPRM, monoCOM, monoMAT, n_get_eigen, ths, maxiter, val, vec)
  end subroutine monolis_precond_ROM_setup

  subroutine monolis_precond_ROM_apply(monoPRM, monoCOM, monoMAT, X, Y)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kint) :: i, j
    real(kdouble) :: X(:), Y(:), phi
    real(kdouble), allocatable :: coef(:)

    allocate(coef(n_get_eigen), source = 0.0d0)

    do i = 1, n_get_eigen
      phi = 0.0d0
      do j = 1, monoMAT%NP*monoMAT%NDOF
        phi = phi + vec(j,i)*X(i)
      enddo
      coef(i) = phi*val(i)
    enddo

    Y = 0.0d0
    do i = 1, n_get_eigen
      do j = 1, monoMAT%NP*monoMAT%NDOF
        Y(j) = Y(j) + vec(j,i)*coef(i)
      enddo
    enddo
  end subroutine monolis_precond_ROM_apply

  subroutine monolis_precond_ROM_clear(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT

  end subroutine monolis_precond_ROM_clear
end module mod_monolis_precond_ROM