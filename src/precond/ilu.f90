module mod_monolis_precond_ilu
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_fact_LU
  use mod_monolis_matrix_fillin

  implicit none

  type(monolis_mat), save :: monoTREE

contains

  subroutine  monolis_precond_ilu_setup(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT

    call monolis_precond_ilu_init(monoPRM, monoCOM, monoMAT, monoTREE)
    call monolis_init_LU_inner(monoPRM, monoCOM, monoTREE)
    call monolis_fact_LU_inner(monoPRM, monoCOM, monoTREE)
  end subroutine monolis_precond_ilu_setup

  subroutine monolis_precond_ilu_apply(monoPRM, monoCOM, monoMAT, X, Y)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    real(kind=kdouble) :: X(:), Y(:)
    integer(kind=kint) :: i, N, NDOF

    N = monoMAT%N
    NDOF = monoMAT%NDOF
    do i = 1, N*NDOF
      monoTREE%B(i) = X(i)
    enddo
    call monolis_solv_LU_inner(monoPRM, monoCOM, monoTREE)
    do i = 1, N*NDOF
      Y(i) = monoTREE%X(i)
    enddo
  end subroutine monolis_precond_ilu_apply

  subroutine monolis_precond_ilu_clear(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT

    call monolis_clear_LU_inner(monoPRM, monoCOM, monoTREE)
  end subroutine monolis_precond_ilu_clear

  subroutine monolis_precond_ilu_init(monoPRM, monoCOM, monoMAT, monoTREE)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    type(monolis_mat) :: monoTREE
    integer(kind=kint) :: N, NDOF
    logical :: is_fillin = .false.
    logical :: is_asym = .false.

    N = monoMAT%N
    NDOF = monoMAT%NDOF
    monoTREE%N = monoMAT%N
    monoTREE%NP = monoMAT%NP
    monoTREE%NDOF = monoMAT%NDOF
    allocate(monoTREE%B(N*NDOF))
    allocate(monoTREE%X(N*NDOF))
    monoTREE%B = 0.0d0
    monoTREE%X = 0.0d0
    call monolis_matrix_get_fillin(monoPRM, monoCOM, monoMAT, monoTREE, is_fillin, is_asym)
    call monolis_matrix_copy_with_fillin(monoPRM, monoCOM, monoMAT, monoTREE, is_asym)
  end subroutine monolis_precond_ilu_init
end module mod_monolis_precond_ilu