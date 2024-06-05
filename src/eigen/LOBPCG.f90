module mod_monolis_eigen_LOBPCG
  use mod_monolis_utils
  use mod_monolis_def_solver
  use mod_monolis_def_mat
  use mod_monolis_solve
  use mod_monolis_converge
  use mod_monolis_eigen_lanczos_util
  use blopex_fortran_hold_vars

  implicit none

contains

  subroutine monolis_eigen_standard_lobpcg(monolis, monoCOM, n_get_eigen, ths, maxiter, &
    & val, vec, is_bc)
    implicit none
    type(monolis_structure) :: monolis
    type(monolis_com) :: monoCOM
    integer(kint) :: n_get_eigen, maxiter
    real(kdouble) :: ths, val(:), vec(:,:)
    logical :: is_bc(:)

    call monolis_eigen_standard_lobpcg_main(monolis%PRM, monoCOM, monolis%MAT, &
      & n_get_eigen, ths, maxiter, val, vec, is_bc)
  end subroutine monolis_eigen_standard_lobpcg

  subroutine monolis_eigen_standard_lobpcg_main(monoPRM, monoCOM, monoMAT, n_get_eigen, ths, maxiter, &
    & val, vec, is_bc)
    implicit none
    type(monolis_prm), target :: monoPRM
    type(monolis_com), target :: monoCOM
    type(monolis_mat), target :: monoMAT
    type(monolis_mat) :: monoPREC
    integer(kint) :: maxiter, n_get_eigen, n, loglevel
    real(kdouble) :: ths
    real(kdouble) :: val(:), vec(:,:)
    logical :: is_bc(:)

    n = monoMAT%NP*monoMAT%NDOF
    N_hold = n
    M_hold = n_get_eigen
    monoPRM_hold => monoPRM
    monoCOM_hold => monoCOM
    monoMAT_hold => monoMAT

    call monolis_precond_setup(monoPRM, monoCOM, monoMAT, monoPREC)

    loglevel = 0
    call blopex_lobpcg_solve(n, n_get_eigen, maxiter, ths, loglevel, val, vec)
  end subroutine monolis_eigen_standard_lobpcg_main
end module mod_monolis_eigen_LOBPCG
