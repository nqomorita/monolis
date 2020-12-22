module mod_monolis_eigen_lanczos
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_solve

  implicit none

contains

  subroutine monolis_eigen_standard_lanczos(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kint) :: N, NP, NDOF, NNDOF
    integer(kint) :: i, iter, maxiter
    real(kdouble) :: alpha, beta
    real(kdouble), allocatable :: p(:), q(:)

    if(monoPRM%is_debug) call monolis_debug_header("monolis_eigen_standard_lanczos")

    N     = monoMAT%N
    NP    = monoMAT%NP
    NDOF  = monoMAT%NDOF
    NNDOF = N*NDOF
    beta = 0.0d0

    do iter = 1, maxiter
      !call monolis_set_RHS(q)
      call monolis_solve_(monoPRM, monoCOM, monoMAT)

      call monolis_vec_AXPY(N, NDOF, -beta, q, monoMAT%X, p)

      call monolis_inner_product_R(monoCOM, N, NDOF, p, q, alpha, monoPRM%tdotp, monoPRM%tcomm_dotp)

      call monolis_vec_AXPY(N, NDOF, -alpha, q, p, p)

      call monolis_inner_product_R(monoCOM, N, NDOF, p, p, beta, monoPRM%tdotp, monoPRM%tcomm_dotp)
      beta = dsqrt(beta)

      do i = 1, NNDOF
        q(i) = p(i)/beta
      enddo
    enddo
  end subroutine monolis_eigen_standard_lanczos

end module mod_monolis_eigen_lanczos
