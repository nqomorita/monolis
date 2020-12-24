module mod_monolis_eigen_lanczos
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_solve
  use mod_monolis_eigen_lanczos_util

  implicit none

contains

  subroutine monolis_eigen_standard_lanczos(monolis)
    implicit none
    type(monolis_structure) :: monolis
    call monolis_eigen_standard_lanczos_(monolis%PRM, monolis%COM, monolis%MAT)
  end subroutine monolis_eigen_standard_lanczos

  subroutine monolis_eigen_standard_lanczos_(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kint) :: N, NP, NDOF, total_dof
    integer(kint) :: i, iter, maxiter
    real(kdouble) :: beta_t
    real(kdouble), allocatable :: p(:), q(:,:), alpha(:), beta(:), eigen_value(:), eigen_mode(:,:)

    if(monoPRM%is_debug) call monolis_debug_header("monolis_eigen_standard_lanczos_")

    N     = monoMAT%N
    NP    = monoMAT%NP
    NDOF  = monoMAT%NDOF

    total_dof = N*NDOF
    call monolis_allreduce_I1(total_dof, monolis_sum, monoCOM%comm)

    maxiter = 100

    allocate(alpha(maxiter), source = 0.0d0)
    allocate(beta(maxiter), source = 0.0d0)
    allocate(eigen_value(maxiter), source = 0.0d0)
    allocate(q(NP*NDOF,maxiter), source = 0.0d0)
    allocate(p(NP*NDOF), source = 0.0d0)
    allocate(eigen_mode(NP*NDOF,maxiter), source = 0.0d0)

    call lanczos_initialze(NP*NDOF, q(:,1))

    do iter = 1, maxiter
      write(*,*)"iter", iter

      call monolis_set_RHS(monoMAT, q(:,iter))

      call monolis_solve_(monoPRM, monoCOM, monoMAT)

      call monolis_vec_AXPY(N, NDOF, -beta(iter), q(:,iter), monoMAT%X, p)

      call monolis_inner_product_R(monoCOM, N, NDOF, p, q(:,iter), alpha(iter), monoPRM%tdotp, monoPRM%tcomm_dotp)

      call monolis_vec_AXPY(N, NDOF, -alpha(iter), q(:,iter), p, p)

      !call monolis_gram_schmidt(p)

      call monolis_inner_product_R(monoCOM, N, NDOF, p, p, beta_t, monoPRM%tdotp, monoPRM%tcomm_dotp)
      beta(iter+1) = dsqrt(beta_t)

      do i = 1, NP*NDOF
        q(i,iter+1) = p(i)/beta(iter+1)
      enddo

      call monolis_get_eigen_pair_from_tridiag(iter, alpha, beta, q, eigen_value, eigen_mode)

      !call check_converge_Lanczos()
      if(iter >= total_dof) exit
    enddo
  end subroutine monolis_eigen_standard_lanczos_

end module mod_monolis_eigen_lanczos
