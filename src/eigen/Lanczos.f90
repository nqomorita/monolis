module mod_monolis_eigen_lanczos
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_solve
  use mod_monolis_eigen_lanczos_util

  implicit none

contains

  subroutine monolis_eigen_inverted_standard_lanczos &
    & (monolis, n_get_eigen, ths, maxiter, val, vec)
    implicit none
    type(monolis_structure) :: monolis
    integer(kint) :: n_get_eigen, maxiter
    real(kdouble) :: ths, val(:), vec(:,:)

    call monolis_eigen_inverted_standard_lanczos_(monolis%PRM, monolis%COM, monolis%MAT, &
      & n_get_eigen, ths, maxiter, val, vec)
  end subroutine monolis_eigen_inverted_standard_lanczos

  subroutine monolis_eigen_inverted_standard_lanczos_ &
    & (monoPRM, monoCOM, monoMAT, n_get_eigen, ths, maxiter, val, vec)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kint) :: N, NP, NDOF, total_dof
    integer(kint) :: i, iter, maxiter, n_get_eigen
    real(kdouble) :: beta_t, ths
    real(kdouble) :: vec(:,:), val(:)
    real(kdouble), allocatable :: p(:), q(:,:), alpha(:), beta(:), eigen_value(:), eigen_mode(:,:)

    if(monoPRM%is_debug) call monolis_debug_header("monolis_eigen_inverted_standard_lanczos_")

    N     = monoMAT%N
    NP    = monoMAT%NP
    NDOF  = monoMAT%NDOF

    total_dof = N*NDOF
    call monolis_allreduce_I1(total_dof, monolis_sum, monoCOM%comm)

    if(n_get_eigen > total_dof) n_get_eigen = total_dof

    allocate(alpha(maxiter), source = 0.0d0)
    allocate(beta(maxiter+1), source = 0.0d0)
    allocate(eigen_value(maxiter), source = 0.0d0)
    allocate(q(NP*NDOF,0:maxiter+1), source = 0.0d0)
    allocate(p(NP*NDOF), source = 0.0d0)
    allocate(eigen_mode(NP*NDOF,maxiter), source = 0.0d0)

    call lanczos_initialze(NP*NDOF, q(:,1))

    do iter = 1, maxiter
      call monolis_set_RHS(monoMAT, q(:,iter))

      call monolis_solve_(monoPRM, monoCOM, monoMAT)

      call monolis_vec_AXPY(N, NDOF, -beta(iter), q(:,iter-1), monoMAT%X, p)

      call monolis_inner_product_R(monoCOM, N, NDOF, p, q(:,iter), alpha(iter), monoPRM%tdotp, monoPRM%tcomm_dotp)

      call monolis_vec_AXPY(N, NDOF, -alpha(iter), q(:,iter), p, p)

      call monolis_gram_schmidt(monoPRM, monoCOM, monoMAT, iter, q, p)

      call monolis_inner_product_R(monoCOM, N, NDOF, p, p, beta_t, monoPRM%tdotp, monoPRM%tcomm_dotp)
      beta(iter+1) = dsqrt(beta_t)

write(*,"(a,i6,a,1p2e12.4)")"iter: ", iter, ", beta: ", beta(iter+1), beta(iter+1)/beta(2)

      beta_t = 1.0d0/beta(iter+1)
      do i = 1, NP*NDOF
        q(i,iter+1) = p(i)*beta_t
      enddo

      if((iter >= n_get_eigen .and. beta(iter+1)/beta(2) < ths) .or. &
       & iter >= total_dof .or. iter == maxiter)then
        call monolis_get_eigen_pair_from_tridiag(iter, alpha, beta, q, eigen_value, eigen_mode)

!write(*,*)"eigen_value"
        do i = 1, n_get_eigen
          !write(*,"(1p2e12.5)")1.0d0/eigen_value(i)
          val(i) = 1.0d0/eigen_value(i)
        enddo
!write(*,*)"e_mode"
!do i = 1, iter
!  write(*,"(1p10e12.5)")eigen_mode(:,i)
!enddo
        vec = eigen_mode(:,1:n_get_eigen)
        exit
      endif
    enddo
  end subroutine monolis_eigen_inverted_standard_lanczos_

end module mod_monolis_eigen_lanczos
