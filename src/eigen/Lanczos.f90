module mod_monolis_eigen_lanczos
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_solve
  use mod_monolis_eigen_lanczos_util
  use mod_monolis_util_debug

  implicit none

contains

  subroutine monolis_eigen_inverted_standard_lanczos &
    & (monolis, n_get_eigen, ths, maxiter, val, vec, is_bc)
    implicit none
    type(monolis_structure) :: monolis
    integer(kint) :: n_get_eigen, maxiter
    real(kdouble) :: ths, val(:), vec(:,:)
    logical, optional :: is_bc(:)

    call monolis_eigen_inverted_standard_lanczos_(monolis%PRM, monolis%COM, monolis%MAT, &
      & n_get_eigen, ths, maxiter, val, vec, is_bc)
  end subroutine monolis_eigen_inverted_standard_lanczos

  subroutine monolis_eigen_inverted_standard_lanczos_ &
    & (monoPRM, monoCOM, monoMAT, n_get_eigen, ths, maxiter, val, vec, is_bc)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kint) :: N, NP, NDOF, total_dof, j, k
    integer(kint) :: i, iter, maxiter, n_get_eigen
    real(kdouble) :: beta_t, ths, norm, tmp
    real(kdouble) :: vec(:,:), val(:)
    real(kdouble), allocatable :: p(:), q(:,:), alpha(:), beta(:), eigen_value(:), eigen_mode(:,:), prev(:)
    logical :: is_converge
    logical, optional :: is_bc(:)

    if(monoPRM%is_debug) call monolis_debug_header("monolis_eigen_inverted_standard_lanczos_")

    N     = monoMAT%N
    NP    = monoMAT%NP
    NDOF  = monoMAT%NDOF
    norm = 0.0d0
    is_converge = .false.

    total_dof = N*NDOF
    call monolis_allreduce_I1(total_dof, monolis_sum, monoCOM%comm)

    if(n_get_eigen > total_dof) n_get_eigen = total_dof

    allocate(alpha(maxiter), source = 0.0d0)
    allocate(beta(maxiter+1), source = 0.0d0)
    allocate(eigen_value(maxiter), source = 0.0d0)
    allocate(prev(maxiter), source = 0.0d0)
    allocate(q(NP*NDOF,0:maxiter+1), source = 0.0d0)
    allocate(p(NP*NDOF), source = 0.0d0)
    allocate(eigen_mode(NP*NDOF,n_get_eigen), source = 0.0d0)

    call lanczos_initialze(NP*NDOF, q(:,1), is_bc)

    do iter = 1, maxiter
      call monolis_set_RHS(monoMAT, q(:,iter))

      call monolis_solve_(monoPRM, monoCOM, monoMAT)

      call monolis_vec_AXPY(N, NDOF, -beta(iter), q(:,iter-1), monoMAT%X, p)

      call monolis_inner_product_R(monoCOM, N, NDOF, p, q(:,iter), alpha(iter))

      call monolis_vec_AXPY(N, NDOF, -alpha(iter), q(:,iter), p, p)

      call monolis_gram_schmidt(monoPRM, monoCOM, monoMAT, iter, q, p)

      call monolis_inner_product_R(monoCOM, N, NDOF, p, p, beta_t)

      beta(iter+1) = dsqrt(beta_t)
      beta_t = 1.0d0/beta(iter+1)
      do i = 1, NP*NDOF
        q(i,iter+1) = p(i)*beta_t
      enddo

      call monolis_get_eigen_pair_from_tridiag(iter, n_get_eigen, alpha, beta, q, eigen_value, eigen_mode)

      if(iter > n_get_eigen)then
        norm = 0.0d0
        do i = 1, n_get_eigen
          tmp = (prev(i) - eigen_value(i))**2/prev(i)**2
          norm = max(norm, tmp)
        enddo
        if(norm < ths) is_converge = .true.
      endif
      prev = eigen_value

write(*,"(a,i6,a,1p2e12.4)")"iter: ", iter, ", ths: ", norm

      if(is_converge .or. iter >= total_dof .or. iter == maxiter)then
        do i = 1, n_get_eigen
          val(i) = eigen_value(i)
          do j = 1, NP*NDOF
            vec(j,i) = eigen_mode(j,i)
          enddo
        enddo
        exit
      endif
    enddo
  end subroutine monolis_eigen_inverted_standard_lanczos_

end module mod_monolis_eigen_lanczos

  subroutine interface_monolis_eigen_inverted_standard_lanczos_ &
    & (monoPRM, monoCOM, monoMAT, NPNDOF, &
    & n_get_eigen, ths, maxiter, val, vec, is_bc)
    use mod_monolis_prm
    use mod_monolis_com
    use mod_monolis_mat
    use mod_monolis_eigen_lanczos
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kint) :: maxiter, n_get_eigen, NPNDOF
    real(kdouble) :: ths
    real(kdouble) :: vec(NPNDOF,n_get_eigen), val(n_get_eigen)
    logical, optional :: is_bc(NPNDOF)

    call monolis_eigen_inverted_standard_lanczos_(monoPRM, monoCOM, monoMAT, &
      & n_get_eigen, ths, maxiter, val, vec, is_bc)
  end subroutine interface_monolis_eigen_inverted_standard_lanczos_
