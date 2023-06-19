!> Lanczos 法モジュール
module mod_monolis_eigen_lanczos
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc
  use mod_monolis_eigen_lanczos_util
  use mod_monolis_inner_product
  use mod_monolis_matvec
  use mod_monolis_vec_util
  use mod_monolis_solve

  implicit none

contains

  !> @ingroup eigen
  !> Lanczos 法（逆反復、最小固有値、実数型、メイン関数）
  subroutine monolis_eigen_inverted_standard_lanczos_R_main( &
    & monoPRM, monoCOM, monoMAT, monoPREC, n_get_eigen, ths, maxiter, val, vec, is_bc)
    implicit none
    !> [in,out] パラメータ構造体
    type(monolis_prm), intent(inout) :: monoPRM
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in,out] 行列構造体
    type(monolis_mat), intent(inout) :: monoMAT
    !> [in,out] 前処理構造体
    type(monolis_mat), intent(inout) :: monoPREC
    !> [in,out] 取得固有値数
    integer(kint), intent(inout) :: n_get_eigen
    !> [in] 収束判定閾値
    real(kdouble), intent(in) :: ths
    !> [in] 最大反復回数
    integer(kint), intent(in) :: maxiter
    !> [out] 固有値
    real(kdouble), intent(out) :: val(:)
    !> [out] 固有ベクトル
    real(kdouble), intent(out) :: vec(:,:)
    !> [in] Dirhchlet 境界条件判定フラグ
    logical, intent(in) :: is_bc(:)
    integer(kint) :: N, NP, NDOF, total_dof, j, k
    integer(kint) :: i, iter, n_bc
    real(kdouble) :: beta_t, norm, tmp
    real(kdouble), allocatable :: p(:), q(:,:), alpha(:), beta(:), eigen_value(:), eigen_mode(:,:), prev(:)
    logical :: is_converge

    call monolis_std_debug_log_header("monolis_eigen_inverted_standard_lanczos_R_main")

    N     = monoMAT%N
    NP    = monoMAT%NP
    NDOF  = monoMAT%NDOF
    norm = 0.0d0
    is_converge = .false.

    total_dof = N*NDOF

    n_bc = 0
    do i = 1, N*NDOF
      if(is_bc(i)) n_bc = n_bc + 1
    enddo

    call monolis_allreduce_I1(total_dof, monolis_mpi_sum, monoCOM%comm)
    call monolis_allreduce_I1(n_bc, monolis_mpi_sum, monoCOM%comm)
    total_dof = total_dof - n_bc

    if(n_get_eigen > total_dof) n_get_eigen = total_dof

    call monolis_alloc_R_1d(alpha, maxiter)
    call monolis_alloc_R_1d(beta, maxiter)
    call monolis_alloc_R_1d(eigen_value, maxiter)
    call monolis_alloc_R_1d(prev, maxiter)
    call monolis_alloc_R_1d(p, NP*NDOF)
    call monolis_alloc_R_2d(q, NP*NDOF, maxiter + 1)
    call monolis_alloc_R_2d(eigen_mode, NP*NDOF, n_get_eigen)

    call lanczos_initialze(monoCOM, N, NDOF, q(:,1), is_bc)

    do iter = 1, maxiter
      call monolis_set_RHS_R(monoMAT, q(:,iter))

      monoPRM%Iarray(monolis_prm_I_is_prec_stored) = monolis_I_true
      call monolis_solve_main_R(monoPRM, monoCOM, monoMAT, monoPREC)

      do i = 1, NP*NDOF
        if(is_bc(i)) monoMAT%R%X(i) = 0.0d0
      enddo

      if(iter > 1)then
        call monolis_vec_AXPBY_R(N, NDOF, -beta(iter-1), q(:,iter-1), 1.0d0, monoMAT%R%X, p)
      else
        p = monoMAT%R%X
      endif

      call monolis_inner_product_main_R(monoCOM, N, NDOF, p, q(:,iter), alpha(iter))

      call monolis_vec_AXPBY_R(N, NDOF, -alpha(iter), q(:,iter), 1.0d0, p, p)

      call monolis_gram_schmidt_R(monoCOM, iter, N, NDOF, p, q)

      call monolis_inner_product_main_R(monoCOM, N, NDOF, p, p, beta_t)

      beta(iter) = dsqrt(beta_t)
      beta_t = 1.0d0/beta(iter)
      do i = 1, NP*NDOF
        q(i,iter+1) = p(i)*beta_t
      enddo

      call monolis_get_inverted_eigen_pair_from_tridiag(iter, n_get_eigen, &
        & alpha, beta, q, eigen_value, eigen_mode, norm)

      if(norm < ths) is_converge = .true.

      if(monolis_mpi_get_local_my_rank(monoCOM%comm) == 0)then
        write(*,"(a,i6,a,1p2e12.4)")"iter: ", iter, ", ths: ", norm
      endif

      if(is_converge .or. iter >= total_dof .or. iter == maxiter)then
        do i = 1, n_get_eigen
          val(i) = eigen_value(i)
          do j = 1, NP*NDOF
            vec(j,i) = eigen_mode(j,i)
          enddo
          call monolis_mpi_update_R(monoCOM, NDOF, vec(:,i), tmp)
        enddo
        exit
      endif
    enddo
  end subroutine monolis_eigen_inverted_standard_lanczos_R_main

  !> @ingroup eigen
  !> Lanczos 法（順反復、最大固有値、実数型、メイン関数）
  subroutine monolis_eigen_standard_lanczos_R_main( &
    & monoPRM, monoCOM, monoMAT, n_get_eigen, ths, maxiter, val, vec, is_bc)
    implicit none
    !> [in] パラメータ構造体
    type(monolis_prm), intent(in) :: monoPRM
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in] 行列構造体
    type(monolis_mat), intent(in) :: monoMAT
    !> [in,out] 取得固有値数
    integer(kint), intent(inout) :: n_get_eigen
    !> [in] 収束判定閾値
    real(kdouble), intent(in) :: ths
    !> [in] 最大反復回数
    integer(kint), intent(in) :: maxiter
    !> [out] 固有値
    real(kdouble), intent(out) :: val(:)
    !> [out] 固有ベクトル
    real(kdouble), intent(out) :: vec(:,:)
    !> [in] Dirhchlet 境界条件判定フラグ
    logical, intent(in) :: is_bc(:)
    integer(kint) :: N, NP, NDOF, total_dof, n_bc, j, k
    integer(kint) :: i, iter
    real(kdouble) :: beta_t, norm, tmp
    real(kdouble) :: tspmv, tcomm_spmv
    real(kdouble), allocatable :: p(:), q(:,:), alpha(:), beta(:), eigen_value(:), eigen_mode(:,:), prev(:)
    logical :: is_converge

    call monolis_std_debug_log_header("monolis_eigen_standard_lanczos_R_main")

    N     = monoMAT%N
    NP    = monoMAT%NP
    NDOF  = monoMAT%NDOF
    norm = 0.0d0
    is_converge = .false.

    total_dof = N*NDOF

    n_bc = 0
    do i = 1, N*NDOF
      if(is_bc(i)) n_bc = n_bc + 1
    enddo

    call monolis_allreduce_I1(total_dof, monolis_mpi_sum, monoCOM%comm)
    call monolis_allreduce_I1(n_bc, monolis_mpi_sum, monoCOM%comm)
    total_dof = total_dof - n_bc

    if(n_get_eigen > total_dof) n_get_eigen = total_dof

    call monolis_alloc_R_1d(alpha, maxiter)
    call monolis_alloc_R_1d(beta, maxiter)
    call monolis_alloc_R_1d(eigen_value, maxiter)
    call monolis_alloc_R_1d(prev, maxiter)
    call monolis_alloc_R_1d(p, NP*NDOF)
    call monolis_alloc_R_2d(q, NP*NDOF, maxiter + 1)
    call monolis_alloc_R_2d(eigen_mode, NP*NDOF, n_get_eigen)

    call lanczos_initialze(monoCOM, N, NDOF, q(:,1), is_bc)

    do iter = 1, maxiter
      call monolis_matvec_product_main_R(monoCOM, monoMAT, q(:,iter), monoMAT%R%X, tspmv, tcomm_spmv)

      call monolis_mpi_update_R(monoCOM, monoMAT%NDOF, monoMAT%R%X, tcomm_spmv)

      do i = 1, N*NDOF
        if(is_bc(i)) monoMAT%R%X(i) = 0.0d0
      enddo

      if(iter > 1)then
        call monolis_vec_AXPBY_R(N, NDOF, -beta(iter-1), q(:,iter-1), 1.0d0, monoMAT%R%X, p)
      else
        p = monoMAT%R%X
      endif

      call monolis_inner_product_main_R(monoCOM, N, NDOF, p, q(:,iter), alpha(iter))

      call monolis_vec_AXPBY_R(N, NDOF, -alpha(iter), q(:,iter), 1.0d0, p, p)

      call monolis_gram_schmidt_R(monoCOM, iter, N, NDOF, p, q)

      call monolis_inner_product_main_R(monoCOM, N, NDOF, p, p, beta_t)

      beta(iter) = dsqrt(beta_t)
      beta_t = 1.0d0/beta(iter)
      do i = 1, NP*NDOF
        q(i,iter+1) = p(i)*beta_t
      enddo

      call monolis_get_inverted_eigen_pair_from_tridiag(iter, n_get_eigen, &
        & alpha, beta, q, eigen_value, eigen_mode, norm)

      do i = 1, min(iter, n_get_eigen)
        eigen_value(i) = 1.0d0/eigen_value(i)
      enddo

      if(norm < ths) is_converge = .true.

      if(monolis_mpi_get_local_my_rank(monoCOM%comm) == 0)then
        write(*,"(a,i6,a,1p2e12.4)")"iter: ", iter, ", ths: ", norm
      endif

      if(is_converge .or. iter >= total_dof .or. iter == maxiter)then
        do i = 1, n_get_eigen
          val(i) = eigen_value(i)
          do j = 1, NP*NDOF
            vec(j,i) = eigen_mode(j,i)
          enddo
          call monolis_mpi_update_R(monoCOM, NDOF, vec(:,i), tmp)
        enddo
        exit
      endif
    enddo
  end subroutine monolis_eigen_standard_lanczos_R_main
end module mod_monolis_eigen_lanczos

