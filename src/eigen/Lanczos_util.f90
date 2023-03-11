!> Lanczos 法 util モジュール
module mod_monolis_eigen_lanczos_util
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc

  implicit none

contains

  !> Lanczos 法の初期ベクトル生成
  subroutine lanczos_initialze(monoCOM, N, NDOF, q, is_bc, beta)
    implicit none
    !> 通信テーブル構造体
    type(monolis_com) :: monoCOM
    !> 計算点数
    integer(kint) :: N
    !> 計算点上の自由度
    integer(kint) :: NDOF
    !> 入力ベクトル
    real(kdouble) :: q(:)
    !> Dirhchlet 境界条件判定フラグ
    logical :: is_bc(:)
    !> パラメータ beta
    real(kdouble) :: beta
    integer(kint) :: i, comm_size
    real(kdouble) :: norm, t1, t2

    comm_size = monolis_mpi_local_comm_size(monoCOM%comm)

    call get_rundom_number(N*NDOF, q, comm_size)

    do i = 1, N*NDOF
      if(is_bc(i)) q(i) = 0.0d0
    enddo

    call monolis_update_R(monoCOM, NDOF, q, t1)

    call monolis_inner_product_R(monoCOM, N, NDOF, q, q, norm, t1, t2)

    beta = dsqrt(norm)
    norm = 1.0d0/dsqrt(norm)
    do i = 1, N*NDOF
      q(i) = q(i)*norm
    enddo
  end subroutine lanczos_initialze

  subroutine monolis_get_eigen_pair_from_tridiag(iter, n_get_eigen, &
    & alpha_t, beta_t, q, e_value, e_mode, norm)
    implicit none
    integer(kint), intent(in) :: iter
    integer(kint) :: i, n, ldz, info, n_get_eigen
    real(kdouble) :: alpha_t(:), beta_t(:), q(:,0:), e_value(:), e_mode(:,:), norm
    integer(kint), allocatable :: isuppz(:), idum(:)
    real(kdouble), allocatable :: alpha(:), beta(:), rdum(:), e_mode_t(:,:)

    !> DSTEVR
    allocate(alpha(iter), source = 0.0d0)
    allocate(beta (max(1,iter-1)), source = 0.0d0)
    allocate(rdum(20*iter), source = 0.0d0)
    allocate(e_mode_t(iter,iter), source = 0.0d0)

    alpha = alpha_t(1:iter)
    beta = beta_t(2:max(1,iter-1)+1)
    if(iter == 1) beta = 0.0d0

    n = iter
    ldz = iter

    call dstev("V", n, alpha, beta, e_mode_t, ldz, rdum, info)

    norm = 0.0d0
    do i = 1, min(iter, n_get_eigen)
      e_value(i) = 1.0d0/alpha(iter - i +1)
      e_mode(:,i) = matmul(q(:,1:iter), e_mode_t(1:iter,iter - i + 1))
      if(norm < sqrt(e_mode_t(iter,iter - i + 1)**2)*beta_t(iter+1))then
        norm = sqrt(e_mode_t(iter,iter - i + 1)**2)*beta_t(iter+1)
      endif
    enddo

    if(info /= 0) stop "monolis_get_eigen_pair_from_tridiag"

    deallocate(alpha)
    deallocate(beta)
    deallocate(e_mode_t)
    deallocate(rdum)
  end subroutine monolis_get_eigen_pair_from_tridiag

end module mod_monolis_eigen_lanczos_util
