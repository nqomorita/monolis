!> Lanczos 法 util モジュール
module mod_monolis_eigen_lanczos_util
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc
  use mod_monolis_inner_product
  use mod_monolis_lapack

  implicit none

contains

  !> Lanczos 法の初期ベクトル生成
  subroutine lanczos_initialze(monoCOM, N, NDOF, q, is_bc)
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
    integer(kint) :: i, comm_size
    real(kdouble) :: norm, t1, t2

    comm_size = monolis_mpi_local_comm_size(monoCOM%comm)

    call monolis_get_rundom_number_R(N*NDOF, q, comm_size)

    do i = 1, N*NDOF
      if(is_bc(i)) q(i) = 0.0d0
    enddo

    call monolis_update_R(monoCOM, NDOF, q, t1)

    call monolis_inner_product_main_R(monoCOM, N, NDOF, q, q, norm, t1, t2)

    norm = 1.0d0/dsqrt(norm)
    do i = 1, N*NDOF
      q(i) = q(i)*norm
    enddo
  end subroutine lanczos_initialze

  !> Lanczos 法における三重対角行列の固有値分解
  subroutine monolis_get_inverted_eigen_pair_from_tridiag(iter, n_get_eigen, &
    & alpha, beta, q, eig_val, eig_mode, norm)
    implicit none
    !> 反復回数
    integer(kint), intent(in) :: iter
    !> 取得したい固有値数
    integer(kint), intent(in) :: n_get_eigen
    !> 対角成分
    real(kdouble) :: alpha(:)
    !> 副対角成分
    real(kdouble) :: beta(:)
    !> Lanczos 法から得られるユニタリ行列
    real(kdouble) :: q(:,:)
    !> 固有値
    real(kdouble) :: eig_val(:)
    !> 固有ベクトル
    real(kdouble) :: eig_mode(:,:)
    !> 固有方程式の残差
    real(kdouble) :: norm
    integer(kint) :: i
    real(kdouble) :: r
    real(kdouble), allocatable :: eig_val_tri(:)
    real(kdouble), allocatable :: eig_mode_tri(:,:)

    call monolis_alloc_R_1d(eig_val_tri, iter)
    call monolis_alloc_R_2d(eig_mode_tri, iter, iter)

    call monolis_lapack_dstev(iter, alpha, beta, eig_val_tri, eig_mode_tri)

    norm = 0.0d0
    do i = 1, min(iter, n_get_eigen)
      eig_val(i) = 1.0d0/eig_val_tri(iter - i +1)
      eig_mode(:,i) = matmul(q(:,1:iter), eig_mode_tri(1:iter,iter - i + 1))

      r = sqrt(eig_mode_tri(iter,iter - i + 1)**2)*beta(iter)
      if(norm < r)then
         norm = r
      endif
    enddo
  end subroutine monolis_get_inverted_eigen_pair_from_tridiag

end module mod_monolis_eigen_lanczos_util
