!> Lanczos 法 util モジュール
module mod_monolis_eigen_lanczos_util
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc
  use mod_monolis_inner_product
  use mod_monolis_lapack

  implicit none

contains

  !> @ingroup eigen
  !> Lanczos 法の初期ベクトル生成
  subroutine lanczos_initialze(monoCOM, N, NDOF, q, is_bc)
    implicit none
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in] 計算点数
    integer(kint), intent(in) :: N
    !> [in] 計算点上の自由度
    integer(kint), intent(in) :: NDOF
    !> [out] 入力ベクトル
    real(kdouble), intent(out) :: q(:)
    !> [in] Dirhchlet 境界条件判定フラグ
    logical, intent(in) :: is_bc(:)
    integer(kint) :: i
    real(kdouble) :: norm, t1, t2
    integer(kint), allocatable :: vtxdist(:)

    call monolis_com_n_vertex_list(N*NDOF, monoCOM%comm, vtxdist)

    call monolis_get_rundom_number_R(N*NDOF, q, vtxdist(monoCOM%my_rank + 1))

    do i = 1, N*NDOF
      if(is_bc(i)) q(i) = 0.0d0
    enddo

    call monolis_mpi_update_R(monoCOM, NDOF, q, t1)

    call monolis_inner_product_main_R(monoCOM, N, NDOF, q, q, norm, t1, t2)

    norm = 1.0d0/dsqrt(norm)
    do i = 1, N*NDOF
      q(i) = q(i)*norm
    enddo

    call monolis_mpi_update_R(monoCOM, NDOF, q, t1)
  end subroutine lanczos_initialze

  !> @ingroup eigen
  !> Lanczos 法における三重対角行列の固有値分解
  subroutine monolis_get_inverted_eigen_pair_from_tridiag(iter, n_get_eigen, &
    & alpha, beta, q, eig_val, eig_mode, norm)
    implicit none
    !> [in] 反復回数
    integer(kint), intent(in) :: iter
    !> [in] 取得したい固有値数
    integer(kint), intent(in) :: n_get_eigen
    !> [in] 対角成分
    real(kdouble), intent(in) :: alpha(:)
    !> [in] 副対角成分
    real(kdouble), intent(in) :: beta(:)
    !> [in] Lanczos 法から得られるユニタリ行列
    real(kdouble), intent(in) :: q(:,:)
    !> [out] 固有値
    real(kdouble), intent(out) :: eig_val(:)
    !> [out] 固有ベクトル
    real(kdouble), intent(out) :: eig_mode(:,:)
    !> [out] 固有方程式の残差
    real(kdouble), intent(out) :: norm
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
