module mod_monolis_matmat
  use mod_monolis_prm
  use mod_monolis_mat
  implicit none

contains

  !> @ingroup linalg
  !> ベクトル内積
  subroutine monolis_matmat_local(monoCOM, A, B, C, tspmv, tcomm)
    implicit none
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in] 行列構造体
    type(monolis_mat), intent(in) :: A
    !> [in] 行列構造体
    type(monolis_mat), intent(in) :: B
    !> [in] 行列構造体
    type(monolis_mat), intent(inout) :: C
    !> [in,out] 計算時間
    real(kdouble), intent(inout) :: tspmv
    !> [in,out] 通信時間
    real(kdouble), intent(inout) :: tcomm
    real(kdouble) :: t1, t2

#ifdef DEBUG
    call monolis_std_debug_log_header("monolis_matmat_local")
#endif
    t1 = monolis_get_time()

    !call monolis_update_pre_R(monoCOM, monoMAT%NDOF, X, tcomm)

    call monolis_matmat_nn(monoCOM, A, B, C, monoMAT%NDOF)

    t2 = monolis_get_time()
    tspmv = tspmv + t2 - t1
  end subroutine monolis_matmat_local

  !> @ingroup linalg
  !> ベクトル内積
  subroutine monolis_matmat_nn(monoCOM, A, B, C, NDOF)
    implicit none
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in] 行列構造体
    type(monolis_mat), intent(in) :: A
    !> [in] 行列構造体
    type(monolis_mat), intent(in) :: B
    !> [in] 行列構造体
    type(monolis_mat), intent(inout) :: C
    !> [in] 計算点が持つ自由度
    integer(kint), intent(in) :: NDOF

    call monolis_matmat_symbolic(monoCOM, A, B, C, NDOF)
    call monolis_matmat_value(monoCOM, A, B, C, NDOF)
  end subroutine monolis_matmat_nn

  subroutine monolis_matmat_symbolic(monoCOM, A, B, C, NDOF)
    implicit none
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in] 行列構造体
    type(monolis_mat), intent(in) :: A
    !> [in] 行列構造体
    type(monolis_mat), intent(in) :: B
    !> [in] 行列構造体
    type(monolis_mat), intent(inout) :: C
    !> [in] 計算点が持つ自由度
    integer(kint), intent(in) :: NDOF

  end subroutine monolis_matmat_symbolic

  subroutine monolis_matmat_value(monoCOM, A, B, C, NDOF)
    implicit none
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in] 行列構造体
    type(monolis_mat), intent(in) :: A
    !> [in] 行列構造体
    type(monolis_mat), intent(in) :: B
    !> [in] 行列構造体
    type(monolis_mat), intent(inout) :: C
    !> [in] 計算点が持つ自由度
    integer(kint), intent(in) :: NDOF

  end subroutine monolis_matmat_value
end module mod_monolis_matmat