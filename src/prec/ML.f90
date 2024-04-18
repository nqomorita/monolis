!> AMG 前処理関連モジュール
module mod_monolis_precond_ML
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc

  implicit none

  private

  public :: monolis_precond_ml_setup_R
  public :: monolis_precond_ml_apply_R
  public :: monolis_precond_ml_clear_R

  type(monolis_mat), pointer, save :: monoMAT_save
  type(monolis_com), pointer, save :: monoCOM_save

contains

  !> @ingroup prec
  !> 前処理生成：ML 前処理（実数型）
  subroutine monolis_precond_ml_setup_R(monoPRM, monoCOM, monoMAT)
    implicit none
    !> [in,out] パラメータ構造体
    type(monolis_prm), intent(inout) :: monoPRM
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in), target :: monoCOM
    !> [in] 行列構造体
    type(monolis_mat), intent(in), target :: monoMAT
    integer(kint) :: sym, ndof, ierr

    monoCOM_save => monoCOM
    monoMAT_save => monoMAT

    call monolis_std_debug_log_header("monolis_precond_ml_setup_R")
    sym = 0
    ndof = monoMAT%ndof
    call monolis_precond_ml_setup(sym, ndof, ierr);
  end subroutine monolis_precond_ml_setup_R

  !> @ingroup prec
  !> 前処理適用：SOR 前処理（実数型）
  subroutine monolis_precond_ml_apply_R(monoPRM, monoCOM, monoMAT, X, Y)
    implicit none
    !> [in,out] パラメータ構造体
    type(monolis_prm), intent(inout) :: monoPRM
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in] 行列構造体
    type(monolis_mat), intent(in) :: monoMAT
    real(kdouble) :: X(:), Y(:)
    integer(kint) :: ierr

    call monolis_std_debug_log_header("monolis_precond_ml_apply_R")
    call monolis_precond_ml_apply(X, ierr)
    Y = X
  end subroutine monolis_precond_ml_apply_R

  !> @ingroup prec
  !> 前処理初期化：SOR 前処理（実数型）
  subroutine monolis_precond_ml_clear_R(monoPRM, monoCOM, monoMAT)
    implicit none
    !> [in,out] パラメータ構造体
    type(monolis_prm), intent(inout) :: monoPRM
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in] 行列構造体
    type(monolis_mat), intent(in) :: monoMAT
    integer(kint) :: ierr

    call monolis_std_debug_log_header("monolis_precond_ml_clear_R")

    call monolis_precond_ml_clear(ierr)
  end subroutine monolis_precond_ml_clear_R

  !> wrapper section
  subroutine monolis_ML_get_nlocal(nlocal, nlocal_allcolumns, ierr)
    use mod_monolis_utils
    implicit none
    integer(kint), intent(out) :: nlocal
    integer(kint), intent(out) :: nlocal_allcolumns
    integer(kint), intent(out) :: ierr
    !type(monolisST_matrix), pointer :: hecMAT
    
    !nlocal = hecMAT%N * hecMAT%NDOF
    !nlocal_allcolumns = hecMAT%NP * hecMAT%NDOF
    !ierr = 0
  end subroutine monolis_ML_get_nlocal

  subroutine monolis_ML_get_opt(opt, ierr)
    use mod_monolis_utils
    implicit none
    integer(kint), intent(out) :: opt(*)
    integer(kint), intent(out) :: ierr
    !type(monolisST_matrix), pointer :: hecMAT
    integer(kint) :: iopt(10)

    !call monolis_mat_get_solver_opt(hecMAT, iopt)
    !opt(1:10) = iopt(1:10)
    !ierr = 0
  end subroutine monolis_ML_get_opt

  subroutine monolis_ML_set_opt(opt, ierr)
    use mod_monolis_utils
    implicit none
    integer(kint), intent(in) :: opt(*)
    integer(kint), intent(out) :: ierr
    integer(kint) :: iopt(10)
    
    !iopt(1:10) = opt(1:10)
    !call monolis_mat_set_solver_opt(hecMAT, iopt)
    !ierr = 0
  end subroutine monolis_ML_set_opt

  subroutine monolis_ML_getrow_nn(n_requested_rows, requested_rows, &
      allocated_space, cols, values, row_lengths, ierr)
    use mod_monolis_utils
    implicit none
    integer(kint), intent(in) :: n_requested_rows
    integer(kint), intent(in) :: requested_rows(n_requested_rows)
    integer(kint), intent(in) :: allocated_space
    integer(kint), intent(out) :: cols(allocated_space)
    real(kdouble), intent(out) :: values(allocated_space)
    integer(kint), intent(out) :: row_lengths(n_requested_rows)
    integer(kint), intent(out) :: ierr
    integer(kint) :: m, i, row, inod, idof, nl, nd, nu, js, je, j, jj, jdof, start, ndof

  !  !call monolis_mat_id_get(id, hecMAT, hecMESH)
  !  ndof = hecMAT%NDOF
  !  m = 1
  !  do i = 1, n_requested_rows
  !    row = requested_rows(i) + 1 ! '+1' for Fortran-numbering
  !    inod = (row-1)/ndof + 1
  !    idof = row - (inod-1)*ndof
  !    nl = (hecMAT%indexL(inod) - hecMAT%indexL(inod-1)) * ndof
  !    nd = ndof
  !    nu = (hecMAT%indexU(inod) - hecMAT%indexU(inod-1)) * ndof
  !    if (allocated_space < m + nl + nd + nu) return
  !    start = m
  !    js = hecMAT%indexL(inod-1)+1
  !    je = hecMAT%indexL(inod)
  !    do j = js, je
  !      jj = hecMAT%itemL(j)
  !      do jdof = 1, ndof
  !        cols(m) = (jj-1)*ndof + jdof - 1 ! '-1' for C-numbering
  !        values(m) = hecMAT%AL((j-1)*ndof*ndof + (idof-1)*ndof + jdof)
  !        m = m+1
  !      enddo
  !    enddo
  !    do jdof = 1, ndof
  !      cols(m) = (inod-1)*ndof + jdof - 1 ! '-1' for C-numbering
  !      values(m) = hecMAT%D((inod-1)*ndof*ndof + (idof-1)*ndof + jdof)
  !      m = m+1
  !    enddo
  !    js = hecMAT%indexU(inod-1)+1
  !    je = hecMAT%indexU(inod)
  !    do j = js, je
  !      jj = hecMAT%itemU(j)
  !      do jdof = 1, ndof
  !        cols(m) = (jj-1)*ndof + jdof - 1 ! '-1' for C-numbering
  !        values(m) = hecMAT%AU((j-1)*ndof*ndof + (idof-1)*ndof + jdof)
  !        m = m+1
  !      enddo
  !    enddo
  !    row_lengths(i) = m - start
  !  enddo
  !  ierr = 1
  end subroutine monolis_ML_getrow_nn

  subroutine monolis_ML_matvec_nn(in_length, p, out_length, ap, ierr)
    use mod_monolis_utils
    implicit none
    integer(kint), intent(in) :: in_length
    real(kdouble), intent(in) :: p(in_length)
    integer(kint), intent(in) :: out_length
    real(kdouble), intent(out) :: ap(out_length)
    integer(kint), intent(out) :: ierr
    real(kdouble), allocatable :: w(:)
    integer(kint) :: i

    !call monolis_mat_id_get(id, hecMAT, hecMESH)
    !allocate(w(hecMAT%NP*hecMAT%NDOF))
    !do i = 1, hecMAT%N*hecMAT%NDOF
    !  w(i) = p(i)
    !enddo
    !call monolis_matvec(hecMESH, hecMAT, w, ap)
    !deallocate(w)
    !ierr = 0
  end subroutine monolis_ML_matvec_nn

  subroutine monolis_ML_comm_nn(x, ierr)
    use mod_monolis_utils
    implicit none
    real(kdouble), intent(inout) :: x(*)
    integer(kint), intent(out) :: ierr

    !type(monolisST_matrix), pointer :: hecMAT
    !call monolis_mat_id_get(id, hecMAT, hecMESH)
    !call monolis_update_R (hecMESH, x, hecMAT%NP, hecMAT%NDOF)
    ierr = 0
  end subroutine monolis_ML_comm_nn

end module mod_monolis_precond_ML