!> 多重フロント法 Cholesky 数値分解
!>
!> 解析フェーズ（[[mod_monolis_fact_analysis]]）の結果を再利用し、
!> 対称正定値行列 A = L L^T の下三角因子 L のみを計算する。
!> 上三角ブロック（upper_update）は対称性により保持しない。
module mod_monolis_fact_factorize_cholesky
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc
  use mod_monolis_lapack
  use mod_monolis_fact_factorize, only: front_block_size

  implicit none

  private
  public :: monolis_fact_factorize_cholesky

contains

  !> @ingroup fact
  !> Cholesky 数値分解フェーズ
  subroutine monolis_fact_factorize_cholesky(monoMAT, lu)
    implicit none
    !> [in] 行列構造体（値配列を参照）
    type(monolis_mat),    intent(in)    :: monoMAT
    !> [in,out] 分解構造体（解析フェーズ済み）
    type(monolis_mat_lu), intent(inout) :: lu

    integer(kint) :: nfronts, order_pos, front, child
    integer(kint) :: fs, npiv, nupd, entry_pos, ierr_front

    call monolis_std_debug_log_header("monolis_fact_factorize_cholesky")

    if (.not. lu%analyzed) then
      call monolis_std_error_string("monolis_fact_factorize_cholesky: analysis not done")
      call monolis_std_error_stop()
    end if

    nfronts = lu%nfronts
    if (allocated(lu%factors)) deallocate(lu%factors)
    allocate(lu%factors(max(1, nfronts)))

    do order_pos = 1, nfronts
      front = lu%front_postorder(order_pos)
      fs   = lu%front_size(front)
      npiv = lu%front_pivot_size(front)
      nupd = lu%front_update_size(front)

      call allocate_front_storage_chol(lu%factors(front), fs, npiv, nupd)

      child = lu%front_first_child(front)
      do while (child /= 0)
        call assemble_child_contribution_chol(lu, child, front, lu%factors(front))
        child = lu%front_next_sibling(child)
      end do

      do entry_pos = lu%orig_ptr(front), lu%orig_ptr(front + 1) - 1
        call add_front_entry_chol(lu%factors(front), &
            lu%orig_row_pos(entry_pos), lu%orig_col_pos(entry_pos), &
            monoMAT%R%A(lu%orig_entry(entry_pos)))
      end do

      call factor_one_front_cholesky(lu%factors(front), front, ierr_front)
      if (ierr_front /= 0) then
        call monolis_std_error_string("monolis_fact_factorize_cholesky: matrix is not positive definite")
        call monolis_std_error_stop()
      end if
    end do

    !> contribution は親に展開済みなので解放
    do front = 1, nfronts
      if (allocated(lu%factors(front)%contribution)) then
        deallocate(lu%factors(front)%contribution)
      end if
    end do

    lu%factorized = .true.
  end subroutine monolis_fact_factorize_cholesky

  !> 子フロントの contribution を親フロントに加算
  subroutine assemble_child_contribution_chol(lu, child, parent, parent_data)
    implicit none
    type(monolis_mat_lu),     intent(inout) :: lu
    integer(kint),            intent(in)    :: child, parent
    type(monolis_mat_frontal), intent(inout) :: parent_data

    integer(kint) :: child_nupd, j, parent_j, run, pos_base, run_first, run_last

    child_nupd = lu%front_update_size(child)
    if (child_nupd <= 0) return
    if (lu%front_parent(child) /= parent) then
      call monolis_std_error_string("assemble_child_contribution_chol: parent mismatch")
      call monolis_std_error_stop()
    end if
    if (.not. allocated(lu%factors(child)%contribution)) then
      call monolis_std_error_string("assemble_child_contribution_chol: missing contribution")
      call monolis_std_error_stop()
    end if

    pos_base  = lu%contrib_pos_ptr(child) - 1
    run_first = lu%contrib_run_ptr(child)
    run_last  = lu%contrib_run_ptr(child + 1) - 1
    do j = 1, child_nupd
      parent_j = lu%contrib_parent_pos(pos_base + j)
      do run = run_first, run_last
        call add_front_column_run_chol(parent_data, &
            lu%contrib_run_parent_first(run), parent_j, &
            lu%contrib_run_len(run), &
            lu%factors(child)%contribution(lu%contrib_run_first(run), j))
      end do
    end do

    deallocate(lu%factors(child)%contribution)
  end subroutine assemble_child_contribution_chol

  !> フロントの数値格納領域を確保しゼロ初期化（下三角のみ、upper_update なし）
  subroutine allocate_front_storage_chol(front_data, fs, npiv, nupd)
    implicit none
    type(monolis_mat_frontal), intent(inout) :: front_data
    integer(kint), intent(in) :: fs, npiv, nupd

    front_data%front_size  = fs
    front_data%pivot_size  = npiv
    front_data%update_size = nupd

    allocate(front_data%factor(max(1, fs), max(1, npiv)))
    call dlaset('A', max(1, fs), max(1, npiv), 0.0d0, 0.0d0, &
        front_data%factor, max(1, fs))

    if (nupd > 0) then
      allocate(front_data%contribution(nupd, nupd))
      call dlaset('A', nupd, nupd, 0.0d0, 0.0d0, &
          front_data%contribution, max(1, nupd))
    end if
  end subroutine allocate_front_storage_chol

  !> 単一の値をフロントの該当ブロックに加算（上三角ブロックは対称性により無視）
  subroutine add_front_entry_chol(front_data, row_pos, col_pos, value)
    implicit none
    type(monolis_mat_frontal), intent(inout) :: front_data
    integer(kint),  intent(in) :: row_pos, col_pos
    real(kdouble),  intent(in) :: value

    integer(kint) :: npiv

    npiv = front_data%pivot_size
    if (col_pos <= npiv) then
      front_data%factor(row_pos, col_pos) = &
          front_data%factor(row_pos, col_pos) + value
    else if (row_pos > npiv) then
      front_data%contribution(row_pos - npiv, col_pos - npiv) = &
          front_data%contribution(row_pos - npiv, col_pos - npiv) + value
    end if
  end subroutine add_front_entry_chol

  !> 連続した行の値（source(1:len)）を、親フロントの該当列に加算
  !> （上三角ブロックへの寄与は対称性により無視）
  subroutine add_front_column_run_chol(front_data, row_first, col_pos, len, source)
    implicit none
    type(monolis_mat_frontal), intent(inout) :: front_data
    integer(kint),  intent(in) :: row_first, col_pos, len
    real(kdouble),  intent(in) :: source(*)

    integer(kint) :: npiv, pivot_len, update_len, update_row, update_col

    if (len <= 0) return
    npiv = front_data%pivot_size
    if (col_pos <= npiv) then
      call daxpy(len, 1.0d0, source, 1, &
          front_data%factor(row_first, col_pos), 1)
      return
    end if

    update_col = col_pos - npiv
    pivot_len = 0
    if (row_first <= npiv) pivot_len = min(len, npiv - row_first + 1)

    update_len = len - pivot_len
    if (update_len > 0) then
      update_row = row_first + pivot_len - npiv
      call daxpy(update_len, 1.0d0, source(1 + pivot_len), 1, &
          front_data%contribution(update_row, update_col), 1)
    end if
  end subroutine add_front_column_run_chol

  !> 単一フロントの Cholesky 分解（パネル化、対角ピボット）
  subroutine factor_one_front_cholesky(front_data, front, ierr)
    implicit none
    type(monolis_mat_frontal), intent(inout) :: front_data
    integer(kint), intent(in)  :: front
    integer(kint), intent(out) :: ierr

    integer(kint) :: fs, npiv, nupd, ldf, k, j
    integer(kint) :: panel_end, block_cols, info

    ierr = 0
    fs   = front_data%front_size
    npiv = front_data%pivot_size
    nupd = front_data%update_size
    ldf  = max(1, fs)
    if (npiv <= 0) return

    k = 1
    do while (k <= npiv)
      panel_end  = min(npiv, k + front_block_size - 1)
      block_cols = panel_end - k + 1

      !> パネル対角ブロックの Cholesky（LAPACK）
      call dpotrf('L', block_cols, front_data%factor(k, k), ldf, info)
      if (info /= 0) then
        ierr = front
        return
      end if
      !> 微小ピボット検出（元実装の pivot <= 100*eps 相当）
      do j = k, panel_end
        if (front_data%factor(j, j) <= sqrt(100.0d0 * epsilon(1.0d0))) then
          ierr = front
          return
        end if
      end do

      !> パネル下部の L: A21 L11^{-T}（対角ブロックの下の行を一括更新）
      if (fs > panel_end) then
        call dtrsm('R', 'L', 'T', 'N', fs - panel_end, block_cols, 1.0d0, &
            front_data%factor(k, k), ldf, front_data%factor(panel_end + 1, k), ldf)
      end if

      !> 対称性を利用し、正方ブロックは dsyrk（下三角のみ）、
      !> それ以外の長方形部分は dgemm で更新する
      if (npiv > panel_end) then
        call dsyrk('L', 'N', npiv - panel_end, block_cols, &
            -1.0d0, front_data%factor(panel_end + 1, k), ldf, 1.0d0, &
            front_data%factor(panel_end + 1, panel_end + 1), ldf)
        if (fs > npiv) then
          call dgemm('N', 'T', fs - npiv, npiv - panel_end, block_cols, &
              -1.0d0, front_data%factor(npiv + 1, k), ldf, &
              front_data%factor(panel_end + 1, k), ldf, 1.0d0, &
              front_data%factor(npiv + 1, panel_end + 1), ldf)
        end if
      end if
      if (nupd > 0) then
        call dsyrk('L', 'N', nupd, block_cols, -1.0d0, &
            front_data%factor(npiv + 1, k), ldf, &
            1.0d0, front_data%contribution, max(1, nupd))
      end if

      k = panel_end + 1
    end do

    !> contribution は親フロントへの assembly が全列 run を読むため、
    !> dsyrk で未計算の狭義上三角を下三角の鏡像で補完する
    if (nupd > 0) then
      do k = 2, nupd
        do j = 1, k - 1
          front_data%contribution(j, k) = front_data%contribution(k, j)
        end do
      end do
    end if
  end subroutine factor_one_front_cholesky

end module mod_monolis_fact_factorize_cholesky
