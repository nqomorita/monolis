!> 疎行列操作関数群（メイン関数）
!# monolis_set_scalar_to_sparse_matrix_main_R(index, item, A, ndof, ci, cj, csub_i, csub_j, val)
!# monolis_set_scalar_to_sparse_matrix_main_C(index, item, A, ndof, ci, cj, csub_i, csub_j, val)
!# monolis_set_block_to_sparse_matrix_main_R(index, item, A, ndof, ci, cj, val)
!# monolis_set_block_to_sparse_matrix_main_C(index, item, A, ndof, ci, cj, val)
!# monolis_get_scalar_from_sparse_matrix_main_R(index, item, A, ndof, ci, cj, csub_i, csub_j, val, is_find)
!# monolis_get_scalar_from_sparse_matrix_main_C(index, item, A, ndof, ci, cj, csub_i, csub_j, val, is_find)
!# monolis_add_scalar_to_sparse_matrix_main_R(index, item, A, ndof, ci, cj, csub_i, csub_j, val)
!# monolis_add_scalar_to_sparse_matrix_main_C(index, item, A, ndof, ci, cj, csub_i, csub_j, val)
!# monolis_add_matrix_to_sparse_matrix_main_R(index, item, A, n1, n2, ndof, e1, e2, val)
!# monolis_add_matrix_to_sparse_matrix_main_C(index, item, A, n1, n2, ndof, e1, e2, val)
!# monolis_set_Dirichlet_bc_main_R(index, item, A, B, indexR, itemR, permA, &
!# monolis_set_Dirichlet_bc_main_C(index, item, A, B, indexR, itemR, permA, &
!# monolis_stop_by_matrix_assemble(ci, cj)
!# monolis_stop_by_submatrix_access(ndof, sub_dof)
!# monolis_stop_by_set_DBC(node_id)
!# monolis_stop_by_set_zero_diag_component(node_id, ndof)
!# monolis_get_max_matrix_component_main_R(monoMAT, monoCOM, max_val)
!# monolis_check_diagonal_zero_component_main_R(monoPRM, monoMAT)
module mod_monolis_spmat_handler_util
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc
#ifdef _OPENACC
  use openacc
#endif

  implicit none

  !> デバイス組立の遅延エラー検出:エラー数（GPU 要素ループからの組立はデバイス上で停止できないため、
  !> パターン外参照等をカウンタに記録し、solve 入口・境界条件・形式変換の同期点で検査する）
  integer(kint) :: monolis_assemble_error_count = 0
  !> デバイス組立の遅延エラー検出:エラー発生位置（行番号）
  integer(kint) :: monolis_assemble_error_i = 0
  !> デバイス組立の遅延エラー検出:エラー発生位置（列番号）
  integer(kint) :: monolis_assemble_error_j = 0
!$acc declare create(monolis_assemble_error_count, monolis_assemble_error_i, monolis_assemble_error_j)

contains

  !> @ingroup dev_matrix
  !> 組立エラーの記録（デバイス呼び出し可能。エラーカウンタと発生位置を更新する）
  subroutine monolis_record_matrix_assemble_error(ci, cj)
    implicit none
!$acc routine seq
    !> [in] 行番号
    integer(kint), intent(in) :: ci
    !> [in] 列番号
    integer(kint), intent(in) :: cj

!$acc atomic update
!$omp atomic update
    monolis_assemble_error_count = monolis_assemble_error_count + 1
    !# 発生位置は競合時にいずれかの失敗位置が残れば十分(報告用)
    monolis_assemble_error_i = ci
    monolis_assemble_error_j = cj
  end subroutine monolis_record_matrix_assemble_error

  !> @ingroup matrix
  !> デバイス組立の遅延エラー検査
  !> @details GPU 要素ループからの組立（\ref monolis_add_matrix_to_sparse_matrix_main_R 等の
  !> デバイス直接呼び出し）で記録されたエラーを検査し、エラーがあれば停止する。
  !> solve 入口・境界条件・形式変換の同期点から自動的に呼ばれるため、通常ユーザが呼ぶ必要はない。
  subroutine monolis_check_matrix_assemble_error()
    implicit none
#ifdef _OPENACC
    !$acc wait(1)
    !$acc update self(monolis_assemble_error_count)
    if(monolis_assemble_error_count > 0)then
      !$acc update self(monolis_assemble_error_i, monolis_assemble_error_j)
      call monolis_stop_by_matrix_assemble(monolis_assemble_error_i, monolis_assemble_error_j)
    endif
#endif
  end subroutine monolis_check_matrix_assemble_error

  !> @ingroup dev_matrix
  !> デバイス組立エラーカウンタの初期化（常駐開始時に呼ばれる）
  subroutine monolis_init_matrix_assemble_error()
    implicit none
    monolis_assemble_error_count = 0
    monolis_assemble_error_i = 0
    monolis_assemble_error_j = 0
#ifdef _OPENACC
    !$acc update device(monolis_assemble_error_count, monolis_assemble_error_i, monolis_assemble_error_j)
#endif
  end subroutine monolis_init_matrix_assemble_error

  !# setter
  !> @ingroup matrix
  !> スカラ値を疎行列に設定（メイン関数、実数型）
  !> @details GPU 要素ループ（OpenACC compute region）内から直接呼び出し可能（acc routine seq）。
  !> デバイス呼び出し時は常駐配列（index, item, A, n_dof_index, n_dof_index2）を present で渡すこと。
  !> パターン外参照等のエラーはカウンタに記録され、\ref monolis_check_matrix_assemble_error で検査される。
  subroutine monolis_set_scalar_to_sparse_matrix_main_R(index, item, A, n_dof_index, n_dof_index2, &
    ci, cj, csub_i, csub_j, val)
    implicit none
!$acc routine seq
    !> [in] index 配列
    integer(kint), intent(in) :: index(:)
    !> [in] item 配列
    integer(kint), intent(in) :: item(:)
    !> [in,out] 係数行列
    real(kdouble), intent(inout) :: A(:)
    !> [in] ブロック自由度リスト
    integer(kint), intent(in) :: n_dof_index(:)
    !> [in] ブロック自由度リスト
    integer(kint), intent(in) :: n_dof_index2(:)
    !> [in] 行番号
    integer(kint), intent(in) :: ci
    !> [in] 列番号
    integer(kint), intent(in) :: cj
    !> [in] ブロック中の行番号
    integer(kint), intent(in) :: csub_i
    !> [in] ブロック中の列番号
    integer(kint), intent(in) :: csub_j
    !> [in] 設定値
    real(kdouble), intent(in) :: val
    integer(kint) :: j, im, jS, jE, ndof

    ndof = n_dof_index(ci + 1) - n_dof_index(ci)
    if(ndof < csub_i)then
      call monolis_record_matrix_assemble_error(ci, cj)
      return
    endif
    ndof = n_dof_index(cj + 1) - n_dof_index(cj)
    if(ndof < csub_j)then
      call monolis_record_matrix_assemble_error(ci, cj)
      return
    endif

    jS = index(ci) + 1
    jE = index(ci + 1)
    !# 行内の item は昇順ソート済みのため二分探索
    call monolis_bsearch_I(item, jS, jE, cj, j)
    if(j == -1)then
      call monolis_record_matrix_assemble_error(ci, cj)
      return
    endif

    im = n_dof_index2(j) + ndof*(csub_i-1) + csub_j
    A(im) = val
  end subroutine monolis_set_scalar_to_sparse_matrix_main_R

  !> @ingroup dev_matrix
  !> スカラ値を疎行列に設定（メイン関数、複素数型）
  subroutine monolis_set_scalar_to_sparse_matrix_main_C(index, item, A, n_dof_index, n_dof_index2, &
    ci, cj, csub_i, csub_j, val)
    implicit none
    !> [in] index 配列
    integer(kint), intent(in) :: index(:)
    !> [in] item 配列
    integer(kint), intent(in) :: item(:)
    !> [in,out] 係数行列
    complex(kdouble), intent(inout) :: A(:)
    !> [in] ブロック自由度リスト
    integer(kint), intent(in) :: n_dof_index(:)
    !> [in] ブロック自由度リスト
    integer(kint), intent(in) :: n_dof_index2(:)
    !> [in] 行番号
    integer(kint), intent(in) :: ci
    !> [in] 列番号
    integer(kint), intent(in) :: cj
    !> [in] ブロック中の行番号
    integer(kint), intent(in) :: csub_i
    !> [in] ブロック中の列番号
    integer(kint), intent(in) :: csub_j
    !> [in] 設定値
    complex(kdouble), intent(in) :: val
    integer(kint) :: j, jn, im, jS, jE, ndof

    ndof = n_dof_index(ci + 1) - n_dof_index(ci)
    if(ndof < csub_i) call monolis_stop_by_submatrix_access(ndof, csub_i)
    ndof = n_dof_index(cj + 1) - n_dof_index(cj)
    if(ndof < csub_j) call monolis_stop_by_submatrix_access(ndof, csub_j)

    jS = index(ci) + 1
    jE = index(ci + 1)
    !# 行内の item は昇順ソート済みのため二分探索
    call monolis_bsearch_I(item, jS, jE, cj, j)
    if(j /= -1)then
      ndof = n_dof_index(cj + 1) - n_dof_index(cj)
      im = n_dof_index2(j) + ndof*(csub_i-1) + csub_j
      A(im) = val
      return
    endif

    call monolis_stop_by_matrix_assemble(ci, cj)
  end subroutine monolis_set_scalar_to_sparse_matrix_main_C

  !> @ingroup matrix
  !> 小行列を疎行列に設定（メイン関数、実数型）
  !> @details GPU 要素ループ（OpenACC compute region）内から直接呼び出し可能（acc routine seq）。
  !> デバイス呼び出し時は常駐配列（index, item, A, n_dof_index, n_dof_index2）を present で渡すこと。
  !> パターン外参照等のエラーはカウンタに記録され、\ref monolis_check_matrix_assemble_error で検査される。
  subroutine monolis_set_block_to_sparse_matrix_main_R(index, item, A, n_dof_index, n_dof_index2, &
    ci, cj, val)
    implicit none
!$acc routine seq
    !> [in] index 配列
    integer(kint), intent(in) :: index(:)
    !> [in] item 配列
    integer(kint), intent(in) :: item(:)
    !> [in,out] 係数行列
    real(kdouble), intent(inout) :: A(:)
    !> [in] ブロック自由度リスト
    integer(kint), intent(in) :: n_dof_index(:)
    !> [in] ブロック自由度リスト
    integer(kint), intent(in) :: n_dof_index2(:)
    !> [in] 行番号
    integer(kint), intent(in) :: ci
    !> [in] 列番号
    integer(kint), intent(in) :: cj
    !> [in] 設定値
    real(kdouble), intent(in) :: val(:,:)
    integer(kint) :: j, im, jS, jE, i1, i2, n1, n2

    jS = index(ci) + 1
    jE = index(ci + 1)
    !# 行内の item は昇順ソート済みのため二分探索
    call monolis_bsearch_I(item, jS, jE, cj, j)
    if(j == -1)then
      call monolis_record_matrix_assemble_error(ci, cj)
      return
    endif

    n1 = n_dof_index(ci + 1) - n_dof_index(ci)
    n2 = n_dof_index(cj + 1) - n_dof_index(cj)
    im = n_dof_index2(j)
    do i1 = 1, n1
    do i2 = 1, n2
      A(im + n2*(i1-1) + i2) = val(i1, i2)
    enddo
    enddo
  end subroutine monolis_set_block_to_sparse_matrix_main_R

  !> @ingroup dev_matrix
  !> 小行列を疎行列に設定（メイン関数、複素数型）
  subroutine monolis_set_block_to_sparse_matrix_main_C(index, item, A, n_dof_index, n_dof_index2, &
    ci, cj, val)
    implicit none
    !> [in] index 配列
    integer(kint), intent(in) :: index(:)
    !> [in] item 配列
    integer(kint), intent(in) :: item(:)
    !> [in,out] 係数行列
    complex(kdouble), intent(inout) :: A(:)
    !> [in] ブロック自由度リスト
    integer(kint), intent(in) :: n_dof_index(:)
    !> [in] ブロック自由度リスト
    integer(kint), intent(in) :: n_dof_index2(:)
    !> [in] 行番号
    integer(kint), intent(in) :: ci
    !> [in] 列番号
    integer(kint), intent(in) :: cj
    !> [in] 設定値
    complex(kdouble), intent(in) :: val(:,:)
    integer(kint) :: j, jn, im, jS, jE, i1, i2, n1, n2

    jS = index(ci) + 1
    jE = index(ci + 1)
    !# 行内の item は昇順ソート済みのため二分探索
    call monolis_bsearch_I(item, jS, jE, cj, j)
    if(j /= -1)then
      n1 = n_dof_index(ci + 1) - n_dof_index(ci)
      n2 = n_dof_index(cj + 1) - n_dof_index(cj)
      do i1 = 1, n1
      do i2 = 1, n2
        im = n_dof_index2(j) + n2*(i1-1) + i2
        A(im) = val(i1, i2)
      enddo
      enddo
      return
    endif
    call monolis_stop_by_matrix_assemble(ci, cj)
  end subroutine monolis_set_block_to_sparse_matrix_main_C

  !# getter
  !> @ingroup dev_matrix
  !> スカラ値を疎行列から取得（メイン関数、実数型）
  subroutine monolis_get_scalar_from_sparse_matrix_main_R(index, item, A, n_dof_index, n_dof_index2, &
    ci, cj, csub_i, csub_j, val, is_find)
    implicit none
    !> [in] index 配列
    integer(kint), intent(in) :: index(:)
    !> [in] item 配列
    integer(kint), intent(in) :: item(:)
    !> [in] 係数行列
    real(kdouble), intent(in) :: A(:)
    !> [in] ブロック自由度リスト
    integer(kint), intent(in) :: n_dof_index(:)
    !> [in] ブロック自由度リスト
    integer(kint), intent(in) :: n_dof_index2(:)
    !> [in] 行番号
    integer(kint), intent(in) :: ci
    !> [in] 列番号
    integer(kint), intent(in) :: cj
    !> [in] ブロック中の行番号
    integer(kint), intent(in) :: csub_i
    !> [in] ブロック中の列番号
    integer(kint), intent(in) :: csub_j
    !> [out] 設定値
    real(kdouble), intent(out) :: val
    !> [out] 取得判定フラグ
    logical, intent(out) :: is_find
    integer(kint) :: j, jn, im, jS, jE, ndof

    val = 0.0d0
    is_find = .false.

    ndof = n_dof_index(ci + 1) - n_dof_index(ci)
    if(ndof < csub_i) return
    ndof = n_dof_index(cj + 1) - n_dof_index(cj)
    if(ndof < csub_j) return

    jS = index(ci) + 1
    jE = index(ci + 1)
    !# 行内の item は昇順ソート済みのため二分探索
    call monolis_bsearch_I(item, jS, jE, cj, j)
    if(j == -1) return

    ndof = n_dof_index(cj + 1) - n_dof_index(cj)
    im = n_dof_index2(j) + ndof*(csub_i-1) + csub_j
    !# GPU 正本モード:ライトスルー方式によりホスト鏡像は常に最新のため、ホスト側の値を返す
    val = A(im)
    is_find = .true.
  end subroutine monolis_get_scalar_from_sparse_matrix_main_R

  !> @ingroup dev_matrix
  !> スカラ値を疎行列から取得（メイン関数、複素数型）
  subroutine monolis_get_scalar_from_sparse_matrix_main_C(index, item, A, n_dof_index, n_dof_index2, &
    ci, cj, csub_i, csub_j, val, is_find)
    implicit none
    !> [in] index 配列
    integer(kint), intent(in) :: index(:)
    !> [in] item 配列
    integer(kint), intent(in) :: item(:)
    !> [in] 係数行列
    complex(kdouble), intent(in) :: A(:)
    !> [in] ブロック自由度リスト
    integer(kint), intent(in) :: n_dof_index(:)
    !> [in] ブロック自由度リスト
    integer(kint), intent(in) :: n_dof_index2(:)
    !> [in] 行番号
    integer(kint), intent(in) :: ci
    !> [in] 列番号
    integer(kint), intent(in) :: cj
    !> [in] ブロック中の行番号
    integer(kint), intent(in) :: csub_i
    !> [in] ブロック中の列番号
    integer(kint), intent(in) :: csub_j
    !> [out] 設定値
    complex(kdouble), intent(out) :: val
    !> [out] 取得判定フラグ
    logical, intent(out) :: is_find
    integer(kint) :: j, jn, im, jS, jE, ndof

    val = 0.0d0
    is_find = .false.

    ndof = n_dof_index(ci + 1) - n_dof_index(ci)
    if(ndof < csub_i) return
    ndof = n_dof_index(cj + 1) - n_dof_index(cj)
    if(ndof < csub_j) return

    jS = index(ci) + 1
    jE = index(ci + 1)
    !# 行内の item は昇順ソート済みのため二分探索
    call monolis_bsearch_I(item, jS, jE, cj, j)
    if(j == -1) return

    ndof = n_dof_index(cj + 1) - n_dof_index(cj)
    im = n_dof_index2(j) + ndof*(csub_i-1) + csub_j
    val = A(im)
    is_find = .true.
  end subroutine monolis_get_scalar_from_sparse_matrix_main_C

  !# adder
  !> @ingroup matrix
  !> スカラ値を疎行列に足込（メイン関数、実数型）
  !> @details GPU 要素ループ（OpenACC compute region）内から直接呼び出し可能（acc routine seq）。
  !> デバイス呼び出し時は常駐配列（index, item, A, n_dof_index, n_dof_index2）を present で渡すこと。
  !> 足し込みは atomic で行われる。パターン外参照等のエラーはカウンタに記録され、
  !> \ref monolis_check_matrix_assemble_error で検査される。
  subroutine monolis_add_scalar_to_sparse_matrix_main_R(index, item, A, n_dof_index, n_dof_index2, &
    ci, cj, csub_i, csub_j, val)
    implicit none
!$acc routine seq
    !> [in] index 配列
    integer(kint), intent(in) :: index(:)
    !> [in] item 配列
    integer(kint), intent(in) :: item(:)
    !> [in,out] 係数行列
    real(kdouble), intent(inout) :: A(:)
    !> [in] ブロック自由度リスト
    integer(kint), intent(in) :: n_dof_index(:)
    !> [in] ブロック自由度リスト
    integer(kint), intent(in) :: n_dof_index2(:)
    !> [in] 行番号
    integer(kint), intent(in) :: ci
    !> [in] 列番号
    integer(kint), intent(in) :: cj
    !> [in] ブロック中の行番号
    integer(kint), intent(in) :: csub_i
    !> [in] ブロック中の列番号
    integer(kint), intent(in) :: csub_j
    !> [in] 設定値
    real(kdouble), intent(in) :: val
    integer(kint) :: j, im, jS, jE, ndof

    ndof = n_dof_index(ci + 1) - n_dof_index(ci)
    if(ndof < csub_i)then
      call monolis_record_matrix_assemble_error(ci, cj)
      return
    endif
    ndof = n_dof_index(cj + 1) - n_dof_index(cj)
    if(ndof < csub_j)then
      call monolis_record_matrix_assemble_error(ci, cj)
      return
    endif

    jS = index(ci) + 1
    jE = index(ci + 1)
    !# 行内の item は昇順ソート済みのため二分探索
    call monolis_bsearch_I(item, jS, jE, cj, j)
    if(j == -1)then
      call monolis_record_matrix_assemble_error(ci, cj)
      return
    endif

    im = n_dof_index2(j) + ndof*(csub_i-1) + csub_j
    !# 複数スレッドから同じ成分に足し込まれる場合があるため atomic
!$acc atomic update
!$omp atomic update
    A(im) = A(im) + val
  end subroutine monolis_add_scalar_to_sparse_matrix_main_R

  !> @ingroup dev_matrix
  !> スカラ値を疎行列に足込（メイン関数、複素数型）
  subroutine monolis_add_scalar_to_sparse_matrix_main_C(index, item, A, n_dof_index, n_dof_index2, &
    ci, cj, csub_i, csub_j, val)
    implicit none
    !> [in] index 配列
    integer(kint), intent(in) :: index(:)
    !> [in] item 配列
    integer(kint), intent(in) :: item(:)
    !> [in,out] 係数行列
    complex(kdouble), intent(inout) :: A(:)
    !> [in] ブロック自由度リスト
    integer(kint), intent(in) :: n_dof_index(:)
    !> [in] ブロック自由度リスト
    integer(kint), intent(in) :: n_dof_index2(:)
    !> [in] 行番号
    integer(kint), intent(in) :: ci
    !> [in] 列番号
    integer(kint), intent(in) :: cj
    !> [in] ブロック中の行番号
    integer(kint), intent(in) :: csub_i
    !> [in] ブロック中の列番号
    integer(kint), intent(in) :: csub_j
    !> [in] 設定値
    complex(kdouble), intent(in) :: val
    integer(kint) :: j, jn, im, jS, jE, ndof

    ndof = n_dof_index(ci + 1) - n_dof_index(ci)
    if(ndof < csub_i) call monolis_stop_by_submatrix_access(ndof, csub_i)
    ndof = n_dof_index(cj + 1) - n_dof_index(cj)
    if(ndof < csub_j) call monolis_stop_by_submatrix_access(ndof, csub_j)

    jS = index(ci) + 1
    jE = index(ci + 1)
    !# 行内の item は昇順ソート済みのため二分探索
    call monolis_bsearch_I(item, jS, jE, cj, j)
    if(j /= -1)then
      ndof = n_dof_index(cj + 1) - n_dof_index(cj)
      im = n_dof_index2(j) + ndof*(csub_i-1) + csub_j
      !# 複数スレッドから同じ成分に足し込まれる場合があるため atomic
!$omp atomic update
      A(im) = A(im) + val
      return
    endif

    call monolis_stop_by_matrix_assemble(ci, cj)
  end subroutine monolis_add_scalar_to_sparse_matrix_main_C

  !> @ingroup matrix
  !> 行列値を疎行列に足込（メイン関数、実数型）
  !> @details 要素剛性行列を全体行列へ足し込む。GPU 要素ループ（OpenACC compute region）内から
  !> 直接呼び出し可能（acc routine seq）。デバイス呼び出し時は常駐配列
  !> （index, item, A, n_dof_index, n_dof_index2）を present で渡し、e1/e2/val は
  !> スレッドローカル（private）配列を渡すこと。各 (行, 列) ブロックは行内二分探索で特定され、
  !> 足し込みは atomic で行われる。パターン外参照のエラーはカウンタに記録され、
  !> \ref monolis_check_matrix_assemble_error で検査される。
  subroutine monolis_add_matrix_to_sparse_matrix_main_R(index, item, A, n1, n2, &
    n_dof_index, n_dof_index2, e1, e2, val)
    implicit none
!$acc routine seq
    !> [in] index 配列
    integer(kint), intent(in) :: index(:)
    !> [in] item 配列
    integer(kint), intent(in) :: item(:)
    !> [in,out] 係数行列
    real(kdouble), intent(inout) :: A(:)
    !> [in] 入力行列の行サイズ
    integer(kint), intent(in) :: n1
    !> [in] 入力行列の列サイズ
    integer(kint), intent(in) :: n2
    !> [in] ブロック自由度リスト
    integer(kint), intent(in) :: n_dof_index(:)
    !> [in] ブロック自由度リスト
    integer(kint), intent(in) :: n_dof_index2(:)
    !> [in] 行番号リスト
    integer(kint), intent(in) :: e1(n1)
    !> [in] 列番号リスト
    integer(kint), intent(in) :: e2(n2)
    !> [in] 設定値
    real(kdouble), intent(in) :: val(:,:)
    integer(kint) :: i1, i2, j, jS, jE, in, jn, kn
    integer(kint) :: k1, k2, im, ioff, joff, ndof1, ndof2

    ioff = 0
    do i1 = 1, n1
      in = e1(i1)
      ndof1 = n_dof_index(in + 1) - n_dof_index(in)
      jS = index(in) + 1
      jE = index(in + 1)
      joff = 0
      do i2 = 1, n2
        jn = e2(i2)
        ndof2 = n_dof_index(jn + 1) - n_dof_index(jn)
        !# 行内の item は昇順ソート済みのため二分探索
        call monolis_bsearch_I(item, jS, jE, jn, j)
        if(j == -1)then
          call monolis_record_matrix_assemble_error(in, jn)
        else
          kn = n_dof_index2(j)
          do k1 = 1, ndof1
          do k2 = 1, ndof2
            im = kn + ndof2*(k1-1) + k2
            !# 複数スレッドから同じ成分に足し込まれる場合があるため atomic
!$acc atomic update
!$omp atomic update
            A(im) = A(im) + val(ioff+k1, joff+k2)
          enddo
          enddo
        endif
        joff = joff + ndof2
      enddo
      ioff = ioff + ndof1
    enddo
  end subroutine monolis_add_matrix_to_sparse_matrix_main_R

  !> @ingroup dev_matrix
  !> 行列値を疎行列に足込（メイン関数、複素数型）
  subroutine monolis_add_matrix_to_sparse_matrix_main_C(index, item, A, n1, n2, &
    n_dof_index, n_dof_index2, e1, e2, val)
    implicit none
    !> [in] index 配列
    integer(kint), intent(in) :: index(:)
    !> [in] item 配列
    integer(kint), intent(in) :: item(:)
    !> [in,out] 係数行列
    complex(kdouble), intent(inout) :: A(:)
    !> [in] 入力行列の行サイズ
    integer(kint), intent(in) :: n1
    !> [in] 入力行列の列サイズ
    integer(kint), intent(in) :: n2
    !> [in] ブロック自由度リスト
    integer(kint), intent(in) :: n_dof_index(:)
    !> [in] ブロック自由度リスト
    integer(kint), intent(in) :: n_dof_index2(:)
    !> [in] 行番号リスト
    integer(kint), intent(in) :: e1(n1)
    !> [in] 列番号リスト
    integer(kint), intent(in) :: e2(n2)
    !> [in] 設定値
    complex(kdouble), intent(in) :: val(:,:)
    integer(kint) :: e1t(n1), e2t(n2)
    integer(kint) :: i, in, jn, im, kn, jS, jE, i1, i2, j1, j2, k1, k2
    integer(kint) :: eperm1(n1), eperm2(n2)
    integer(kint), allocatable :: idx1(:), idx2(:)
    integer(kint), allocatable :: idx1t(:), idx2t(:)
    complex(kdouble), allocatable :: temp(:,:)

    e1t = e1
    call monolis_get_sequence_array_I(eperm1, n1, 1, 1)
    call monolis_qsort_I_2d(e1t, eperm1, 1, n1)

    e2t = e2
    call monolis_get_sequence_array_I(eperm2, n2, 1, 1)
    call monolis_qsort_I_2d(e2t, eperm2, 1, n2)

    call monolis_alloc_I_1d(idx1, n1 + 1)
    call monolis_alloc_I_1d(idx2, n2 + 1)
    call monolis_alloc_I_1d(idx1t, n1 + 1)
    call monolis_alloc_I_1d(idx2t, n2 + 1)

    do i = 1, n1
      in = e1(i)
      idx1(i + 1) = idx1(i) + n_dof_index(in + 1) - n_dof_index(in)
    enddo

    do i = 1, n2
      in = e2(i)
      idx2(i + 1) = idx2(i) + n_dof_index(in + 1) - n_dof_index(in)
    enddo

    do i = 1, n1
      in = e1t(i)
      idx1t(i + 1) = idx1t(i) + n_dof_index(in + 1) - n_dof_index(in)
    enddo

    do i = 1, n2
      in = e2t(i)
      idx2t(i + 1) = idx2t(i) + n_dof_index(in + 1) - n_dof_index(in)
    enddo

    call monolis_alloc_C_2d(temp, idx1t(n1 + 1), idx2t(n2 + 1))

    do i2 = 1, n2
      j2 = eperm2(i2)
      do i1 = 1, n1
        j1 = eperm1(i1)
        do k2 = 1, idx2t(i2 + 1) - idx2t(i2)
        do k1 = 1, idx1t(i1 + 1) - idx1t(i1)
          temp(idx1t(i1)+k1, idx2t(i2)+k2) = val(idx1(j1)+k1, idx2(j2)+k2)
        enddo
        enddo
      enddo
    enddo

    do i1 = 1, n1
      in = e1t(i1)
      jS = index(in) + 1
      jE = index(in + 1)
      aa:do i2 = 1, n2
        do j1 = jS, jE
          jn = item(j1)
          if(jn == e2t(i2))then
            kn = idx2t(i2 + 1) - idx2t(i2)
            do k1 = 1, idx1t(i1 + 1) - idx1t(i1)
            do k2 = 1, kn
              im = n_dof_index2(j1) + kn*(k1-1) + k2
              !# 複数スレッドから同じ成分に足し込まれる場合があるため atomic
!$omp atomic update
              A(im) = A(im) + temp(idx1t(i1)+k1, idx2t(i2)+k2)
            enddo
            enddo
            jS = j1
            cycle aa
          endif
        enddo
      call monolis_stop_by_matrix_assemble(e1t(i1), e2t(i2))
      enddo aa
    enddo
  end subroutine monolis_add_matrix_to_sparse_matrix_main_C

  !> @ingroup dev_matrix
  !> 境界条件処理（実数型、メイン関数）
  subroutine monolis_set_Dirichlet_bc_main_R(index, item, A, B, indexR, itemR, permA, &
    & n_dof_index, n_dof_index2, node_id, ndof_bc, val)
    implicit none
    !> [in] index 配列
    integer(kint), intent(in) :: index(:)
    !> [in] item 配列
    integer(kint), intent(in) :: item(:)
    !> [in,out] 係数行列
    real(kdouble), intent(inout) :: A(:)
    !> [in,out] 右辺ベクトル
    real(kdouble), intent(inout) :: B(:)
    !> [in] index 配列（CSC 形式）
    integer(kint), intent(in) :: indexR(:)
    !> [in] item 配列（CSC 形式）
    integer(kint), intent(in) :: itemR(:)
    !> [in] 行列成分の CSC 形式との置換ベクトル
    integer(kint), intent(in) :: permA(:)
    !> [in] ブロック自由度リスト
    integer(kint), intent(in) :: n_dof_index(:)
    !> [in] ブロック自由度リスト
    integer(kint), intent(in) :: n_dof_index2(:)
    !> [in] 自由度番号
    integer(kint), intent(in) :: node_id
    !> [in] ブロック番号
    integer(kint), intent(in) :: ndof_bc
    !> [in] 境界条件の設定値
    real(kdouble), intent(in) :: val
    integer(kint) :: j, k, jn, kn, m, n, jS, jE, ndof1, ndof2
    logical :: is_add
#ifdef _OPENACC
    logical :: is_dev, is_dev_B
#endif

    ndof1 = n_dof_index(node_id + 1) - n_dof_index(node_id)
    if(ndof1 < ndof_bc) call monolis_stop_by_submatrix_access(ndof1, ndof_bc)

    is_add = .false.

#ifdef _OPENACC
    !# GPU 正本モード:組み立て中の非同期カーネルを待ち、デバイス組立の遅延エラーを検査する
    is_dev = acc_is_present(A)
    is_dev_B = .false.
    if(is_dev) is_dev_B = acc_is_present(B)
    call monolis_check_matrix_assemble_error()

    !# デバイス完結版:右辺 B もデバイス常駐の場合(GPU 要素ループ組立のデバイス専有フロー)、
    !# 列・行の消去、右辺の更新、対角 1 の設定を全てデバイス上で実行する(ホスト鏡像は更新しない)
    if(is_dev_B)then
      jS = indexR(node_id) + 1
      jE = indexR(node_id + 1)
      !< ndof1: 横、ndof2: 縦
!$acc parallel loop private(jn, kn, ndof2, k, m, n) &
!$acc & present(A, B, indexR, itemR, permA, n_dof_index, n_dof_index2)
      do j = jS, jE
        jn = itemR(j)
        kn = permA(j)
        ndof2 = n_dof_index(jn + 1) - n_dof_index(jn)
        do k = 1, ndof2
          m = n_dof_index2(kn) + ndof1*(k-1) + ndof_bc
          n = n_dof_index(jn) + k
          B(n) = B(n) - val*A(m)
          A(m) = 0.0d0
        enddo
      enddo
!$acc end parallel loop

      jS = index(node_id) + 1
      jE = index(node_id + 1)
!$acc parallel loop private(jn, ndof2, k, m) reduction(.or.:is_add) &
!$acc & present(A, index, item, n_dof_index, n_dof_index2)
      do j = jS, jE
        jn = item(j)
        !< ndof2: 横、ndof1: 縦
        ndof2 = n_dof_index(jn + 1) - n_dof_index(jn)
        do k = 1, ndof2
          m = n_dof_index2(j) + ndof2*(ndof_bc-1)
          A(m + k) = 0.0d0
        enddo

        if(jn == node_id)then
          m = n_dof_index2(j) + (ndof2+1)*(ndof_bc-1) + 1
          A(m) = 1.0d0
          is_add = .true.
        endif
      enddo
!$acc end parallel loop

      if(.not. is_add)then
        call monolis_stop_by_set_DBC(node_id)
      endif

!$acc serial present(B, n_dof_index)
      B(n_dof_index(node_id) + ndof_bc) = val
!$acc end serial
      return
    endif
#endif

    jS = indexR(node_id) + 1
    jE = indexR(node_id + 1)
    !< ndof1: 横、ndof2: 縦
    do j = jS, jE
      jn = itemR(j)
      kn = permA(j)
      ndof2 = n_dof_index(jn + 1) - n_dof_index(jn)
      do k = 1, ndof2
        m = n_dof_index2(kn) + ndof1*(k-1) + ndof_bc
        n = n_dof_index(jn) + k
        B(n) = B(n) - val*A(m)
        A(m) = 0.0d0
      enddo
    enddo

#ifdef _OPENACC
    if(is_dev)then
!$acc parallel loop private(jn, kn, ndof2, k, m) async(1) &
!$acc & present(A, indexR, itemR, permA, n_dof_index, n_dof_index2)
    do j = jS, jE
      jn = itemR(j)
      kn = permA(j)
      ndof2 = n_dof_index(jn + 1) - n_dof_index(jn)
      do k = 1, ndof2
        m = n_dof_index2(kn) + ndof1*(k-1) + ndof_bc
        A(m) = 0.0d0
      enddo
    enddo
!$acc end parallel loop
    endif
#endif

    jS = index(node_id) + 1
    jE = index(node_id + 1)
    do j = jS, jE
      jn = item(j)
      !< ndof2: 横、ndof1: 縦
      ndof2 = n_dof_index(jn + 1) - n_dof_index(jn)
      do k = 1, ndof2
        m = n_dof_index2(j) + ndof2*(ndof_bc-1)
        A(m + k) = 0.0d0
      enddo

      if(jn == node_id)then
        m = n_dof_index2(j) + (ndof2+1)*(ndof_bc-1) + 1
        A(m) = 1.0d0
        is_add = .true.
      endif
    enddo

#ifdef _OPENACC
    if(is_dev)then
!$acc parallel loop private(jn, ndof2, k, m) async(1) &
!$acc & present(A, index, item, n_dof_index, n_dof_index2)
    do j = jS, jE
      jn = item(j)
      ndof2 = n_dof_index(jn + 1) - n_dof_index(jn)
      do k = 1, ndof2
        m = n_dof_index2(j) + ndof2*(ndof_bc-1)
        A(m + k) = 0.0d0
      enddo

      if(jn == node_id)then
        m = n_dof_index2(j) + (ndof2+1)*(ndof_bc-1) + 1
        A(m) = 1.0d0
      endif
    enddo
!$acc end parallel loop
    !$acc wait(1)
    endif
#endif

    if(.not. is_add)then
      call monolis_stop_by_set_DBC(node_id)
    endif

    B(n_dof_index(node_id) + ndof_bc) = val
  end subroutine monolis_set_Dirichlet_bc_main_R

  !> @ingroup matrix
  !> Dirichlet 境界条件の一括付与（実数型、メイン関数）
  !> @details 全自由度の固定フラグと設定値を受け取り、全行を 1 パスで処理する
  !> （行オーナーシップ方式。1 スレッドが 1 ブロック行の A・B を専有するため atomic・CSC 不要）。
  !> 既存の \ref monolis_set_Dirichlet_bc_main_R を固定自由度へ昇順に適用した結果と数学的に同値。
  !> GPU 正本モードでは、行列 A と右辺 B が共にデバイス常駐の場合は全てデバイス上で実行し、
  !> A のみ常駐の場合はホストで実行して A の消去をデバイスへ反映する（ライトスルー）。
  subroutine monolis_set_Dirichlet_bc_all_main_R(index, item, A, B, &
    & n_dof_index, n_dof_index2, NP, is_bc, bc_val)
    implicit none
    !> [in] index 配列
    integer(kint), intent(in) :: index(:)
    !> [in] item 配列
    integer(kint), intent(in) :: item(:)
    !> [in,out] 係数行列
    real(kdouble), intent(inout) :: A(:)
    !> [in,out] 右辺ベクトル
    real(kdouble), intent(inout) :: B(:)
    !> [in] ブロック自由度リスト
    integer(kint), intent(in) :: n_dof_index(:)
    !> [in] ブロック自由度リスト
    integer(kint), intent(in) :: n_dof_index2(:)
    !> [in] 全ブロック行数
    integer(kint), intent(in) :: NP
    !> [in] 自由度ごとの固定フラグ（長さ n_dof_index(NP+1)）
    logical, intent(in) :: is_bc(:)
    !> [in] 固定自由度に課す値（長さ n_dof_index(NP+1)。非固定自由度の値は参照されない）
    real(kdouble), intent(in) :: bc_val(:)
    integer(kint) :: i, j, jS, jE, jn, k1, k2, r, c, m, ndof1, ndof2
    integer(kint) :: NPNDOF, missing_row
    logical :: found
#ifdef _OPENACC
    logical :: is_dev, is_dev_B
#endif

    NPNDOF = n_dof_index(NP + 1)
    missing_row = NP + 1

#ifdef _OPENACC
    !# GPU 正本モード:組み立て中の非同期カーネルを待ち、デバイス組立の遅延エラーを検査する
    is_dev = acc_is_present(A)
    is_dev_B = .false.
    if(is_dev) is_dev_B = acc_is_present(B)
    call monolis_check_matrix_assemble_error()

    !# デバイス完結版:右辺 B もデバイス常駐の場合(デバイス専有フロー)、全行の処理を
    !# 1 カーネルで実行する(ホスト鏡像は更新しない)。固定フラグ・設定値はホスト正本とし、
    !# 既にデバイス常駐の場合も最新のホスト値で上書きする
    if(is_dev_B)then
!$acc data copyin(is_bc(1:NPNDOF), bc_val(1:NPNDOF))
!$acc update device(is_bc(1:NPNDOF), bc_val(1:NPNDOF))
!$acc parallel loop private(j, jS, jE, jn, k1, k2, r, c, m, ndof1, ndof2, found) &
!$acc & present(A, B, index, item, n_dof_index, n_dof_index2, is_bc, bc_val) &
!$acc & reduction(min:missing_row)
      do i = 1, NP
        ndof1 = n_dof_index(i + 1) - n_dof_index(i)
        jS = index(i) + 1
        jE = index(i + 1)
        do k1 = 1, ndof1
          r = n_dof_index(i) + k1
          if(is_bc(r))then
            !# 固定行:行の全成分をゼロ化し、対角成分に 1 を設定する
            found = .false.
            do j = jS, jE
              jn = item(j)
              ndof2 = n_dof_index(jn + 1) - n_dof_index(jn)
              do k2 = 1, ndof2
                A(n_dof_index2(j) + ndof2*(k1-1) + k2) = 0.0d0
              enddo
              if(jn == i)then
                A(n_dof_index2(j) + ndof2*(k1-1) + k1) = 1.0d0
                found = .true.
              endif
            enddo
            B(r) = bc_val(r)
            if(.not. found .and. i < missing_row) missing_row = i
          else
            !# 自由行:固定列の成分を右辺へ移項してゼロ化する(列昇順)
            do j = jS, jE
              jn = item(j)
              ndof2 = n_dof_index(jn + 1) - n_dof_index(jn)
              do k2 = 1, ndof2
                c = n_dof_index(jn) + k2
                if(is_bc(c))then
                  m = n_dof_index2(j) + ndof2*(k1-1) + k2
                  B(r) = B(r) - bc_val(c)*A(m)
                  A(m) = 0.0d0
                endif
              enddo
            enddo
          endif
        enddo
      enddo
!$acc end parallel loop
!$acc end data

      if(missing_row <= NP)then
        call monolis_stop_by_set_DBC(missing_row)
      endif
      return
    endif
#endif

    !# ホスト実行(デバイス専有フロー以外)
!$omp parallel default(none) &
!$omp & shared(index, item, A, B, n_dof_index, n_dof_index2, is_bc, bc_val) &
!$omp & firstprivate(NP) &
!$omp & private(i, j, jS, jE, jn, k1, k2, r, c, m, ndof1, ndof2, found) &
!$omp & reduction(min:missing_row)
!$omp do
    do i = 1, NP
      ndof1 = n_dof_index(i + 1) - n_dof_index(i)
      jS = index(i) + 1
      jE = index(i + 1)
      do k1 = 1, ndof1
        r = n_dof_index(i) + k1
        if(is_bc(r))then
          !# 固定行:行の全成分をゼロ化し、対角成分に 1 を設定する
          found = .false.
          do j = jS, jE
            jn = item(j)
            ndof2 = n_dof_index(jn + 1) - n_dof_index(jn)
            do k2 = 1, ndof2
              A(n_dof_index2(j) + ndof2*(k1-1) + k2) = 0.0d0
            enddo
            if(jn == i)then
              A(n_dof_index2(j) + ndof2*(k1-1) + k1) = 1.0d0
              found = .true.
            endif
          enddo
          B(r) = bc_val(r)
          if(.not. found .and. i < missing_row) missing_row = i
        else
          !# 自由行:固定列の成分を右辺へ移項してゼロ化する(列昇順)
          do j = jS, jE
            jn = item(j)
            ndof2 = n_dof_index(jn + 1) - n_dof_index(jn)
            do k2 = 1, ndof2
              c = n_dof_index(jn) + k2
              if(is_bc(c))then
                m = n_dof_index2(j) + ndof2*(k1-1) + k2
                B(r) = B(r) - bc_val(c)*A(m)
                A(m) = 0.0d0
              endif
            enddo
          enddo
        endif
      enddo
    enddo
!$omp end do
!$omp end parallel

    if(missing_row <= NP)then
      call monolis_stop_by_set_DBC(missing_row)
    endif

#ifdef _OPENACC
    !# ライトスルー:デバイス常駐の行列 A に消去と対角 1 の設定を反映する(B はホストのみ)
    if(is_dev)then
!$acc data copyin(is_bc(1:NPNDOF))
!$acc update device(is_bc(1:NPNDOF))
!$acc parallel loop private(j, jS, jE, jn, k1, k2, c, ndof1, ndof2) &
!$acc & present(A, index, item, n_dof_index, n_dof_index2, is_bc)
      do i = 1, NP
        ndof1 = n_dof_index(i + 1) - n_dof_index(i)
        jS = index(i) + 1
        jE = index(i + 1)
        do k1 = 1, ndof1
          if(is_bc(n_dof_index(i) + k1))then
            do j = jS, jE
              jn = item(j)
              ndof2 = n_dof_index(jn + 1) - n_dof_index(jn)
              do k2 = 1, ndof2
                A(n_dof_index2(j) + ndof2*(k1-1) + k2) = 0.0d0
              enddo
              if(jn == i)then
                A(n_dof_index2(j) + ndof2*(k1-1) + k1) = 1.0d0
              endif
            enddo
          else
            do j = jS, jE
              jn = item(j)
              ndof2 = n_dof_index(jn + 1) - n_dof_index(jn)
              do k2 = 1, ndof2
                c = n_dof_index(jn) + k2
                if(is_bc(c))then
                  A(n_dof_index2(j) + ndof2*(k1-1) + k2) = 0.0d0
                endif
              enddo
            enddo
          endif
        enddo
      enddo
!$acc end parallel loop
!$acc end data
    endif
#endif
  end subroutine monolis_set_Dirichlet_bc_all_main_R

  !> @ingroup dev_matrix
  !> 境界条件処理（複素数型、メイン関数）
  subroutine monolis_set_Dirichlet_bc_main_C(index, item, A, B, indexR, itemR, permA, &
    & n_dof_index, n_dof_index2, node_id, ndof_bc, val)
    implicit none
    !> [in] index 配列
    integer(kint), intent(in) :: index(:)
    !> [in] item 配列
    integer(kint), intent(in) :: item(:)
    !> [in,out] 係数行列
    complex(kdouble), intent(inout) :: A(:)
    !> [in,out] 右辺ベクトル
    complex(kdouble), intent(inout) :: B(:)
    !> [in] index 配列（CSC 形式）
    integer(kint), intent(in) :: indexR(:)
    !> [in] item 配列（CSC 形式）
    integer(kint), intent(in) :: itemR(:)
    !> [in] 行列成分の CSC 形式との置換ベクトル
    integer(kint), intent(in) :: permA(:)
    !> [in] ブロック自由度リスト
    integer(kint), intent(in) :: n_dof_index(:)
    !> [in] ブロック自由度リスト
    integer(kint), intent(in) :: n_dof_index2(:)
    !> [in] 自由度番号
    integer(kint), intent(in) :: node_id
    !> [in] ブロック番号
    integer(kint), intent(in) :: ndof_bc
    !> [in] 境界条件の設定値
    complex(kdouble), intent(in) :: val
    integer(kint) :: j, k, jn, kn, jS, jE, ndof1, ndof2
    logical :: is_add

    ndof1 = n_dof_index(node_id + 1) - n_dof_index(node_id)
    if(ndof1 < ndof_bc) call monolis_stop_by_submatrix_access(ndof1, ndof_bc)

    is_add = .false.

    jS = indexR(node_id) + 1
    jE = indexR(node_id + 1)
    do j = jS, jE
      jn = itemR(j)
      kn = permA(j)
      ndof2 = n_dof_index(jn + 1) - n_dof_index(jn)
      do k = 1, ndof2
        B(n_dof_index(jn)+k) = B(n_dof_index(jn)+k) - val*A(n_dof_index2(kn) + ndof1*(k-1) + ndof_bc)
        A(n_dof_index2(kn) + ndof1*(k-1) + ndof_bc) = 0.0d0
      enddo
    enddo

    jS = index(node_id) + 1
    jE = index(node_id + 1)
    do j = jS, jE
      jn = item(j)
      ndof2 = n_dof_index(jn + 1) - n_dof_index(jn)
      do k = 1, ndof2
        A(n_dof_index2(j) + ndof1*(ndof_bc-1) + k) = 0.0d0
      enddo

      jn = item(j)
      if(jn == node_id)then
        A(n_dof_index2(j) + (ndof1+1)*(ndof_bc-1) + 1) = 1.0d0
        is_add = .true.
      endif
    enddo

    if(.not. is_add)then
      call monolis_stop_by_set_DBC(node_id)
    endif

    B(n_dof_index(node_id) + ndof_bc) = val
  end subroutine monolis_set_Dirichlet_bc_main_C

  !> @ingroup dev_matrix
  !> 非零要素がメモリ確保されていない場合にエラーストップ
  subroutine monolis_stop_by_matrix_assemble(ci, cj)
    !> [in] 行番号
    integer(kint), intent(in) :: ci
    !> [in] 列番号
    integer(kint), intent(in) :: cj

    call monolis_std_error_string("")
    write(*,"(a,i0,a,i0,a)") "The non-zero element at (", ci, ", ", cj, &
      & ") is not allocated. The value is not accessible."
    call monolis_std_error_stop()
  end subroutine monolis_stop_by_matrix_assemble

  !> @ingroup dev_matrix
  !> ブロックサイズより大きい部分へのアクサス時にエラーストップ
  subroutine monolis_stop_by_submatrix_access(ndof, sub_dof)
    !> [in] 行番号
    integer(kint), intent(in) :: ndof
    !> [in] ブロック番号
    integer(kint), intent(in) :: sub_dof

    call monolis_std_error_string("")
    write(*,"(a)")   "input value is greater than the DoF of submatrix."
    write(*,"(a,i8)")"       the DoF of submatrix: ", ndof
    write(*,"(a,i8)")"       value:                ", sub_dof
    call monolis_std_error_stop()
  end subroutine monolis_stop_by_submatrix_access

  !> @ingroup dev_matrix
  !> 対角成分が見つからない場合にエラーストップ
  subroutine monolis_stop_by_set_DBC(node_id)
    !> [in] 行番号
    integer(kint), intent(in) :: node_id

    call monolis_std_error_string("")
    write(*,"(a)")   "The diagonal component is not found."
    write(*,"(a,i8)")"       tdiagonal number     : ", node_id
    call monolis_std_error_stop()
  end subroutine monolis_stop_by_set_DBC

  !> @ingroup dev_matrix
  !> 対角成分が零の場合にエラーストップ
  subroutine monolis_stop_by_set_zero_diag_component(node_id, ndof)
    !> [in] 行番号
    integer(kint), intent(in) :: node_id
    !> [in] ブロック自由度
    integer(kint), intent(in) :: ndof

    call monolis_std_error_string("")
    write(*,"(a)")   "The diagonal component is zero."
    write(*,"(a,i8)")"       diagonal number     : ", node_id
    write(*,"(a,i8)")"       the DoF of submatrix: ", ndof
    call monolis_std_error_stop()
  end subroutine monolis_stop_by_set_zero_diag_component

  !> @ingroup dev_matrix
  !> 行列の最大値の取得（実数型、メイン関数）
  subroutine monolis_get_max_matrix_component_main_R(monoMAT, monoCOM, max_val)
    implicit none
    !> [in] monolis MAT 構造体
    type(monolis_mat), intent(in) :: monoMAT
    !> [in] monolis COM 構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [out] 最大値
    real(kdouble), intent(out) :: max_val

    max_val = maxval(monoMAT%R%A)
    call monolis_allreduce_R1(max_val, monolis_mpi_max, monoCOM%comm)
  end subroutine monolis_get_max_matrix_component_main_R

  !> @ingroup dev_matrix
  !> 行列の対角成分の 0 チェック（実数型、メイン関数）
  subroutine monolis_check_diagonal_zero_component_main_R(monoPRM, monoMAT)
    implicit none
    !> [in] monolis PRM 構造体
    type(monolis_prm), intent(in) :: monoPRM
    !> [in] monolis MAT 構造体
    type(monolis_mat), intent(in) :: monoMAT
    integer(kint) :: i, j, k, jS, jE, in, kn, N, n_dof

    if(monoPRM%Iarray(monolis_prm_I_is_check_diag) == monolis_I_false) return

    N = monoMAT%N

!$omp parallel default(none) &
!$omp & shared(monoMAT) &
!$omp & firstprivate(N) &
!$omp & private(i, j, k, jS, jE, in, kn, n_dof)
!$omp do
    do i = 1, N
      jS = monoMAT%CSR%index(i) + 1
      jE = monoMAT%CSR%index(i + 1)
      do j = jS, jE
        in = monoMAT%CSR%item(j)
        if(i == in)then
          n_dof = monoMAT%n_dof_index(i + 1) - monoMAT%n_dof_index(i)
          do k = 1, n_dof
            kn = monoMAT%n_dof_index2(j) + (n_dof+1)*(k-1) + 1
            if(monoMAT%R%A(kn) == 0.0d0)then
              call monolis_stop_by_set_zero_diag_component(i, k)
              stop
            endif
          enddo
        endif
      enddo
    enddo
!$omp end do
!$omp end parallel
  end subroutine monolis_check_diagonal_zero_component_main_R

  !> @ingroup dev_matrix
  !> 疎行列から密行列に変換（実数型、メイン関数）
  subroutine monolis_convert_sparse_matrix_to_dense_matrix_R(monoMAT, monoCOM, dense)
    implicit none
    !> [in] monolis MAT 構造体
    type(monolis_mat), intent(in) :: monoMAT
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [out] 疎行列を変換した密行列
    real(kdouble), intent(inout), allocatable :: dense(:,:)
    integer(kint) :: N, NT
    integer(kint) :: i, j, k1, k2, jS, jE, jn, kn, n1, n2, nd1, nd2
    integer(kint) :: in, my_rank
    integer(kint), allocatable :: vertex_id(:)
    integer(kint), allocatable :: vtxdist(:), gindex(:)
    real(kdouble) :: t1

    N = monoMAT%N
    if(monolis_mpi_get_local_comm_size(monoCOM%comm) > 1) N = monoCOM%n_internal_vertex

    call monolis_alloc_I_1d(vertex_id, monoMAT%NP)
    call monolis_generate_global_vertex_id(N, monoMAT%NP, vertex_id, monoCOM)

    NT = monoMAT%n_dof_index(N + 1)
    call monolis_allreduce_I1(NT, monolis_mpi_sum, monoCOM%comm)
    call monolis_alloc_R_2d(dense, monoMAT%n_dof_index(N + 1), NT)

    in = monolis_mpi_get_local_comm_size(monoCOM%comm)
    call monolis_alloc_I_1d(vtxdist, in + 1)
    call monolis_com_n_vertex_list(monoMAT%n_dof_index(N + 1), monoCOM%comm, vtxdist)

    my_rank = monolis_mpi_get_local_my_rank(monoCOM%comm)
    call monolis_alloc_I_1d(gindex, monoMAT%NP)
    gindex = vtxdist(my_rank + 1)
    do i = 1, monoMAT%N - 1
      gindex(i + 1) = gindex(i) + monoMAT%n_dof_list(i)
    enddo
    call monolis_mpi_update_I(monoCOM, 1, gindex, t1)

!$omp parallel default(none) &
!$omp & shared(monoMAT, dense, gindex) &
!$omp & firstprivate(N) &
!$omp & private(i, j, k1, k2, jS, jE, jn, kn, n1, n2, nd1, nd2)
!$omp do
    do i = 1, N
      jS = monoMAT%CSR%index(i) + 1
      jE = monoMAT%CSR%index(i + 1)
      do j = jS, jE
        jn = monoMAT%CSR%item(j)
        kn = monoMAT%n_dof_index2(j)
        nd1 = monoMAT%n_dof_list(i)
        nd2 = monoMAT%n_dof_list(jn)
        do k1 = 1, nd1
        do k2 = 1, nd2
          n1 = monoMAT%n_dof_index(i) + k1
          n2 = gindex(jn) + k2
          dense(n1, n2) = monoMAT%R%A(kn + nd2*(k1 - 1) + k2)
        enddo
        enddo
      enddo
    enddo
!$omp end do
!$omp end parallel
  end subroutine monolis_convert_sparse_matrix_to_dense_matrix_R

  !> @ingroup dev_matrix
  !> CSR 形式の疎行列を CSC 形式に変換（実数型）
  subroutine monolis_convert_CSR_matrix_to_CSC_matrix_R(monoCOM, A_in, A_out)
    implicit none
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in] monolis MAT 構造体
    type(monolis_mat), intent(in) :: A_in
    !> [in] monolis MAT 構造体
    type(monolis_mat), intent(inout) :: A_out
    integer(kint) :: i, j, k, jS, jE, col, ndof, nnz, N, NP
    integer(kint), allocatable :: col_count(:), col_offset(:), perm(:)

    N = A_in%N
    NP = A_in%NP
    nnz = A_in%CSR%index(NP + 1)
    ndof = A_in%NDOF

    A_out%N = N
    A_out%NP = NP
    A_out%NDOF = ndof

    call monolis_palloc_I_1d(A_out%CSR%index, NP + 1)
    call monolis_palloc_I_1d(A_out%CSR%item, nnz)
    call monolis_palloc_R_1d(A_out%R%A, ndof*ndof*nnz)
    call monolis_alloc_I_1d(col_count, NP)
    call monolis_alloc_I_1d(perm, nnz)

    do i = 1, N
      jS = A_in%CSR%index(i) + 1
      jE = A_in%CSR%index(i + 1)
      do j = jS, jE
        col = A_in%CSR%item(j)
        col_count(col) = col_count(col) + 1
      enddo
    enddo

    do i = 1, NP
      A_out%CSR%index(i + 1) = A_out%CSR%index(i) + col_count(i)
    enddo

    call monolis_alloc_I_1d(col_offset, NP)

    !# CSR を走査して CSC の item を埋める
    do i = 1, N
      jS = A_in%CSR%index(i) + 1
      jE = A_in%CSR%index(i + 1)
      do j = jS, jE
        col = A_in%CSR%item(j)
        k = A_out%CSR%index(col) + col_offset(col) + 1
        A_out%CSR%item(k) = i
        perm(k) = j
        col_offset(col) = col_offset(col) + 1
      enddo
    enddo

    !# 行列値をコピー（perm を使って並び替え）
!$omp parallel default(none) &
!$omp & shared(A_in, A_out, perm) &
!$omp & firstprivate(nnz) &
!$omp & private(i)
!$omp do
    do i = 1, nnz
      A_out%R%A(i) = A_in%R%A(perm(i))
    enddo
!$omp end do
!$omp end parallel
  end subroutine monolis_convert_CSR_matrix_to_CSC_matrix_R

end module mod_monolis_spmat_handler_util
