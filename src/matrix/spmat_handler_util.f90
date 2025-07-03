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

  implicit none

contains

  !# setter
  !> @ingroup dev_matrix
  !> スカラ値を疎行列に設定（メイン関数、実数型）
  subroutine monolis_set_scalar_to_sparse_matrix_main_R(index, item, A, n_dof_index, n_dof_index2, &
    ci, cj, csub_i, csub_j, val)
    implicit none
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
    integer(kint) :: j, jn, im, jS, jE, ndof

    ndof = n_dof_index(ci + 1) - n_dof_index(ci)
    if(ndof < csub_i) call monolis_stop_by_submatrix_access(ndof, csub_i)
    ndof = n_dof_index(cj + 1) - n_dof_index(cj)
    if(ndof < csub_j) call monolis_stop_by_submatrix_access(ndof, csub_j)

    jS = index(ci) + 1
    jE = index(ci + 1)
    do j = jS, jE
      jn = item(j)
      if(jn == cj)then
        ndof = n_dof_index(cj + 1) - n_dof_index(cj)
        im = n_dof_index2(j) + ndof*(csub_i-1) + csub_j
        A(im) = val
        return
      endif
    enddo

    call monolis_stop_by_matrix_assemble(ci, cj)
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
    do j = jS, jE
      jn = item(j)
      if(jn == cj)then
        ndof = n_dof_index(cj + 1) - n_dof_index(cj)
        im = n_dof_index2(j) + ndof*(csub_i-1) + csub_j
        A(im) = val
        return
      endif
    enddo

    call monolis_stop_by_matrix_assemble(ci, cj)
  end subroutine monolis_set_scalar_to_sparse_matrix_main_C

  !> @ingroup dev_matrix
  !> 小行列を疎行列に設定（メイン関数、実数型）
  subroutine monolis_set_block_to_sparse_matrix_main_R(index, item, A, n_dof_index, n_dof_index2, &
    ci, cj, val)
    implicit none
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
    integer(kint) :: j, jn, im, jS, jE, i1, i2, n1, n2

    jS = index(ci) + 1
    jE = index(ci + 1)
    do j = jS, jE
      jn = item(j)
      if(jn == cj)then
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
    enddo
    call monolis_stop_by_matrix_assemble(ci, cj)
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
    do j = jS, jE
      jn = item(j)
      if(jn == cj)then
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
    enddo
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
    do j = jS, jE
      jn = item(j)
      if(jn == cj)then
        ndof = n_dof_index(cj + 1) - n_dof_index(cj)
        im = n_dof_index2(j) + ndof*(csub_i-1) + csub_j
        val = A(im)
        is_find = .true.
        return
      endif
    enddo
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
    do j = jS, jE
      jn = item(j)
      if(jn == cj)then
        ndof = n_dof_index(cj + 1) - n_dof_index(cj)
        im = n_dof_index2(j) + ndof*(csub_i-1) + csub_j
        val = A(im)
        is_find = .true.
        return
      endif
    enddo
  end subroutine monolis_get_scalar_from_sparse_matrix_main_C

  !# adder
  !> @ingroup dev_matrix
  !> スカラ値を疎行列に足込（メイン関数、実数型）
  subroutine monolis_add_scalar_to_sparse_matrix_main_R(index, item, A, n_dof_index, n_dof_index2, &
    ci, cj, csub_i, csub_j, val)
    implicit none
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
    integer(kint) :: j, jn, im, jS, jE, ndof

    ndof = n_dof_index(ci + 1) - n_dof_index(ci)
    if(ndof < csub_i) call monolis_stop_by_submatrix_access(ndof, csub_i)
    ndof = n_dof_index(cj + 1) - n_dof_index(cj)
    if(ndof < csub_j) call monolis_stop_by_submatrix_access(ndof, csub_j)

    jS = index(ci) + 1
    jE = index(ci + 1)
    do j = jS, jE
      jn = item(j)
      if(jn == cj)then
        ndof = n_dof_index(cj + 1) - n_dof_index(cj)
        im = n_dof_index2(j) + ndof*(csub_i-1) + csub_j
        A(im) = A(im) + val
        return
      endif
    enddo

    call monolis_stop_by_matrix_assemble(ci, cj)
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
    do j = jS, jE
      jn = item(j)
      if(jn == cj)then
        ndof = n_dof_index(cj + 1) - n_dof_index(cj)
        im = n_dof_index2(j) + ndof*(csub_i-1) + csub_j
        A(im) = A(im) + val
        return
      endif
    enddo

    call monolis_stop_by_matrix_assemble(ci, cj)
  end subroutine monolis_add_scalar_to_sparse_matrix_main_C

  !> @ingroup dev_matrix
  !> 行列値を疎行列に足込（メイン関数、実数型）
  subroutine monolis_add_matrix_to_sparse_matrix_main_R(index, item, A, n1, n2, &
    n_dof_index, n_dof_index2, e1, e2, val)
    implicit none
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
    integer(kint) :: e1t(n1), e2t(n2)
    integer(kint) :: i, in, jn, im, kn, jS, jE, i1, i2, j1, j2, k1, k2
    integer(kint) :: eperm1(n1), eperm2(n2)
    integer(kint), allocatable :: idx1(:), idx2(:)
    integer(kint), allocatable :: idx1t(:), idx2t(:)
    real(kdouble), allocatable :: temp(:,:)

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

    call monolis_alloc_R_2d(temp, idx1t(n1 + 1), idx2t(n2 + 1))

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

    call monolis_alloc_C_2d(temp, idx1(n1 + 1), idx2(n2 + 1))

    do i2 = 1, n2
      j2 = eperm2(i2)
      do i1 = 1, n1
        j1 = eperm1(i1)
        do k2 = 1, idx2t(i2 + 1) - idx2t(i2)
        do k1 = 1, idx1t(i1 + 1) - idx1t(i1)
          temp(idx1(i1)+k1, idx2(i2)+k2) = val(idx1t(j1)+k1, idx2t(j2)+k2)
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
              A(im) = A(im) + temp(idx1(i1)+k1, idx2(i2)+k2)
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
  end subroutine monolis_set_Dirichlet_bc_main_R

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

    N =  monoMAT%N

    do i = 1, N
      jS = monoMAT%CSR%index(i) + 1
      jE = monoMAT%CSR%index(i + 1)
      do j = jS, jE
        in = monoMAT%CSR%item(j)
        if(i == in)then
          n_dof = monoMAT%n_dof_index(in + 1) - monoMAT%n_dof_index(in)
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
  end subroutine monolis_convert_sparse_matrix_to_dense_matrix_R
end module mod_monolis_spmat_handler_util
