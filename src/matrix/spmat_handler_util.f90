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
  !> スカラ値を疎行列に設定（メイン関数、実数型）
  subroutine monolis_set_scalar_to_sparse_matrix_main_R(index, item, A, ndof, ci, cj, csub_i, csub_j, val)
    implicit none
    !> index 配列
    integer(kint), intent(in) :: index(:)
    !> item 配列
    integer(kint), intent(in) :: item(:)
    !> 係数行列
    real(kdouble), intent(inout) :: A(:)
    !> ブロック自由度
    integer(kint), intent(in) :: ndof
    !> 行番号
    integer(kint), intent(in) :: ci
    !> 列番号
    integer(kint), intent(in) :: cj
    !> ブロック中の行番号
    integer(kint), intent(in) :: csub_i
    !> ブロック中の列番号
    integer(kint), intent(in) :: csub_j
    !> 設定値
    real(kdouble), intent(in) :: val
    integer(kint) :: j, jn, im, jS, jE, NDOF2

    if(ndof < csub_i) call monolis_stop_by_submatrix_access(ndof, csub_i)
    if(ndof < csub_j) call monolis_stop_by_submatrix_access(ndof, csub_j)

    NDOF2 = ndof*ndof

    jS = index(ci) + 1
    jE = index(ci + 1)
    do j = jS, jE
      jn = item(j)
      if(jn == cj)then
        im = NDOF2*(j-1) + ndof*(csub_i-1) + csub_j
        A(im) = val
        return
      endif
    enddo

    call monolis_stop_by_matrix_assemble(ci, cj)
  end subroutine monolis_set_scalar_to_sparse_matrix_main_R

  !> スカラ値を疎行列に設定（メイン関数、複素数型）
  subroutine monolis_set_scalar_to_sparse_matrix_main_C(index, item, A, ndof, ci, cj, csub_i, csub_j, val)
    implicit none
    !> index 配列
    integer(kint), intent(in) :: index(:)
    !> item 配列
    integer(kint), intent(in) :: item(:)
    !> 係数行列
    complex(kdouble), intent(inout) :: A(:)
    !> ブロック自由度
    integer(kint), intent(in) :: ndof
    !> 行番号
    integer(kint), intent(in) :: ci
    !> 列番号
    integer(kint), intent(in) :: cj
    !> ブロック中の行番号
    integer(kint), intent(in) :: csub_i
    !> ブロック中の列番号
    integer(kint), intent(in) :: csub_j
    !> 設定値
    complex(kdouble), intent(in) :: val
    integer(kint) :: j, jn, im, jS, jE, NDOF2

    if(ndof < csub_i) call monolis_stop_by_submatrix_access(ndof, csub_i)
    if(ndof < csub_j) call monolis_stop_by_submatrix_access(ndof, csub_j)

    NDOF2 = ndof*ndof

    jS = index(ci) + 1
    jE = index(ci + 1)
    do j = jS, jE
      jn = item(j)
      if(jn == cj)then
        im = NDOF2*(j-1) + ndof*(csub_i-1) + csub_j
        A(im) = val
        return
      endif
    enddo

    call monolis_stop_by_matrix_assemble(ci, cj)
  end subroutine monolis_set_scalar_to_sparse_matrix_main_C

  !> 小行列を疎行列に設定（メイン関数、実数型）
  subroutine monolis_set_block_to_sparse_matrix_main_R(index, item, A, ndof, ci, cj, val)
    implicit none
    !> index 配列
    integer(kint), intent(in) :: index(:)
    !> item 配列
    integer(kint), intent(in) :: item(:)
    !> 係数行列
    real(kdouble), intent(inout) :: A(:)
    !> ブロック自由度
    integer(kint), intent(in) :: ndof
    !> 行番号
    integer(kint), intent(in) :: ci
    !> 列番号
    integer(kint), intent(in) :: cj
    !> 設定値
    real(kdouble), intent(in) :: val(ndof,ndof)
    integer(kint) :: k, jn, im, jS, jE, i1, i2, NDOF2

    NDOF2 = ndof*ndof

    jS = index(ci) + 1
    jE = index(ci + 1)
    do k = jS, jE
      jn = item(k)
      if(jn == cj)then
        do i1 = 1, ndof
        do i2 = 1, ndof
          im = NDOF2*(k-1) + ndof*(i1-1) + i2
          A(im) = val(i1, i2)
        enddo
        enddo
        return
      endif
    enddo
    call monolis_stop_by_matrix_assemble(ci, cj)
  end subroutine monolis_set_block_to_sparse_matrix_main_R

  !> 小行列を疎行列に設定（メイン関数、複素数型）
  subroutine monolis_set_block_to_sparse_matrix_main_C(index, item, A, ndof, ci, cj, val)
    implicit none
    !> index 配列
    integer(kint), intent(in) :: index(:)
    !> item 配列
    integer(kint), intent(in) :: item(:)
    !> 係数行列
    complex(kdouble), intent(inout) :: A(:)
    !> ブロック自由度
    integer(kint), intent(in) :: ndof
    !> 行番号
    integer(kint), intent(in) :: ci
    !> 列番号
    integer(kint), intent(in) :: cj
    !> 設定値
    complex(kdouble), intent(in) :: val(ndof,ndof)
    integer(kint) :: k, jn, im, jS, jE, i1, i2, NDOF2

    NDOF2 = ndof*ndof

    jS = index(ci) + 1
    jE = index(ci + 1)
    do k = jS, jE
      jn = item(k)
      if(jn == cj)then
        do i1 = 1, ndof
        do i2 = 1, ndof
          im = NDOF2*(k-1) + ndof*(i1-1) + i2
          A(im) = val(i1, i2)
        enddo
        enddo
        return
      endif
    enddo
    call monolis_stop_by_matrix_assemble(ci, cj)
  end subroutine monolis_set_block_to_sparse_matrix_main_C

  !# getter
  !> スカラ値を疎行列から取得（メイン関数、実数型）
  subroutine monolis_get_scalar_from_sparse_matrix_main_R(index, item, A, ndof, ci, cj, csub_i, csub_j, val, is_find)
    implicit none
    !> index 配列
    integer(kint), intent(in) :: index(:)
    !> item 配列
    integer(kint), intent(in) :: item(:)
    !> 係数行列
    real(kdouble), intent(in) :: A(:)
    !> ブロック自由度
    integer(kint), intent(in) :: ndof
    !> 行番号
    integer(kint), intent(in) :: ci
    !> 列番号
    integer(kint), intent(in) :: cj
    !> ブロック中の行番号
    integer(kint), intent(in) :: csub_i
    !> ブロック中の列番号
    integer(kint), intent(in) :: csub_j
    !> 設定値
    real(kdouble), intent(out) :: val
    !> 取得判定フラグ
    logical :: is_find
    integer(kint) :: j, jn, im, jS, jE, NDOF2

    val = 0.0d0
    is_find = .false.

    NDOF2 = ndof*ndof
    if(ndof < csub_i) return
    if(ndof < csub_j) return

    jS = index(ci) + 1
    jE = index(ci + 1)
    do j = jS, jE
      jn = item(j)
      if(jn == cj)then
        im = NDOF2*(j-1) + ndof*(csub_i-1) + csub_j
        val = A(im)
        is_find = .true.
        return
      endif
    enddo
  end subroutine monolis_get_scalar_from_sparse_matrix_main_R

  !> スカラ値を疎行列から取得（メイン関数、複素数型）
  subroutine monolis_get_scalar_from_sparse_matrix_main_C(index, item, A, ndof, ci, cj, csub_i, csub_j, val, is_find)
    implicit none
    !> index 配列
    integer(kint), intent(in) :: index(:)
    !> item 配列
    integer(kint), intent(in) :: item(:)
    !> 係数行列
    complex(kdouble), intent(in) :: A(:)
    !> ブロック自由度
    integer(kint), intent(in) :: ndof
    !> 行番号
    integer(kint), intent(in) :: ci
    !> 列番号
    integer(kint), intent(in) :: cj
    !> ブロック中の行番号
    integer(kint), intent(in) :: csub_i
    !> ブロック中の列番号
    integer(kint), intent(in) :: csub_j
    !> 設定値
    complex(kdouble), intent(out) :: val
    !> 取得判定フラグ
    logical :: is_find
    integer(kint) :: j, jn, im, jS, jE, NDOF2

    val = 0.0d0
    is_find = .false.

    NDOF2 = ndof*ndof
    if(ndof < csub_i) return
    if(ndof < csub_j) return

    jS = index(ci) + 1
    jE = index(ci + 1)
    do j = jS, jE
      jn = item(j)
      if(jn == cj)then
        im = NDOF2*(j-1) + ndof*(csub_i-1) + csub_j
        val = A(im)
        is_find = .true.
        return
      endif
    enddo
  end subroutine monolis_get_scalar_from_sparse_matrix_main_C

  !# adder
  !> スカラ値を疎行列に足込（メイン関数、実数型）
  subroutine monolis_add_scalar_to_sparse_matrix_main_R(index, item, A, ndof, ci, cj, csub_i, csub_j, val)
    implicit none
    !> index 配列
    integer(kint), intent(in) :: index(:)
    !> item 配列
    integer(kint), intent(in) :: item(:)
    !> 係数行列
    real(kdouble), intent(inout) :: A(:)
    !> ブロック自由度
    integer(kint), intent(in) :: ndof
    !> 行番号
    integer(kint), intent(in) :: ci
    !> 列番号
    integer(kint), intent(in) :: cj
    !> ブロック中の行番号
    integer(kint), intent(in) :: csub_i
    !> ブロック中の列番号
    integer(kint), intent(in) :: csub_j
    !> 設定値
    real(kdouble), intent(in) :: val
    integer(kint) :: j, jn, im, jS, jE, NDOF2

    if(ndof < csub_i) call monolis_stop_by_submatrix_access(ndof, csub_i)
    if(ndof < csub_j) call monolis_stop_by_submatrix_access(ndof, csub_j)

    NDOF2 = ndof*ndof

    jS = index(ci) + 1
    jE = index(ci + 1)
    do j = jS, jE
      jn = item(j)
      if(jn == cj)then
        im = NDOF2*(j-1) + ndof*(csub_i-1) + csub_j
        A(im) = A(im) + val
        return
      endif
    enddo

    call monolis_stop_by_matrix_assemble(ci, cj)
  end subroutine monolis_add_scalar_to_sparse_matrix_main_R

  !> スカラ値を疎行列に足込（メイン関数、複素数型）
  subroutine monolis_add_scalar_to_sparse_matrix_main_C(index, item, A, ndof, ci, cj, csub_i, csub_j, val)
    implicit none
    !> index 配列
    integer(kint), intent(in) :: index(:)
    !> item 配列
    integer(kint), intent(in) :: item(:)
    !> 係数行列
    complex(kdouble), intent(inout) :: A(:)
    !> ブロック自由度
    integer(kint), intent(in) :: ndof
    !> 行番号
    integer(kint), intent(in) :: ci
    !> 列番号
    integer(kint), intent(in) :: cj
    !> ブロック中の行番号
    integer(kint), intent(in) :: csub_i
    !> ブロック中の列番号
    integer(kint), intent(in) :: csub_j
    !> 設定値
    complex(kdouble), intent(in) :: val
    integer(kint) :: j, jn, im, jS, jE, NDOF2

    if(ndof < csub_i) call monolis_stop_by_submatrix_access(ndof, csub_i)
    if(ndof < csub_j) call monolis_stop_by_submatrix_access(ndof, csub_j)

    NDOF2 = ndof*ndof

    jS = index(ci) + 1
    jE = index(ci + 1)
    do j = jS, jE
      jn = item(j)
      if(jn == cj)then
        im = NDOF2*(j-1) + ndof*(csub_i-1) + csub_j
        A(im) = A(im) + val
        return
      endif
    enddo

    call monolis_stop_by_matrix_assemble(ci, cj)
  end subroutine monolis_add_scalar_to_sparse_matrix_main_C

  !> 行列値を疎行列に足込（メイン関数、実数型）
  subroutine monolis_add_matrix_to_sparse_matrix_main_R(index, item, A, n1, n2, ndof, e1, e2, val)
    implicit none
    !> index 配列
    integer(kint), intent(in) :: index(:)
    !> item 配列
    integer(kint), intent(in) :: item(:)
    !> 係数行列
    real(kdouble), intent(inout) :: A(:)
    !> 入力行列の行サイズ
    integer(kint), intent(in) :: n1
    !> 入力行列の列サイズ
    integer(kint), intent(in) :: n2
    !> ブロック自由度
    integer(kint), intent(in) :: ndof
    !> 行番号リスト
    integer(kint), intent(in) :: e1(n1)
    !> 列番号リスト
    integer(kint), intent(in) :: e2(n2)
    !> 設定値
    real(kdouble), intent(in) :: val(:,:)
    integer(kint) :: e1t(n1), e2t(n2)
    integer(kint) :: i, j, k, in, jn, im, jS, jE, i2, j2, i1, j1, NDOF2
    integer(kint) :: eperm1(n1), eperm2(n2)
    real(kdouble) :: temp(n1*ndof,n2*ndof)

    NDOF2 = ndof*ndof

    e1t = e1
    call monolis_get_sequence_array_I(eperm1, n1, 1, 1)
    call monolis_qsort_I_2d(e1t, eperm1, 1, n1)

    e2t = e2
    call monolis_get_sequence_array_I(eperm2, n2, 1, 1)
    call monolis_qsort_I_2d(e2t, eperm2, 1, n2)

    temp = 0.0d0
    do i = 1, n2
      i1 = eperm2(i)
      do j = 1, n1
        j1 = eperm1(j)
        do i2 = 1, ndof
          do j2 = 1, ndof
            temp(ndof*(i-1)+i2, ndof*(j-1)+j2) = val(ndof*(i1-1)+i2, ndof*(j1-1)+j2)
          enddo
        enddo
      enddo
    enddo

    do i = 1, n1
      in = e1t(i)
      jS = index(in) + 1
      jE = index(in + 1)
      aa:do j = 1, n2
        do k = jS, jE
          jn = item(k)
          if(jn == e2t(j))then
            do i1 = 1, ndof
            do i2 = 1, ndof
              im = NDOF2*(k-1) + ndof*(i1-1) + i2
              A(im) = A(im) + temp(ndof*(i-1)+i1, ndof*(j-1)+i2)
            enddo
            enddo
            jS = k + 1
            cycle aa
          endif
        enddo
      call monolis_stop_by_matrix_assemble(e1t(i), e2t(j))
      enddo aa
    enddo
  end subroutine monolis_add_matrix_to_sparse_matrix_main_R

  !> 行列値を疎行列に足込（メイン関数、複素数型）
  subroutine monolis_add_matrix_to_sparse_matrix_main_C(index, item, A, n1, n2, ndof, e1, e2, val)
    implicit none
    !> index 配列
    integer(kint), intent(in) :: index(:)
    !> item 配列
    integer(kint), intent(in) :: item(:)
    !> 係数行列
    complex(kdouble), intent(inout) :: A(:)
    !> 入力行列の行サイズ
    integer(kint), intent(in) :: n1
    !> 入力行列の列サイズ
    integer(kint), intent(in) :: n2
    !> ブロック自由度
    integer(kint), intent(in) :: ndof
    !> 行番号リスト
    integer(kint), intent(in) :: e1(n1)
    !> 列番号リスト
    integer(kint), intent(in) :: e2(n2)
    !> 設定値
    complex(kdouble), intent(in) :: val(:,:)
    integer(kint) :: e1t(n1), e2t(n2)
    integer(kint) :: i, j, k, in, jn, im, jS, jE, i2, j2, i1, j1, NDOF2
    integer(kint) :: eperm1(n1), eperm2(n2)
    complex(kdouble) :: temp(n1*ndof,n2*ndof)

    NDOF2 = ndof*ndof

    e1t = e1
    call monolis_get_sequence_array_I(eperm1, n1, 1, 1)
    call monolis_qsort_I_2d(e1t, eperm1, 1, n1)

    e2t = e2
    call monolis_get_sequence_array_I(eperm2, n2, 1, 1)
    call monolis_qsort_I_2d(e2t, eperm2, 1, n2)

    temp = 0.0d0
    do i = 1, n2
      i1 = eperm2(i)
      do j = 1, n1
        j1 = eperm1(j)
        do i2 = 1, ndof
          do j2 = 1, ndof
            temp(ndof*(i-1)+i2, ndof*(j-1)+j2) = val(ndof*(i1-1)+i2, ndof*(j1-1)+j2)
          enddo
        enddo
      enddo
    enddo

    do i = 1, n1
      in = e1t(i)
      jS = index(in) + 1
      jE = index(in + 1)
      aa:do j = 1, n2
        do k = jS, jE
          jn = item(k)
          if(jn == e2t(j))then
            do i1 = 1, ndof
            do i2 = 1, ndof
              im = NDOF2*(k-1) + ndof*(i1-1) + i2
              A(im) = A(im) + temp(ndof*(i-1)+i1, ndof*(j-1)+i2)
            enddo
            enddo
            jS = k + 1
            cycle aa
          endif
        enddo
      call monolis_stop_by_matrix_assemble(e1t(i), e2t(j))
      enddo aa
    enddo
  end subroutine monolis_add_matrix_to_sparse_matrix_main_C

  !> 境界条件処理（実数型、メイン関数）
  subroutine monolis_set_Dirichlet_bc_main_R(index, item, A, B, indexR, itemR, permA, &
    & ndof, node_id, ndof_bc, val)
    implicit none
    !> index 配列
    integer(kint), intent(in) :: index(:)
    !> item 配列
    integer(kint), intent(in) :: item(:)
    !> 係数行列
    real(kdouble), intent(inout) :: A(:)
    !> 右辺ベクトル
    real(kdouble), intent(inout) :: B(:)
    !> index 配列（CSC 形式）
    integer(kint), intent(in) :: indexR(:)
    !> item 配列（CSC 形式）
    integer(kint), intent(in) :: itemR(:)
    !> 行列成分の CSC 形式との置換ベクトル
    integer(kint), intent(in) :: permA(:)
    !> ブロック自由度
    integer(kint), intent(in) :: ndof
    !> 自由度番号
    integer(kint), intent(in) :: node_id
    !> ブロック番号
    integer(kint), intent(in) :: ndof_bc
    !> 境界条件の設定値
    real(kdouble), intent(in) :: val
    integer(kint) :: j, k, jn, kn, jS, jE, NDOF2
    logical :: is_add

    if(ndof < ndof_bc) call monolis_stop_by_submatrix_access(ndof, ndof_bc)

    is_add = .false.
    NDOF2 = ndof*ndof

    jS = indexR(node_id) + 1
    jE = indexR(node_id + 1)
    do j = jS, jE
      jn = itemR(j)
      kn = permA(j)
      do k = 1, ndof
        B(ndof*(jn-1)+k) = B(ndof*(jn-1)+k) - val*A(NDOF2*(kn-1) + ndof*(k-1) + ndof_bc)
        A(NDOF2*(kn-1) + ndof*(k-1) + ndof_bc) = 0.0d0
      enddo
    enddo

    jS = index(node_id) + 1
    jE = index(node_id + 1)
    do j = jS, jE
      do k = 1, ndof
        A(NDOF2*(j-1) + ndof*(ndof_bc-1) + k) = 0.0d0
      enddo

      jn = item(j)
      if(jn == node_id)then
        A(NDOF2*(j-1) + (ndof+1)*(ndof_bc-1) + 1) = 1.0d0
        is_add = .true.
      endif
    enddo

    if(.not. is_add)then
      call monolis_stop_by_set_DBC(node_id)
    endif

    B(ndof*node_id - ndof + ndof_bc) = val
  end subroutine monolis_set_Dirichlet_bc_main_R

  !> 境界条件処理（複素数型、メイン関数）
  subroutine monolis_set_Dirichlet_bc_main_C(index, item, A, B, indexR, itemR, permA, &
    & ndof, node_id, ndof_bc, val)
    implicit none
    !> index 配列
    integer(kint), intent(in) :: index(:)
    !> item 配列
    integer(kint), intent(in) :: item(:)
    !> 係数行列
    complex(kdouble), intent(inout) :: A(:)
    !> 右辺ベクトル
    complex(kdouble), intent(inout) :: B(:)
    !> index 配列（CSC 形式）
    integer(kint), intent(in) :: indexR(:)
    !> item 配列（CSC 形式）
    integer(kint), intent(in) :: itemR(:)
    !> 行列成分の CSC 形式との置換ベクトル
    integer(kint), intent(in) :: permA(:)
    !> ブロック自由度
    integer(kint), intent(in) :: ndof
    !> 自由度番号
    integer(kint), intent(in) :: node_id
    !> ブロック番号
    integer(kint), intent(in) :: ndof_bc
    !> 境界条件の設定値
    complex(kdouble), intent(in) :: val
    integer(kint) :: j, k, jn, kn, jS, jE, NDOF2
    logical :: is_add

    if(ndof < ndof_bc) call monolis_stop_by_submatrix_access(ndof, ndof_bc)

    is_add = .false.
    NDOF2 = ndof*ndof

    jS = indexR(node_id) + 1
    jE = indexR(node_id + 1)
    do j = jS, jE
      jn = itemR(j)
      kn = permA(j)
      do k = 1, ndof
        B(ndof*(jn-1)+k) = B(ndof*(jn-1)+k) - val*A(NDOF2*(kn-1) + ndof*(k-1) + ndof_bc)
        A(NDOF2*(kn-1) + ndof*(k-1) + ndof_bc) = 0.0d0
      enddo
    enddo

    jS = index(node_id) + 1
    jE = index(node_id + 1)
    do j = jS, jE
      do k = 1, ndof
        A(NDOF2*(j-1) + ndof*(ndof_bc-1) + k) = 0.0d0
      enddo

      jn = item(j)
      if(jn == node_id)then
        A(NDOF2*(j-1) + (ndof+1)*(ndof_bc-1) + 1) = 1.0d0
        is_add = .true.
      endif
    enddo

    if(.not. is_add)then
      call monolis_stop_by_set_DBC(node_id)
    endif

    B(ndof*node_id - ndof + ndof_bc) = val
  end subroutine monolis_set_Dirichlet_bc_main_C

  !> 非零要素がメモリ確保されていない場合にエラーストップ
  subroutine monolis_stop_by_matrix_assemble(ci, cj)
    !> 行番号
    integer(kint), intent(in) :: ci
    !> 列番号
    integer(kint), intent(in) :: cj

    call monolis_std_error_string("")
    write(*,"(a,i0,a,i0,a)") "The non-zero element at (", ci, ", ", cj, &
      & ") is not allocated. The value is not accessible."
    call monolis_std_error_stop()
  end subroutine monolis_stop_by_matrix_assemble

  !> ブロックサイズより大きい部分へのアクサス時にエラーストップ
  subroutine monolis_stop_by_submatrix_access(ndof, sub_dof)
    !> 行番号
    integer(kint), intent(in) :: ndof
    !> ブロック番号
    integer(kint), intent(in) :: sub_dof

    call monolis_std_error_string("")
    write(*,"(a)")   "input value is greater than the DoF of submatrix."
    write(*,"(a,i8)")"       the DoF of submatrix: ", ndof
    write(*,"(a,i8)")"       value:                ", sub_dof
    call monolis_std_error_stop()
  end subroutine monolis_stop_by_submatrix_access

  !> 対角成分が見つからない場合にエラーストップ
  subroutine monolis_stop_by_set_DBC(node_id)
    !> 行番号
    integer(kint), intent(in) :: node_id

    call monolis_std_error_string("")
    write(*,"(a)")   "The diagonal component is not found."
    write(*,"(a,i8)")"       tdiagonal number     : ", node_id
    call monolis_std_error_stop()
  end subroutine monolis_stop_by_set_DBC

  !> 対角成分が零の場合にエラーストップ
  subroutine monolis_stop_by_set_zero_diag_component(node_id, ndof)
    !> 行番号
    integer(kint), intent(in) :: node_id
    !> ブロック自由度
    integer(kint), intent(in) :: ndof

    call monolis_std_error_string("")
    write(*,"(a)")   "The diagonal component is zero."
    write(*,"(a,i8)")"       diagonal number     : ", node_id
    write(*,"(a,i8)")"       the DoF of submatrix: ", ndof
    call monolis_std_error_stop()
  end subroutine monolis_stop_by_set_zero_diag_component

  !> 行列の最大値の取得（実数型、メイン関数）
  subroutine monolis_get_max_matrix_component_main_R(monoMAT, monoCOM, max_val)
    implicit none
    !> monolis MAT 構造体
    type(monolis_mat) :: monoMAT
    !> monolis COM 構造体
    type(monolis_com) :: monoCOM
    !> 最大値
    real(kdouble) :: max_val

    max_val = maxval(monoMAT%R%A)
    call monolis_allreduce_R1(max_val, monolis_mpi_max, monoCOM%comm)
  end subroutine monolis_get_max_matrix_component_main_R

  !> 行列の対角成分の 0 チェック（実数型、メイン関数）
  subroutine monolis_check_diagonal_zero_component_main_R(monoPRM, monoMAT)
    implicit none
    !> monolis PRM 構造体
    type(monolis_prm) :: monoPRM
    !> monolis MAT 構造体
    type(monolis_mat) :: monoMAT
    integer(kint) :: i, j, k, jS, jE, in, kn, N, NDOF, NDOF2

    if(monoPRM%Iarray(monolis_prm_I_is_check_diag) == monolis_I_false) return

    N =  monoMAT%N
    NDOF  = monoMAT%NDOF
    NDOF2 = NDOF*NDOF

    do i = 1, N
      jS = monoMAT%CSR%index(i) + 1
      jE = monoMAT%CSR%index(i + 1)
      do j = jS, jE
        in = monoMAT%CSR%item(j)
        if(i == in)then
          do k = 1, NDOF
            kn = NDOF2*(j-1) + (NDOF+1)*(k-1) + 1
            if(monoMAT%R%A(kn) == 0.0d0)then
              call monolis_stop_by_set_zero_diag_component(i, k)
              stop
            endif
          enddo
        endif
      enddo
    enddo
  end subroutine monolis_check_diagonal_zero_component_main_R
end module mod_monolis_spmat_handler_util
