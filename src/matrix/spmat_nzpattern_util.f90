!> 非零構造決定モジュール（メイン関数群）
module mod_monolis_spmat_nonzero_pattern_util
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc

  implicit none

contains

  !> @ingroup dev_matrix
  !> 節点グラフから疎行列パターンを決定（メイン関数）
  subroutine monolis_get_nonzero_pattern_by_nodal_graph_main(MAT, n_node, ndof, index, item)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_mat), intent(inout) :: MAT
    !> [in] 節点数
    integer(kint), intent(in) :: n_node
    !> [in] 計算点が持つ自由度
    integer(kint), intent(in) :: ndof
    !> [in] index 配列
    integer(kint), intent(in) :: index(:)
    !> [in] item 配列
    integer(kint), intent(in) :: item(:)
    integer(kint) :: i, j, nz, jS, jE

    MAT%N = n_node
    MAT%NP = n_node
    MAT%NDOF = ndof

    call monolis_palloc_I_1d(MAT%CSR%index, n_node + 1)

    do i = 1, n_node
      MAT%CSR%index(i + 1) = index(i + 1) + i
    enddo

    nz = MAT%CSR%index(n_node + 1)

    call monolis_palloc_I_1d(MAT%CSR%item, nz)

    do i = 1, n_node
      jS = MAT%CSR%index(i) + 1
      jE = MAT%CSR%index(i + 1)
      MAT%CSR%item(jS) = i
      do j = jS + 1, jE
        MAT%CSR%item(j) = item(j - i)
      enddo
      call monolis_qsort_I_1d(MAT%CSR%item(jS:jE), 1, jE - jS + 1)
    enddo

    call monolis_palloc_I_1d(MAT%CSC%index, n_node + 1)
    call monolis_palloc_I_1d(MAT%CSC%item, nz)
    call monolis_palloc_I_1d(MAT%CSC%perm, nz)

    call monolis_get_CSC_format(MAT%N, MAT%N, nz, &
      & MAT%CSR%index, MAT%CSR%item, &
      & MAT%CSC%index, MAT%CSC%item, MAT%CSC%perm)
  end subroutine monolis_get_nonzero_pattern_by_nodal_graph_main

  !> @ingroup dev_matrix
  !> 節点グラフから疎行列パターンを決定（任意節点自由度、メイン関数）
  subroutine monolis_get_nonzero_pattern_by_nodal_graph_with_arbit_main( &
    MAT, n_node, n_dof_list, index, item)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_mat), intent(inout) :: MAT
    !> [in] 節点数
    integer(kint), intent(in) :: n_node
    !> [in] 自由度リスト
    integer(kint), intent(in) :: n_dof_list(:)
    !> [in] index 配列
    integer(kint), intent(in) :: index(:)
    !> [in] item 配列
    integer(kint), intent(in) :: item(:)
    integer(kint) :: i, in, j, jn, nz, jS, jE, k, kE, kS, l
    integer(kint) :: total_dof, ncol, nrow
    integer(kint), allocatable :: n_dof_index(:)

    total_dof = 0
    do i = 1, n_node
      total_dof = total_dof + n_dof_list(i)
    enddo

    call monolis_alloc_I_1d(n_dof_index, n_node)

    call monolis_get_n_dof_index(n_node, n_dof_list, n_dof_index)

    MAT%N = total_dof
    MAT%NP = total_dof
    MAT%NDOF = 1

    call monolis_palloc_I_1d(MAT%CSR%index, total_dof + 1)

    !# count nz
    nz = 0
    do i = 1, n_node
      jS = index(i) + 1
      jE = index(i + 1)
      nz = nz + n_dof_list(i)*n_dof_list(i)
      do j = jS, jE
        jn = item(j)
        nz = nz +  n_dof_list(i)*n_dof_list(jn)
      enddo
    enddo

    call monolis_palloc_I_1d(MAT%CSR%item, nz)

    !# construct index and item
    in = 0
    ncol = 0
    do i = 1, n_node
      jS = index(i) + 1
      jE = index(i+1)
      do k = 1, n_dof_list(i)
        ncol = ncol + 1
        nrow = n_dof_list(i)
        do j = 1, n_dof_list(i)
          in = in + 1
          MAT%CSR%item(in) = n_dof_index(i) + j
        enddo
        do j = jS, jE
          jn = item(j)
          nrow = nrow + n_dof_list(jn)
          do l = 1, n_dof_list(jn)
            in = in + 1
            MAT%CSR%item(in) = n_dof_index(jn) + l
          enddo
        enddo
        MAT%CSR%index(ncol + 1) = MAT%CSR%index(ncol) + nrow

        kS = MAT%CSR%index(ncol) + 1
        kE = MAT%CSR%index(ncol + 1)
        call monolis_qsort_I_1d(MAT%CSR%item(kS:kE), 1, kE-kS+1)
      enddo
    enddo

    call monolis_palloc_I_1d(MAT%CSC%index, total_dof + 1)
    call monolis_palloc_I_1d(MAT%CSC%item, nz)
    call monolis_palloc_I_1d(MAT%CSC%perm, nz)

    call monolis_get_CSC_format(MAT%N, MAT%N, nz, &
      & MAT%CSR%index, MAT%CSR%item, &
      & MAT%CSC%index, MAT%CSC%item, MAT%CSC%perm)
  end subroutine monolis_get_nonzero_pattern_by_nodal_graph_with_arbit_main

  !> @ingroup dev_matrix
  !> 節点数 index の作成
  subroutine monolis_get_n_dof_index(n_node, n_dof_list, n_dof_index)
    implicit none
    !> [in] 節点数
    integer(kint), intent(in) :: n_node
    !> [in] 節点自由度リスト
    integer(kint), intent(in) :: n_dof_list(:)
    !> [in,out] 節点自由度 index
    integer(kint), intent(inout) :: n_dof_index(:)
    integer(kint) :: i

    do i = 1, n_node - 1
      n_dof_index(i+1) = n_dof_index(i) + n_dof_list(i)
    enddo
  end subroutine monolis_get_n_dof_index

  !> @ingroup dev_matrix
  !> 疎行列の行列成分のメモリ確保（実数型）
  subroutine monolis_alloc_nonzero_pattern_mat_val_R(MAT)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_mat), intent(inout) :: MAT
    integer(kint) :: NP, NDOF, NZ

    NDOF = MAT%NDOF
    NP = MAT%NP
    NZ = MAT%CSR%index(NP + 1)

    call monolis_palloc_R_1d(MAT%R%A, NZ*NDOF*NDOF)
    call monolis_palloc_R_1d(MAT%R%B, NP*NDOF)
    call monolis_palloc_R_1d(MAT%R%X, NP*NDOF)
  end subroutine monolis_alloc_nonzero_pattern_mat_val_R

  !> @ingroup dev_matrix
  !> 疎行列の行列成分のメモリ確保（複素数型）
  subroutine monolis_alloc_nonzero_pattern_mat_val_C(MAT)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_mat), intent(inout) :: MAT
    integer(kint) :: NDOF, NP, NZ

    NDOF = MAT%NDOF
    NP = MAT%NP
    NZ = MAT%CSR%index(NP + 1)

    call monolis_palloc_C_1d(MAT%C%A, NZ*NDOF*NDOF)
    call monolis_palloc_C_1d(MAT%C%B, NP*NDOF)
    call monolis_palloc_C_1d(MAT%C%X, NP*NDOF)
  end subroutine monolis_alloc_nonzero_pattern_mat_val_C

  !> @ingroup dev_matrix
  !> CSR 形式から CSC 形式のデータを生成
  subroutine monolis_get_CSC_format(NC, NR, NZ, index, item, indexR, itemR, permR)
    implicit none
    !> [in] 行数
    integer(kint), intent(in) :: NC
    !> [in] 列数
    integer(kint), intent(in) :: NR
    !> [in] 非零要素数
    integer(kint), intent(in) :: NZ
    !> [in] index 配列（CSR 形式）
    integer(kint), intent(in) :: index(:)
    !> [in] item 配列（CSR 形式）
    integer(kint), intent(in) :: item(:)
    !> [out] index 配列（CSC 形式）
    integer(kint), intent(out) :: indexR(:)
    !> [in,out] index 配列（CSC 形式）
    integer(kint), intent(inout) :: itemR(:)
    !> [in,out] index 配列（CSC 形式）
    integer(kint), intent(inout) :: permR(:)
    integer(kint), allocatable :: temp(:)
    integer(kint) :: i, j, in, jS, jE, m, p

    call monolis_alloc_I_1d(temp, NR)

    do i = 1, NC
      jS = index(i) + 1
      jE = index(i + 1)
      do j = jS, jE
        in = item(j)
        temp(in) = temp(in) + 1
      enddo
    enddo

    indexR = 0
    do i = 1, NR
      indexR(i + 1) = indexR(i) + temp(i)
    enddo

    temp = 0
    do i = 1, NC
      jS = index(i) + 1
      jE = index(i + 1)
      do j = jS, jE
        in = item(j)
        m = indexR(in)
        temp(in) = temp(in) + 1
        p = temp(in)
        itemR(m + p) = i
        permR(m + p) = j
      enddo
    enddo

    do i = 1, NR
      jS = indexR(i) + 1
      jE = indexR(i + 1)
      call monolis_qsort_I_1d(itemR(jS:jE), 1, jE - jS + 1)
    enddo
  end subroutine monolis_get_CSC_format
end module mod_monolis_spmat_nonzero_pattern_util
