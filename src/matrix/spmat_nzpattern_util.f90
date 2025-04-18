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
  !> 疎行列構造体に自由度リストを登録
  subroutine monolis_set_n_dof_index(MAT, n_dof_list)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_mat), intent(inout) :: MAT
    !> [in] 節点自由度リスト
    integer(kint), intent(in) :: n_dof_list(:)
    integer(kint) :: i, j, jS, jE, in, NZ

    NZ = MAT%CSR%index(MAT%NP + 1)

    call monolis_palloc_I_1d(MAT%n_dof_list, MAT%NP)
    call monolis_palloc_I_1d(MAT%n_dof_index, MAT%NP + 1)
    call monolis_palloc_I_1d(MAT%n_dof_index2, NZ)

    MAT%n_dof_list = n_dof_list

    do i = 1, MAT%NP
      MAT%n_dof_index(i + 1) = MAT%n_dof_index(i) + n_dof_list(i)
    enddo

    do i = 1, MAT%NP
      jS = MAT%CSR%index(i) + 1
      jE = MAT%CSR%index(i + 1)
      do j = jS, jE
        in = MAT%CSR%item(j)
        MAT%n_dof_index2(j + 1) = MAT%n_dof_index2(j) + n_dof_list(i)*n_dof_list(in)
      enddo
    enddo
  end subroutine monolis_set_n_dof_index

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
  !> 疎行列の行列成分のメモリ確保（実数型）
  subroutine monolis_alloc_nonzero_pattern_mat_val_V_R(MAT)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_mat), intent(inout) :: MAT
    integer(kint) :: N, NZ, in

    N  = MAT%n_dof_index (MAT%NP + 1)
    in = MAT%CSR%index(MAT%NP + 1)
    NZ = MAT%n_dof_index2(in + 1)

    call monolis_palloc_R_1d(MAT%R%A, NZ)
    call monolis_palloc_R_1d(MAT%R%B, N)
    call monolis_palloc_R_1d(MAT%R%X, N)
  end subroutine monolis_alloc_nonzero_pattern_mat_val_V_R

  !> @ingroup dev_matrix
  !> 疎行列の行列成分のメモリ確保（複素数型）
  subroutine monolis_alloc_nonzero_pattern_mat_val_V_C(MAT)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_mat), intent(inout) :: MAT
    integer(kint) :: N, NZ, in

    N  = MAT%n_dof_index (MAT%NP + 1)
    in = MAT%CSR%index(MAT%NP + 1)
    NZ = MAT%n_dof_index2(in + 1)

    call monolis_palloc_C_1d(MAT%C%A, NZ)
    call monolis_palloc_C_1d(MAT%C%B, N)
    call monolis_palloc_C_1d(MAT%C%X, N)
  end subroutine monolis_alloc_nonzero_pattern_mat_val_V_C

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
