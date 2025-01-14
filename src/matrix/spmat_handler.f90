!> 疎行列操作関数群
module mod_monolis_spmat_handler
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc
  use mod_monolis_spmat_handler_util
  use mod_monolis_spmat_nonzero_pattern_util
  use mod_gedatsu

  implicit none

contains

  !# setter
  !> @ingroup matrix
  !> スカラ値を疎行列に設定（実数型）
  subroutine monolis_set_scalar_to_sparse_matrix_R(monolis, i, j, sub_i, sub_j, val)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_structure), intent(inout) :: monolis
    !> [in] 行番号
    integer(kint), intent(in) :: i
    !> [in] 列番号
    integer(kint), intent(in) :: j
    !> [in] ブロック中の行番号
    integer(kint), intent(in) :: sub_i
    !> [in] ブロック中の列番号
    integer(kint), intent(in) :: sub_j
    !> [in] 行列値
    real(kdouble), intent(in) :: val

    call monolis_set_scalar_to_sparse_matrix_main_R(monolis%MAT%CSR%index, monolis%MAT%CSR%item, &
      & monolis%MAT%R%A, monolis%MAT%ndof, i, j, sub_i, sub_j, val)
  end subroutine monolis_set_scalar_to_sparse_matrix_R

  !> @ingroup matrix
  !> スカラ値を疎行列に設定（複素数型）
  subroutine monolis_set_scalar_to_sparse_matrix_C(monolis, i, j, sub_i, sub_j, val)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_structure), intent(inout) :: monolis
    !> [in] 行番号
    integer(kint), intent(in) :: i
    !> [in] 列番号
    integer(kint), intent(in) :: j
    !> [in] ブロック中の行番号
    integer(kint), intent(in) :: sub_i
    !> [in] ブロック中の列番号
    integer(kint), intent(in) :: sub_j
    !> [in] 行列値
    complex(kdouble), intent(in) :: val

    call monolis_set_scalar_to_sparse_matrix_main_C(monolis%MAT%CSR%index, monolis%MAT%CSR%item, &
      & monolis%MAT%C%A, monolis%MAT%ndof, i, j, sub_i, sub_j, val)
  end subroutine monolis_set_scalar_to_sparse_matrix_C

  !> @ingroup matrix
  !> ブロック行列を疎行列に設定（実数型）
  subroutine monolis_set_block_to_sparse_matrix_R(monolis, i, j, val)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_structure), intent(inout) :: monolis
    !> [in] 行番号
    integer(kint), intent(in) :: i
    !> [in] 列番号
    integer(kint), intent(in) :: j
    !> [in] 行列値（サイズ：[ndof, ndof]）
    real(kdouble), intent(in) :: val(:,:)

    call monolis_set_block_to_sparse_matrix_main_R(monolis%MAT%CSR%index, monolis%MAT%CSR%item, &
      & monolis%MAT%R%A, monolis%MAT%ndof, i, j, val)
  end subroutine monolis_set_block_to_sparse_matrix_R

  !> @ingroup matrix
  !> ブロック行列を疎行列に設定（複素数型）
  subroutine monolis_set_block_to_sparse_matrix_C(monolis, i, j, val)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_structure), intent(inout) :: monolis
    !> [in] 行番号
    integer(kint), intent(in) :: i
    !> [in] 列番号
    integer(kint), intent(in) :: j
    !> [in] 行列値（サイズ：[ndof, ndof]）
    complex(kdouble), intent(in) :: val(:,:)

    call monolis_set_block_to_sparse_matrix_main_C(monolis%MAT%CSR%index, monolis%MAT%CSR%item, &
      & monolis%MAT%C%A, monolis%MAT%ndof, i, j, val)
  end subroutine monolis_set_block_to_sparse_matrix_C

  !# getter
  !> @ingroup matrix
  !> スカラ値を疎行列から取得（実数型）
  subroutine monolis_get_scalar_from_sparse_matrix_R(monolis, i, j, sub_i, sub_j, val, is_find)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_structure), intent(inout) :: monolis
    !> 行番号
    integer(kint), intent(in) :: i
    !> 列番号
    integer(kint), intent(in) :: j
    !> ブロック中の行番号
    integer(kint), intent(in) :: sub_i
    !> [in] ブロック中の列番号
    integer(kint), intent(in) :: sub_j
    !> [out] 行列値
    real(kdouble), intent(out) :: val
    !> [out] 取得判定フラグ
    logical, intent(out) :: is_find

    call monolis_get_scalar_from_sparse_matrix_main_R(monolis%MAT%CSR%index, monolis%MAT%CSR%item, &
      & monolis%MAT%R%A, monolis%MAT%ndof, i, j, sub_i, sub_j, val, is_find)
  end subroutine monolis_get_scalar_from_sparse_matrix_R

  !> @ingroup matrix
  !> スカラ値を疎行列から取得（複素数型）
  subroutine monolis_get_scalar_from_sparse_matrix_C(monolis, i, j, sub_i, sub_j, val, is_find)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_structure), intent(inout) :: monolis
    !> [in] 行番号
    integer(kint), intent(in) :: i
    !> [in] 列番号
    integer(kint), intent(in) :: j
    !> [in] ブロック中の行番号
    integer(kint), intent(in) :: sub_i
    !> [in] ブロック中の列番号
    integer(kint), intent(in) :: sub_j
    !> [out] 行列値
    complex(kdouble), intent(out) :: val
    !> [out] 取得判定フラグ
    logical, intent(out) :: is_find

    call monolis_get_scalar_from_sparse_matrix_main_C(monolis%MAT%CSR%index, monolis%MAT%CSR%item, &
      & monolis%MAT%C%A, monolis%MAT%ndof, i, j, sub_i, sub_j, val, is_find)
  end subroutine monolis_get_scalar_from_sparse_matrix_C

  !# adder
  !> @ingroup matrix
  !> スカラ値を疎行列に足込（実数型）
  subroutine monolis_add_scalar_to_sparse_matrix_R(monolis, i, j, sub_i, sub_j, val)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_structure), intent(inout) :: monolis
    !> [in] 行番号
    integer(kint), intent(in) :: i
    !> [in] 列番号
    integer(kint), intent(in) :: j
    !> [in] ブロック中の行番号
    integer(kint), intent(in) :: sub_i
    !> [in] ブロック中の列番号
    integer(kint), intent(in) :: sub_j
    !> [in] 行列値
    real(kdouble), intent(in) :: val

    call monolis_add_scalar_to_sparse_matrix_main_R(monolis%MAT%CSR%index, monolis%MAT%CSR%item, &
      & monolis%MAT%R%A, monolis%MAT%ndof, i, j, sub_i, sub_j, val)
  end subroutine monolis_add_scalar_to_sparse_matrix_R

  !> @ingroup matrix
  !> スカラ値を疎行列に足込（複素数型）
  subroutine monolis_add_scalar_to_sparse_matrix_C(monolis, i, j, sub_i, sub_j, val)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_structure), intent(inout) :: monolis
    !> [in] 行番号
    integer(kint), intent(in) :: i
    !> [in] 列番号
    integer(kint), intent(in) :: j
    !> [in] ブロック中の行番号
    integer(kint), intent(in) :: sub_i
    !> [in] ブロック中の列番号
    integer(kint), intent(in) :: sub_j
    !> [in] 行列値
    complex(kdouble), intent(in) :: val

    call monolis_add_scalar_to_sparse_matrix_main_C(monolis%MAT%CSR%index, monolis%MAT%CSR%item, &
      & monolis%MAT%C%A, monolis%MAT%ndof, i, j, sub_i, sub_j, val)
  end subroutine monolis_add_scalar_to_sparse_matrix_C

  !> @ingroup matrix
  !> 行列を疎行列に足込（実数型）
  subroutine monolis_add_matrix_to_sparse_matrix_R(monolis, n_base, connectivity, mat)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_structure), intent(inout) :: monolis
    !> [in] 節点数
    integer(kint), intent(in) :: n_base
    !> [in] 要素コネクティビティ
    integer(kint), intent(in) :: connectivity(n_base)
    !> [in] 行列値
    real(kdouble), intent(in) :: mat(:,:)

    call monolis_add_matrix_to_sparse_matrix_main_R(monolis%MAT%CSR%index, monolis%MAT%CSR%item, &
      & monolis%MAT%R%A, n_base, n_base, monolis%MAT%ndof, connectivity, connectivity, mat)
  end subroutine monolis_add_matrix_to_sparse_matrix_R

  !> @ingroup matrix
  !> 行列を疎行列に足込（複素数型）
  subroutine monolis_add_matrix_to_sparse_matrix_C(monolis, n_base, connectivity, mat)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_structure), intent(inout) :: monolis
    !> [in] 節点数
    integer(kint), intent(in) :: n_base
    !> [in] 要素コネクティビティ
    integer(kint), intent(in) :: connectivity(n_base)
    !> [in] 行列値
    complex(kdouble), intent(in) :: mat(:,:)

    call monolis_add_matrix_to_sparse_matrix_main_C(monolis%MAT%CSR%index, monolis%MAT%CSR%item, &
      & monolis%MAT%C%A, n_base, n_base, monolis%MAT%ndof, connectivity, connectivity, mat)
  end subroutine monolis_add_matrix_to_sparse_matrix_C

  !> @ingroup matrix
  !> 行列を疎行列の非対角部分に足込（実数型）
  subroutine monolis_add_matrix_to_sparse_matrix_offdiag_R(monolis, n_base_i, n_base_j, &
    & connectivity_i, connectivity_j, mat)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_structure), intent(inout) :: monolis
    !> [in] 節点数
    integer(kint), intent(in) :: n_base_i
    !> [in] 節点数
    integer(kint), intent(in) :: n_base_j
    !> [in] 要素コネクティビティ
    integer(kint), intent(in) :: connectivity_i(n_base_i)
    !> [in] 要素コネクティビティ
    integer(kint), intent(in) :: connectivity_j(n_base_j)
    !> [in] 行列値
    real(kdouble), intent(in) :: mat(:,:)

    call monolis_add_matrix_to_sparse_matrix_main_R(monolis%MAT%CSR%index, monolis%MAT%CSR%item, &
      & monolis%MAT%R%A, n_base_i, n_base_j, monolis%MAT%ndof, connectivity_i, connectivity_j, mat)
  end subroutine monolis_add_matrix_to_sparse_matrix_offdiag_R

  !> @ingroup matrix
  !> 行列を疎行列の非対角部分に足込（複素数型）
  subroutine monolis_add_matrix_to_sparse_matrix_offdiag_C(monolis, n_base_i, n_base_j, &
    & connectivity_i, connectivity_j, mat)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_structure), intent(inout) :: monolis
    !> [in] 節点数
    integer(kint), intent(in) :: n_base_i
    !> [in] 節点数
    integer(kint), intent(in) :: n_base_j
    !> [in] 要素コネクティビティ
    integer(kint), intent(in) :: connectivity_i(n_base_i)
    !> [in] 要素コネクティビティ
    integer(kint), intent(in) :: connectivity_j(n_base_j)
    !> [in] 行列値
    complex(kdouble), intent(in) :: mat(:,:)

    call monolis_add_matrix_to_sparse_matrix_main_C(monolis%MAT%CSR%index, monolis%MAT%CSR%item, &
      & monolis%MAT%C%A, n_base_i, n_base_j, monolis%MAT%ndof, connectivity_i, connectivity_j, mat)
  end subroutine monolis_add_matrix_to_sparse_matrix_offdiag_C

  !# CSR data setter
  !> @ingroup matrix
  !> BCSR 形式の疎行列を直接設定（実数型）
  subroutine monolis_set_matrix_BCSR_R(monolis, N, NP, NDOF, NZ, A, index, item)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_structure), intent(inout) :: monolis
    !> [in] 内部自由度数
    integer(kint), intent(in) :: N
    !> [in] 全自由度数
    integer(kint), intent(in) :: NP
    !> [in] ブロックサイズ
    integer(kint), intent(in) :: NDOF
    !> [in] 非零要素数
    integer(kint), intent(in) :: NZ
    !> [in] index 配列
    integer(kint), intent(in) :: index(:)
    !> [in] item 配列
    integer(kint), intent(in) :: item(:)
    !> [in] 行列値配列
    real(kdouble), intent(in) :: A(:)
    integer(kint) :: i

    call monolis_pdealloc_I_1d(monolis%MAT%CSR%index)
    call monolis_pdealloc_I_1d(monolis%MAT%CSR%item)
    call monolis_pdealloc_I_1d(monolis%MAT%CSC%index)
    call monolis_pdealloc_I_1d(monolis%MAT%CSC%item)
    call monolis_pdealloc_I_1d(monolis%MAT%CSC%perm)
    call monolis_pdealloc_R_1d(monolis%MAT%R%A)
    call monolis_pdealloc_R_1d(monolis%MAT%R%B)
    call monolis_pdealloc_R_1d(monolis%MAT%R%X)

    monolis%MAT%N = N
    monolis%MAT%NP = NP
    monolis%MAT%NDOF = NDOF

    call monolis_palloc_I_1d(monolis%MAT%CSR%index, NP + 1)
    call monolis_palloc_I_1d(monolis%MAT%CSR%item, NZ)
    call monolis_palloc_R_1d(monolis%MAT%R%A, NDOF*NDOF*NZ)
    call monolis_palloc_R_1d(monolis%MAT%R%B, NDOF*NP)
    call monolis_palloc_R_1d(monolis%MAT%R%X, NDOF*NP)

    do i = 1, NP
      monolis%MAT%CSR%index(i) = index(i)
    enddo

    do i = 1, NZ
      monolis%MAT%CSR%item(i) = item(i)
    enddo

    do i = 1, NDOF*NDOF*NZ
      monolis%MAT%R%A(i) = A(i)
    enddo

    call monolis_palloc_I_1d(monolis%MAT%CSC%index, NP + 1)
    call monolis_palloc_I_1d(monolis%MAT%CSC%item, NZ)
    call monolis_palloc_I_1d(monolis%MAT%CSC%perm, NZ)

    call monolis_get_CSC_format(N, NP, NZ, &
      & monolis%MAT%CSR%index, monolis%MAT%CSR%item, &
      & monolis%MAT%CSC%index, monolis%MAT%CSC%item, monolis%MAT%CSC%perm)
  end subroutine monolis_set_matrix_BCSR_R

  !> @ingroup matrix
  !> BCSR 形式の疎行列の行列値を直接設定（実数型）
  subroutine monolis_set_matrix_BCSR_mat_val_R(monolis, NDOF, NZ, A)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_structure), intent(inout) :: monolis
    !> [in] ブロックサイズ
    integer(kint), intent(in) :: NDOF
    !> [in] 非零要素数
    integer(kint), intent(in) :: NZ
    !> [in] 行列値配列
    real(kdouble), intent(in) :: A(:)
    integer(kint) :: i

    do i = 1, NDOF*NDOF*NZ
      monolis%MAT%R%A(i) = A(i)
    enddo
  end subroutine monolis_set_matrix_BCSR_mat_val_R

  !# boundary condition
  !> @ingroup matrix
  !> 境界条件処理（実数型）
  subroutine monolis_set_Dirichlet_bc_R(monolis, B, node_id, ndof_bc, val)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_structure), intent(inout) :: monolis
    !> [in,out] 右辺ベクトル
    real(kdouble), intent(inout) :: B(:)
    !> [in] 自由度番号
    integer(kint), intent(in) :: node_id
    !> [in] ブロック番号
    integer(kint), intent(in) :: ndof_bc
    !> [in] 境界条件の設定値
    real(kdouble), intent(in) :: val

    call monolis_set_Dirichlet_bc_main_R(monolis%MAT%CSR%index, monolis%MAT%CSR%item, monolis%MAT%R%A, B, &
      & monolis%MAT%CSC%index, monolis%MAT%CSC%item, monolis%MAT%CSC%perm, &
      & monolis%MAT%ndof, node_id, ndof_bc, val)
  end subroutine monolis_set_Dirichlet_bc_R

  !> @ingroup matrix
  !> 境界条件処理（複素数型）
  subroutine monolis_set_Dirichlet_bc_C(monolis, B, node_id, ndof_bc, val)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_structure), intent(inout) :: monolis
    !> [in,out] 右辺ベクトル
    complex(kdouble), intent(inout) :: B(:)
    !> [in] 自由度番号
    integer(kint), intent(in) :: node_id
    !> [in] ブロック番号
    integer(kint), intent(in) :: ndof_bc
    !> [in] 境界条件の設定値
    complex(kdouble), intent(in) :: val

    call monolis_set_Dirichlet_bc_main_C(monolis%MAT%CSR%index, monolis%MAT%CSR%item, monolis%MAT%C%A, B, &
      & monolis%MAT%CSC%index, monolis%MAT%CSC%item, monolis%MAT%CSC%perm, &
      & monolis%MAT%ndof, node_id, ndof_bc, val)
  end subroutine monolis_set_Dirichlet_bc_C

  !# 密行列から疎行列の取得
  !> @ingroup matrix
  !> 密行列から疎行列情報の取得（実数型）
  subroutine monolis_get_sparse_matrix_from_dense_matrix_R(MAT, n, m, n_dof, dense)
    implicit none
    !> [in] monolis MAT 構造体
    type(monolis_mat), intent(inout) :: MAT
    integer(kint), intent(in) :: n
    integer(kint), intent(in) :: m
    integer(kint), intent(in) :: n_dof
    real(kdouble), intent(in) :: dense(:,:)
    !logical, intent(in) :: is_symmetric
    integer(kint) :: n_elem, n_base
    integer(kint), allocatable :: elem(:,:), index(:), item(:)
    integer(kint), allocatable :: nodal_index(:), nodal_item(:)

    !# get nonzero pattern
    call monolis_get_nonzero_pattern_from_dense(n, n_dof, dense, n_elem, n_base, elem)

    call gedatsu_convert_simple_mesh_to_connectivity_graph(n_elem, n_base, elem, index, item)

    call gedatsu_convert_connectivity_graph_to_nodal_graph(n, n_elem, index, item, nodal_index, nodal_item)

    call monolis_get_nonzero_pattern_by_nodal_graph_main(MAT, n, n_dof, nodal_index, nodal_item)

    call monolis_alloc_nonzero_pattern_mat_val_R(MAT)

    !# set nonzero value
    call monolis_set_value_from_dense(MAT, n, m, n_dof, dense)
  end subroutine monolis_get_sparse_matrix_from_dense_matrix_R

  subroutine monolis_set_value_from_dense(MAT, n, m, n_dof, dense)
    implicit none
    !> [in] monolis MAT 構造体
    type(monolis_mat), intent(inout) :: MAT
    integer(kint), intent(in) :: n
    integer(kint), intent(in) :: m
    integer(kint), intent(in) :: n_dof
    real(kdouble), intent(in) :: dense(:,:)
    real(kdouble) :: block(n_dof,n_dof)
    integer(kint) :: i, j, in, jn, im, jm, ip, jp, iq, jq

    do i = 1, m
    aa:do j = 1, n
      do in = 1, n_dof
      do jn = 1, n_dof
        im = n_dof*(i - 1) + in
        jm = n_dof*(j - 1) + jn
        if(dense(jm,im) /= 0.0d0)then
          do ip = 1, n_dof
          do jp = 1, n_dof
            iq = n_dof*(i - 1) + ip
            jq = n_dof*(j - 1) + jp
            block(jp,ip) = dense(jq,iq)
          enddo
          enddo
          call monolis_set_block_to_sparse_matrix_main_R(MAT%CSR%index, MAT%CSR%item, MAT%R%A, MAT%ndof, j, i, block)
          cycle aa
        endif
      enddo
      enddo
    enddo aa
    enddo
  end subroutine monolis_set_value_from_dense
end module mod_monolis_spmat_handler
