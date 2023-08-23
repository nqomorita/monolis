!> 非零構造決定モジュール
module mod_monolis_spmat_nonzero_pattern
  use mod_monolis_utils
  use mod_gedatsu
  use mod_monolis_def_mat
  use mod_monolis_def_struc
  use mod_monolis_spmat_nonzero_pattern_util

  implicit none

contains

  !> @ingroup nzpattern
  !> 単一メッシュデータから疎行列パターンを決定（実数型）
  subroutine monolis_get_nonzero_pattern_by_simple_mesh_R(monolis, n_node, n_base, ndof, n_elem, elem)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_structure), intent(inout) :: monolis
    !> [in] 節点数
    integer(kint), intent(in) :: n_node
    !> [in] 要素を構成する形状関数の数
    integer(kint), intent(in) :: n_base
    !> [in] 計算点が持つ自由度
    integer(kint), intent(in) :: ndof
    !> [in] 要素数
    integer(kint), intent(in) :: n_elem
    !> [in] 要素コネクティビティ
    integer(kint), intent(in) :: elem(:,:)
    integer(kint), allocatable :: index(:)
    integer(kint), allocatable :: item(:)
    integer(kint), allocatable :: nodal_index(:)
    integer(kint), allocatable :: nodal_item(:)

    call gedatsu_convert_simple_mesh_to_connectivity_graph &
      & (n_elem, n_base, elem, index, item)

    call gedatsu_convert_connectivity_graph_to_nodal_graph &
      & (n_node, n_elem, index, item, nodal_index, nodal_item)

    call monolis_get_nonzero_pattern_by_nodal_graph_main &
      & (monolis%MAT, n_node, ndof, nodal_index, nodal_item)

    call monolis_alloc_nonzero_pattern_mat_val_R(monolis%MAT)
  end subroutine monolis_get_nonzero_pattern_by_simple_mesh_R

  !> @ingroup nzpattern
  !> 単一メッシュデータから疎行列パターンを決定（複素数型）
  subroutine monolis_get_nonzero_pattern_by_simple_mesh_C(monolis, n_node, n_base, ndof, n_elem, elem)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_structure), intent(inout) :: monolis
    !> [in] 節点数
    integer(kint), intent(in) :: n_node
    !> [in] 要素を構成する形状関数の数
    integer(kint), intent(in) :: n_base
    !> [in] 計算点が持つ自由度
    integer(kint), intent(in) :: ndof
    !> [in] 要素数
    integer(kint), intent(in) :: n_elem
    !> [in] 要素コネクティビティ
    integer(kint), intent(in) :: elem(:,:)
    integer(kint), allocatable :: index(:)
    integer(kint), allocatable :: item(:)
    integer(kint), allocatable :: nodal_index(:)
    integer(kint), allocatable :: nodal_item(:)

    call gedatsu_convert_simple_mesh_to_connectivity_graph &
      & (n_elem, n_base, elem, index, item)

    call gedatsu_convert_connectivity_graph_to_nodal_graph &
      & (n_node, n_elem, index, item, nodal_index, nodal_item)

    call monolis_get_nonzero_pattern_by_nodal_graph_main &
      & (monolis%MAT, n_node, ndof, nodal_index, nodal_item)

    call monolis_alloc_nonzero_pattern_mat_val_C(monolis%MAT)
  end subroutine monolis_get_nonzero_pattern_by_simple_mesh_C

  !> @ingroup nzpattern
  !> 単一メッシュデータから疎行列パターンを決定（任意節点自由度、実数型）
  subroutine monolis_get_nonzero_pattern_by_simple_mesh_with_arbitrary_dof_R( &
    monolis, n_node, n_base, n_dof_list, n_elem, elem)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_structure), intent(inout) :: monolis
    !> [in] 節点数
    integer(kint), intent(in) :: n_node
    !> [in] 要素を構成する形状関数の数
    integer(kint), intent(in) :: n_base
    !> [in] 自由度数リスト
    integer(kint), intent(in) :: n_dof_list(:)
    !> [in] 要素数
    integer(kint), intent(in) :: n_elem
    !> [in] 要素コネクティビティ
    integer(kint), intent(in) :: elem(:,:)
    integer(kint), allocatable :: index(:)
    integer(kint), allocatable :: item(:)
    integer(kint), allocatable :: nodal_index(:)
    integer(kint), allocatable :: nodal_item(:)

    call gedatsu_convert_simple_mesh_to_connectivity_graph &
      & (n_elem, n_base, elem, index, item)

    call gedatsu_convert_connectivity_graph_to_nodal_graph &
      & (n_node, n_elem, index, item, nodal_index, nodal_item)

    call monolis_get_nonzero_pattern_by_nodal_graph_with_arbit_main &
      & (monolis%MAT, n_node, n_dof_list, nodal_index, nodal_item)

    call monolis_alloc_nonzero_pattern_mat_val_R(monolis%MAT)
  end subroutine monolis_get_nonzero_pattern_by_simple_mesh_with_arbitrary_dof_R

  !> @ingroup nzpattern
  !> 単一メッシュデータから疎行列パターンを決定（任意節点自由度、複素数型）
  subroutine monolis_get_nonzero_pattern_by_simple_mesh_with_arbitrary_dof_C( &
    monolis, n_node, n_base, n_dof_list, n_elem, elem)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_structure), intent(inout) :: monolis
    !> [in] 節点数
    integer(kint), intent(in) :: n_node
    !> [in] 要素を構成する形状関数の数
    integer(kint), intent(in) :: n_base
    !> [in] 自由度数リスト
    integer(kint), intent(in) :: n_dof_list(:)
    !> [in] 要素数
    integer(kint), intent(in) :: n_elem
    !> [in] 要素コネクティビティ
    integer(kint), intent(in) :: elem(:,:)
    integer(kint), allocatable :: index(:)
    integer(kint), allocatable :: item(:)
    integer(kint), allocatable :: nodal_index(:)
    integer(kint), allocatable :: nodal_item(:)

    call gedatsu_convert_simple_mesh_to_connectivity_graph &
      & (n_elem, n_base, elem, index, item)

    call gedatsu_convert_connectivity_graph_to_nodal_graph &
      & (n_node, n_elem, index, item, nodal_index, nodal_item)

    call monolis_get_nonzero_pattern_by_nodal_graph_with_arbit_main &
      & (monolis%MAT, n_node, n_dof_list, nodal_index, nodal_item)

    call monolis_alloc_nonzero_pattern_mat_val_C(monolis%MAT)
  end subroutine monolis_get_nonzero_pattern_by_simple_mesh_with_arbitrary_dof_C

  !> @ingroup nzpattern
  !> コネクティビティグラフから疎行列パターンを決定（実数型）
  subroutine monolis_get_nonzero_pattern_by_connectivity_R( &
      & monolis, n_node, ndof, n_elem, index, item)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_structure), intent(inout) :: monolis
    !> [in] 節点数
    integer(kint), intent(in) :: n_node
    !> [in] 計算点が持つ自由度
    integer(kint), intent(in) :: ndof
    !> [in] 要素数
    integer(kint), intent(in) :: n_elem
    !> [out] コネクティビティ index 配列
    integer(kint), intent(out) :: index(:)
    !> [out] コネクティビティ item 配列
    integer(kint), intent(out) :: item(:)
    integer(kint), allocatable :: nodal_index(:)
    integer(kint), allocatable :: nodal_item(:)

    call gedatsu_convert_connectivity_graph_to_nodal_graph &
      & (n_node, n_elem, index, item, nodal_index, nodal_item)

    call monolis_get_nonzero_pattern_by_nodal_graph_main &
      & (monolis%MAT, n_node, ndof, nodal_index, nodal_item)

    call monolis_alloc_nonzero_pattern_mat_val_R(monolis%MAT)
  end subroutine monolis_get_nonzero_pattern_by_connectivity_R

  !> @ingroup nzpattern
  !> コネクティビティグラフから疎行列パターンを決定（複素数型）
  subroutine monolis_get_nonzero_pattern_by_connectivity_C( &
      & monolis, n_node, ndof, n_elem, index, item)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_structure), intent(inout) :: monolis
    !> [in] 節点数
    integer(kint), intent(in) :: n_node
    !> [in] 計算点が持つ自由度
    integer(kint), intent(in) :: ndof
    !> [in] 要素数
    integer(kint), intent(in) :: n_elem
    !> [out] コネクティビティ index 配列
    integer(kint), intent(out) :: index(:)
    !> [out] コネクティビティ item 配列
    integer(kint), intent(out) :: item(:)
    integer(kint), allocatable :: nodal_index(:)
    integer(kint), allocatable :: nodal_item(:)

    call gedatsu_convert_connectivity_graph_to_nodal_graph &
      & (n_node, n_elem, index, item, nodal_index, nodal_item)

    call monolis_get_nonzero_pattern_by_nodal_graph_main &
      & (monolis%MAT, n_node, ndof, nodal_index, nodal_item)

    call monolis_alloc_nonzero_pattern_mat_val_C(monolis%MAT)
  end subroutine monolis_get_nonzero_pattern_by_connectivity_C

  !> @ingroup nzpattern
  !> 節点グラフから疎行列パターンを決定（実数型）
  subroutine monolis_get_nonzero_pattern_by_nodal_graph_R(monolis, n_node, ndof, index, item)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_structure), intent(inout) :: monolis
    !> [in] 節点数
    integer(kint), intent(in) :: n_node
    !> [in] 計算点が持つ自由度
    integer(kint), intent(in) :: ndof
    !> [in] index 配列
    integer(kint), intent(in) :: index(:)
    !> [in] item 配列
    integer(kint), intent(in) :: item(:)

    call monolis_get_nonzero_pattern_by_nodal_graph_main &
      & (monolis%MAT, n_node, ndof, index, item)

    call monolis_alloc_nonzero_pattern_mat_val_R(monolis%MAT)
  end subroutine monolis_get_nonzero_pattern_by_nodal_graph_R

  !> @ingroup nzpattern
  !> 節点グラフから疎行列パターンを決定（複素数型）
  subroutine monolis_get_nonzero_pattern_by_nodal_graph_C(monolis, n_node, ndof, index, item)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_structure), intent(inout) :: monolis
    !> [in] 節点数
    integer(kint), intent(in) :: n_node
    !> [in] 計算点が持つ自由度
    integer(kint), intent(in) :: ndof
    !> [in] index 配列
    integer(kint), intent(in) :: index(:)
    !> [in] item 配列
    integer(kint), intent(in) :: item(:)

    call monolis_get_nonzero_pattern_by_nodal_graph_main &
      & (monolis%MAT, n_node, ndof, index, item)

    call monolis_alloc_nonzero_pattern_mat_val_C(monolis%MAT)
  end subroutine monolis_get_nonzero_pattern_by_nodal_graph_C

  !> @ingroup nzpattern
  !> 節点グラフから疎行列パターンを決定（任意節点自由度、実数型）
  subroutine monolis_get_nonzero_pattern_by_nodal_graph_with_arbitrary_dof_R( &
    monolis, n_node, n_dof_list, index, item)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_structure), intent(inout) :: monolis
    !> [in] 節点数
    integer(kint), intent(in) :: n_node
    !> [in] 自由度リスト
    integer(kint), intent(in) :: n_dof_list(:)
    !> [in] index 配列
    integer(kint), intent(in) :: index(:)
    !> [in] item 配列
    integer(kint), intent(in) :: item(:)

    call monolis_get_nonzero_pattern_by_nodal_graph_with_arbit_main &
      & (monolis%MAT, n_node, n_dof_list, index, item)

    call monolis_alloc_nonzero_pattern_mat_val_R(monolis%MAT)
  end subroutine monolis_get_nonzero_pattern_by_nodal_graph_with_arbitrary_dof_R

  !> @ingroup nzpattern
  !> 節点グラフから疎行列パターンを決定（任意節点自由度、複素数型）
  subroutine monolis_get_nonzero_pattern_by_nodal_graph_with_arbitrary_dof_C( &
    monolis, n_node, n_dof_list, index, item)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_structure), intent(inout) :: monolis
    !> [in] 節点数
    integer(kint), intent(in) :: n_node
    !> [in] 自由度リスト
    integer(kint), intent(in) :: n_dof_list(:)
    !> [in] index 配列
    integer(kint), intent(in) :: index(:)
    !> [in] item 配列
    integer(kint), intent(in) :: item(:)

    call monolis_get_nonzero_pattern_by_nodal_graph_with_arbit_main &
      & (monolis%MAT, n_node, n_dof_list, index, item)

    call monolis_alloc_nonzero_pattern_mat_val_C(monolis%MAT)
  end subroutine monolis_get_nonzero_pattern_by_nodal_graph_with_arbitrary_dof_C
end module mod_monolis_spmat_nonzero_pattern
