!> 非零構造決定モジュール
module mod_monolis_spmat_nonzero_pattern
  use mod_monolis_utils
  use mod_gedatsu
  use mod_monolis_def_mat
  use mod_monolis_def_struc
  use mod_monolis_spmat_nonzero_pattern_util

  implicit none

contains

  !> 単一メッシュデータから疎行列パターンを決定（実数型）
  subroutine monolis_get_nonzero_pattern_by_simple_mesh_R(monolis, n_node, n_base, ndof, n_elem, elem)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> 節点数
    integer(kint) :: n_node
    !> 要素を構成する節点数
    integer(kint) :: n_base
    !> 自由度数
    integer(kint) :: ndof
    !> 要素数
    integer(kint) :: n_elem
    !> 要素
    integer(kint) :: elem(:,:)
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

  !> 単一メッシュデータから疎行列パターンを決定（複素数型）
  subroutine monolis_get_nonzero_pattern_by_simple_mesh_C(monolis, n_node, n_base, ndof, n_elem, elem)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> 節点数
    integer(kint) :: n_node
    !> 要素を構成する節点数
    integer(kint) :: n_base
    !> 自由度数
    integer(kint) :: ndof
    !> 要素数
    integer(kint) :: n_elem
    !> 要素
    integer(kint) :: elem(:,:)
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

  !> 単一メッシュデータから疎行列パターンを決定（任意節点自由度、実数型）
  subroutine monolis_get_nonzero_pattern_by_simple_mesh_with_arbitrary_dof_R &
    (monolis, n_node, n_base, n_dof_list, n_elem, elem)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> 節点数
    integer(kint) :: n_node
    !> 要素を構成する節点数
    integer(kint) :: n_base
    !> 自由度数リスト
    integer(kint) :: n_dof_list(:)
    !> 要素数
    integer(kint) :: n_elem
    !> 要素
    integer(kint) :: elem(:,:)
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

  !> 単一メッシュデータから疎行列パターンを決定（任意節点自由度、複素数型）
  subroutine monolis_get_nonzero_pattern_by_simple_mesh_with_arbitrary_dof_C &
    (monolis, n_node, n_base, n_dof_list, n_elem, elem)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> 節点数
    integer(kint) :: n_node
    !> 要素を構成する節点数
    integer(kint) :: n_base
    !> 自由度数リスト
    integer(kint) :: n_dof_list(:)
    !> 要素数
    integer(kint) :: n_elem
    !> 要素
    integer(kint) :: elem(:,:)
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

  !> コネクティビティグラフから疎行列パターンを決定（実数型）
  subroutine monolis_get_nonzero_pattern_by_connectivity_R &
      & (monolis, n_node, n_base, ndof, n_elem, index, item)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> 節点数
    integer(kint) :: n_node
    !> 要素を構成する節点数
    integer(kint) :: n_base
    !> 自由度数
    integer(kint) :: ndof
    !> 要素数
    integer(kint) :: n_elem
    !> コネクティビティ index 配列
    integer(kint) :: index(:)
    !> コネクティビティ item 配列
    integer(kint) :: item(:)
    integer(kint), allocatable :: nodal_index(:)
    integer(kint), allocatable :: nodal_item(:)

    call gedatsu_convert_connectivity_graph_to_nodal_graph &
      & (n_node, n_elem, index, item, nodal_index, nodal_item)

    call monolis_get_nonzero_pattern_by_nodal_graph_main &
      & (monolis%MAT, n_node, ndof, nodal_index, nodal_item)

    call monolis_alloc_nonzero_pattern_mat_val_R(monolis%MAT)
  end subroutine monolis_get_nonzero_pattern_by_connectivity_R

  !> コネクティビティグラフから疎行列パターンを決定（複素数型）
  subroutine monolis_get_nonzero_pattern_by_connectivity_C &
      & (monolis, n_node, n_base, ndof, n_elem, index, item)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> 節点数
    integer(kint) :: n_node
    !> 要素を構成する節点数
    integer(kint) :: n_base
    !> 自由度数
    integer(kint) :: ndof
    !> 要素数
    integer(kint) :: n_elem
    !> コネクティビティ index 配列
    integer(kint) :: index(:)
    !> コネクティビティ item 配列
    integer(kint) :: item(:)
    integer(kint), allocatable :: nodal_index(:)
    integer(kint), allocatable :: nodal_item(:)

    call gedatsu_convert_connectivity_graph_to_nodal_graph &
      & (n_node, n_elem, index, item, nodal_index, nodal_item)

    call monolis_get_nonzero_pattern_by_nodal_graph_main &
      & (monolis%MAT, n_node, ndof, nodal_index, nodal_item)

    call monolis_alloc_nonzero_pattern_mat_val_C(monolis%MAT)
  end subroutine monolis_get_nonzero_pattern_by_connectivity_C

  !> 節点グラフから疎行列パターンを決定（実数型）
  subroutine monolis_get_nonzero_pattern_by_nodal_graph_R(monolis, n_node, ndof, index, item)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> 節点数
    integer(kint) :: n_node
    !> 自由度数
    integer(kint) :: ndof
    !> index 配列
    integer(kint) :: index(:)
    !> item 配列
    integer(kint) :: item(:)

    call monolis_get_nonzero_pattern_by_nodal_graph_main &
      & (monolis%MAT, n_node, ndof, index, item)

    call monolis_alloc_nonzero_pattern_mat_val_R(monolis%MAT)
  end subroutine monolis_get_nonzero_pattern_by_nodal_graph_R

  !> 節点グラフから疎行列パターンを決定（複素数型）
  subroutine monolis_get_nonzero_pattern_by_nodal_graph_C(monolis, n_node, ndof, index, item)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> 節点数
    integer(kint) :: n_node
    !> 自由度数
    integer(kint) :: ndof
    !> index 配列
    integer(kint) :: index(:)
    !> item 配列
    integer(kint) :: item(:)

    call monolis_get_nonzero_pattern_by_nodal_graph_main &
      & (monolis%MAT, n_node, ndof, index, item)

    call monolis_alloc_nonzero_pattern_mat_val_C(monolis%MAT)
  end subroutine monolis_get_nonzero_pattern_by_nodal_graph_C

  !> 節点グラフから疎行列パターンを決定（任意節点自由度、実数型）
  subroutine monolis_get_nonzero_pattern_by_nodal_graph_with_arbitrary_dof_R &
    (monolis, n_node, n_dof_list, index, item)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> 節点数
    integer(kint) :: n_node
    !> 自由度リスト
    integer(kint) :: n_dof_list(:)
    !> index 配列
    integer(kint) :: index(:)
    !> item 配列
    integer(kint) :: item(:)

    call monolis_get_nonzero_pattern_by_nodal_graph_with_arbit_main &
      & (monolis%MAT, n_node, n_dof_list, index, item)

    call monolis_alloc_nonzero_pattern_mat_val_R(monolis%MAT)
  end subroutine monolis_get_nonzero_pattern_by_nodal_graph_with_arbitrary_dof_R

  !> 節点グラフから疎行列パターンを決定（任意節点自由度、複素数型）
  subroutine monolis_get_nonzero_pattern_by_nodal_graph_with_arbitrary_dof_C &
    (monolis, n_node, n_dof_list, index, item)
    implicit none
    !> monolis 構造体
    type(monolis_structure) :: monolis
    !> 節点数
    integer(kint) :: n_node
    !> 自由度リスト
    integer(kint) :: n_dof_list(:)
    !> index 配列
    integer(kint) :: index(:)
    !> item 配列
    integer(kint) :: item(:)

    call monolis_get_nonzero_pattern_by_nodal_graph_with_arbit_main &
      & (monolis%MAT, n_node, n_dof_list, index, item)

    call monolis_alloc_nonzero_pattern_mat_val_C(monolis%MAT)
  end subroutine monolis_get_nonzero_pattern_by_nodal_graph_with_arbitrary_dof_C
end module mod_monolis_spmat_nonzero_pattern
