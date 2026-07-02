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
    integer(kint), allocatable :: n_dof_list(:)

    call gedatsu_convert_simple_mesh_to_connectivity_graph &
      & (n_elem, n_base, elem, index, item)

    call gedatsu_convert_connectivity_graph_to_nodal_graph &
      & (n_node, n_elem, index, item, nodal_index, nodal_item)

    call monolis_get_nonzero_pattern_by_nodal_graph_main &
      & (monolis%MAT, n_node, ndof, nodal_index, nodal_item)

    call monolis_alloc_I_1d(n_dof_list, n_node)
    n_dof_list = ndof
    call monolis_set_n_dof_index(monolis%MAT, n_dof_list)

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
    integer(kint), allocatable :: n_dof_list(:)

    call gedatsu_convert_simple_mesh_to_connectivity_graph &
      & (n_elem, n_base, elem, index, item)

    call gedatsu_convert_connectivity_graph_to_nodal_graph &
      & (n_node, n_elem, index, item, nodal_index, nodal_item)

    call monolis_get_nonzero_pattern_by_nodal_graph_main &
      & (monolis%MAT, n_node, ndof, nodal_index, nodal_item)

    call monolis_alloc_I_1d(n_dof_list, n_node)
    n_dof_list = ndof
    call monolis_set_n_dof_index(monolis%MAT, n_dof_list)

    call monolis_alloc_nonzero_pattern_mat_val_C(monolis%MAT)
  end subroutine monolis_get_nonzero_pattern_by_simple_mesh_C

  !> @ingroup nzpattern
  !> 単一メッシュデータから疎行列パターンを決定（任意節点自由度、実数型）
  subroutine monolis_get_nonzero_pattern_by_simple_mesh_V_R( &
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

    call monolis_get_nonzero_pattern_by_nodal_graph_main &
      & (monolis%MAT, n_node, -1, nodal_index, nodal_item)

    call monolis_set_n_dof_index(monolis%MAT, n_dof_list)

    call monolis_alloc_nonzero_pattern_mat_val_V_R(monolis%MAT)
  end subroutine monolis_get_nonzero_pattern_by_simple_mesh_V_R

  !> @ingroup nzpattern
  !> 単一メッシュデータから疎行列パターンを決定（任意節点自由度、複素数型）
  subroutine monolis_get_nonzero_pattern_by_simple_mesh_V_C( &
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

    call monolis_get_nonzero_pattern_by_nodal_graph_main &
      & (monolis%MAT, n_node, -1, nodal_index, nodal_item)

    call monolis_set_n_dof_index(monolis%MAT, n_dof_list)

    call monolis_alloc_nonzero_pattern_mat_val_V_C(monolis%MAT)
  end subroutine monolis_get_nonzero_pattern_by_simple_mesh_V_C

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
    integer(kint), allocatable :: n_dof_list(:)

    call gedatsu_convert_connectivity_graph_to_nodal_graph &
      & (n_node, n_elem, index, item, nodal_index, nodal_item)

    call monolis_get_nonzero_pattern_by_nodal_graph_main &
      & (monolis%MAT, n_node, ndof, nodal_index, nodal_item)

    call monolis_alloc_I_1d(n_dof_list, n_node)
    n_dof_list = ndof
    call monolis_set_n_dof_index(monolis%MAT, n_dof_list)

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
    integer(kint), allocatable :: n_dof_list(:)

    call gedatsu_convert_connectivity_graph_to_nodal_graph &
      & (n_node, n_elem, index, item, nodal_index, nodal_item)

    call monolis_get_nonzero_pattern_by_nodal_graph_main &
      & (monolis%MAT, n_node, ndof, nodal_index, nodal_item)

    call monolis_alloc_I_1d(n_dof_list, n_node)
    n_dof_list = ndof
    call monolis_set_n_dof_index(monolis%MAT, n_dof_list)

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
    integer(kint), allocatable :: n_dof_list(:)

    call monolis_get_nonzero_pattern_by_nodal_graph_main &
      & (monolis%MAT, n_node, ndof, index, item)

    call monolis_alloc_I_1d(n_dof_list, n_node)
    n_dof_list = ndof
    call monolis_set_n_dof_index(monolis%MAT, n_dof_list)

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
    integer(kint), allocatable :: n_dof_list(:)

    call monolis_get_nonzero_pattern_by_nodal_graph_main &
      & (monolis%MAT, n_node, ndof, index, item)

    call monolis_alloc_I_1d(n_dof_list, n_node)
    n_dof_list = ndof
    call monolis_set_n_dof_index(monolis%MAT, n_dof_list)

    call monolis_alloc_nonzero_pattern_mat_val_C(monolis%MAT)
  end subroutine monolis_get_nonzero_pattern_by_nodal_graph_C

  !> @ingroup nzpattern
  !> 節点グラフから疎行列パターンを決定（任意節点自由度、実数型）
  subroutine monolis_get_nonzero_pattern_by_nodal_graph_V_R( &
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

    call monolis_get_nonzero_pattern_by_nodal_graph_main &
      & (monolis%MAT, n_node, -1, index, item)

    call monolis_set_n_dof_index(monolis%MAT, n_dof_list)

    call monolis_alloc_nonzero_pattern_mat_val_V_R(monolis%MAT)
  end subroutine monolis_get_nonzero_pattern_by_nodal_graph_V_R

  !> @ingroup nzpattern
  !> 節点グラフから疎行列パターンを決定（任意節点自由度、複素数型）
  subroutine monolis_get_nonzero_pattern_by_nodal_graph_V_C( &
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

    call monolis_get_nonzero_pattern_by_nodal_graph_main &
      & (monolis%MAT, n_node, -1, index, item)

    call monolis_set_n_dof_index(monolis%MAT, n_dof_list)

    call monolis_alloc_nonzero_pattern_mat_val_V_C(monolis%MAT)
  end subroutine monolis_get_nonzero_pattern_by_nodal_graph_V_C

  !> @ingroup nzpattern
  !> マージン（固定バンド幅）付きで疎行列パターンを確保（実数型）
  !> @details 粒子法のように毎ステップ近傍が変化する問題向けに、各行を n_bandwidth スロットの
  !> 固定幅 BCSR として一度だけ確保する。以降は \ref monolis_update_nonzero_pattern_with_margin_by_nodal_graph_R
  !> で構造のみを更新でき、疎構造の再構築コストを回避できる。
  !> n_bandwidth は対角成分を含む 1 行あたりの最大非零ブロック数（許容近傍数は n_bandwidth-1）。
  subroutine monolis_alloc_nonzero_pattern_with_margin_R(monolis, n_node, ndof, n_bandwidth)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_structure), intent(inout) :: monolis
    !> [in] 節点数
    integer(kint), intent(in) :: n_node
    !> [in] 計算点が持つ自由度
    integer(kint), intent(in) :: ndof
    !> [in] 1 行あたりの固定バンド幅（対角を含む）
    integer(kint), intent(in) :: n_bandwidth

    call monolis_alloc_nonzero_pattern_with_margin_main(monolis%MAT, n_node, ndof, n_bandwidth)

    call monolis_pdealloc_R_1d(monolis%MAT%R%A)
    call monolis_pdealloc_R_1d(monolis%MAT%R%B)
    call monolis_pdealloc_R_1d(monolis%MAT%R%X)
    call monolis_palloc_R_1d(monolis%MAT%R%A, n_node*n_bandwidth*ndof*ndof)
    call monolis_palloc_R_1d(monolis%MAT%R%B, n_node*ndof)
    call monolis_palloc_R_1d(monolis%MAT%R%X, n_node*ndof)
  end subroutine monolis_alloc_nonzero_pattern_with_margin_R

  !> @ingroup nzpattern
  !> マージン（固定バンド幅）付きで疎行列パターンを確保（複素数型）
  !> @details 実数版 \ref monolis_alloc_nonzero_pattern_with_margin_R の複素数版。
  subroutine monolis_alloc_nonzero_pattern_with_margin_C(monolis, n_node, ndof, n_bandwidth)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_structure), intent(inout) :: monolis
    !> [in] 節点数
    integer(kint), intent(in) :: n_node
    !> [in] 計算点が持つ自由度
    integer(kint), intent(in) :: ndof
    !> [in] 1 行あたりの固定バンド幅（対角を含む）
    integer(kint), intent(in) :: n_bandwidth

    call monolis_alloc_nonzero_pattern_with_margin_main(monolis%MAT, n_node, ndof, n_bandwidth)

    call monolis_pdealloc_C_1d(monolis%MAT%C%A)
    call monolis_pdealloc_C_1d(monolis%MAT%C%B)
    call monolis_pdealloc_C_1d(monolis%MAT%C%X)
    call monolis_palloc_C_1d(monolis%MAT%C%A, n_node*n_bandwidth*ndof*ndof)
    call monolis_palloc_C_1d(monolis%MAT%C%B, n_node*ndof)
    call monolis_palloc_C_1d(monolis%MAT%C%X, n_node*ndof)
  end subroutine monolis_alloc_nonzero_pattern_with_margin_C

  !> @ingroup nzpattern
  !> 固定バンド幅パターンの非零構造を節点グラフで更新（実数型）
  !> @details \ref monolis_alloc_nonzero_pattern_with_margin_R で確保済みの固定幅構造に対し、
  !> 各行の対角成分と近傍（並び順を保持）を再設定し、残りをパディングで埋める。
  !> いずれかの行で必要バンド幅が確保済みバンド幅を超える、もしくは節点数が変化した場合は、
  !> 必要幅で全配列を再確保する（値配列は 0 初期化される）。値配列の中身は変更しない。
  subroutine monolis_update_nonzero_pattern_with_margin_by_nodal_graph_R(monolis, n_node, index, item)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_structure), intent(inout) :: monolis
    !> [in] 節点数
    integer(kint), intent(in) :: n_node
    !> [in] 節点グラフの index 配列（0 始まりプレフィックス）
    integer(kint), intent(in) :: index(:)
    !> [in] 節点グラフの item 配列（近傍列番号、1 始まり）
    integer(kint), intent(in) :: item(:)
    integer(kint) :: i, deg, max_deg, needed_bw, cap_bw

    cap_bw = 0
    if(associated(monolis%MAT%CSR%item) .and. monolis%MAT%NP > 0) &
      & cap_bw = size(monolis%MAT%CSR%item)/monolis%MAT%NP

    max_deg = 0
    do i = 1, n_node
      deg = index(i + 1) - index(i)
      if(deg > max_deg) max_deg = deg
    enddo
    needed_bw = max_deg + 1

    if(needed_bw > cap_bw .or. n_node /= monolis%MAT%NP)then
      call monolis_std_warning_string &
        & ("monolis_update_nonzero_pattern_with_margin_by_nodal_graph_R: bandwidth margin exceeded, reallocating")
      call monolis_alloc_nonzero_pattern_with_margin_main(monolis%MAT, n_node, monolis%MAT%NDOF, needed_bw)
      call monolis_pdealloc_R_1d(monolis%MAT%R%A)
      call monolis_pdealloc_R_1d(monolis%MAT%R%B)
      call monolis_pdealloc_R_1d(monolis%MAT%R%X)
      call monolis_palloc_R_1d(monolis%MAT%R%A, n_node*needed_bw*monolis%MAT%NDOF*monolis%MAT%NDOF)
      call monolis_palloc_R_1d(monolis%MAT%R%B, n_node*monolis%MAT%NDOF)
      call monolis_palloc_R_1d(monolis%MAT%R%X, n_node*monolis%MAT%NDOF)
    endif

    call monolis_update_nonzero_pattern_fixed_width_main(monolis%MAT, n_node, index, item)
  end subroutine monolis_update_nonzero_pattern_with_margin_by_nodal_graph_R

  !> @ingroup nzpattern
  !> 固定バンド幅パターンの非零構造を節点グラフで更新（複素数型）
  !> @details 実数版 \ref monolis_update_nonzero_pattern_with_margin_by_nodal_graph_R の複素数版。
  subroutine monolis_update_nonzero_pattern_with_margin_by_nodal_graph_C(monolis, n_node, index, item)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_structure), intent(inout) :: monolis
    !> [in] 節点数
    integer(kint), intent(in) :: n_node
    !> [in] 節点グラフの index 配列（0 始まりプレフィックス）
    integer(kint), intent(in) :: index(:)
    !> [in] 節点グラフの item 配列（近傍列番号、1 始まり）
    integer(kint), intent(in) :: item(:)
    integer(kint) :: i, deg, max_deg, needed_bw, cap_bw

    cap_bw = 0
    if(associated(monolis%MAT%CSR%item) .and. monolis%MAT%NP > 0) &
      & cap_bw = size(monolis%MAT%CSR%item)/monolis%MAT%NP

    max_deg = 0
    do i = 1, n_node
      deg = index(i + 1) - index(i)
      if(deg > max_deg) max_deg = deg
    enddo
    needed_bw = max_deg + 1

    if(needed_bw > cap_bw .or. n_node /= monolis%MAT%NP)then
      call monolis_std_warning_string &
        & ("monolis_update_nonzero_pattern_with_margin_by_nodal_graph_C: bandwidth margin exceeded, reallocating")
      call monolis_alloc_nonzero_pattern_with_margin_main(monolis%MAT, n_node, monolis%MAT%NDOF, needed_bw)
      call monolis_pdealloc_C_1d(monolis%MAT%C%A)
      call monolis_pdealloc_C_1d(monolis%MAT%C%B)
      call monolis_pdealloc_C_1d(monolis%MAT%C%X)
      call monolis_palloc_C_1d(monolis%MAT%C%A, n_node*needed_bw*monolis%MAT%NDOF*monolis%MAT%NDOF)
      call monolis_palloc_C_1d(monolis%MAT%C%B, n_node*monolis%MAT%NDOF)
      call monolis_palloc_C_1d(monolis%MAT%C%X, n_node*monolis%MAT%NDOF)
    endif

    call monolis_update_nonzero_pattern_fixed_width_main(monolis%MAT, n_node, index, item)
  end subroutine monolis_update_nonzero_pattern_with_margin_by_nodal_graph_C
end module mod_monolis_spmat_nonzero_pattern
