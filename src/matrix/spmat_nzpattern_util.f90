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
    call monolis_palloc_I_1d(MAT%n_dof_index2, NZ + 1)

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
  
  !> @ingroup dev_matrix
  !> 疎行列パターンから節点グラフを決定（メイン関数）
  subroutine monolis_get_nodal_graph_by_nonzero_pattern(MAT, N, index, item)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_mat), intent(in) :: MAT
    !> [in] 節点数
    integer(kint), intent(out) :: N
    !> [in] index 配列
    integer(kint), allocatable :: index(:)
    !> [in] item 配列
    integer(kint), allocatable :: item(:)
    integer(kint) :: i, j, NZ, jS, jE, in, jn

    N = MAT%N
    NZ = MAT%CSR%index(N + 1) - N

    call monolis_alloc_I_1d(index, N + 1)
    call monolis_alloc_I_1d(item, NZ)

    in = 0
    do i = 1, N
      index(i + 1) = MAT%CSR%index(i + 1) - i
      jS = MAT%CSR%index(i) + 1
      jE = MAT%CSR%index(i + 1)
      do j = jS, jE
        jn = MAT%CSR%item(j)
        if(i == jn) cycle
        in = in + 1
        item(in) = jn
      enddo
    enddo
  end subroutine monolis_get_nodal_graph_by_nonzero_pattern

  !> @ingroup dev_matrix
  !> マージン（固定バンド幅）付きの疎行列パターンのメモリを確保（メイン関数）
  !> @details 各行を n_bandwidth スロットの固定幅として CSR/CSC 配列のメモリのみを確保する。
  !> CSR%index・CSR%item の具体値と CSC 形式は \ref monolis_update_nonzero_pattern_fixed_width_main で構築する。
  !> 自由度情報（一様自由度）のみ閉形式で設定する。
  !> n_bandwidth は対角成分を含む 1 行あたりの最大非零ブロック数（許容する近傍数は n_bandwidth-1）。
  subroutine monolis_alloc_nonzero_pattern_with_margin_main(MAT, n_node, ndof, n_bandwidth)
    implicit none
    !> [in,out] monolis 行列構造体
    type(monolis_mat), intent(inout) :: MAT
    !> [in] 節点数
    integer(kint), intent(in) :: n_node
    !> [in] 計算点が持つ自由度
    integer(kint), intent(in) :: ndof
    !> [in] 1 行あたりの固定バンド幅（対角を含む）
    integer(kint), intent(in) :: n_bandwidth
    integer(kint) :: i, j, NZ

    call monolis_std_debug_log_header("monolis_alloc_nonzero_pattern_with_margin_main")

    if(n_bandwidth < 1)then
      call monolis_std_error_string("monolis_alloc_nonzero_pattern_with_margin_main: n_bandwidth must be >= 1")
      call monolis_std_error_stop()
    endif

    MAT%N = n_node
    MAT%NP = n_node
    MAT%NDOF = ndof
    NZ = n_node*n_bandwidth

    !> CSR index
    call monolis_pdealloc_I_1d(MAT%CSR%index)
    call monolis_palloc_I_1d(MAT%CSR%index, n_node + 1)

    !> CSR item
    call monolis_pdealloc_I_1d(MAT%CSR%item)
    call monolis_palloc_I_1d(MAT%CSR%item, NZ)

    !> CSC 形式
    call monolis_pdealloc_I_1d(MAT%CSC%index)
    call monolis_pdealloc_I_1d(MAT%CSC%item)
    call monolis_pdealloc_I_1d(MAT%CSC%perm)
    call monolis_palloc_I_1d(MAT%CSC%index, n_node + 1)
    call monolis_palloc_I_1d(MAT%CSC%item, NZ)
    call monolis_palloc_I_1d(MAT%CSC%perm, NZ)

    !> 自由度情報（一様自由度の閉形式で設定）
    call monolis_pdealloc_I_1d(MAT%n_dof_list)
    call monolis_pdealloc_I_1d(MAT%n_dof_index)
    call monolis_pdealloc_I_1d(MAT%n_dof_index2)
    call monolis_palloc_I_1d(MAT%n_dof_list, n_node)
    call monolis_palloc_I_1d(MAT%n_dof_index, n_node + 1)
    call monolis_palloc_I_1d(MAT%n_dof_index2, NZ + 1)
  end subroutine monolis_alloc_nonzero_pattern_with_margin_main

  !> @ingroup dev_matrix
  !> 固定バンド幅パターンの非零構造を節点グラフで更新（メイン関数）
  !> @details 確保済みの固定幅 BCSR 構造に対し、CSR%index（固定ストライド）を設定し、
  !> 各行に対角成分と近傍、および対角・近傍と重複しない連続した列番号のパディングを格納する。
  !> メモリアクセス最適化のため各行を昇順ソートし、CSC 形式を再構築する。値配列は変更しない。
  subroutine monolis_update_nonzero_pattern_fixed_width_main(MAT, n_node, index, item)
    implicit none
    !> [in,out] monolis 行列構造体
    type(monolis_mat), intent(inout) :: MAT
    !> [in] 節点数
    integer(kint), intent(in) :: n_node
    !> [in] 節点グラフの index 配列（0 始まりプレフィックス）
    integer(kint), intent(in) :: index(:)
    !> [in] 節点グラフの item 配列（近傍列番号、1 始まり）
    integer(kint), intent(in) :: item(:)
    integer(kint) :: i, k, jS, jE, deg, bw, NZ, c, slot, ntouch
    integer(kint), allocatable :: used(:), touched(:)

    call monolis_std_debug_log_header("monolis_update_nonzero_pattern_fixed_width_main")

    !> 固定バンド幅（確保済み item 長から算出）
    bw = size(MAT%CSR%item)/MAT%NP

    if(bw > n_node)then
      call monolis_std_error_string &
        & ("monolis_update_nonzero_pattern_fixed_width_main: n_bandwidth must be <= n_node")
      call monolis_std_error_stop()
    endif

    !> CSR index（固定ストライド）
    do i = 1, n_node
      MAT%CSR%index(i + 1) = i*bw
    enddo

    call monolis_alloc_I_1d(used, n_node)
    call monolis_alloc_I_1d(touched, bw)

    do i = 1, n_node
      jS = MAT%CSR%index(i) + 1
      jE = MAT%CSR%index(i + 1)
      deg = index(i + 1) - index(i)
      if(deg + 1 > bw)then
        call monolis_std_error_string("monolis_update_nonzero_pattern_fixed_width_main: bandwidth overflow")
        call monolis_std_error_stop()
      endif

      ntouch = 0
      !> 対角成分
      MAT%CSR%item(jS) = i
      ntouch = ntouch + 1; touched(ntouch) = i; used(i) = 1
      !> 近傍成分
      do k = 1, deg
        c = item(index(i) + k)
        MAT%CSR%item(jS + k) = c
        ntouch = ntouch + 1; touched(ntouch) = c; used(c) = 1
      enddo
      !> パディング（対角・近傍と重複しない連続した列番号、値 0）
      c = i
      do slot = jS + deg + 1, jE
        do
          c = mod(c, n_node) + 1
          if(used(c) == 0) exit
        enddo
        MAT%CSR%item(slot) = c
        ntouch = ntouch + 1; touched(ntouch) = c; used(c) = 1
      enddo
      !> used 配列をリセット
      do k = 1, ntouch
        used(touched(k)) = 0
      enddo
      !> 行内の列番号を昇順ソート
      call monolis_qsort_I_1d(MAT%CSR%item(jS:jE), 1, jE - jS + 1)
    enddo

    !> CSC 形式を再構築
    NZ = MAT%CSR%index(n_node + 1)
    call monolis_get_CSC_format(MAT%N, MAT%N, NZ, &
      & MAT%CSR%index, MAT%CSR%item, MAT%CSC%index, MAT%CSC%item, MAT%CSC%perm)
  end subroutine monolis_update_nonzero_pattern_fixed_width_main
end module mod_monolis_spmat_nonzero_pattern_util
