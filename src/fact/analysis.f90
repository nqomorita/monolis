module mod_monolis_fact_analysis
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc
  use mod_monolis_pord_const
  use mod_monolis_pord_types
  use mod_monolis_pord_graph
  use mod_monolis_pord_tree
  use mod_monolis_pord_ordering

  implicit none

contains

  !> @ingroup fact
  !> 解析フェーズ：PORD でフィル削減順序を計算し、monoLU に格納する
  subroutine monolis_fact_analysis(monoMAT, lu)
    implicit none
    !> [in] 行列構造体
    type(monolis_mat),    intent(in)    :: monoMAT
    !> [in,out] LU 分解構造体（perm / iperm が出力される）
    type(monolis_mat_lu), intent(inout) :: lu

    type(graph_t)    :: G
    type(elimtree_t) :: T
    integer(kint) :: n, nz_off, i, j, k, jS, jE, cnt, pos
    integer(ip)   :: nvtx_ip
    integer(ip), allocatable :: perm_ip(:)
    integer(ip)   :: opts_dummy(1)

    n = monoMAT%N
    lu%N = n

    if (n <= 0) return

    ! -----------------------------------------------------------------------
    ! 1. 非対角の非ゼロ要素数をカウント
    !    （monolis の CSR は構造的対称性を仮定）
    ! -----------------------------------------------------------------------
    nz_off = 0
    do i = 1, n
      jS = monoMAT%CSR%index(i) + 1
      jE = monoMAT%CSR%index(i + 1)
      do k = jS, jE
        if (monoMAT%CSR%item(k) /= i) nz_off = nz_off + 1
      end do
    end do

    ! -----------------------------------------------------------------------
    ! 2. PORD グラフ構築（0-based の xadj/adjncy）
    ! -----------------------------------------------------------------------
    nvtx_ip = int(n, ip)
    call newGraph(G, nvtx_ip, int(nz_off, ip))

    G%xadj(0) = 0
    do i = 1, n
      jS = monoMAT%CSR%index(i) + 1
      jE = monoMAT%CSR%index(i + 1)
      cnt = 0
      do k = jS, jE
        if (monoMAT%CSR%item(k) /= i) cnt = cnt + 1
      end do
      G%xadj(i) = G%xadj(i - 1) + int(cnt, ip)
    end do

    do i = 1, n
      jS = monoMAT%CSR%index(i) + 1
      jE = monoMAT%CSR%index(i + 1)
      pos = int(G%xadj(i - 1), kint)
      do k = jS, jE
        j = monoMAT%CSR%item(k)
        if (j /= i) then
          G%adjncy(pos) = int(j - 1, ip)   ! 0-based
          pos = pos + 1
        end if
      end do
    end do

    ! -----------------------------------------------------------------------
    ! 3. PORD によるフィル削減順序計算（既定オプション）
    ! -----------------------------------------------------------------------
    opts_dummy = 0
    call monolis_pord_ordering(G, opts_dummy, .true., T)

    ! -----------------------------------------------------------------------
    ! 4. perm / iperm 抽出（1-based）
    ! -----------------------------------------------------------------------
    call monolis_dealloc_I_1d(lu%perm)
    call monolis_dealloc_I_1d(lu%iperm)
    call monolis_alloc_I_1d(lu%perm,  n)
    call monolis_alloc_I_1d(lu%iperm, n)

    allocate(perm_ip(n))
    call monolis_pord_perm_from_elimtree(T, nvtx_ip, perm_ip)
    do i = 1, n
      lu%perm(i) = int(perm_ip(i), kint)
    end do
    deallocate(perm_ip)

    do i = 1, n
      lu%iperm(lu%perm(i)) = i
    end do

    ! -----------------------------------------------------------------------
    ! 5. 後片付け
    ! -----------------------------------------------------------------------
    call freeElimTree(T)
    call freeGraph(G)
  end subroutine monolis_fact_analysis

end module mod_monolis_fact_analysis
