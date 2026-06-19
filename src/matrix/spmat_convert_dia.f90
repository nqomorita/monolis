!> CSR 形式から DIA 形式への変換モジュール
module mod_monolis_spmat_convert_dia
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc

  implicit none
  private

  public :: monolis_convert_CSR_to_DIA_R
  public :: monolis_convert_CSR_to_DIA_V_R
  public :: monolis_dealloc_DIA_R

contains

  !> @ingroup matrix
  !> CSR 形式から DIA 形式への変換（実数型、nxn ブロック）
  !> @details 各非ゼロブロック (i, item(j)) のブロック列オフセット item(j) - i に存在フラグを
  !>          立て、フラグ配列を昇順に走査して対角線本数 Ndiag と offset(:) を確定する。
  !>          値は Adia((d-1)*NDOF2*N + ((k-1)*NDOF + (l-1))*N + i) の配置で詰める。
  !>          行 i を最内（ストライド1）にし、GPU スレッド間のメモリアクセスをコアレッシングさせる。
  subroutine monolis_convert_CSR_to_DIA_R(monoMAT)
    implicit none
    !> [in,out] 行列構造体
    type(monolis_mat), intent(inout) :: monoMAT
    integer(kint) :: N, NP, NDOF, NDOF2, Ndiag
    integer(kint) :: i, j, k, l, in, jS, jE, d, idx
    integer(kint), allocatable :: diag_of_offset(:)

    call monolis_std_debug_log_header("monolis_convert_CSR_to_DIA_R")

    N    = monoMAT%N
    NP   = monoMAT%NP
    NDOF = monoMAT%NDOF
    NDOF2 = NDOF*NDOF

    !# 各オフセットの存在フラグを立てる（オフセット値域は -(NP-1)..(NP-1)）
    call monolis_alloc_I_1d(diag_of_offset, 2*NP - 1)
    do i = 1, N
      jS = monoMAT%CSR%index(i) + 1
      jE = monoMAT%CSR%index(i + 1)
      do j = jS, jE
        in = monoMAT%CSR%item(j)
        diag_of_offset(in - i + NP) = 1
      enddo
    enddo

    !# フラグ配列を昇順に走査して対角線本数と offset を確定
    Ndiag = 0
    do idx = 1, 2*NP - 1
      if(diag_of_offset(idx) /= 0) Ndiag = Ndiag + 1
    enddo

    monoMAT%DIA%Ndiag = Ndiag
    call monolis_pdealloc_I_1d(monoMAT%DIA%offset)
    call monolis_palloc_I_1d(monoMAT%DIA%offset, Ndiag)

    d = 0
    do idx = 1, 2*NP - 1
      if(diag_of_offset(idx) /= 0)then
        d = d + 1
        monoMAT%DIA%offset(d) = idx - NP
        diag_of_offset(idx) = d
      endif
    enddo

    !# 値配列の確保（パディング込みで 0 初期化）
    call monolis_pdealloc_R_1d(monoMAT%R%Adia)
    call monolis_palloc_R_1d(monoMAT%R%Adia, NDOF2*N*Ndiag)

    !# CSR 値を DIA 配置に詰める
    do i = 1, N
      jS = monoMAT%CSR%index(i) + 1
      jE = monoMAT%CSR%index(i + 1)
      do j = jS, jE
        in = monoMAT%CSR%item(j)
        d = diag_of_offset(in - i + NP)
        do k = 1, NDOF
          do l = 1, NDOF
            monoMAT%R%Adia((d-1)*NDOF2*N + ((k-1)*NDOF + (l-1))*N + i) = &
              monoMAT%R%A(NDOF2*(j-1) + NDOF*(k-1) + l)
          enddo
        enddo
      enddo
    enddo

    call monolis_dealloc_I_1d(diag_of_offset)
  end subroutine monolis_convert_CSR_to_DIA_R

  !> @ingroup matrix
  !> CSR 形式から DIA 形式への変換（実数型、可変ブロック）
  !> @details ノードごとにブロックサイズが異なる可変ブロック（NDOF == -1）版。
  !>          対角本数 Ndiag と offset(:) の確定は均一ブロック版と同一。
  !>          ブロック (対角 d, 行 i) のサイズは n_dof_list(i) × n_dof_list(i+offset(d)) と行ごとに変わるため、
  !>          値配列 Adia は各ブロックを行優先で連続配置し、ブロック先頭の Adia 内オフセット（0-based）を
  !>          prefix sum 配列 Vptr((d-1)*N + i) に保持する。範囲内だが構造的に非ゼロでない (d, i) には
  !>          0 埋めの領域を確保しておき、カーネル側は範囲内の全対角を一律に処理する。
  subroutine monolis_convert_CSR_to_DIA_V_R(monoMAT)
    implicit none
    !> [in,out] 行列構造体
    type(monolis_mat), intent(inout) :: monoMAT
    integer(kint) :: N, NP, Ndiag
    integer(kint) :: i, j, m, in, jS, jE, d, idx, p, jcol, n1, n2, kn, aoff, total
    integer(kint), allocatable :: diag_of_offset(:)

    call monolis_std_debug_log_header("monolis_convert_CSR_to_DIA_V_R")

    N  = monoMAT%N
    NP = monoMAT%NP

    !# 各オフセットの存在フラグを立てる（オフセット値域は -(NP-1)..(NP-1)）
    call monolis_alloc_I_1d(diag_of_offset, 2*NP - 1)
    do i = 1, N
      jS = monoMAT%CSR%index(i) + 1
      jE = monoMAT%CSR%index(i + 1)
      do j = jS, jE
        in = monoMAT%CSR%item(j)
        diag_of_offset(in - i + NP) = 1
      enddo
    enddo

    !# フラグ配列を昇順に走査して対角線本数と offset を確定
    Ndiag = 0
    do idx = 1, 2*NP - 1
      if(diag_of_offset(idx) /= 0) Ndiag = Ndiag + 1
    enddo

    monoMAT%DIA%Ndiag = Ndiag
    call monolis_pdealloc_I_1d(monoMAT%DIA%offset)
    call monolis_palloc_I_1d(monoMAT%DIA%offset, Ndiag)

    !# offset を確定し、同じ配列を対角番号の逆引き表として再利用
    d = 0
    do idx = 1, 2*NP - 1
      if(diag_of_offset(idx) /= 0)then
        d = d + 1
        monoMAT%DIA%offset(d) = idx - NP
        diag_of_offset(idx) = d
      endif
    enddo

    !# 各 (対角, 行) ブロックの値配列内開始オフセットを prefix sum で確定
    call monolis_pdealloc_I_1d(monoMAT%DIA%Vptr)
    call monolis_palloc_I_1d(monoMAT%DIA%Vptr, N*Ndiag + 1)

    do d = 1, Ndiag
      do i = 1, N
        p = (d-1)*N + i
        jcol = i + monoMAT%DIA%offset(d)
        if(jcol >= 1 .and. jcol <= NP)then
          monoMAT%DIA%Vptr(p + 1) = monoMAT%n_dof_list(i) * monoMAT%n_dof_list(jcol)
        endif
      enddo
    enddo
    do p = 1, N*Ndiag
      monoMAT%DIA%Vptr(p + 1) = monoMAT%DIA%Vptr(p + 1) + monoMAT%DIA%Vptr(p)
    enddo
    total = monoMAT%DIA%Vptr(N*Ndiag + 1)

    !# 値配列の確保（パディング込みで 0 初期化）
    call monolis_pdealloc_R_1d(monoMAT%R%Adia)
    call monolis_palloc_R_1d(monoMAT%R%Adia, total)

    !# CSR 値を DIA 配置（ブロック連続）に詰める
    do i = 1, N
      n1 = monoMAT%n_dof_list(i)
      jS = monoMAT%CSR%index(i) + 1
      jE = monoMAT%CSR%index(i + 1)
      do j = jS, jE
        in = monoMAT%CSR%item(j)
        d = diag_of_offset(in - i + NP)
        p = (d-1)*N + i
        n2 = monoMAT%n_dof_list(in)
        aoff = monoMAT%DIA%Vptr(p)
        kn = monoMAT%n_dof_index2(j)
        do m = 1, n1*n2
          monoMAT%R%Adia(aoff + m) = monoMAT%R%A(kn + m)
        enddo
      enddo
    enddo

    call monolis_dealloc_I_1d(diag_of_offset)
  end subroutine monolis_convert_CSR_to_DIA_V_R

  !> @ingroup matrix
  !> DIA 形式の配列を解放（実数型）
  subroutine monolis_dealloc_DIA_R(monoMAT)
    implicit none
    !> [in,out] 行列構造体
    type(monolis_mat), intent(inout) :: monoMAT

    call monolis_std_debug_log_header("monolis_dealloc_DIA_R")

    monoMAT%DIA%Ndiag = 0
    call monolis_pdealloc_I_1d(monoMAT%DIA%offset)
    call monolis_pdealloc_I_1d(monoMAT%DIA%Vptr)
    call monolis_pdealloc_R_1d(monoMAT%R%Adia)
  end subroutine monolis_dealloc_DIA_R
end module mod_monolis_spmat_convert_dia
