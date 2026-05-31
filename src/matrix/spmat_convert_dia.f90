!> CSR 形式から DIA 形式への変換モジュール
module mod_monolis_spmat_convert_dia
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc

  implicit none
  private

  public :: monolis_convert_CSR_to_DIA_R
  public :: monolis_dealloc_DIA_R

contains

  !> @ingroup matrix
  !> CSR 形式から DIA 形式への変換（実数型、nxn ブロック）
  !> @details 各非ゼロブロック (i, item(j)) のブロック列オフセット item(j) - i を
  !>          全行から収集・ソート・一意化して対角線本数 Ndiag と offset(:) を確定する。
  !>          値は Adia(NDOF2*((d-1)*N + (i-1)) + NDOF*(k-1) + l) の column-major 配置で詰める。
  subroutine monolis_convert_CSR_to_DIA_R(monoMAT)
    implicit none
    !> [in,out] 行列構造体
    type(monolis_mat), intent(inout) :: monoMAT
    integer(kint) :: N, NP, NDOF, NDOF2, nz, Ndiag
    integer(kint) :: i, j, k, l, in, jS, jE, d, idx, newlen
    integer(kint), allocatable :: off_all(:)
    integer(kint), allocatable :: diag_of_offset(:)

    call monolis_std_debug_log_header("monolis_convert_CSR_to_DIA_R")

    N    = monoMAT%N
    NP   = monoMAT%NP
    NDOF = monoMAT%NDOF
    NDOF2 = NDOF*NDOF

    nz = monoMAT%CSR%index(N + 1)

    !# 全非ゼロのブロック列オフセットを収集
    call monolis_alloc_I_1d(off_all, nz)
    idx = 0
    do i = 1, N
      jS = monoMAT%CSR%index(i) + 1
      jE = monoMAT%CSR%index(i + 1)
      do j = jS, jE
        in = monoMAT%CSR%item(j)
        idx = idx + 1
        off_all(idx) = in - i
      enddo
    enddo

    !# ソート・一意化して対角線本数と offset を確定
    if(nz > 0)then
      call monolis_qsort_I_1d(off_all, 1, nz)
      call monolis_get_uniq_array_I(off_all, nz, newlen)
    else
      newlen = 0
    endif
    Ndiag = newlen

    monoMAT%DIA%Ndiag = Ndiag
    call monolis_pdealloc_I_1d(monoMAT%DIA%offset)
    call monolis_palloc_I_1d(monoMAT%DIA%offset, Ndiag)
    do d = 1, Ndiag
      monoMAT%DIA%offset(d) = off_all(d)
    enddo

    !# 各オフセットの対角番号を引く逆引き表（オフセット値域は -(NP-1)..(NP-1)）
    call monolis_alloc_I_1d(diag_of_offset, 2*NP - 1)
    do d = 1, Ndiag
      diag_of_offset(monoMAT%DIA%offset(d) + NP) = d
    enddo

    !# 値配列の確保（パディング込みで 0 初期化）
    call monolis_pdealloc_R_1d(monoMAT%R%Adia)
    call monolis_palloc_R_1d(monoMAT%R%Adia, NDOF2*N*Ndiag)
    monoMAT%R%Adia = 0.0d0

    !# CSR 値を DIA 配置に詰める
    do i = 1, N
      jS = monoMAT%CSR%index(i) + 1
      jE = monoMAT%CSR%index(i + 1)
      do j = jS, jE
        in = monoMAT%CSR%item(j)
        d = diag_of_offset(in - i + NP)
        do k = 1, NDOF
          do l = 1, NDOF
            monoMAT%R%Adia(NDOF2*((d-1)*N + (i-1)) + NDOF*(k-1) + l) = &
              monoMAT%R%A(NDOF2*(j-1) + NDOF*(k-1) + l)
          enddo
        enddo
      enddo
    enddo

    call monolis_dealloc_I_1d(off_all)
    call monolis_dealloc_I_1d(diag_of_offset)
  end subroutine monolis_convert_CSR_to_DIA_R

  !> @ingroup matrix
  !> DIA 形式の配列を解放（実数型）
  subroutine monolis_dealloc_DIA_R(monoMAT)
    implicit none
    !> [in,out] 行列構造体
    type(monolis_mat), intent(inout) :: monoMAT

    call monolis_std_debug_log_header("monolis_dealloc_DIA_R")

    monoMAT%DIA%Ndiag = 0
    call monolis_pdealloc_I_1d(monoMAT%DIA%offset)
    call monolis_pdealloc_R_1d(monoMAT%R%Adia)
  end subroutine monolis_dealloc_DIA_R
end module mod_monolis_spmat_convert_dia
