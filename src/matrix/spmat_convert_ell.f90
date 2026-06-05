!> CSR 形式から ELL 形式への変換モジュール
module mod_monolis_spmat_convert_ell
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc

  implicit none
  private

  public :: monolis_convert_CSR_to_ELL_R
  public :: monolis_dealloc_ELL_R

contains

  !> @ingroup matrix
  !> CSR 形式から ELL 形式への変換（実数型、nxn ブロック）
  !> @details 各行の非ゼロブロック数の最大値 Nmaxcol を求め、行ごとに最大 Nmaxcol 個の
  !>          非ゼロブロックをスロットへ詰める。スロット j のブロック列番号を col((j-1)*N + i)、
  !>          値を Aell((j-1)*NDOF2*N + ((k-1)*NDOF + (l-1))*N + i) の配置で格納する。
  !>          行 i を最内（ストライド1）にし、GPU スレッド間のメモリアクセスをコアレッシングさせる。
  !>          非ゼロブロックが Nmaxcol 未満の行のパディングスロットは col = 0、値 = 0 とする。
  subroutine monolis_convert_CSR_to_ELL_R(monoMAT)
    implicit none
    !> [in,out] 行列構造体
    type(monolis_mat), intent(inout) :: monoMAT
    integer(kint) :: N, NDOF, NDOF2, Nmaxcol
    integer(kint) :: i, j, k, l, in, jS, jE, slot, nrow

    call monolis_std_debug_log_header("monolis_convert_CSR_to_ELL_R")

    N    = monoMAT%N
    NDOF = monoMAT%NDOF
    NDOF2 = NDOF*NDOF

    !# 1 行あたりの最大非ゼロブロック数を求める
    Nmaxcol = 0
    do i = 1, N
      nrow = monoMAT%CSR%index(i + 1) - monoMAT%CSR%index(i)
      if(nrow > Nmaxcol) Nmaxcol = nrow
    enddo

    monoMAT%ELL%Nmaxcol = Nmaxcol

    !# 列番号配列の確保（パディング込みで 0 初期化）
    call monolis_pdealloc_I_1d(monoMAT%ELL%col)
    call monolis_palloc_I_1d(monoMAT%ELL%col, N*Nmaxcol)
    monoMAT%ELL%col = 0

    !# 値配列の確保（パディング込みで 0 初期化）
    call monolis_pdealloc_R_1d(monoMAT%R%Aell)
    call monolis_palloc_R_1d(monoMAT%R%Aell, NDOF2*N*Nmaxcol)
    monoMAT%R%Aell = 0.0d0

    !# CSR 値を ELL 配置に詰める
    do i = 1, N
      jS = monoMAT%CSR%index(i) + 1
      jE = monoMAT%CSR%index(i + 1)
      slot = 0
      do j = jS, jE
        in = monoMAT%CSR%item(j)
        slot = slot + 1
        monoMAT%ELL%col((slot-1)*N + i) = in
        do k = 1, NDOF
          do l = 1, NDOF
            monoMAT%R%Aell((slot-1)*NDOF2*N + ((k-1)*NDOF + (l-1))*N + i) = &
              monoMAT%R%A(NDOF2*(j-1) + NDOF*(k-1) + l)
          enddo
        enddo
      enddo
    enddo
  end subroutine monolis_convert_CSR_to_ELL_R

  !> @ingroup matrix
  !> ELL 形式の配列を解放（実数型）
  subroutine monolis_dealloc_ELL_R(monoMAT)
    implicit none
    !> [in,out] 行列構造体
    type(monolis_mat), intent(inout) :: monoMAT

    call monolis_std_debug_log_header("monolis_dealloc_ELL_R")

    monoMAT%ELL%Nmaxcol = 0
    call monolis_pdealloc_I_1d(monoMAT%ELL%col)
    call monolis_pdealloc_R_1d(monoMAT%R%Aell)
  end subroutine monolis_dealloc_ELL_R
end module mod_monolis_spmat_convert_ell
