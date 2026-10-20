!> CSR 形式から ELL 形式への変換モジュール
module mod_monolis_spmat_convert_ell
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc
#ifdef _OPENACC
  use openacc
#endif

  implicit none
  private

  public :: monolis_convert_CSR_to_ELL_R
  public :: monolis_convert_CSR_to_ELL_V_R
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
    integer(kint) :: i, j, k, l, in, jS, jE, slot
    integer(kint), pointer, contiguous :: csrIndex(:), csrItem(:), ellCol(:)
    real(kdouble), pointer, contiguous :: csrA(:), ellA(:)
#ifdef _OPENACC
    logical :: is_dev
#endif

    call monolis_std_debug_log_header("monolis_convert_CSR_to_ELL_R")

    N    = monoMAT%N
    NDOF = monoMAT%NDOF
    NDOF2 = NDOF*NDOF
    csrIndex => monoMAT%CSR%index
    csrItem  => monoMAT%CSR%item
    csrA     => monoMAT%R%A

    !# 1 行あたりの最大非ゼロブロック数を求める
    Nmaxcol = 0
!$omp parallel default(none) &
!$omp & shared(csrIndex, Nmaxcol) &
!$omp & firstprivate(N) &
!$omp & private(i)
!$omp do reduction(max:Nmaxcol)
    do i = 1, N
      if(csrIndex(i + 1) - csrIndex(i) > Nmaxcol) Nmaxcol = csrIndex(i + 1) - csrIndex(i)
    enddo
!$omp end do
!$omp end parallel

    monoMAT%ELL%Nmaxcol = Nmaxcol

    !# 列番号配列の確保（パディング込みで 0 初期化）
    call monolis_pdealloc_I_1d(monoMAT%ELL%col)
    call monolis_palloc_I_1d(monoMAT%ELL%col, N*Nmaxcol)

    !# 値配列の確保（パディング込みで 0 初期化）
    call monolis_pdealloc_R_1d(monoMAT%R%Aell)
    call monolis_palloc_R_1d(monoMAT%R%Aell, NDOF2*N*Nmaxcol)

    ellCol => monoMAT%ELL%col
    ellA   => monoMAT%R%Aell

#ifdef _OPENACC
    !# GPU 正本モード:行列値がデバイス常駐の場合、ELL 配列をデバイス上にも確保し、
    !# 詰め替えカーネルをデバイスでも実行する(ライトスルー)
    is_dev = acc_is_present(csrA)
    if(is_dev)then
      !$acc enter data create(ellCol, ellA)
      !# 組み立て中の非同期カーネルを待つ(詰め替えは組み立て完了後の値を参照する)
      !$acc wait(1)
    endif
#endif

    !# CSR 値を ELL 配置に詰める(パディングスロットは col = 0、値 = 0 とする)
!$omp parallel default(none) &
!$omp & shared(csrIndex, csrItem, csrA, ellCol, ellA) &
!$omp & firstprivate(N, NDOF, NDOF2, Nmaxcol) &
!$omp & private(i, j, k, l, in, jS, jE, slot)
!$omp do
    do i = 1, N
      jS = csrIndex(i) + 1
      jE = csrIndex(i + 1)
      do slot = 1, Nmaxcol
        j = jS + slot - 1
        if(j <= jE)then
          in = csrItem(j)
          ellCol((slot-1)*N + i) = in
          do k = 1, NDOF
            do l = 1, NDOF
              ellA((slot-1)*NDOF2*N + ((k-1)*NDOF + (l-1))*N + i) = &
                csrA(NDOF2*(j-1) + NDOF*(k-1) + l)
            enddo
          enddo
        else
          ellCol((slot-1)*N + i) = 0
          do k = 1, NDOF
            do l = 1, NDOF
              ellA((slot-1)*NDOF2*N + ((k-1)*NDOF + (l-1))*N + i) = 0.0d0
            enddo
          enddo
        endif
      enddo
    enddo
!$omp end do
!$omp end parallel

#ifdef _OPENACC
    if(is_dev)then
!$acc parallel loop gang vector private(j, k, l, in, jS, jE, slot) &
!$acc & present(csrIndex, csrItem, csrA, ellCol, ellA)
    do i = 1, N
      jS = csrIndex(i) + 1
      jE = csrIndex(i + 1)
!$acc loop seq
      do slot = 1, Nmaxcol
        j = jS + slot - 1
        if(j <= jE)then
          in = csrItem(j)
          ellCol((slot-1)*N + i) = in
          do k = 1, NDOF
            do l = 1, NDOF
              ellA((slot-1)*NDOF2*N + ((k-1)*NDOF + (l-1))*N + i) = &
                csrA(NDOF2*(j-1) + NDOF*(k-1) + l)
            enddo
          enddo
        else
          ellCol((slot-1)*N + i) = 0
          do k = 1, NDOF
            do l = 1, NDOF
              ellA((slot-1)*NDOF2*N + ((k-1)*NDOF + (l-1))*N + i) = 0.0d0
            enddo
          enddo
        endif
      enddo
    enddo
!$acc end parallel loop
    endif
#endif
  end subroutine monolis_convert_CSR_to_ELL_R

  !> @ingroup matrix
  !> CSR 形式から ELL 形式への変換（実数型、可変ブロック）
  !> @details ノードごとにブロックサイズが異なる可変ブロック（NDOF == -1）版。
  !>          列番号配列 col は均一ブロック版と同一で、スロット j のブロック列番号を col((j-1)*N + i) に格納する。
  !>          ブロック (スロット j, 行 i) のサイズは n_dof_list(i) × n_dof_list(col) と行ごとに変わるため、
  !>          値配列 Aell は各ブロックを行優先で連続配置し、ブロック先頭の Aell 内オフセット（0-based）を
  !>          prefix sum 配列 Vptr((j-1)*N + i) に保持する。均一版のような行 i 最内のコアレッシング配置は
  !>          ブロックサイズが行ごとに異なるため適用しない。
  subroutine monolis_convert_CSR_to_ELL_V_R(monoMAT)
    implicit none
    !> [in,out] 行列構造体
    type(monolis_mat), intent(inout) :: monoMAT
    integer(kint) :: N, Nmaxcol
    integer(kint) :: i, j, m, in, jS, jE, slot, nrow, p, n1, n2, kn, aoff, total
    integer(kint), pointer, contiguous :: csrIndex(:), csrItem(:), vptr(:)
    integer(kint), pointer, contiguous :: ndofList(:), ndofIndex2(:)
    real(kdouble), pointer, contiguous :: csrA(:), ellA(:)
#ifdef _OPENACC
    logical :: is_dev
#endif

    call monolis_std_debug_log_header("monolis_convert_CSR_to_ELL_V_R")

    N = monoMAT%N
    csrIndex   => monoMAT%CSR%index
    csrItem    => monoMAT%CSR%item
    csrA       => monoMAT%R%A
    ndofList   => monoMAT%n_dof_list
    ndofIndex2 => monoMAT%n_dof_index2

    !# 1 行あたりの最大非ゼロブロック数を求める
    Nmaxcol = 0
!$omp parallel default(none) &
!$omp & shared(monoMAT, Nmaxcol) &
!$omp & firstprivate(N) &
!$omp & private(i, nrow)
!$omp do reduction(max:Nmaxcol)
    do i = 1, N
      nrow = monoMAT%CSR%index(i + 1) - monoMAT%CSR%index(i)
      if(nrow > Nmaxcol) Nmaxcol = nrow
    enddo
!$omp end do
!$omp end parallel

    monoMAT%ELL%Nmaxcol = Nmaxcol

    !# 列番号配列の確保（パディング込みで 0 初期化）
    call monolis_pdealloc_I_1d(monoMAT%ELL%col)
    call monolis_palloc_I_1d(monoMAT%ELL%col, N*Nmaxcol)

    !# 列番号を ELL 配置に詰める
!$omp parallel default(none) &
!$omp & shared(monoMAT) &
!$omp & firstprivate(N) &
!$omp & private(i, j, in, jS, jE, slot)
!$omp do
    do i = 1, N
      jS = monoMAT%CSR%index(i) + 1
      jE = monoMAT%CSR%index(i + 1)
      slot = 0
      do j = jS, jE
        in = monoMAT%CSR%item(j)
        slot = slot + 1
        monoMAT%ELL%col((slot-1)*N + i) = in
      enddo
    enddo
!$omp end do
!$omp end parallel

    !# 各 (スロット, 行) ブロックの値配列内開始オフセットを prefix sum で確定
    call monolis_pdealloc_I_1d(monoMAT%ELL%Vptr)
    call monolis_palloc_I_1d(monoMAT%ELL%Vptr, N*Nmaxcol + 1)

!$omp parallel default(none) &
!$omp & shared(monoMAT) &
!$omp & firstprivate(N, Nmaxcol) &
!$omp & private(slot, i, p, in)
!$omp do collapse(2)
    do slot = 1, Nmaxcol
      do i = 1, N
        p = (slot-1)*N + i
        in = monoMAT%ELL%col(p)
        if(in >= 1)then
          monoMAT%ELL%Vptr(p + 1) = monoMAT%n_dof_list(i) * monoMAT%n_dof_list(in)
        endif
      enddo
    enddo
!$omp end do
!$omp end parallel
    do p = 1, N*Nmaxcol
      monoMAT%ELL%Vptr(p + 1) = monoMAT%ELL%Vptr(p + 1) + monoMAT%ELL%Vptr(p)
    enddo
    total = monoMAT%ELL%Vptr(N*Nmaxcol + 1)

    !# 値配列の確保（パディング込みで 0 初期化）
    call monolis_pdealloc_R_1d(monoMAT%R%Aell)
    call monolis_palloc_R_1d(monoMAT%R%Aell, total)

    vptr => monoMAT%ELL%Vptr
    ellA => monoMAT%R%Aell

#ifdef _OPENACC
    !# GPU 正本モード:行列値がデバイス常駐の場合、Vptr をデバイスへ転送し、
    !# 値配列をデバイス上にも確保して詰め替えカーネルをデバイスでも実行する(ライトスルー)
    is_dev = acc_is_present(csrA)
    if(is_dev)then
      !$acc enter data copyin(vptr)
      !$acc enter data create(ellA)
      !# 組み立て中の非同期カーネルを待つ(詰め替えは組み立て完了後の値を参照する)
      !$acc wait(1)
    endif
#endif

    !# CSR 値を ELL 配置（ブロック連続）に詰める
!$omp parallel default(none) &
!$omp & shared(csrIndex, csrItem, csrA, ellA, vptr, ndofList, ndofIndex2) &
!$omp & firstprivate(N) &
!$omp & private(i, j, m, in, jS, jE, slot, p, n1, n2, kn, aoff)
!$omp do
    do i = 1, N
      n1 = ndofList(i)
      jS = csrIndex(i) + 1
      jE = csrIndex(i + 1)
      slot = 0
      do j = jS, jE
        in = csrItem(j)
        slot = slot + 1
        p = (slot-1)*N + i
        n2 = ndofList(in)
        aoff = vptr(p)
        kn = ndofIndex2(j)
        do m = 1, n1*n2
          ellA(aoff + m) = csrA(kn + m)
        enddo
      enddo
    enddo
!$omp end do
!$omp end parallel

#ifdef _OPENACC
    if(is_dev)then
!$acc parallel loop gang vector private(j, m, in, jS, jE, slot, p, n1, n2, kn, aoff) &
!$acc & present(csrIndex, csrItem, csrA, ellA, vptr, ndofList, ndofIndex2)
    do i = 1, N
      n1 = ndofList(i)
      jS = csrIndex(i) + 1
      jE = csrIndex(i + 1)
      slot = 0
!$acc loop seq
      do j = jS, jE
        in = csrItem(j)
        slot = slot + 1
        p = (slot-1)*N + i
        n2 = ndofList(in)
        aoff = vptr(p)
        kn = ndofIndex2(j)
        do m = 1, n1*n2
          ellA(aoff + m) = csrA(kn + m)
        enddo
      enddo
    enddo
!$acc end parallel loop
    endif
#endif
  end subroutine monolis_convert_CSR_to_ELL_V_R

  !> @ingroup matrix
  !> ELL 形式の配列を解放（実数型）
  subroutine monolis_dealloc_ELL_R(monoMAT)
    implicit none
    !> [in,out] 行列構造体
    type(monolis_mat), intent(inout) :: monoMAT

    call monolis_std_debug_log_header("monolis_dealloc_ELL_R")

#ifdef _OPENACC
    !# GPU 正本モード:デバイス常駐領域があれば解放する
    if(associated(monoMAT%ELL%col))then
      !$acc exit data delete(monoMAT%ELL%col) if(acc_is_present(monoMAT%ELL%col))
    endif
    if(associated(monoMAT%ELL%Vptr))then
      !$acc exit data delete(monoMAT%ELL%Vptr) if(acc_is_present(monoMAT%ELL%Vptr))
    endif
    if(associated(monoMAT%R%Aell))then
      !$acc exit data delete(monoMAT%R%Aell) if(acc_is_present(monoMAT%R%Aell))
    endif
#endif

    monoMAT%ELL%Nmaxcol = 0
    call monolis_pdealloc_I_1d(monoMAT%ELL%col)
    call monolis_pdealloc_I_1d(monoMAT%ELL%Vptr)
    call monolis_pdealloc_R_1d(monoMAT%R%Aell)
  end subroutine monolis_dealloc_ELL_R
end module mod_monolis_spmat_convert_ell
