module mod_monolis_matmat
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc
  implicit none

contains

  !> @ingroup linalg
  !> ベクトル内積
  subroutine monolis_matmat_local(monoCOM, A, B, C, is_CSC_A, is_CSC_B, tspmv, tcomm)
    implicit none
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in] 行列構造体
    type(monolis_mat), intent(in) :: A
    !> [in] 行列構造体
    type(monolis_mat), intent(in) :: B
    !> [in] 行列構造体
    type(monolis_mat), intent(inout) :: C
    !>
    logical :: is_CSC_A
    !>
    logical :: is_CSC_B
    !> [in,out] 計算時間
    real(kdouble), intent(inout) :: tspmv
    !> [in,out] 通信時間
    real(kdouble), intent(inout) :: tcomm
    real(kdouble) :: t1, t2

#ifdef DEBUG
    call monolis_std_debug_log_header("monolis_matmat_local")
#endif
    t1 = monolis_get_time()

    !call monolis_update_pre_R(monoCOM, monoMAT%NDOF, X, tcomm)

    call monolis_matmat_nn(monoCOM, A, B, C, A%NDOF)

    t2 = monolis_get_time()
    tspmv = tspmv + t2 - t1
  end subroutine monolis_matmat_local

  !> @ingroup linalg
  !> ベクトル内積
  subroutine monolis_matmat_nn(monoCOM, A, B, C, NDOF)
    implicit none
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in] 行列構造体
    type(monolis_mat), intent(in) :: A
    !> [in] 行列構造体
    type(monolis_mat), intent(in) :: B
    !> [in] 行列構造体
    type(monolis_mat), intent(inout) :: C
    !> [in] 計算点が持つ自由度
    integer(kint), intent(in) :: NDOF

    call monolis_matmat_symbolic(monoCOM, A, B, C, NDOF)
    call monolis_matmat_value_nn(monoCOM, A, B, C, NDOF)
  end subroutine monolis_matmat_nn

  subroutine monolis_matmat_symbolic(monoCOM, A, B, C, NDOF)
    implicit none
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in] 行列構造体
    type(monolis_mat), intent(in) :: A
    !> [in] 行列構造体
    type(monolis_mat), intent(in) :: B
    !> [in] 行列構造体
    type(monolis_mat), intent(inout) :: C
    !> [in] 計算点が持つ自由度
    integer(kint), intent(in) :: NDOF
    integer(kint) :: i, j, k, jA, jB, jS, jE, col_A, col_B, nz, NP
    integer(kint), allocatable :: row_C(:)
    logical, allocatable :: is_nonzero(:)

#ifdef DEBUG
    call monolis_std_debug_log_header("monolis_matmat_symbolic")
#endif

    NP = A%NP

    !# 行列 C の構造を初期化
    C%N = A%N
    C%NP = A%NP
    C%NDOF = NDOF

    call monolis_alloc_L_1d(is_nonzero, B%NP)
    call monolis_alloc_I_1d(row_C, B%NP)
    call monolis_palloc_I_1d(C%CSR%index, NP + 1)

    !# 各行について非零要素数をカウント
    do i = 1, NP
      nz = 0

      !# A の i 行の各非零要素について
      jS = A%CSR%index(i) + 1
      jE = A%CSR%index(i + 1)
      do jA = jS, jE
        col_A = A%CSR%item(jA)

        !# B の col_A 行の各非零要素について
        jS = B%CSR%index(col_A) + 1
        jE = B%CSR%index(col_A + 1)
        do jB = jS, jE
          col_B = B%CSR%item(jB)

          !# C(i, col_B) が非零になる
          if(.not. is_nonzero(col_B))then
            is_nonzero(col_B) = .true.
            nz = nz + 1
            row_C(nz) = col_B
          endif
        enddo
      enddo

      !# 訪問した列のみフラグを解除（全リセット O(NP) の回避）
      do j = 1, nz
        is_nonzero(row_C(j)) = .false.
      enddo

      C%CSR%index(i + 1) = C%CSR%index(i) + nz
    enddo

    !# C%CSR%item の確保と設定
    nz = C%CSR%index(NP + 1)
    call monolis_palloc_I_1d(C%CSR%item, nz)

    !# 各行の列インデックスを設定
    do i = 1, NP
      nz = 0

      !# A の i 行の各非零要素について
      jS = A%CSR%index(i) + 1
      jE = A%CSR%index(i + 1)
      do jA = jS, jE
        col_A = A%CSR%item(jA)

        !# B の col_A 行の各非零要素について
        jS = B%CSR%index(col_A) + 1
        jE = B%CSR%index(col_A + 1)
        do jB = jS, jE
          col_B = B%CSR%item(jB)

          !# C(i, col_B) が非零になる
          if(.not. is_nonzero(col_B))then
            is_nonzero(col_B) = .true.
            nz = nz + 1
            row_C(nz) = col_B
          endif
        enddo
      enddo

      !# 訪問した列のみフラグを解除（全リセット O(NP) の回避）
      do j = 1, nz
        is_nonzero(row_C(j)) = .false.
      enddo

      ! この行の列インデックスをソートして格納
      call monolis_qsort_I_1d(row_C(1:nz), 1, nz)
      jS = C%CSR%index(i) + 1
      jE = C%CSR%index(i + 1)
      do j = jS, jE
        C%CSR%item(j) = row_C(j - jS + 1)
      enddo
    enddo

    call monolis_dealloc_L_1d(is_nonzero)
    call monolis_dealloc_I_1d(row_C)
  end subroutine monolis_matmat_symbolic

  subroutine monolis_matmat_value_nn(monoCOM, A, B, C, NDOF)
    implicit none
    !> [in] 通信テーブル構造体
    type(monolis_com), intent(in) :: monoCOM
    !> [in] 行列構造体
    type(monolis_mat), intent(in) :: A
    !> [in] 行列構造体
    type(monolis_mat), intent(in) :: B
    !> [in] 行列構造体
    type(monolis_mat), intent(inout) :: C
    !> [in] 計算点が持つ自由度
    integer(kint), intent(in) :: NDOF
    integer(kint) :: i, jC, jA, jB, jS, jE, col_A, col_B, NP, nz, NDOF2
    integer(kint) :: k, l, m, kn_A, kn_B, kn_C
    integer(kint), allocatable :: col_to_pos(:)

#ifdef DEBUG
    call monolis_std_debug_log_header("monolis_matmat_value_nn")
#endif

    NP = A%NP
    NDOF2 = NDOF * NDOF
    nz = C%CSR%index(NP + 1)
    call monolis_palloc_R_1d(C%R%A, NDOF2 * nz)

    !# Gustavson 法：行ごとに C の列番号 → 格納位置マップを構築し、
    !# A(i,:) × B(col_A,:) を C へ直接足し込む（探索なし、O(flops)）
!$omp parallel default(none) &
!$omp & shared(A, B, C) &
!$omp & firstprivate(NP, NDOF, NDOF2) &
!$omp & private(i, jC, jA, jB, jS, jE, col_A, col_B, k, l, m, kn_A, kn_B, kn_C, col_to_pos)
    call monolis_alloc_I_1d(col_to_pos, B%NP)
!$omp do
    do i = 1, NP
      !# C の i 行の列番号 → 格納位置マップを構築
      jS = C%CSR%index(i) + 1
      jE = C%CSR%index(i + 1)
      do jC = jS, jE
        col_to_pos(C%CSR%item(jC)) = jC
      enddo

      !# A の i 行の各非零要素
      jS = A%CSR%index(i) + 1
      jE = A%CSR%index(i + 1)
      do jA = jS, jE
        col_A = A%CSR%item(jA)
        kn_A = NDOF2 * (jA - 1)

        !# B の col_A 行の各非零要素を C(i, col_B) に足し込む
        jS = B%CSR%index(col_A) + 1
        jE = B%CSR%index(col_A + 1)
        do jB = jS, jE
          col_B = B%CSR%item(jB)
          jC = col_to_pos(col_B)
          kn_B = NDOF2 * (jB - 1)
          kn_C = NDOF2 * (jC - 1)

          !# ブロック行列の積: C(i,col_B) += A(i,col_A) * B(col_A,col_B)
          do k = 1, NDOF
            do l = 1, NDOF
              do m = 1, NDOF
                C%R%A(kn_C + NDOF*(k-1) + l) = C%R%A(kn_C + NDOF*(k-1) + l) &
                  & + A%R%A(kn_A + NDOF*(k-1) + m) * B%R%A(kn_B + NDOF*(m-1) + l)
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
!$omp end do
    call monolis_dealloc_I_1d(col_to_pos)
!$omp end parallel
  end subroutine monolis_matmat_value_nn
end module mod_monolis_matmat