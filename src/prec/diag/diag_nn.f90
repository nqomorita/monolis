!> 対角スケーリング前処理（nxn ブロック）
module mod_monolis_precond_diag_nn
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc

  implicit none

contains

  !> @ingroup prec
  !> 前処理生成：対角スケーリング前処理（nxn ブロック、実数型）
  subroutine monolis_precond_diag_nn_setup_R(monoMAT, monoPREC)
    implicit none
    !> [in] 行列構造体
    type(monolis_mat), target, intent(in) :: monoMAT
    !> [in,out] 前処理構造体
    type(monolis_mat), target, intent(inout) :: monoPREC
    integer(kint) :: i, ii, j, jS, jE, in, k, l, N, NDOF, NDOF2
    real(kdouble), allocatable :: T(:), LU(:,:)
    integer(kint), pointer :: index(:), item(:)
    real(kdouble), pointer :: A(:), ALU(:)

    N =  monoMAT%N
    NDOF =  monoMAT%NDOF
    NDOF2 = NDOF*NDOF
    A => monoMAT%R%A
    index => monoMAT%CSR%index
    item => monoMAT%CSR%item

    call monolis_alloc_R_1d(T, NDOF)
    call monolis_alloc_R_2d(LU, NDOF, NDOF)
    call monolis_palloc_R_1d(monoPREC%R%D, NDOF2*N)
    ALU => monoPREC%R%D
    monoPREC%N = monoMAT%N

!$omp parallel default(none) &
!$omp & shared(A, ALU, index, item) &
!$omp & firstprivate(N, NDOF, NDOF2) &
!$omp & private(T, LU, i, j, k, jS, jE, in)
!$omp do
    do i = 1, N
      jS = index(i) + 1
      jE = index(i + 1)
      do ii = jS, jE
        in = item(ii)
        if(i == in)then
          do j = 1, NDOF
            do k = 1, NDOF
              LU(j,k) = A(NDOF2*(ii-1) + NDOF*(j-1) + k)
            enddo
          enddo

          do k = 1, NDOF
            if(LU(k,k) == 0.0d0) stop "** monolis error: zero diag in monolis_precond_diag_nn_setup"
            LU(k,k) = 1.0d0/LU(k,k)
            do l = k+1, NDOF
              LU(l,k) = LU(l,k)*LU(k,k)
              do j = k+1, NDOF
                T(j) = LU(l,j) - LU(l,k)*LU(k,j)
              enddo
              do j = k+1, NDOF
                LU(l,j) = T(j)
              enddo
            enddo
          enddo
          do j = 1, NDOF
            do k = 1, NDOF
              ALU(NDOF2*(i-1) + NDOF*(j-1) + k) = LU(j,k)
            enddo
          enddo
        endif
      enddo
    enddo
!$omp end do
!$omp end parallel

    deallocate(T)
    deallocate(LU)
  end subroutine monolis_precond_diag_nn_setup_R

  !> @ingroup prec
  !> 前処理生成：対角スケーリング前処理（nxn ブロック、複素数型）
  subroutine monolis_precond_diag_nn_setup_C(monoMAT, monoPREC)
    implicit none
    !> [in] 行列構造体
    type(monolis_mat), target, intent(in) :: monoMAT
    !> [in,out] 前処理構造体
    type(monolis_mat), target, intent(inout) :: monoPREC
    integer(kint) :: i, ii, j, jS, jE, in, k, l, N, NDOF, NDOF2
    complex(kdouble), allocatable :: T(:), LU(:,:)
    integer(kint), pointer :: index(:), item(:)
    complex(kdouble), pointer :: A(:), ALU(:)

    N =  monoMAT%N
    NDOF =  monoMAT%NDOF
    NDOF2 = NDOF*NDOF
    A => monoMAT%C%A
    index => monoMAT%CSR%index
    item => monoMAT%CSR%item

    call monolis_alloc_C_1d(T, NDOF)
    call monolis_alloc_C_2d(LU, NDOF, NDOF)
    call monolis_palloc_C_1d(monoPREC%C%D, NDOF2*N)
    ALU => monoPREC%C%D
    monoPREC%N = monoMAT%N

!$omp parallel default(none) &
!$omp & shared(A, ALU, index, item) &
!$omp & firstprivate(N, NDOF, NDOF2) &
!$omp & private(T, LU, i, j, k, jS, jE, in)
!$omp do
    do i = 1, N
      jS = index(i) + 1
      jE = index(i + 1)
      do ii = jS, jE
        in = item(ii)
        if(i == in)then
          do j = 1, NDOF
            do k = 1, NDOF
              LU(j,k) = A(NDOF2*(ii-1) + NDOF*(j-1) + k)
            enddo
          enddo

          do k = 1, NDOF
            if(LU(k,k) == 0.0d0) stop "** monolis error: zero diag in monolis_precond_diag_nn_setup"
            LU(k,k) = 1.0d0/LU(k,k)
            do l = k+1, NDOF
              LU(l,k) = LU(l,k)*LU(k,k)
              do j = k+1, NDOF
                T(j) = LU(l,j) - LU(l,k)*LU(k,j)
              enddo
              do j = k+1, NDOF
                LU(l,j) = T(j)
              enddo
            enddo
          enddo
          do j = 1, NDOF
            do k = 1, NDOF
              ALU(NDOF2*(i-1) + NDOF*(j-1) + k) = LU(j,k)
            enddo
          enddo
        endif
      enddo
    enddo
!$omp end do
!$omp end parallel

    deallocate(T)
    deallocate(LU)
  end subroutine monolis_precond_diag_nn_setup_C

  !> @ingroup prec
  !> 前処理適用：対角スケーリング前処理（nxn ブロック、実数型）
  subroutine monolis_precond_diag_nn_apply_R(monoMAT, monoPREC, X, Y)
    implicit none
    !> [in] 行列構造体
    type(monolis_mat), target, intent(in) :: monoMAT
    !> [in] 前処理構造体
    type(monolis_mat), target, intent(in) :: monoPREC
    integer(kint) :: i, j, k, N, NDOF, NDOF2
    real(kdouble) :: X(:), Y(:)
    real(kdouble), pointer :: ALU(:)
    real(kdouble), allocatable :: T(:)

    N =  monoMAT%N
    NDOF  = monoMAT%NDOF
    NDOF2 = NDOF*NDOF
    ALU => monoPREC%R%D

    call monolis_alloc_R_1d(T, NDOF)

!$omp parallel default(none) &
!$omp & shared(ALU, X, Y) &
!$omp & firstprivate(N, NDOF, NDOF2) &
!$omp & private(i, j, k, T)
!$omp do
    do i = 1, N
      do j = 1, NDOF
        T(j) = X(NDOF*(i-1) + j)
      enddo
      do j = 2, NDOF
        do k = 1, j-1
          T(j) = T(j) - ALU(NDOF2*(i-1) + NDOF*(j-1) + k)*T(k)
        enddo
      enddo
      do j = NDOF, 1, -1
        do k = NDOF, j+1, -1
          T(j) = T(j) - ALU(NDOF2*(i-1) + NDOF*(j-1) + k)*T(k)
        enddo
        T(j) = ALU(NDOF2*(i-1) + (NDOF+1)*(j-1) + 1)*T(j)
      enddo
      do k = 1, NDOF
        Y(NDOF*(i-1) + k) = T(k)
      enddo
    enddo
!$omp end do
!$omp end parallel

    deallocate(T)
  end subroutine monolis_precond_diag_nn_apply_R

  !> @ingroup prec
  !> 前処理適用：対角スケーリング前処理（nxn ブロック、複素数型）
  subroutine monolis_precond_diag_nn_apply_C(monoMAT, monoPREC, X, Y)
    implicit none
    !> [in] 行列構造体
    type(monolis_mat), target, intent(in) :: monoMAT
    !> [in] 前処理構造体
    type(monolis_mat), target, intent(in) :: monoPREC
    integer(kint) :: i, j, k, N, NDOF, NDOF2
    complex(kdouble) :: X(:), Y(:)
    complex(kdouble), pointer :: ALU(:)
    complex(kdouble), allocatable :: T(:)

    N =  monoMAT%N
    NDOF  = monoMAT%NDOF
    NDOF2 = NDOF*NDOF
    ALU => monoPREC%C%D

    call monolis_alloc_C_1d(T, NDOF)

!$omp parallel default(none) &
!$omp & shared(ALU, X, Y) &
!$omp & firstprivate(N, NDOF, NDOF2) &
!$omp & private(i, j, k, T)
!$omp do
    do i = 1, N
      do j = 1, NDOF
        T(j) = X(NDOF*(i-1) + j)
      enddo
      do j = 2, NDOF
        do k = 1, j-1
          T(j) = T(j) - ALU(NDOF2*(i-1) + NDOF*(j-1) + k)*T(k)
        enddo
      enddo
      do j = NDOF, 1, -1
        do k = NDOF, j+1, -1
          T(j) = T(j) - ALU(NDOF2*(i-1) + NDOF*(j-1) + k)*T(k)
        enddo
        T(j) = ALU(NDOF2*(i-1) + (NDOF+1)*(j-1) + 1)*T(j)
      enddo
      do k = 1, NDOF
        Y(NDOF*(i-1) + k) = T(k)
      enddo
    enddo
!$omp end do
!$omp end parallel

    deallocate(T)
  end subroutine monolis_precond_diag_nn_apply_C

  !> @ingroup prec
  !> 前処理初期化：対角スケーリング前処理（nxn ブロック、実数型）
  subroutine monolis_precond_diag_nn_clear_R(monoPREC)
    implicit none
    !> [in,out] 前処理構造体
    type(monolis_mat), intent(inout) :: monoPREC

    call monolis_pdealloc_R_1d(monoPREC%R%D)
  end subroutine monolis_precond_diag_nn_clear_R

  !> @ingroup prec
  !> 前処理初期化：対角スケーリング前処理（nxn ブロック、複素数型）
  subroutine monolis_precond_diag_nn_clear_C(monoPREC)
    implicit none
    !> [in,out] 前処理構造体
    type(monolis_mat), intent(inout) :: monoPREC

    call monolis_pdealloc_C_1d(monoPREC%C%D)
  end subroutine monolis_precond_diag_nn_clear_C
end module mod_monolis_precond_diag_nn
