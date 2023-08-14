!> 対角スケーリング前処理（3x3 ブロック）
module mod_monolis_precond_diag_33
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc

  implicit none

contains

  !> @ingroup prec
  !> 前処理生成：対角スケーリング前処理（3x3 ブロック）
  subroutine monolis_precond_diag_33_setup_R(monoMAT, monoPREC)
    implicit none
    !> [in] 行列構造体
    type(monolis_mat), target, intent(in) :: monoMAT
    !> [in,out] 前処理構造体
    type(monolis_mat), target, intent(inout) :: monoPREC
    integer(kint) :: i, j, jS, jE, in, k, l, N
    real(kdouble) :: T(3,3), P(3)
    integer(kint), pointer :: index(:), item(:)
    real(kdouble), pointer :: A(:), ALU(:)

    N =  monoMAT%N
    A => monoMAT%R%A
    index => monoMAT%CSR%index
    item => monoMAT%CSR%item

    call monolis_palloc_R_1d(monoPREC%R%D, 9*N)
    ALU => monoPREC%R%D

!$omp parallel default(none) &
!$omp & shared(A, ALU, index, item) &
!$omp & firstprivate(N) &
!$omp & private(T, P, i, j, k, jS, jE, in)
!$omp do
    do i = 1, N
      jS = index(i) + 1
      jE = index(i + 1)
      do j = jS, jE
        in = item(j)
        if(i == in)then
          ALU(9*i-8) = A(9*j-8)
          ALU(9*i-7) = A(9*j-7)
          ALU(9*i-6) = A(9*j-6)
          ALU(9*i-5) = A(9*j-5)
          ALU(9*i-4) = A(9*j-4)
          ALU(9*i-3) = A(9*j-3)
          ALU(9*i-2) = A(9*j-2)
          ALU(9*i-1) = A(9*j-1)
          ALU(9*i  ) = A(9*j  )
        endif
      enddo
    enddo
!$omp end do
!$omp do
    do l = 1, N
      T(1,1) = ALU(9*l-8)
      T(1,2) = ALU(9*l-7)
      T(1,3) = ALU(9*l-6)
      T(2,1) = ALU(9*l-5)
      T(2,2) = ALU(9*l-4)
      T(2,3) = ALU(9*l-3)
      T(3,1) = ALU(9*l-2)
      T(3,2) = ALU(9*l-1)
      T(3,3) = ALU(9*l  )
      do k = 1, 3
        T(k,k) = 1.0d0/T(k,k)
        do i = k+1, 3
          T(i,k) = T(i,k) * T(k,k)
          do j = k+1, 3
            P(j) = T(i,j) - T(i,k)*T(k,j)
          enddo
          do j = k+1, 3
            T(i,j) = P(j)
          enddo
        enddo
      enddo
      ALU(9*l-8) = T(1,1)
      ALU(9*l-7) = T(1,2)
      ALU(9*l-6) = T(1,3)
      ALU(9*l-5) = T(2,1)
      ALU(9*l-4) = T(2,2)
      ALU(9*l-3) = T(2,3)
      ALU(9*l-2) = T(3,1)
      ALU(9*l-1) = T(3,2)
      ALU(9*l  ) = T(3,3)
    enddo
!$omp end do
!$omp end parallel
  end subroutine monolis_precond_diag_33_setup_R

  !> @ingroup prec
  !> 前処理適用：対角スケーリング前処理（3x3 ブロック）
  subroutine monolis_precond_diag_33_apply_R(monoMAT, monoPREC, X, Y)
    implicit none
    !> [in] 行列構造体
    type(monolis_mat), intent(in) :: monoMAT
    !> [in,out] 前処理構造体
    type(monolis_mat), target, intent(inout) :: monoPREC
    integer(kint) :: i, N
    real(kdouble) :: X1, X2, X3
    real(kdouble) :: X(:), Y(:)
    real(kdouble), pointer :: ALU(:)

    N =  monoMAT%N
    ALU => monoPREC%R%D

!$omp parallel default(none) &
!$omp & shared(monoMAT, ALU, X, Y) &
!$omp & private(i, X1, X2, X3)
!$omp do
    do i = 1, N
      X1 = X(3*i-2)
      X2 = X(3*i-1)
      X3 = X(3*i  )
      X2 = X2 - ALU(9*i-5)*X1
      X3 = X3 - ALU(9*i-2)*X1 - ALU(9*i-1)*X2
      X3 = ALU(9*i  )* X3
      X2 = ALU(9*i-4)*(X2 - ALU(9*i-3)*X3)
      X1 = ALU(9*i-8)*(X1 - ALU(9*i-6)*X3 - ALU(9*i-7)*X2)
      Y(3*i-2) = X1
      Y(3*i-1) = X2
      Y(3*i  ) = X3
    enddo
!$omp end do
!$omp end parallel
  end subroutine monolis_precond_diag_33_apply_R

  !> @ingroup prec
  !> 前処理初期化：対角スケーリング前処理（3x3 ブロック）
  subroutine monolis_precond_diag_33_clear_R(monoPREC)
    implicit none
    !> [in,out] 前処理構造体
    type(monolis_mat), intent(inout) :: monoPREC

    call monolis_pdealloc_R_1d(monoPREC%R%D)
  end subroutine monolis_precond_diag_33_clear_R

end module mod_monolis_precond_diag_33