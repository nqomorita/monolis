!> 対角スケーリング前処理（任意ブロック）
module mod_monolis_precond_diag_V
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc

  implicit none

contains

  !> @ingroup prec
  !> 前処理生成：対角スケーリング前処理（任意ブロック、実数型）
  subroutine monolis_precond_diag_V_setup_R(monoMAT, monoPREC)
    implicit none
    !> [in] 行列構造体
    type(monolis_mat), target, intent(in) :: monoMAT
    !> [in,out] 前処理構造体
    type(monolis_mat), target, intent(inout) :: monoPREC
    integer(kint) :: i, ii, j, jS, jE, in, jn, k, l, N, NZ, NDOF_MAX
    integer(kint) :: NDOF, kn
    real(kdouble), allocatable :: T(:), LU(:,:)
    integer(kint), pointer :: index(:), item(:)
    real(kdouble), pointer :: A(:), ALU(:)

    N =  monoMAT%N
    A => monoMAT%R%A
    index => monoMAT%CSR%index
    item => monoMAT%CSR%item

    NZ = 0
    do i = 1, N
      NZ = NZ + monoMAT%n_dof_list(i)*monoMAT%n_dof_list(i)
    enddo

    NDOF_MAX = maxval(monoMAT%n_dof_list)

    call monolis_alloc_R_1d(T, NDOF_MAX)
    call monolis_alloc_R_2d(LU, NDOF_MAX, NDOF_MAX)
    call monolis_palloc_R_1d(monoPREC%R%D, NZ)
    ALU => monoPREC%R%D
    monoPREC%N = monoMAT%N

    kn = 0
    do i = 1, N
      jS = index(i) + 1
      jE = index(i + 1)
      do ii = jS, jE
        in = item(ii)
        if(i == in)then
          NDOF = monoMAT%n_dof_list(i)
          jn = monoMAT%n_dof_index2(ii)
          do j = 1, NDOF
            do k = 1, NDOF
              LU(j,k) = A(jn + NDOF*(j-1) + k)
            enddo
          enddo

          do k = 1, NDOF
            if(LU(k,k) == 0.0d0) stop "** monolis error: zero diag in monolis_precond_diag_V_setup_R"
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
              ALU(kn + NDOF*(j-1) + k) = LU(j,k)
            enddo
          enddo
          kn = kn + NDOF*NDOF
        endif
      enddo
    enddo

    deallocate(T)
    deallocate(LU)
  end subroutine monolis_precond_diag_V_setup_R

  !> @ingroup prec
  !> 前処理生成：対角スケーリング前処理（任意ブロック、複素数型）
  subroutine monolis_precond_diag_V_setup_C(monoMAT, monoPREC)
    implicit none
    !> [in] 行列構造体
    type(monolis_mat), target, intent(in) :: monoMAT
    !> [in,out] 前処理構造体
    type(monolis_mat), target, intent(inout) :: monoPREC
    integer(kint) :: i, ii, j, jS, jE, in, jn, k, l, N, NZ, NDOF_MAX
    integer(kint) :: NDOF, kn
    complex(kdouble), allocatable :: T(:), LU(:,:)
    integer(kint), pointer :: index(:), item(:)
    complex(kdouble), pointer :: A(:), ALU(:)

    N =  monoMAT%N
    A => monoMAT%C%A
    index => monoMAT%CSR%index
    item => monoMAT%CSR%item

    NZ = 0
    do i = 1, N
      NZ = NZ + monoMAT%n_dof_list(i)*monoMAT%n_dof_list(i)
    enddo

    NDOF_MAX = maxval(monoMAT%n_dof_list)

    call monolis_alloc_C_1d(T, NDOF_MAX)
    call monolis_alloc_C_2d(LU, NDOF_MAX, NDOF_MAX)
    call monolis_palloc_C_1d(monoPREC%C%D, NZ)
    ALU => monoPREC%C%D
    monoPREC%N = monoMAT%N

    kn = 0
    do i = 1, N
      jS = index(i) + 1
      jE = index(i + 1)
      do ii = jS, jE
        in = item(ii)
        if(i == in)then
          NDOF = monoMAT%n_dof_list(i)
          jn = monoMAT%n_dof_index2(ii)
          do j = 1, NDOF
            do k = 1, NDOF
              LU(j,k) = A(jn + NDOF*(j-1) + k)
            enddo
          enddo

          do k = 1, NDOF
            if(LU(k,k) == 0.0d0) stop "** monolis error: zero diag in monolis_precond_diag_V_setup_C"
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
              ALU(kn + NDOF*(j-1) + k) = LU(j,k)
            enddo
          enddo
          kn = kn + NDOF*NDOF
        endif
      enddo
    enddo

    deallocate(T)
    deallocate(LU)
  end subroutine monolis_precond_diag_V_setup_C

  !> @ingroup prec
  !> 前処理適用：対角スケーリング前処理（任意ブロック、実数型）
  subroutine monolis_precond_diag_V_apply_R(monoMAT, monoPREC, X, Y)
    implicit none
    !> [in] 行列構造体
    type(monolis_mat), target, intent(in) :: monoMAT
    !> [in] 前処理構造体
    type(monolis_mat), target, intent(in) :: monoPREC
    integer(kint) :: i, j, k, in, kn, N, NDOF, NDOF_MAX
    real(kdouble) :: X(:), Y(:)
    real(kdouble), pointer :: ALU(:)
    real(kdouble), allocatable :: T(:)

    N =  monoMAT%N
    ALU => monoPREC%R%D
    NDOF_MAX = maxval(monoMAT%n_dof_list)

    call monolis_alloc_R_1d(T, NDOF_MAX)

    kn = 0
    do i = 1, N
      NDOF = monoMAT%n_dof_list(i)
      in = monoMAT%n_dof_index(i)
      do j = 1, NDOF
        T(j) = X(in + j)
      enddo
      do j = 2, NDOF
        do k = 1, j-1
          T(j) = T(j) - ALU(kn + NDOF*(j-1) + k)*T(k)
        enddo
      enddo
      do j = NDOF, 1, -1
        do k = NDOF, j+1, -1
          T(j) = T(j) - ALU(kn + NDOF*(j-1) + k)*T(k)
        enddo
        T(j) = ALU(kn + NDOF*(j-1) + j)*T(j)
      enddo
      kn = kn + NDOF*NDOF
      do j = 1, NDOF
        Y(in + j) = T(j)
      enddo
    enddo

    deallocate(T)
  end subroutine monolis_precond_diag_V_apply_R

  !> @ingroup prec
  !> 前処理適用：対角スケーリング前処理（任意ブロック、複素数型）
  subroutine monolis_precond_diag_V_apply_C(monoMAT, monoPREC, X, Y)
    implicit none
    !> [in] 行列構造体
    type(monolis_mat), target, intent(in) :: monoMAT
    !> [in] 前処理構造体
    type(monolis_mat), target, intent(in) :: monoPREC
    integer(kint) :: i, j, k, in, kn, N, NDOF, NDOF_MAX
    complex(kdouble) :: X(:), Y(:)
    complex(kdouble), pointer :: ALU(:)
    complex(kdouble), allocatable :: T(:)

    N =  monoMAT%N
    ALU => monoPREC%C%D
    NDOF_MAX = maxval(monoMAT%n_dof_list)

    call monolis_alloc_C_1d(T, NDOF_MAX)

    kn = 0
    do i = 1, N
      NDOF = monoMAT%n_dof_list(i)
      in = monoMAT%n_dof_index(i)
      do j = 1, NDOF
        T(j) = X(in + j)
      enddo
      do j = 2, NDOF
        do k = 1, j-1
          T(j) = T(j) - ALU(kn + NDOF*(j-1) + k)*T(k)
        enddo
      enddo
      do j = NDOF, 1, -1
        do k = NDOF, j+1, -1
          T(j) = T(j) - ALU(kn + NDOF*(j-1) + k)*T(k)
        enddo
        T(j) = ALU(kn + (NDOF+1)*(j-1) + 1)*T(j)
      enddo
      kn = kn + NDOF*NDOF
      do j = 1, NDOF
        Y(in + j) = T(j)
      enddo
    enddo

    deallocate(T)
  end subroutine monolis_precond_diag_V_apply_C

  !> @ingroup prec
  !> 前処理初期化：対角スケーリング前処理（任意ブロック、実数型）
  subroutine monolis_precond_diag_V_clear_R(monoPREC)
    implicit none
    !> [in,out] 前処理構造体
    type(monolis_mat), intent(inout) :: monoPREC

    call monolis_pdealloc_R_1d(monoPREC%R%D)
  end subroutine monolis_precond_diag_V_clear_R

  !> @ingroup prec
  !> 前処理初期化：対角スケーリング前処理（任意ブロック、複素数型）
  subroutine monolis_precond_diag_V_clear_C(monoPREC)
    implicit none
    !> [in,out] 前処理構造体
    type(monolis_mat), intent(inout) :: monoPREC

    call monolis_pdealloc_C_1d(monoPREC%C%D)
  end subroutine monolis_precond_diag_V_clear_C
end module mod_monolis_precond_diag_V
