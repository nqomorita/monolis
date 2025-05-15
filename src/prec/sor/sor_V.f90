!> SOR 前処理（nxn ブロック）
module mod_monolis_precond_sor_V
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc

  implicit none

contains

  !> @ingroup prec
  !> 前処理生成：SOR 前処理（任意ブロック、実数型）
  subroutine monolis_precond_sor_V_setup_R(monoMAT, monoPREC)
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
            if(LU(k,k) == 0.0d0) stop "** monolis error: zero diag in monolis_precond_sor_V_setup_R"
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
  end subroutine monolis_precond_sor_V_setup_R

  !> @ingroup prec
  !> 前処理生成：SOR 前処理（任意ブロック、複素数型）
  subroutine monolis_precond_sor_V_setup_C(monoMAT, monoPREC)
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
            if(LU(k,k) == 0.0d0) stop "** monolis error: zero diag in monolis_precond_sor_V_setup_C"
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
  end subroutine monolis_precond_sor_V_setup_C

  !> @ingroup prec
  !> 前処理適用：SOR 前処理（任意ブロック、実数型）
  subroutine monolis_precond_sor_V_apply_R(monoMAT, monoPREC, X, Y)
    implicit none
    !> [in] 行列構造体
    type(monolis_mat), target, intent(in) :: monoMAT
    !> [in] 前処理構造体
    type(monolis_mat), target, intent(in) :: monoPREC
    integer(kint) :: i, j, jE, jS, in, jn, kn, ln, k, l, NDOF
    integer(kint) :: n1, n2, nz
    integer(kint) :: NNDOF, NPNDOF, NDOF_MAX
    integer(kint), pointer :: index(:)
    integer(kint), pointer :: item(:)
    real(kdouble) :: X(:), Y(:)
    real(kdouble), pointer :: A(:), ALU(:)
    real(kdouble), allocatable :: XT(:), YT(:), ST(:)

    A => monoMAT%R%A
    ALU => monoPREC%R%D
    index => monoMAT%CSR%index
    item => monoMAT%CSR%item

    NDOF_MAX = maxval(monoMAT%n_dof_list)

    call monolis_get_vec_size(monoMAT%N, monoMAT%NP, monoMAT%NDOF, &
      monoMAT%n_dof_index, NNDOF, NPNDOF)

    do i = 1, NNDOF
      Y(i) = X(i)
    enddo

    call monolis_alloc_R_1d(XT, NDOF_MAX)
    call monolis_alloc_R_1d(YT, NDOF_MAX)
    call monolis_alloc_R_1d(ST, NDOF_MAX)

    nz = 0
    do i = 1, monoMAT%N
      n1 = monoMAT%n_dof_list(i)
      in = monoMAT%n_dof_index(i)
      do j = 1, n1
        ST(j) = Y(in + j)
      enddo
      jS = index(i) + 1
      jE = index(i + 1)
      do j = jS, jE
        jn = item(j)
        n2 = monoMAT%n_dof_list(jn)
        if(jn < i)then
          kn = monoMAT%n_dof_index(jn)
          do k = 1, n2
            XT(k) = Y(kn + k)
          enddo
          ln = monoMAT%n_dof_index2(j)
          do k = 1, n1
            do l = 1, n2
              ST(k) = ST(k) - A(ln + n2*(k-1) + l)*XT(l)
            enddo
          enddo
        endif
      enddo

      do j = 1, n1
        XT(j) = ST(j)
      enddo
      do j = 2, n1
        do k = 1, j-1
          XT(j) = XT(j) - ALU(nz + n1*(j-1) + k)*XT(k)
        enddo
      enddo
      do j = n1, 1, -1
        do k = n1, j+1, -1
          XT(j) = XT(j) - ALU(nz + n1*(j-1) + k)*XT(k)
        enddo
        XT(j) = ALU(nz + (n1+1)*(j-1) + 1)*XT(j)
      enddo
      do k = 1, n1
        Y(in + k) = XT(k)
      enddo
      nz = nz + n1*n1
    enddo

    in = monoMAT%n_dof_list(monoMAT%N)
    nz = monoMAT%n_dof_index2(monoMAT%N) - in*in
    do i = monoMAT%N, 1, -1
      n1 = monoMAT%n_dof_list(i)
      do j = 1, n1
        ST(j) = 0.0d0
      enddo
      jS = index(i) + 1
      jE = index(i + 1)
      do j = jE, jS, -1
        jn = item(j)
        n2 = monoMAT%n_dof_list(jn)
        if(i < jn)then
          kn = monoMAT%n_dof_index(jn)
          do k = 1, n2
            XT(k) = Y(kn + k)
          enddo
          ln = monoMAT%n_dof_index2(j)
          do k = 1, n1
            do l = 1, n2
              ST(k) = ST(k) + A(ln + n2*(k-1) + l)*XT(l)
            enddo
          enddo
        endif
      enddo

      do j = 1, n1
        XT(j) = ST(j)
      enddo
      do j = 2, n1
        do k = 1, j-1
          XT(j) = XT(j) - ALU(nz + n1*(j-1) + k)*XT(k)
        enddo
      enddo
      do j = n1, 1, -1
        do k = n1, j+1, -1
          XT(j) = XT(j) - ALU(nz + n1*(j-1) + k)*XT(k)
        enddo
        XT(j) = ALU(nz + (n1+1)*(j-1) + 1)*XT(j)
      enddo
      in = monoMAT%n_dof_index(i)
      do k = 1, n1
        Y(in + k) = Y(in + k) - XT(k)
      enddo
      in = monoMAT%n_dof_list(i - 1)
      nz = nz - in*in
    enddo

    deallocate(XT)
    deallocate(YT)
    deallocate(ST)
  end subroutine monolis_precond_sor_V_apply_R

  !> 前処理適用：SOR 前処理（任意ブロック、複素数型）
  subroutine monolis_precond_sor_V_apply_C(monoMAT, monoPREC, X, Y)
    implicit none
    !> [in] 行列構造体
    type(monolis_mat), target, intent(in) :: monoMAT
    !> [in,out] 前処理構造体
    type(monolis_mat), target, intent(inout) :: monoPREC
    integer(kint) :: i, j, jE, jS, jn, k, l, N, NP, NDOF, NDOF2
    integer(kint), pointer :: index(:)
    integer(kint), pointer :: item(:)
    complex(kdouble) :: X(:), Y(:)
    complex(kdouble), pointer :: A(:), ALU(:)
    complex(kdouble), allocatable :: XT(:), YT(:), ST(:)

    N =  monoMAT%N
    NP = monoMAT%NP
    NDOF = monoMAT%NDOF
    NDOF2 = NDOF*NDOF
    A => monoMAT%C%A
    ALU => monoPREC%C%D
    index => monoMAT%CSR%index
    item => monoMAT%CSR%item

    do i = 1, NP*NDOF
      Y(i) = X(i)
    enddo

    call monolis_alloc_C_1d(XT, NDOF)
    call monolis_alloc_C_1d(YT, NDOF)
    call monolis_alloc_C_1d(ST, NDOF)

    do i = 1, N
      do j = 1, NDOF
        ST(j) = Y(NDOF*(i-1)+j)
      enddo
      jS = index(i) + 1
      jE = index(i + 1)
      do j = jS, jE
        jn = item(j)
        if(jn < i)then
          do k = 1, NDOF
            XT(k) = Y(NDOF*(jn-1)+k)
          enddo
          do k = 1, NDOF
            do l = 1, NDOF
              ST(k) = ST(k) - A(NDOF2*(j-1)+NDOF*(k-1)+l) * XT(l)
            enddo
          enddo
        endif
      enddo

      do j = 1, NDOF
        XT(j) = ST(j)
      enddo
      do j = 2, NDOF
        do k = 1, j-1
          XT(j) = XT(j) - ALU(NDOF2*(i-1) + NDOF*(j-1) + k)*XT(k)
        enddo
      enddo
      do j = NDOF, 1, -1
        do k = NDOF, j+1, -1
          XT(j) = XT(j) - ALU(NDOF2*(i-1) + NDOF*(j-1) + k)*XT(k)
        enddo
        XT(j) = ALU(NDOF2*(i-1) + (NDOF+1)*(j-1) + 1)*XT(j)
      enddo
      do k = 1, NDOF
        Y(NDOF*(i-1)+k) = XT(k)
      enddo
    enddo

    do i = N, 1, -1
      do j = 1, NDOF
        ST(j) = 0.0d0
      enddo
      jS = index(i) + 1
      jE = index(i + 1)
      do j = jE, jS, -1
        jn = item(j)
        if(i < jn)then
          do k = 1, NDOF
            XT(k) = Y(NDOF*(jn-1)+k)
          enddo
          do k = 1, NDOF
            do l = 1, NDOF
              ST(k) = ST(k) + A(NDOF2*(j-1)+(k-1)*NDOF+l)*XT(l)
            enddo
          enddo
        endif
      enddo

      do j = 1, NDOF
        XT(j) = ST(j)
      enddo
      do j = 2, NDOF
        do k = 1, j-1
          XT(j) = XT(j) - ALU(NDOF2*(i-1) + NDOF*(j-1) + k)*XT(k)
        enddo
      enddo
      do j = NDOF, 1, -1
        do k = NDOF, j+1, -1
          XT(j) = XT(j) - ALU(NDOF2*(i-1) + NDOF*(j-1) + k)*XT(k)
        enddo
        XT(j) = ALU(NDOF2*(i-1) + (NDOF+1)*(j-1) + 1)*XT(j)
      enddo
      do k = 1, NDOF
        Y(NDOF*(i-1)+k) = Y(NDOF*(i-1)+k) - XT(k)
      enddo
    enddo

    deallocate(XT)
    deallocate(YT)
    deallocate(ST)
  end subroutine monolis_precond_sor_V_apply_C

  !> @ingroup prec
  !> 前処理初期化：SOR 前処理（任意ブロック、実数型）
  subroutine monolis_precond_sor_V_clear_R(monoPREC)
    implicit none
    !> [in,out] 前処理構造体
    type(monolis_mat), intent(inout) :: monoPREC

    call monolis_pdealloc_R_1d(monoPREC%R%D)
  end subroutine monolis_precond_sor_V_clear_R

  !> @ingroup prec
  !> 前処理初期化：SOR 前処理（任意ブロック、複素数型）
  subroutine monolis_precond_sor_V_clear_C(monoPREC)
    implicit none
    !> [in,out] 前処理構造体
    type(monolis_mat), intent(inout) :: monoPREC

    call monolis_pdealloc_C_1d(monoPREC%C%D)
  end subroutine monolis_precond_sor_V_clear_C
end module mod_monolis_precond_sor_V
