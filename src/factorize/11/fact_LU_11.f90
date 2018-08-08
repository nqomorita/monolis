module mod_monolis_fact_LU_11
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat

  implicit none
  private
  public :: monolis_init_LU_inner_11
  public :: monolis_fact_LU_inner_11
  public :: monolis_solv_LU_inner_11
  public :: monolis_clear_LU_inner_11

  real(kind=kdouble), save, allocatable :: AD(:)
  real(kind=kdouble), save, pointer :: AU(:)

contains

  subroutine monolis_init_LU_inner_11(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat_LDU) :: monoMAT

    AU => monoMAT%U
  end subroutine monolis_init_LU_inner_11

  subroutine monolis_clear_LU_inner_11(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat_LDU) :: monoMAT

    deallocate(AD)
  end subroutine monolis_clear_LU_inner_11

  subroutine monolis_fact_LU_inner_11(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat_LDU) :: monoMAT
    integer(kind=kint) :: N, NP
    integer(kind=kint) :: i, j, k, jS, jE, in, jn, kn, nn
    integer(kind=kint) :: imax
    integer(kind=kint), allocatable :: is_fill(:)
    integer(kind=kint), pointer :: idxU(:), itemU(:)
    real(kind=kdouble) :: d1, u1, l1
    real(kind=kdouble) :: t1, t2
    real(kind=kdouble) :: L(1)
    real(kind=kdouble), allocatable :: U(:)

    N  = monoMAT%N
    NP = monoMAT%NP
    idxU => monoMAT%indexU
    itemU => monoMAT%itemU

    imax = 0
    do i = 1, N
      in = idxU(i) - idxU(i-1) - 1
      imax = max(imax, in)
    enddo

    allocate(U(imax))
    allocate(AD(N))
    allocate(is_fill(NP))
    AD = 0.0d0

    !factorization
    do i = 1, N
      in = idxU(i-1)+1
      if(AU(in) == 0.0d0) stop "*** zero diagonal entry"

      d1 = 1.0d0/AU(in)
      AU(in) = d1

      !U multiple section
      do j = i, NP
        is_fill(j) = 0
      enddo
      in = 1
      do j = idxU(i-1)+2, idxU(i)
        is_fill(itemU(j)) = in
        U(in) = AU(j)
        in = in + 1
      enddo

      nn = 1
      do j = idxU(i-1)+2, idxU(i)
        L(1) = U(nn)
        nn = nn + 1
        jn = itemU(j)
        do k = idxU(jn-1)+1, idxU(jn)
          !outer product section
          kn = itemU(k)
          if(is_fill(kn) /= 0)then
            in = is_fill(kn)
            AU(k) = AU(k) - L(1)*U(in)*d1
          endif
        enddo
      enddo
    enddo

    do i = 1, N
      in = idxU(i-1)+1
      AD(i) = 1.0d0 / AU(in)
    enddo

    deallocate(is_fill)
    deallocate(U)
  end subroutine monolis_fact_LU_inner_11

  subroutine monolis_solv_LU_inner_11(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat_LDU) :: monoMAT
    integer(kind=kint) :: N, NP
    integer(kind=kint) :: i, j, k, in, jS, jE, kn
    integer(kind=kint), pointer :: idxU(:), itemU(:)
    real(kind=kdouble) :: t1, t2
    real(kind=kdouble) :: X1
    real(kind=kdouble) :: A1
    real(kind=kdouble), allocatable :: S(:)
    real(kind=kdouble), pointer :: X(:)

    N  = monoMAT%N
    NP = monoMAT%NP
    allocate(S(NP))
    S = 0.0d0
    monoMAT%X = monoMAT%B
    X => monoMAT%X
    idxU => monoMAT%indexU
    itemU => monoMAT%itemU

    !L
    do i = 1, N
      A1 = S(i)
      in = idxU(i-1)+1
      X(i) = AU(in)*(X(i) - A1)
      X1 = X(i)
      jS = idxU(i-1)+2
      jE = idxU(i)
      do j = jS, jE
        in = itemU(j)
        S(in) = S(in) + AU(j)*X1
      enddo
    enddo
    !D
    do i = 1, N
      X(i) = X(i)*AD(i)
    enddo
    !U
    do i = N, 1, -1
      A1 = 0.0d0
      jS = idxU(i-1)+2
      jE = idxU(i)
      do j = jE, jS, -1
        in = itemU(j)
        X1 = X(in)
        A1 = A1 + AU(j)*X1
      enddo
      in = idxU(i-1)+1
      X(i) = AU(in)*(X(i) - A1)
    enddo

    deallocate(S)
  end subroutine monolis_solv_LU_inner_11
end module mod_monolis_fact_LU_11
