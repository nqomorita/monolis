module mod_monolis_fact_LU_33
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat

  implicit none
  private
  public :: monolis_init_LU_inner_33
  public :: monolis_fact_LU_inner_33
  public :: monolis_solv_LU_inner_33
  public :: monolis_clear_LU_inner_33

  real(kind=kdouble), save, allocatable :: AD(:)
  real(kind=kdouble), save, pointer :: AU(:)

contains

  subroutine monolis_init_LU_inner_33(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat_LDU) :: monoMAT

    AU => monoMAT%U
  end subroutine monolis_init_LU_inner_33

  subroutine monolis_clear_LU_inner_33(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat_LDU) :: monoMAT

    deallocate(AD)
  end subroutine monolis_clear_LU_inner_33

  subroutine monolis_fact_LU_inner_33(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat_LDU) :: monoMAT
    integer(kind=kint) :: N, NP
    integer(kind=kint) :: i, j, k, in, jn, kn, nn
    integer(kind=kint) :: imax
    integer(kind=kint), allocatable :: is_fill(:)
    integer(kind=kint), pointer :: idxU(:), itemU(:)
    real(kind=kdouble) :: d1, d2, d3
    real(kind=kdouble) :: u1, u2, u3
    real(kind=kdouble) :: t1, t2
    real(kind=kdouble) :: L(9)
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

    allocate(U(9*imax))
    allocate(AD(3*N))
    allocate(is_fill(NP))
    AD = 0.0d0

    !factorization
    do i = 1,N
      in = idxU(i-1)+1
      if(AU(9*in-8) == 0.0d0) stop "*** zero diagonal entry"
      if(AU(9*in-4) == 0.0d0) stop "*** zero diagonal entry"
      if(AU(9*in  ) == 0.0d0) stop "*** zero diagonal entry"

      d1 = 1.0d0/AU(9*in-8)
      AU(9*in-8) = d1
      AU(9*in-4) = AU(9*in-4) - AU(9*in-7)*AU(9*in-7)*d1
      AU(9*in-3) = AU(9*in-3) - AU(9*in-7)*AU(9*in-6)*d1
      AU(9*in  ) = AU(9*in  ) - AU(9*in-6)*AU(9*in-6)*d1

      d2 = 1.0d0/AU(9*in-4)
      AU(9*in-4) = d2
      AU(9*in  ) = AU(9*in  ) - AU(9*in-3)*AU(9*in-3)*d2

      d3 = 1.0d0/AU(9*in)
      AU(9*in) = d3

      u1 = AU(9*in-7)
      u2 = AU(9*in-6)
      u3 = AU(9*in-3)

      !U multiple section
      do j = i, NP
        is_fill(j) = 0
      enddo
      in = 1
      do j = idxU(i-1)+2, idxU(i)
        is_fill(itemU(j)) = in
        AU(9*j-5) = AU(9*j-5) - AU(9*j-8)*u1*d1
        AU(9*j-4) = AU(9*j-4) - AU(9*j-7)*u1*d1
        AU(9*j-3) = AU(9*j-3) - AU(9*j-6)*u1*d1
        AU(9*j-2) = AU(9*j-2) - AU(9*j-8)*u2*d1 - AU(9*j-5)*u3*d2
        AU(9*j-1) = AU(9*j-1) - AU(9*j-7)*u2*d1 - AU(9*j-4)*u3*d2
        AU(9*j  ) = AU(9*j  ) - AU(9*j-6)*u2*d1 - AU(9*j-3)*u3*d2
        U(9*in-8) = AU(9*j-8)
        U(9*in-7) = AU(9*j-7)
        U(9*in-6) = AU(9*j-6)
        U(9*in-5) = AU(9*j-5)
        U(9*in-4) = AU(9*j-4)
        U(9*in-3) = AU(9*j-3)
        U(9*in-2) = AU(9*j-2)
        U(9*in-1) = AU(9*j-1)
        U(9*in  ) = AU(9*j  )
        in = in + 1
      enddo

      nn = 1
      do j = idxU(i-1)+2, idxU(i)
        L(1) = U(9*nn-8)
        L(4) = U(9*nn-7)
        L(7) = U(9*nn-6)
        L(2) = U(9*nn-5)
        L(5) = U(9*nn-4)
        L(8) = U(9*nn-3)
        L(3) = U(9*nn-2)
        L(6) = U(9*nn-1)
        L(9) = U(9*nn  )
        nn = nn + 1
        jn = itemU(j)
        do k = idxU(jn-1)+1, idxU(jn)
          !outer product section
          kn = itemU(k)
          if(is_fill(kn) /= 0)then
            in = is_fill(kn)
            AU(9*k-8) = AU(9*k-8) - L(1)*U(9*in-8)*d1 - L(2)*U(9*in-5)*d2 - L(3)*U(9*in-2)*d3
            AU(9*k-7) = AU(9*k-7) - L(1)*U(9*in-7)*d1 - L(2)*U(9*in-4)*d2 - L(3)*U(9*in-1)*d3
            AU(9*k-6) = AU(9*k-6) - L(1)*U(9*in-6)*d1 - L(2)*U(9*in-3)*d2 - L(3)*U(9*in  )*d3
            AU(9*k-5) = AU(9*k-5) - L(4)*U(9*in-8)*d1 - L(5)*U(9*in-5)*d2 - L(6)*U(9*in-2)*d3
            AU(9*k-4) = AU(9*k-4) - L(4)*U(9*in-7)*d1 - L(5)*U(9*in-4)*d2 - L(6)*U(9*in-1)*d3
            AU(9*k-3) = AU(9*k-3) - L(4)*U(9*in-6)*d1 - L(5)*U(9*in-3)*d2 - L(6)*U(9*in  )*d3
            AU(9*k-2) = AU(9*k-2) - L(7)*U(9*in-8)*d1 - L(8)*U(9*in-5)*d2 - L(9)*U(9*in-2)*d3
            AU(9*k-1) = AU(9*k-1) - L(7)*U(9*in-7)*d1 - L(8)*U(9*in-4)*d2 - L(9)*U(9*in-1)*d3
            AU(9*k  ) = AU(9*k  ) - L(7)*U(9*in-6)*d1 - L(8)*U(9*in-3)*d2 - L(9)*U(9*in  )*d3
          endif
        enddo
      enddo
    enddo

    do i = 1, N
      in = idxU(i-1)+1
      AD(3*i-2) = 1.0d0 / AU(9*in-8)
      AD(3*i-1) = 1.0d0 / AU(9*in-4)
      AD(3*i  ) = 1.0d0 / AU(9*in  )
    enddo

    deallocate(is_fill)
    deallocate(U)
  end subroutine monolis_fact_LU_inner_33

  subroutine monolis_solv_LU_inner_33(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat_LDU) :: monoMAT
    integer(kind=kint) :: N, NP
    integer(kind=kint) :: i, j, in, jS, jE
    integer(kind=kint), pointer :: idxU(:), itemU(:)
    real(kind=kdouble) :: t1, t2
    real(kind=kdouble) :: X1, X2, X3
    real(kind=kdouble) :: A1, A2, A3
    real(kind=kdouble), allocatable :: S(:)
    real(kind=kdouble), pointer :: X(:)

    N  = monoMAT%N
    NP = monoMAT%NP
    allocate(S(3*NP))
    S = 0.0d0
    monoMAT%X = monoMAT%B
    X => monoMAT%X
    idxU => monoMAT%indexU
    itemU => monoMAT%itemU

    !L
    do i = 1, N
      A1 = S(3*i-2)
      A2 = S(3*i-1)
      A3 = S(3*i  )
      in = idxU(i-1)+1
      X(3*i-2) = AU(9*in-8)*(X(3*i-2) - A1)
      X(3*i-1) = AU(9*in-4)*(X(3*i-1) - A2 - AU(9*in-7)*X(3*i-2))
      X(3*i  ) = AU(9*in  )*(X(3*i  ) - A3 - AU(9*in-6)*X(3*i-2) - AU(9*in-3)*X(3*i-1))
      X1 = X(3*i-2)
      X2 = X(3*i-1)
      X3 = X(3*i  )
      jS = idxU(i-1)+2
      jE = idxU(i)
      do j = jS, jE
        in = itemU(j)
        S(3*in-2) = S(3*in-2) + AU(9*j-8)*X1 + AU(9*j-5)*X2 + AU(9*j-2)*X3
        S(3*in-1) = S(3*in-1) + AU(9*j-7)*X1 + AU(9*j-4)*X2 + AU(9*j-1)*X3
        S(3*in  ) = S(3*in  ) + AU(9*j-6)*X1 + AU(9*j-3)*X2 + AU(9*j  )*X3
      enddo
    enddo
    !D
    do i = 1, N
      X(3*i-2) = X(3*i-2)*AD(3*i-2)
      X(3*i-1) = X(3*i-1)*AD(3*i-1)
      X(3*i  ) = X(3*i  )*AD(3*i  )
    enddo
    !U
    do i = N, 1, -1
      A1 = 0.0d0
      A2 = 0.0d0
      A3 = 0.0d0
      jS = idxU(i-1)+2
      jE = idxU(i)
      do j = jE, jS, -1
        in = itemU(j)
        X1 = X(3*in-2)
        X2 = X(3*in-1)
        X3 = X(3*in  )
        A1 = A1 + AU(9*j-8)*X1 + AU(9*j-7)*X2 + AU(9*j-6)*X3
        A2 = A2 + AU(9*j-5)*X1 + AU(9*j-4)*X2 + AU(9*j-3)*X3
        A3 = A3 + AU(9*j-2)*X1 + AU(9*j-1)*X2 + AU(9*j  )*X3
      enddo
      in = idxU(i-1)+1
      X(3*i  ) = AU(9*in  )*(X(3*i  ) - A3)
      X(3*i-1) = AU(9*in-4)*(X(3*i-1) - A2 - AU(9*in-3)*X(3*i))
      X(3*i-2) = AU(9*in-8)*(X(3*i-2) - A1 - AU(9*in-7)*X(3*i-1) - AU(9*in-6)*X(3*i))
    enddo

    deallocate(S)
  end subroutine monolis_solv_LU_inner_33
end module mod_monolis_fact_LU_33
