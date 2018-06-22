module mod_monolis_fact_LU_asym_33
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat

  implicit none
  private
  public :: monolis_init_LU_inner_asym_33
  public :: monolis_fact_LU_inner_asym_33
  public :: monolis_solv_LU_inner_asym_33
  public :: monolis_clear_LU_inner_asym_33

  integer(kind=kint), pointer :: idxL(:)  => null()
  integer(kind=kint), pointer :: idxU(:)  => null()
  integer(kind=kint), pointer :: itemL(:) => null()
  integer(kind=kint), pointer :: itemU(:) => null()
  real(kind=kdouble), save, allocatable :: AD(:)
  real(kind=kdouble), save, allocatable :: AL(:)
  real(kind=kdouble), save, allocatable :: AU(:)

contains

  subroutine monolis_init_LU_inner_asym_33(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT

  end subroutine monolis_init_LU_inner_asym_33

  subroutine monolis_clear_LU_inner_asym_33(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT

  end subroutine monolis_clear_LU_inner_asym_33

  subroutine monolis_fact_LU_inner_asym_33(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kind=kint) :: N, NP, NIN, NDOF, NNDOF
    integer(kind=kint) :: i, j, k, jS, jE, in, jn, kn, nn
    integer(kind=kint) :: shift, imax
    integer(kind=kint), allocatable :: isFill(:)
    real(kind=kdouble) :: d1, d2, d3
    real(kind=kdouble) :: u1, u2, u3
    real(kind=kdouble) :: l1, l2, l3
    real(kind=kdouble) :: t1, t2
    real(kind=kdouble) :: D(9), L(9)
    real(kind=kdouble), allocatable :: U(:)

    N  = monoMAT%N
    NP = monoMAT%NP

    imax = 0
    do i = 1,N
      in = idxU(i) - idxU(i-1) - 1
      imax = max(imax, in)
    enddo
    allocate(U(9*imax))
    allocate(isFill(NP)) !first allocation only?

    !factorization
    do i = 1,N
      in = idxU(i-1)+1
      if(AU(9*in-8) == 0.0d0) stop "*** zero diagonal entry"
      if(AU(9*in-4) == 0.0d0) stop "*** zero diagonal entry"
      if(AU(9*in  ) == 0.0d0) stop "*** zero diagonal entry"

      d1 = 1.0d0/AU(9*in-8)
      AU(9*in-8) = d1
      AU(9*in-5) = d1 * AU(9*in-5)
      AU(9*in-2) = d1 * AU(9*in-2)
      AU(9*in-4) = AU(9*in-4) - AU(9*in-5) * AU(9*in-7)
      AU(9*in-3) = AU(9*in-3) - AU(9*in-5) * AU(9*in-6)
      AU(9*in-1) = AU(9*in-1) - AU(9*in-2) * AU(9*in-7)
      AU(9*in  ) = AU(9*in  ) - AU(9*in-2) * AU(9*in-6)

      d2 = 1.0d0/AU(9*in-4)
      AU(9*in-4) = d2
      AU(9*in-1) = d2 * AU(9*in-1)
      AU(9*in  ) = AU(9*in  ) - AU(9*in-1) * AU(9*in-3)

      d3 = 1.0d0/AU(9*in)
      AU(9*in) = d3

      l1 = AU(9*in-5)
      l2 = AU(9*in-2)
      l3 = AU(9*in-1)
      u1 = AU(9*in-7)
      u2 = AU(9*in-6)
      u3 = AU(9*in-3)

      !U multiple section
      !in = idxU(i) - idxU(i-1) -1

      do j=i,NP
        isFill(j) = 0
      enddo

      in = 1
      do j = idxU(i-1)+2, idxU(i)
        isFill(itemU(j)) = in
        AU(9*j-5) = AU(9*j-5) - AU(9*j-8)*l1
        AU(9*j-4) = AU(9*j-4) - AU(9*j-7)*l1
        AU(9*j-3) = AU(9*j-3) - AU(9*j-6)*l1
        AU(9*j-2) = AU(9*j-2) - AU(9*j-8)*l2 - AU(9*j-5)*l3
        AU(9*j-1) = AU(9*j-1) - AU(9*j-7)*l2 - AU(9*j-4)*l3
        AU(9*j  ) = AU(9*j  ) - AU(9*j-6)*l2 - AU(9*j-3)*l3
        !AU(9*j-2) = AU(9*j-2) - AU(9*j-5)*l3
        !AU(9*j-1) = AU(9*j-1) - AU(9*j-4)*l3
        !AU(9*j  ) = AU(9*j  ) - AU(9*j-3)*l3
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

      do j = idxU(i-1)+2, idxU(i)
        jn = itemU(j)
        do k = idxL(jn-1)+1, idxL(jn)-1
          kn = itemL(k)
          if(kn == i)then
            !L inverse section
            AL(9*k-8) = AL(9*k-8) * d1
            AL(9*k-5) = AL(9*k-5) * d1
            AL(9*k-2) = AL(9*k-2) * d1
            L (1) = AL(9*k-8)
            L (4) = AL(9*k-5)
            L (7) = AL(9*k-2)
            AL(9*k-7) = AL(9*k-7) - AL(9*k-8)*u1
            AL(9*k-6) = AL(9*k-6) - AL(9*k-8)*u2
            AL(9*k-4) = AL(9*k-4) - AL(9*k-5)*u1
            AL(9*k-3) = AL(9*k-3) - AL(9*k-5)*u2
            AL(9*k-1) = AL(9*k-1) - AL(9*k-2)*u1
            AL(9*k  ) = AL(9*k  ) - AL(9*k-2)*u2

            AL(9*k-7) = AL(9*k-7) * d2
            AL(9*k-4) = AL(9*k-4) * d2
            AL(9*k-1) = AL(9*k-1) * d2
            L (2) = AL(9*k-7)
            L (5) = AL(9*k-4)
            L (8) = AL(9*k-1)
            AL(9*k-6) = AL(9*k-6) - AL(9*k-7)*u3
            AL(9*k-3) = AL(9*k-3) - AL(9*k-4)*u3
            AL(9*k  ) = AL(9*k  ) - AL(9*k-1)*u3

            AL(9*k-6) = AL(9*k-6) * d3
            AL(9*k-3) = AL(9*k-3) * d3
            AL(9*k  ) = AL(9*k  ) * d3
            L (3) = AL(9*k-6)
            L (6) = AL(9*k-3)
            L (9) = AL(9*k  )
          elseif(i<kn)then
            if(isFill(kn)/=0)then
              !outer product section
              in = isFill(kn)
              AL(9*k-8) = AL(9*k-8) - L(1)*U(9*in-8) - L(2)*U(9*in-5) - L(3)*U(9*in-2)
              AL(9*k-7) = AL(9*k-7) - L(1)*U(9*in-7) - L(2)*U(9*in-4) - L(3)*U(9*in-1)
              AL(9*k-6) = AL(9*k-6) - L(1)*U(9*in-6) - L(2)*U(9*in-3) - L(3)*U(9*in  )
              AL(9*k-5) = AL(9*k-5) - L(4)*U(9*in-8) - L(5)*U(9*in-5) - L(6)*U(9*in-2)
              AL(9*k-4) = AL(9*k-4) - L(4)*U(9*in-7) - L(5)*U(9*in-4) - L(6)*U(9*in-1)
              AL(9*k-3) = AL(9*k-3) - L(4)*U(9*in-6) - L(5)*U(9*in-3) - L(6)*U(9*in  )
              AL(9*k-2) = AL(9*k-2) - L(7)*U(9*in-8) - L(8)*U(9*in-5) - L(9)*U(9*in-2)
              AL(9*k-1) = AL(9*k-1) - L(7)*U(9*in-7) - L(8)*U(9*in-4) - L(9)*U(9*in-1)
              AL(9*k  ) = AL(9*k  ) - L(7)*U(9*in-6) - L(8)*U(9*in-3) - L(9)*U(9*in  )
            endif
          endif
        enddo
        do k = idxU(jn-1)+1, idxU(jn)
          !outer product section
          kn = itemU(k)
          if(isFill(kn)/=0)then
            in = isFill(kn)
            AU(9*k-8) = AU(9*k-8) - L(1)*U(9*in-8) - L(2)*U(9*in-5) - L(3)*U(9*in-2)
            AU(9*k-7) = AU(9*k-7) - L(1)*U(9*in-7) - L(2)*U(9*in-4) - L(3)*U(9*in-1)
            AU(9*k-6) = AU(9*k-6) - L(1)*U(9*in-6) - L(2)*U(9*in-3) - L(3)*U(9*in  )
            AU(9*k-5) = AU(9*k-5) - L(4)*U(9*in-8) - L(5)*U(9*in-5) - L(6)*U(9*in-2)
            AU(9*k-4) = AU(9*k-4) - L(4)*U(9*in-7) - L(5)*U(9*in-4) - L(6)*U(9*in-1)
            AU(9*k-3) = AU(9*k-3) - L(4)*U(9*in-6) - L(5)*U(9*in-3) - L(6)*U(9*in  )
            AU(9*k-2) = AU(9*k-2) - L(7)*U(9*in-8) - L(8)*U(9*in-5) - L(9)*U(9*in-2)
            AU(9*k-1) = AU(9*k-1) - L(7)*U(9*in-7) - L(8)*U(9*in-4) - L(9)*U(9*in-1)
            AU(9*k  ) = AU(9*k  ) - L(7)*U(9*in-6) - L(8)*U(9*in-3) - L(9)*U(9*in  )
          endif
        enddo
      enddo
    enddo

    deallocate(isFill)
    deallocate(U)
  end subroutine monolis_fact_LU_inner_asym_33

  subroutine monolis_solv_LU_inner_asym_33(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kind=kint) :: N, NP, NNDOF
    integer(kind=kint) :: i, j, k, in, jS, jE, kn
    real(kind=kdouble) :: t1, t2
    real(kind=kdouble) :: X1, X2, X3
    real(kind=kdouble) :: A1, A2, A3
    real(kind=kdouble), allocatable :: S(:), X(:)

    N = monoMAT%N
    NP = monoMAT%NP
    NDOF = monoMAT%NDOF
    myrank = hecMESH%my_rank

    !L
    do i = 1,N
      A1 = 0.0d0
      A2 = 0.0d0
      A3 = 0.0d0
      do j = idxL(i-1)+1, idxL(i)-1
        in = itemL(j)
        X1 = X(3*in-2)
        X2 = X(3*in-1)
        X3 = X(3*in  )
        A1 = A1 + AL(9*j-8)*X1 + AL(9*j-7)*X2 + AL(9*j-6)*X3
        A2 = A2 + AL(9*j-5)*X1 + AL(9*j-4)*X2 + AL(9*j-3)*X3
        A3 = A3 + AL(9*j-2)*X1 + AL(9*j-1)*X2 + AL(9*j  )*X3
      enddo
      in = idxU(i-1)+1
      X(3*i-2) = X(3*i-2) - A1
      X(3*i-1) = X(3*i-1) - A2 - AU(9*in-5)*X(3*i-2)
      X(3*i  ) = X(3*i  ) - A3 - AU(9*in-2)*X(3*i-2) - AU(9*in-1)*X(3*i-1)
    enddo
    !U
    do i = N,1,-1
      A1 = 0.0d0
      A2 = 0.0d0
      A3 = 0.0d0
      do j = idxU(i),idxU(i-1)+2,-1
        in = itemU(j)
        X1 = X(3*in-2)
        X2 = X(3*in-1)
        X3 = X(3*in  )
        A1 = A1 + AU(9*j-8)*X1 + AU(9*j-7)*X2 + AU(9*j-6)*X3
        A2 = A2 + AU(9*j-5)*X1 + AU(9*j-4)*X2 + AU(9*j-3)*X3
        A3 = A3 + AU(9*j-2)*X1 + AU(9*j-1)*X2 + AU(9*j  )*X3
      enddo
      in = idxU(i-1)+1
      X(3*i  ) = AU(9*in  ) * ( X(3*i  ) - A3 )
      X(3*i-1) = AU(9*in-4) * ( X(3*i-1) - A2 - AU(9*in-3)*X(3*i) )
      X(3*i-2) = AU(9*in-8) * ( X(3*i-2) - A1 - AU(9*in-7)*X(3*i-1) - AU(9*in-6)*X(3*i) )
    enddo
  end subroutine monolis_solv_LU_inner_asym_33
end module mod_monolis_fact_LU_asym_33
