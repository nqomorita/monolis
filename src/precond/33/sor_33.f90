module mod_monolis_precond_sor_33
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat

  implicit none

  private
  public :: monolis_precond_sor_33_setup
  public :: monolis_precond_sor_33_apply
  public :: monolis_precond_sor_33_clear

contains

  subroutine monolis_precond_sor_33_setup(monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kind=kint) :: i, j, jS, jE, in, k, l, N
    integer(kind=kint), pointer :: index(:), item(:)
    real(kind=kdouble) :: T(3,3), P(3), sigma
    real(kind=kdouble), pointer :: A(:), ALU(:)

    N =  monoMAT%N
    A => monoMAT%A
    index => monoMAT%index
    item => monoMAT%item
    sigma = 1.0d0

    allocate(monoMAT%monoTree%D(9*N))
    ALU => monoMAT%monoTree%D
    ALU = 0.0d0

    do i = 1, N
      jS = index(i-1) + 1
      jE = index(i)
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
  end subroutine monolis_precond_sor_33_setup

  subroutine monolis_precond_sor_33_apply(monoMAT, X, Y)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kind=kint) :: i, j, jS, jE, jn
    integer(kind=kint), pointer :: index(:)
    integer(kind=kint), pointer :: item(:)
    real(kind=kdouble) :: X1, X2, X3, S1, S2, S3
    real(kind=kdouble) :: X(:), Y(:)
    real(kind=kdouble), pointer :: A(:), ALU(:)

    index => monoMAT%index
    item => monoMAT%item
    A => monoMAT%A
    ALU => monoMAT%monoTree%D

    do i = 1, monoMAT%NP*monoMAT%NDOF
      Y(i) = X(i)
    enddo

    do i = 1, monoMAT%N
      S1 = Y(3*i-2)
      S2 = Y(3*i-1)
      S3 = Y(3*i  )
      jS = index(i-1) + 1
      jE = index(i)
      do j = jS, jE
        jn = item(j)
        if(jn < i)then
          X1 = Y(3*jn-2)
          X2 = Y(3*jn-1)
          X3 = Y(3*jn  )
          S1 = S1 - A(9*j-8)*X1 - A(9*j-7)*X2 - A(9*j-6)*X3
          S2 = S2 - A(9*j-5)*X1 - A(9*j-4)*X2 - A(9*j-3)*X3
          S3 = S3 - A(9*j-2)*X1 - A(9*j-1)*X2 - A(9*j  )*X3
        endif
      enddo
      X1 = S1
      X2 = S2
      X3 = S3
      X2 = X2 - ALU(9*i-5)*X1
      X3 = X3 - ALU(9*i-2)*X1 - ALU(9*i-1)*X2
      X3 = ALU(9*i  )* X3
      X2 = ALU(9*i-4)*(X2 - ALU(9*i-3)*X3)
      X1 = ALU(9*i-8)*(X1 - ALU(9*i-6)*X3 - ALU(9*i-7)*X2)
      Y(3*i-2) = X1
      Y(3*i-1) = X2
      Y(3*i  ) = X3
    enddo

    do i = monoMAT%N, 1, -1
      S1 = 0.0d0
      S2 = 0.0d0
      S3 = 0.0d0
      jS = index(i-1) + 1
      jE = index(i)
      do j = jE, jS, -1
        jn = item(j)
        if(i < jn)then
          X1 = Y(3*jn-2)
          X2 = Y(3*jn-1)
          X3 = Y(3*jn  )
          S1 = S1 + A(9*j-8)*X1 + A(9*j-7)*X2 + A(9*j-6)*X3
          S2 = S2 + A(9*j-5)*X1 + A(9*j-4)*X2 + A(9*j-3)*X3
          S3 = S3 + A(9*j-2)*X1 + A(9*j-1)*X2 + A(9*j  )*X3
        endif
      enddo
      X1 = S1
      X2 = S2
      X3 = S3
      X2 = X2 - ALU(9*i-5)*X1
      X3 = X3 - ALU(9*i-2)*X1 - ALU(9*i-1)*X2
      X3 = ALU(9*i  )* X3
      X2 = ALU(9*i-4)*(X2 - ALU(9*i-3)*X3)
      X1 = ALU(9*i-8)*(X1 - ALU(9*i-6)*X3 - ALU(9*i-7)*X2)
      Y(3*i-2) = Y(3*i-2) - X1
      Y(3*i-1) = Y(3*i-1) - X2
      Y(3*i  ) = Y(3*i  ) - X3
    enddo
  end subroutine monolis_precond_sor_33_apply

  subroutine monolis_precond_sor_33_clear(monoMAT)
    implicit none
    type(monolis_mat) :: monoMAT
    real(kind=kdouble), pointer :: ALU(:)
    ALU => monoMAT%monoTree%D
    deallocate(ALU)
  end subroutine monolis_precond_sor_33_clear

end module mod_monolis_precond_sor_33