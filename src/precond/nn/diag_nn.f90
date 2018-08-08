module mod_monolis_precond_diag_nn
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat

  implicit none

  private
  public :: monolis_precond_diag_nn_setup
  public :: monolis_precond_diag_nn_apply
  public :: monolis_precond_diag_nn_clear

contains

  subroutine monolis_precond_diag_nn_setup(monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kind=kint) :: i, ii, j, jS, jE, in, k, l, N, NDOF, NDOF2
    integer(kind=kint), pointer :: index(:), item(:)
    real(kind=kdouble), allocatable :: T(:), LU(:,:)
    real(kind=kdouble), pointer :: A(:), ALU(:)

    N =  monoMAT%N
    NDOF  = monoMAT%NDOF
    NDOF2 = NDOF*NDOF
    A => monoMAT%A
    index => monoMAT%index
    item => monoMAT%item

    allocate(T(NDOF))
    allocate(LU(NDOF,NDOF))
    allocate(monoMAT%monoTree%D(NDOF2*N))
    ALU => monoMAT%monoTree%D
    T   = 0.0d0
    ALU = 0.0d0
    LU  = 0.0d0

    do i = 1, N
      jS = index(i-1) + 1
      jE = index(i)
      do ii = jS, jE
        in = item(ii)
        if(i == in)then
          do j = 1, NDOF
            do k = 1, NDOF
              LU(j,k) = A(NDOF2*(i-1) + NDOF*(j-1) + k)
            enddo
          enddo

          do k = 1, NDOF
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

    deallocate(T)
    deallocate(LU)
  end subroutine monolis_precond_diag_nn_setup

  subroutine monolis_precond_diag_nn_apply(monoMAT, X, Y)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kind=kint) :: i, j, k, N, NDOF, NDOF2
    real(kind=kdouble) :: X(:), Y(:)
    real(kind=kdouble), pointer :: ALU(:)
    real(kind=kdouble), allocatable :: T(:)

    N     = monoMAT%N
    NDOF  = monoMAT%NDOF
    NDOF2 = NDOF*NDOF
    ALU => monoMAT%monoTree%D

    allocate(T(NDOF))
    T = 0.0d0

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

    deallocate(T)
  end subroutine monolis_precond_diag_nn_apply

  subroutine monolis_precond_diag_nn_clear(monoMAT)
    implicit none
    type(monolis_mat) :: monoMAT
    real(kind=kdouble), pointer :: ALU(:)
    ALU => monoMAT%monoTree%D
    deallocate(ALU)
  end subroutine monolis_precond_diag_nn_clear

end module mod_monolis_precond_diag_nn
