module mod_monolis_precond_rom_gs_33
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat

  implicit none

  private
  public :: monolis_precond_rom_gs_33_setup
  public :: monolis_precond_rom_gs_33_apply
  public :: monolis_precond_rom_gs_33_clear

contains

  subroutine monolis_precond_rom_gs_33_setup(monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kint) :: i, j, jS, jE, in, k, l, N
    integer(kint), pointer :: index(:), item(:)
    real(kdouble) :: T(3,3), P(3), sigma
    real(kdouble), pointer :: A(:), ALU(:)

    N =  monoMAT%N
    A => monoMAT%A
    index => monoMAT%index
    item => monoMAT%item
    sigma = 1.0d0

    allocate(monoMAT%monoTree%D(9*N), source = 0.0d0)
    ALU => monoMAT%monoTree%D

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
  end subroutine monolis_precond_rom_gs_33_setup

  subroutine monolis_precond_rom_gs_33_apply(monoMAT, B, X, maxiter)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kint) :: i, j, jS, jE, jn
    integer(kint) :: iter, maxiter
    integer(kint), pointer :: index(:)
    integer(kint), pointer :: item(:)
    real(kdouble) :: X1, X2, X3, S1, S2, S3
    real(kdouble) :: X(:), B(:)
    real(kdouble), pointer :: A(:), ALU(:)

    index => monoMAT%index
    item => monoMAT%item
    A => monoMAT%A
    ALU => monoMAT%monoTree%D

    do iter = 1, maxiter
      do i = 1, monoMAT%N
        S1 = B(3*i-2)
        S2 = B(3*i-1)
        S3 = B(3*i  )
        jS = index(i-1) + 1
        jE = index(i)
        aa:do j = jS, jE
          jn = item(j)
          if(i == jn) cycle aa
          X1 = X(3*jn-2)
          X2 = X(3*jn-1)
          X3 = X(3*jn  )
          S1 = S1 - A(9*j-8)*X1 - A(9*j-7)*X2 - A(9*j-6)*X3
          S2 = S2 - A(9*j-5)*X1 - A(9*j-4)*X2 - A(9*j-3)*X3
          S3 = S3 - A(9*j-2)*X1 - A(9*j-1)*X2 - A(9*j  )*X3
        enddo aa
        X1 = X(3*i-2)
        X2 = X(3*i-1)
        X3 = X(3*i  )
        S1 = S1                 - ALU(9*i-7)*X2 - ALU(9*i-6)*X3
        S2 = S2 - ALU(9*i-5)*X1                 - ALU(9*i-3)*X3
        S3 = S3 - ALU(9*i-2)*X1 - ALU(9*i-1)*X2
        X(3*i-2) = S1/ALU(9*i-8)
        X(3*i-1) = S2/ALU(9*i-4)
        X(3*i  ) = S3/ALU(9*i  )
      enddo
    enddo
  end subroutine monolis_precond_rom_gs_33_apply

  subroutine monolis_precond_rom_gs_33_clear(monoMAT)
    implicit none
    type(monolis_mat) :: monoMAT
    real(kdouble), pointer :: ALU(:)
    ALU => monoMAT%monoTree%D
    deallocate(ALU)
  end subroutine monolis_precond_rom_gs_33_clear

end module mod_monolis_precond_rom_gs_33