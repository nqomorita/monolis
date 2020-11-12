module mod_monolis_solver_SOR
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_precond
  use mod_monolis_matvec
  use mod_monolis_linalg
  use mod_monolis_linalg_util
  use mod_monolis_converge

  implicit none
  private
  public monolis_solver_SOR

  real(kind=kdouble), allocatable :: ALU(:)

contains

  subroutine monolis_solver_SOR(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kind=kint) :: N, NP, NDOF, NDOF2, NNDOF
    integer(kind=kint) :: iter
    real(kind=kdouble) :: tol, resid, R2, B2
    real(kind=kdouble), pointer :: B(:), X(:), ALU(:)
    real(kind=kdouble), allocatable :: R(:)
    logical :: is_converge

    ALU => monoMAT%monoTree%D
    N     = monoMAT%N
    NP    = monoMAT%NP
    NDOF  = monoMAT%NDOF
    NDOF2 = NDOF*NDOF
    NNDOF = N*NDOF
    X => monoMAT%X
    B => monoMAT%B
    tol = monoPRM%tol

    if(monoPRM%is_init_x) X = 0.0d0

    allocate(R(NDOF*NP)); R = 0.0d0

    call monolis_set_converge(monoPRM, monoCOM, monoMAT, B, B2, is_converge, monoPRM%tdotp, monoPRM%tcomm)
    if(is_converge) return
    call monolis_solver_SOR_setup(monoMAT)
    call monolis_inner_product_R(monoCOM, N, NDOF, B, B, B2, monoPRM%tdotp, monoPRM%tcomm)

    do iter = 1, monoPRM%maxiter
      call monolis_solver_SOR_matvec(monoCOM, monoMAT, NDOF, X, B, monoPRM%tspmv, monoPRM%tcomm)
      call monolis_residual(monoCOM, monoMAT, X, B, R, monoPRM%tspmv, monoPRM%tcomm)
      call monolis_inner_product_R(monoCOM, N, NDOF, R, R, R2, monoPRM%tdotp, monoPRM%tcomm)
      resid = dsqrt(R2/B2)

      if(monoCOM%myrank == 0 .and. monoPRM%show_iterlog) write (*,"(i7, 1pe16.6)") iter, resid
      if(resid <= tol) exit
    enddo

    call monolis_update_R(monoCOM, NDOF, X, monoPRM%tcomm)

    deallocate(R)
    deallocate(ALU)
  end subroutine monolis_solver_SOR

  subroutine monolis_solver_SOR_setup(monoMAT)
    implicit none
    type(monolis_mat) :: monoMAT
    integer(kind=kint) :: i, ii, j, jS, jE, in, k, l, N, NP, NDOF, NDOF2
    integer(kind=kint), pointer :: index(:), item(:)
    real(kind=kdouble), pointer :: A(:), ALU(:)
    real(kind=kdouble), allocatable :: T(:), LU(:,:)

    N     = monoMAT%N
    NP    = monoMAT%NP
    NDOF  = monoMAT%NDOF
    NDOF2 = NDOF*NDOF
    A => monoMAT%A
    index => monoMAT%index
    item => monoMAT%item

    allocate(T(NDOF))
    allocate(LU(NDOF,NDOF))
    allocate(monoMAT%monoTree%D(NDOF2*NP))
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
  end subroutine monolis_solver_SOR_setup

  subroutine monolis_solver_SOR_matvec(monoCOM, monoMAT, NDOF, X, B, tspmv, tcomm)
    implicit none
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kind=kint) :: i, j, k, l, in, N, NDOF, NDOF2, jS, jE
    integer(kind=kint), pointer :: index(:), item(:)
    real(kind=kdouble) :: X(:), B(:), XT(NDOF), YT(NDOF), DT(NDOF), WT(NDOF)
    real(kind=kdouble), pointer :: A(:), ALU(:)
    real(kind=kdouble) :: t1, t2, omega
    real(kind=kdouble) :: tspmv, tcomm

    ALU => monoMAT%monoTree%D
    N     = monoMAT%N
    NDOF  = monoMAT%NDOF
    NDOF2 = NDOF*NDOF
    A => monoMAT%A
    index => monoMAT%index
    item  => monoMAT%item
    omega = 1.0d0

    call monolis_update_R(monoCOM, monoMAT%NDOF, X, tcomm)

    do i = 1, N
      DT = 0.0d0
      do k = 1, NDOF
        XT(k) = X(NDOF*(i-1) + k)
      enddo
      do j = 1, NDOF
        do k = 1, NDOF
          DT(j) = DT(j) + A(NDOF2*(i-1) + NDOF*(j-1) + k)*XT(k)
        enddo
      enddo

      YT = 0.0d0
      jS = index(i-1) + 1
      jE = index(i  )
      do j = jS, jE
        in = item(j)
        if(in < i)then
          do k = 1, NDOF
            XT(k) = X(NDOF*(in-1) + k)
          enddo
          do k = 1, NDOF
            do l = 1, NDOF
              YT(k) = YT(k) - A(NDOF2*(j-1) + NDOF*(k-1) + l)*XT(l)
            enddo
          enddo
        endif
      enddo

      jS = index(i-1) + 1
      jE = index(i  )
      do j = jS, jE
        in = item(j)
        if(i < in)then
          do k = 1, NDOF
            XT(k) = X(NDOF*(in-1) + k)
          enddo
          do k = 1, NDOF
            do l = 1, NDOF
              YT(k) = YT(k) - A(NDOF2*(j-1) + NDOF*(k-1) + l)*XT(l)
            enddo
          enddo
        endif
      enddo

      do k = 1, NDOF
        WT(k) = omega*B(NDOF*(i-1) + k) + (1.0d0 - omega)*DT(k) + omega*YT(k)
      enddo

      do j = 2, NDOF
        do k = 1, j-1
          WT(j) = WT(j) - ALU(NDOF2*(i-1) + NDOF*(j-1) + k)*WT(k)
        enddo
      enddo
      do j = NDOF, 1, -1
        do k = NDOF, j+1, -1
          WT(j) = WT(j) - ALU(NDOF2*(i-1) + NDOF*(j-1) + k)*WT(k)
        enddo
        WT(j) = ALU(NDOF2*(i-1) + (NDOF+1)*(j-1) + 1)*WT(j)
      enddo
      do k = 1, NDOF
        X(NDOF*(i-1) + k) = WT(k)
      enddo
    enddo
  end subroutine monolis_solver_SOR_matvec
end module mod_monolis_solver_SOR
