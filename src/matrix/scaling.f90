module mod_monolis_scaling
  use mod_monolis_prm
  use mod_monolis_mat
  use mod_monolis_linalg
  use mod_monolis_util
  implicit none

  private
  public :: monolis_scaling_fw
  public :: monolis_scaling_bk

contains

  subroutine monolis_scaling_fw(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kint) :: N, NP, NDOF, NDOF2
    integer(kint) :: inod
    integer(kint) :: i, j, jS, jE, in, k, l
    real(kdouble) :: tcomm
    real(kdouble), pointer :: A(:), X(:), B(:), diag(:)
    integer(kint), pointer :: index(:), item(:)
    real(kdouble) :: t1, t2

    if(monoPRM%is_debug) call monolis_debug_header("monolis_scaling_fw")
    if(.not. monoPRM%is_scaling) return
    t1 = monolis_get_time()

    N = monoMAT%N
    NP = monoMAT%NP
    NDOF = monoMAT%NDOF
    NDOF2 = NDOF*NDOF
    X => monoMAT%X
    B => monoMAT%B
    A => monoMAT%A
    index => monoMAT%index
    item => monoMAT%item

    allocate(monoMAT%diag(NDOF*NP))
    diag => monoMAT%diag

    do i = 1, NP
      jS = index(i-1) + 1
      jE = index(i)
      do j = jS, jE
        in = item(j)
        if(i == in)then
          do k = 1, NDOF
            diag(NDOF*(i-1)+k) = 1.0d0 / dsqrt(dabs(A(NDOF2*(j-1) + (k-1)*(NDOF+1) + 1)))
          enddo
        endif
      enddo
    enddo

    call monolis_update_R(monoCOM, NDOF, diag, tcomm)

    do i = 1, NP
      jS = index(i-1) + 1
      jE = index(i)
      do j = jS, jE
        in = item(j)
        do k = 1, NDOF
          do l = 1, NDOF
            A(NDOF2*(j-1) + NDOF*(k-1) + l) = A(NDOF2*(j-1) + NDOF*(k-1) + l)*diag(NDOF*(i-1)+k)*diag(NDOF*(in-1)+l)
          enddo
        enddo
      enddo
    enddo

    do i = 1, N
      do k = 1, NDOF
        B(NDOF*(i-1) + k) = B(NDOF*(i-1) + k)*diag(NDOF*(i-1) + k)
      enddo
    enddo

    if(.not. monoPRM%is_init_x)then
      do i = 1, N
        do k = 1, NDOF
          X(NDOF*(i-1) + k) = X(NDOF*(i-1) + k)/diag(NDOF*(i-1) + k)
        enddo
      enddo
    endif

    t2 = monolis_get_time()
    monoPRM%tprep = monoPRM%tprep + t2 - t1
  end subroutine monolis_scaling_fw

  subroutine monolis_scaling_bk(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kint) :: N, NP, NDOF, NDOF2
    integer(kint) :: i, j, k, l, in, jS, jE
    real(kdouble) :: tcomm
    real(kdouble), pointer :: A(:), B(:), X(:), diag(:)
    integer(kint), pointer :: index(:), item(:)
    real(kdouble) :: t1, t2

    if(monoPRM%is_debug) call monolis_debug_header("monolis_scaling_bk")
    if(.not. monoPRM%is_scaling) return
    t1 = monolis_get_time()

    N = monoMAT%N
    NP = monoMAT%NP
    NDOF = monoMAT%NDOF
    NDOF2 = NDOF*NDOF
    B => monoMAT%B
    X => monoMAT%X
    A => monoMAT%A
    index => monoMAT%index
    item => monoMAT%item
    diag => monoMAT%diag

    do i = 1, NP
      jS = index(i-1) + 1
      jE = index(i)
      do j = jS, jE
        in = item(j)
        do k = 1, NDOF
          do l = 1, NDOF
            A(NDOF2*(j-1) + NDOF*(k-1) + l) = A(NDOF2*(j-1) + NDOF*(k-1) + l)/(diag(NDOF*(i-1)+k)*diag(NDOF*(in-1)+l))
          enddo
        enddo
      enddo
    enddo

    do i = 1, N
      do k = 1, NDOF
        X(NDOF*(i-1) + k) = X(NDOF*(i-1) + k)*diag(NDOF*(i-1) + k)
        B(NDOF*(i-1) + k) = B(NDOF*(i-1) + k)/diag(NDOF*(i-1) + k)
      enddo
    enddo

    deallocate(diag)
    t2 = monolis_get_time()
    monoPRM%tprep = monoPRM%tprep + t2 - t1
  end subroutine monolis_scaling_bk
end module mod_monolis_scaling