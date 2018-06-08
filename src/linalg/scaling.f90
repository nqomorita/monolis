module mod_monolis_scaling
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_linalg
  use mod_monolis_linalg_util
  implicit none

  private
  real(kind=kdouble), private, allocatable :: diag(:)
  public :: monolis_scaling_fw
  public :: monolis_scaling_bk

contains

  subroutine monolis_scaling_fw(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kind=kint) :: N, NP, NDOF, NDOF2
    integer(kind=kint) :: isL, ieL, isU, ieU, inod
    integer(kind=kint) :: i, j, k, ii, ij, ip(monoMAT%NDOF), iq(monoMAT%NDOF)
    real(kind=kdouble) :: tcomm
    real(kind=kdouble), pointer :: D(:), AL(:), AU(:), B(:)
    integer(kind=kint), pointer :: indexL(:), itemL(:), indexU(:), itemU(:)

    if(.not. monoPRM%is_scaling) return

    N = monoMAT%N
    NP = monoMAT%NP
    NDOF = monoMAT%NDOF
    NDOF2 = NDOF*NDOF
    B => monoMAT%B
    D => monoMAT%D
    AL => monoMAT%AL
    AU => monoMAT%AU
    indexL => monoMAT%indexL
    itemL => monoMAT%itemL
    indexU => monoMAT%indexU
    itemU => monoMAT%itemU

    allocate(diag(NDOF*NP))

    do i = 1, N
      do k = 1, NDOF
        diag(NDOF*(i-1)+k) = 1.0d0 / dsqrt(dabs(D(NDOF*NDOF*(i-1) + (k-1)*(NDOF+1) + 1)))
      enddo
    enddo

    call monolis_update_R(monoCOM, NDOF, diag, tcomm)

    do i = 1, NP
      do j = 1, NDOF
        ip(j) = NDOF*(i-1) + j
      enddo
      do j = 1, NDOF
        do k = 1, NDOF
          D(NDOF2*(i-1) + NDOF*(j-1) + k) = D(NDOF2*(i-1) + NDOF*(j-1) + k)*diag(ip(j))*diag(ip(k))
        enddo
      enddo

      isL = indexL(i-1) + 1
      ieL = indexL(i  )
      do k = isL, ieL
        inod = itemL(k)
        do ii = 1, NDOF
          iq(ii) = NDOF*(inod-1) + ii
        enddo
        do ii = 1, NDOF
          do ij = 1, NDOF
            AL(NDOF2*(k-1) + NDOF*(ii-1) + ij) = AL(NDOF2*(k-1) + NDOF*(ii-1) + ij)*diag(ip(ii))*diag(iq(ij))
          enddo
        enddo
      enddo

      isU = indexU(i-1) + 1
      ieU = indexU(i  )
      do k = isU, ieU
        inod = itemU(k)
        do ii = 1, NDOF
          iq(ii) = NDOF*(inod-1) + ii
        enddo
        do ii = 1, NDOF
          do ij = 1, NDOF
            AU(NDOF2*(k-1) + NDOF*(ii-1) + ij) = AU(NDOF2*(k-1) + NDOF*(ii-1) + ij)*diag(ip(ii))*diag(iq(ij))
          enddo
        enddo
      enddo
    enddo

    do i = 1, N
      do k = 1, NDOF
        B(NDOF*(i-1) + k) = B(NDOF*(i-1) + k)*diag(NDOF*(i-1) + k)
      enddo
    enddo
  end subroutine monolis_scaling_fw

  subroutine monolis_scaling_bk(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kind=kint) :: N, NP, NDOF, NDOF2
    integer(kind=kint) :: isL, ieL, isU, ieU, inod
    integer(kind=kint) :: i, j, k, ii, ij, ip(monoMAT%NDOF), iq(monoMAT%NDOF)
    real(kind=kdouble) :: tcomm
    real(kind=kdouble), pointer :: D(:), AL(:), AU(:), B(:), X(:)
    integer(kind=kint), pointer :: indexL(:), itemL(:), indexU(:), itemU(:)

    if(.not. monoPRM%is_scaling) return

    N = monoMAT%N
    NP = monoMAT%NP
    NDOF = monoMAT%NDOF
    NDOF2 = NDOF*NDOF
    B => monoMAT%B
    X => monoMAT%X
    D => monoMAT%D
    AL => monoMAT%AL
    AU => monoMAT%AU
    indexL => monoMAT%indexL
    itemL => monoMAT%itemL
    indexU => monoMAT%indexU
    itemU => monoMAT%itemU

    do i = 1, N
      do k = 1, NDOF
        X(NDOF*(i-1) + k) = X(NDOF*(i-1) + k)*diag(NDOF*(i-1) + k)
        B(NDOF*(i-1) + k) = B(NDOF*(i-1) + k)/diag(NDOF*(i-1) + k)
      enddo
    enddo

    do i = 1, NP
      do j = 1, NDOF
        ip(j) = NDOF*(i-1) + j
      enddo
      do j = 1, NDOF
        do k = 1, NDOF
          D(NDOF2*(i-1) + NDOF*(j-1) + k) = D(NDOF2*(i-1) + NDOF*(j-1) + k)/(diag(ip(j))*diag(ip(k)))
        enddo
      enddo

      isL = indexL(i-1) + 1
      ieL = indexL(i  )
      do k = isL, ieL
        inod = itemL(k)
        do ii = 1, NDOF
          iq(ii) = NDOF*(inod-1) + ii
        enddo
        do ii = 1, NDOF
          do ij = 1, NDOF
            AL(NDOF2*(k-1) + NDOF*(ii-1) + ij) = AL(NDOF2*(k-1) + NDOF*(ii-1) + ij)/(diag(ip(ii))*diag(iq(ij)))
          enddo
        enddo
      enddo

      isU = indexU(i-1) + 1
      ieU = indexU(i  )
      do k = isU, ieU
        inod = itemU(k)
        do ii = 1, NDOF
          iq(ii) = NDOF*(inod-1) + ii
        enddo
        do ii = 1, NDOF
          do ij = 1, NDOF
            AU(NDOF2*(k-1) + NDOF*(ii-1) + ij) = AU(NDOF2*(k-1) + NDOF*(ii-1) + ij)/(diag(ip(ii))*diag(iq(ij)))
          enddo
        enddo
      enddo
    enddo

    deallocate(diag)
  end subroutine monolis_scaling_bk
end module mod_monolis_scaling