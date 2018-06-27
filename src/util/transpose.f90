module mod_monolis_transpose
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat

  implicit none

contains

  subroutine monolis_transpose_fw(N, NP, NDOF, NPU, NPL, D, AU, AL, Dt, AUt, ALt)
    implicit none
    integer(kind=kint) :: i, j, jn
    integer(kind=kint) :: N, NP, NDOF, NDOF2, NPU, NPL
    integer(kind=kint), allocatable :: index(:)
    real(kind=kdouble), pointer :: D(:)
    real(kind=kdouble), pointer :: AU(:)
    real(kind=kdouble), pointer :: AL(:)
    real(kind=kdouble), pointer :: Dt(:)
    real(kind=kdouble), pointer :: AUt(:)
    real(kind=kdouble), pointer :: ALt(:)

    NDOF2 = NDOF*NDOF
    allocate(index(NDOF2))
    allocate(Dt (NDOF2*NP ))
    allocate(AUt(NDOF2*NPU))
    allocate(ALt(NDOF2*NPL))
    Dt  = 0.0d0
    AUt = 0.0d0
    ALt = 0.0d0
    index = 0

    do i = 1, NDOF
      do j = 1, NDOF
        index(NDOF*(i-1) + j) = NDOF*(j-1) + i
      enddo
    enddo

    do i = 1, NP
      do j = 1, NDOF2
        jn = index(j)
        Dt(NDOF2*(i-1) + j) = D(NDOF2*(i-1) + jn)
      enddo
    enddo

    do i = 1, NPU
      do j = 1, NDOF2
        jn = index(j)
        AUt(NDOF2*(i-1) + j) = AU(NDOF2*(i-1) + jn)
      enddo
    enddo

    do i = 1, NPL
      do j = 1, NDOF2
        jn = index(j)
        ALt(NDOF2*(i-1) + j) = AL(NDOF2*(i-1) + jn)
      enddo
    enddo
  end subroutine monolis_transpose_fw

  subroutine monolis_transpose_bk(Dt, AUt, ALt)
    implicit none
    integer(kind=kint) :: i
    real(kind=kdouble), pointer :: Dt(:)
    real(kind=kdouble), pointer :: AUt(:)
    real(kind=kdouble), pointer :: ALt(:)

    deallocate(Dt)
    deallocate(AUt)
    deallocate(ALt)
  end subroutine monolis_transpose_bk
end module mod_monolis_transpose