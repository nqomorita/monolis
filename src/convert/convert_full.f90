module mod_monolis_convert_full
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat

  implicit none

contains

  subroutine monolis_convert_full_matrix(Nf, NDOFf, Af, thresh, &
    & N, NDOF, NPU, NPL, &
    & D, AU, AL, indexU, itemU, indexL, itemL)
    implicit none
    real(kind=kdouble), pointer :: Af(:)
    real(kind=kdouble), pointer :: D(:), AU(:), AL(:)
    integer(kind=kint), pointer :: indexU(:)
    integer(kind=kint), pointer :: indexL(:)
    integer(kind=kint), pointer :: itemU(:)
    integer(kind=kint), pointer :: itemL(:)
    integer(kind=kint) :: Nf, NDOFf
    integer(kind=kint) :: N, NDOF, NPU, NPL, NDOF2
    integer(kind=kint) :: i, j, k, jS, jE, iu, il
    real(kind=kdouble) :: thresh, temp

    if(Nf < 1 .or. NDOFf < 1 .or. (.not. associated(Af)))then
      stop "  ** monolis error: monolis_convert_full_matrix_main"
    endif

    if(associated(indexU)) deallocate(indexU)
    if(associated(indexL)) deallocate(indexL)
    allocate(indexU(0:Nf))
    allocate(indexL(0:Nf))
    indexU = 0
    indexL = 0

    N = Nf
    NDOF = NDOFf
    NDOF2 = NDOFf*NDOFf
    do i = 1, Nf
      il = 0
      iu = 0
      do j = 1, Nf
        temp = 0.0d0
        do k = 1, NDOF2
          temp = temp + dabs(Af(NDOF2*Nf*(i-1) + NDOF2*(j-1) + k))
        enddo
        if(thresh < temp)then
          if(i <  j) il = il + 1
          !if(i == j) id = id + 1
          if(j <  i) iu = iu + 1
        endif
      enddo
      indexL(i) = indexL(i-1) + il
      indexU(i) = indexU(i-1) + iu
    enddo

    NPL = indexL(N)
    NPU = indexU(N)
    if(associated(itemL)) deallocate(itemL)
    if(associated(itemU)) deallocate(itemU)
    if(associated(D)) deallocate(D)
    if(associated(AL)) deallocate(AL)
    if(associated(AU)) deallocate(AU)
    allocate(itemL(NPL))
    allocate(itemU(NPU))
    allocate(D(NDOF2*N))
    allocate(AL(NDOF2*NPL))
    allocate(AU(NDOF2*NPU))
    itemL = 0
    itemU = 0
    D = 0.0d0
    AL = 0.0d0
    AU = 0.0d0

    il = 0
    iu = 0
    do i = 1, Nf
      do j = 1, Nf
        temp = 0.0d0
        do k = 1, NDOF2
          temp = temp + dabs(Af(NDOF2*Nf*(i-1) + NDOF2*(j-1) + k))
        enddo
        if(thresh < temp)then
          if(i <  j)then
            il = il + 1
            itemL(il) = j
            do k = 1, NDOF2
              AL(NDOF2*(il-1) + k) = Af(NDOF2*Nf*(i-1) + NDOF2*(j-1) + k)
            enddo
          endif
          if(i == j)then
            do k = 1, NDOF2
              D(NDOF2*(j-1) + k) = Af(NDOF2*Nf*(i-1) + NDOF2*(j-1) + k)
            enddo
          endif
          if(j <  i)then
            iu = iu + 1
            itemU(iu) = j
            do k = 1, NDOF2
              AU(NDOF2*(iu-1) + k) = Af(NDOF2*Nf*(i-1) + NDOF2*(j-1) + k)
            enddo
          endif
        endif
      enddo
    enddo
  end subroutine monolis_convert_full_matrix
end module mod_monolis_convert_full