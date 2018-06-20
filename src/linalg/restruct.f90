module mod_monolis_restruct
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat

  implicit none
  private
  public :: monolis_restruct_matrix

contains

  subroutine monolis_restruct_matrix(monoMAT, monoMAT_reorder, perm, iperm)
    implicit none
    type(monolis_mat) :: monoMAT
    type(monolis_mat) :: monoMAT_reorder
    integer(kind=kint) :: perm(:), iperm(:)
    integer(kind=kint) :: N, NP, NPL, NPU, NDOF, NDOF2

    N = monoMAT%N
    NP = monoMAT%NP
    NPL = monoMAT%indexL(NP)
    NPU = monoMAT%indexU(NP)
    NDOF = monoMAT%NDOF
    NDOF2 = NDOF*NDOF

    monoMAT_reorder%N = N
    monoMAT_reorder%NP = NP
    monoMAT_reorder%NPL = NPL
    monoMAT_reorder%NPU = NPU
    monoMAT_reorder%NDOF = NDOF
    allocate(monoMAT_reorder%indexL(0:NP))
    allocate(monoMAT_reorder%indexU(0:NP))
    allocate(monoMAT_reorder%itemL(NPL))
    allocate(monoMAT_reorder%itemU(NPU))
    call monolis_restruct_matrix_profile(NP, perm, iperm, &
       & monoMAT%indexL, monoMAT%indexU, monoMAT%itemL, monoMAT%itemU, &
       & monoMAT_reorder%indexL, monoMAT_reorder%indexU, monoMAT_reorder%itemL, monoMAT_reorder%itemU)

    allocate(monoMAT_reorder%D (NDOF2*NP ))
    allocate(monoMAT_reorder%AL(NDOF2*NPL))
    allocate(monoMAT_reorder%AU(NDOF2*NPU))
    call monolis_restruct_matrix_values(NP, NDOF, perm, iperm, &
       & monoMAT%indexL, monoMAT%indexU, monoMAT%itemL, monoMAT%itemU, &
       & monoMAT%AL, monoMAT%AU, monoMAT%D, &
       & monoMAT_reorder%indexL, monoMAT_reorder%indexU, monoMAT_reorder%itemL, monoMAT_reorder%itemU, &
       & monoMAT_reorder%AL, monoMAT_reorder%AU, monoMAT_reorder%D)

    allocate(monoMAT_reorder%X(NDOF*NP))
    allocate(monoMAT_reorder%B(NDOF*NP))
  end subroutine monolis_restruct_matrix

  subroutine monolis_restruct_matrix_profile(N, perm, iperm, &
    & indexL, indexU, itemL, itemU, indexLp, indexUp, itemLp, itemUp)
    implicit none
    integer(kind=kint) :: N
    integer(kind=kint) :: perm(:), iperm(:)
    integer(kind=kint) :: indexL(0:), indexU(0:)
    integer(kind=kint) :: itemL(:), itemU(:)
    integer(kind=kint) :: indexLp(0:), indexUp(0:)
    integer(kind=kint) :: itemLp(:), itemUp(:)
    integer(kind=kint) :: cntL, cntU, i, in, j, jo, jn

    cntL = 0
    cntU = 0
    indexLp(0) = 0
    indexUp(0) = 0
    do i = 1, N
      in = perm(i)
      do j = indexL(in-1)+1, indexL(in)
        jo = itemL(j)
        jn = iperm(jo)
        if(jn < i)then
          cntL = cntL + 1
          itemLp(cntL) = jn
        else
          cntU = cntU + 1
          itemUp(cntU) = jn
        endif
      enddo

      do j = indexU(in-1)+1, indexU(in)
        jo = itemU(j)
        if(jo > N) cycle
        jn = iperm(jo)
        if(jn < i)then
          cntL = cntL + 1
          itemLp(cntL) = jn
        else
          cntU = cntU + 1
          itemUp(cntU) = jn
        endif
      enddo
      indexLp(i) = cntL
      indexUp(i) = cntU
      call sort_int_array(itemLp, indexLp(i-1)+1, indexLp(i))
      call sort_int_array(itemUp, indexUp(i-1)+1, indexUp(i))
    enddo
  end subroutine monolis_restruct_matrix_profile

  subroutine monolis_restruct_matrix_values(N, NDOF, perm, iperm, &
      & indexL, indexU, itemL, itemU, AL, AU, D, &
      & indexLp, indexUp, itemLp, itemUp, ALp, AUp, Dp)
    implicit none
    integer(kind=kint) :: N, NDOF
    integer(kind=kint) :: perm(:), iperm(:)
    integer(kind=kint) :: indexL(0:), indexU(0:)
    integer(kind=kint) :: itemL(:), itemU(:)
    real(kind=kdouble) :: AL(:), AU(:), D(:)
    integer(kind=kint) :: indexLp(0:), indexUp(0:)
    integer(kind=kint) :: itemLp(:), itemUp(:)
    real(kind=kdouble) :: ALp(:), AUp(:), Dp(:)
    Dp  = 0.0d0
    ALp = 0.0d0
    AUp = 0.0d0
    call reorder_diag(N, NDOF, iperm, D, Dp)
    call reorder_off_diag(N, NDOF, iperm, indexL, itemL, AL, &
       & indexLp, indexUp, itemLp, itemUp, ALp, AUp)
    call reorder_off_diag(N, NDOF, iperm, indexU, itemU, AU, &
       & indexLp, indexUp, itemLp, itemUp, ALp, AUp)
  end subroutine monolis_restruct_matrix_values

  subroutine reorder_diag(N, NDOF, iperm, D, Dp)
    implicit none
    integer(kind=kint) :: N, NDOF
    integer(kind=kint) :: iperm(:)
    real(kind=kdouble) :: D(:)
    real(kind=kdouble) :: Dp(:)
    integer(kind=kint) :: NDOF2, i, in, jn, jo, j

    NDOF2 = NDOF*NDOF
    do in = 1, N
      i = iperm(in)
      jo = (in-1)*NDOF2
      jn = (i -1)*NDOF2
      do j = 1, NDOF2
        Dp(jn + j) = D(jo + j)
      enddo
    enddo
  end subroutine reorder_diag

  subroutine reorder_off_diag(N, NDOF, iperm, indexX, itemX, AX, &
      & indexLp, indexUp, itemLp, itemUp, ALp, AUp)
    implicit none
    integer(kind=kint) :: N, NDOF
    integer(kind=kint) :: iperm(:)
    integer(kind=kint) :: indexX(0:)
    integer(kind=kint) :: itemX(:)
    real(kind=kdouble) :: AX(:)
    integer(kind=kint) :: indexLp(0:), indexUp(0:)
    integer(kind=kint) :: itemLp(:), itemUp(:)
    real(kind=kdouble) :: ALp(:), AUp(:)
    integer(kind=kint) :: NDOF2, in, i
    integer(kind=kint) :: jsnewL, jenewL, jsnewU, jenewU
    integer(kind=kint) :: jo, ko, kn, jn, lo, ln, l

    NDOF2 = NDOF*NDOF
    do in = 1, N
      i = iperm(in)
      jsnewL = indexLp(i-1)+1
      jenewL = indexLp(i)
      jsnewU = indexUp(i-1)+1
      jenewU = indexUp(i)
      do jo = indexX(in-1)+1, indexX(in)
        ko = itemX(jo)
        if (ko > N) cycle
        kn = iperm(ko)
        if(kn < i)then
          call bsearch_int_array(itemLp, jsnewL, jenewL, kn, jn)
          if(jn < 0)then
            write(*,*) "** monolis error: jn < 0 in reorder_off_diag"
          endif
          lo = (jo-1)*NDOF2
          ln = (jn-1)*NDOF2
          do l = 1, NDOF2
            ALp(ln + l) = AX(lo + l)
          enddo
        else
          call bsearch_int_array(itemUp, jsnewU, jenewU, kn, jn)
          if(jn < 0)then
            write(*,*) "** monolis error: jn < 0 in reorder_off_diag"
          endif
          lo = (jo-1)*NDOF2
          ln = (jn-1)*NDOF2
          do l = 1, NDOF2
            AUp(ln + l) = AX(lo + l)
          enddo
        endif
      enddo
    enddo
  end subroutine reorder_off_diag

  recursive subroutine sort_int_array(array, istart, iend)
    implicit none
    integer(kind=kint) :: array(:)
    integer(kind=kint) :: istart, iend
    integer(kind=kint) :: left, right, center
    integer(kind=kint) :: pivot, tmp

    if(istart >= iend) return
    center = (istart + iend) / 2
    pivot = array(center)
    left = istart
    right = iend
    do
      do while (array(left) < pivot)
        left = left + 1
      enddo
      do while (pivot < array(right))
        right = right - 1
      enddo
      if(left >= right) exit
      tmp = array(left)
      array(left) = array(right)
      array(right) = tmp
      left = left + 1
      right = right - 1
    enddo
    if(istart  < left-1) call sort_int_array(array, istart, left-1)
    if(right+1 < iend  ) call sort_int_array(array, right+1, iend)
  end subroutine sort_int_array

  subroutine bsearch_int_array(array, istart, iend, val, idx)
    implicit none
    integer(kind=kint) :: array(:)
    integer(kind=kint) :: istart, iend
    integer(kind=kint) :: val
    integer(kind=kint) :: idx
    integer(kind=kint) :: center, left, right, pivot
    left = istart
    right = iend
    do
      if(left > right)then
        idx = -1
        exit
      endif
      center = (left + right) / 2
      pivot = array(center)
      if(val < pivot)then
        right = center - 1
        cycle
      elseif(pivot < val)then
        left = center + 1
        cycle
      else
        idx = center
        exit
      endif
    enddo
  end subroutine bsearch_int_array
end module mod_monolis_restruct