module mod_monolis_stdlib
  use mod_monolis_prm
  implicit none

contains

  function monolis_normalize_cross_product_3d(a, b)
    implicit none
    real(kdouble) :: monolis_normalize_cross_product_3d(3), a(3), b(3), norm(3), l2

    norm(1) = a(2)*b(3) - a(3)*b(2)
    norm(2) = a(3)*b(1) - a(1)*b(3)
    norm(3) = a(1)*b(2) - a(2)*b(1)
    l2 = dsqrt(norm(1)*norm(1) + norm(2)*norm(2) + norm(3)*norm(3))
    if(l2 == 0.0d0)then
      monolis_normalize_cross_product_3d = 0.0d0
      return
    endif
    monolis_normalize_cross_product_3d = norm/l2
  end function monolis_normalize_cross_product_3d

  function monolis_normalize_vector(n, a)
    implicit none
    integer(kint) :: i, n
    real(kdouble) :: monolis_normalize_vector(n), a(n), l2

    do i = 1, n
      l2 = a(i)*a(i)
    enddo
    l2 = dsqrt(l2)
    if(l2 == 0.0d0)then
      monolis_normalize_vector = 0.0d0
      return
    endif
    monolis_normalize_vector = a/l2
  end function monolis_normalize_vector

  function monolis_get_l2_norm(n, a)
    implicit none
    integer(kint) :: i, n
    real(kdouble) :: monolis_get_l2_norm(n), a(n), l2

    do i = 1, n
      l2 = a(i)*a(i)
    enddo
    monolis_get_l2_norm = dsqrt(l2)
  end function monolis_get_l2_norm

  subroutine monolis_get_inverse_matrix(n, a, inv)
    implicit none
    integer(kint) :: n, i, j, k
    real(kdouble) :: a(n,n), inv(n,n), tmp

    inv = 0.0d0
    do i = 1, n
      inv(i,i) = 1.0d0
    enddo

    do i = 1, n
      tmp = 1.0d0/a(i,i)
      do j = 1, n
          a(j,i) =   a(j,i) * tmp
        inv(j,i) = inv(j,i) * tmp
      enddo
      do j = 1, n
        if(i /= j) then
          tmp = a(i,j)
          do k = 1, n
              a(k,j) =   a(k,j) -   a(k,i) * tmp
            inv(k,j) = inv(k,j) - inv(k,i) * tmp
          enddo
        endif
      enddo
    enddo
  end subroutine monolis_get_inverse_matrix

  recursive subroutine monolis_qsort_int(array, iS, iE)
    implicit none
    integer(kint), intent(inout) :: array(:)
    integer(kint), intent(in) :: iS, iE
    integer(kint) :: pivot, center, left, right, tmp

    if (iS >= iE) return

    center = (iS + iE) / 2
    pivot = array(center)
    left = iS
    right = iE

    do
      do while (array(left) < pivot)
        left = left + 1
      enddo
      do while (pivot < array(right))
        right = right - 1
      enddo

      if (left >= right) exit

      tmp = array(left)
      array(left) = array(right)
      array(right) = tmp

      left = left + 1
      right = right - 1
    enddo

    if(iS < left-1)  call monolis_qsort_int(array, iS, left-1)
    if(right+1 < iE) call monolis_qsort_int(array, right+1, iE)
  end subroutine monolis_qsort_int

  recursive subroutine monolis_qsort_int_with_perm(array, iS, iE, perm)
    implicit none
    integer(kint), intent(inout) :: array(:), perm(:)
    integer(kint), intent(in) :: iS, iE
    integer(kint) :: pivot, center, left, right, tmp

    if (iS >= iE) return

    center = (iS + iE) / 2
    pivot = array(center)
    left = iS
    right = iE

    do
      do while (array(left) < pivot)
        left = left + 1
      enddo
      do while (pivot < array(right))
        right = right - 1
      enddo

      if (left >= right) exit

      tmp = array(left)
      array(left) = array(right)
      array(right) = tmp

      tmp = perm(left)
      perm(left) = perm(right)
      perm(right) = tmp

      left = left + 1
      right = right - 1
    enddo

    if(iS < left-1)  call monolis_qsort_int_with_perm(array, iS, left-1, perm)
    if(right+1 < iE) call monolis_qsort_int_with_perm(array, right+1, iE, perm)
  end subroutine monolis_qsort_int_with_perm

  subroutine monolis_bsearch_int(array, iS, iE, val, idx)
    implicit none
    integer(kint), intent(in) :: array(:)
    integer(kint), intent(in) :: iS, iE
    integer(kint), intent(in) :: val
    integer(kint), intent(out) :: idx
    integer(kint) :: center, left, right, pivot

    left = iS
    right = iE

    do
      if (left > right) then
        idx = -1
        exit
      endif
      center = (left + right) / 2
      pivot = array(center)
      if (val < pivot) then
        right = center - 1
        cycle
      elseif (pivot < val) then
        left = center + 1
        cycle
      else
        idx = center
        exit
      endif
    enddo
  end subroutine monolis_bsearch_int

  subroutine monolis_bsearch_int_with_position(array, iS, iE, val, idx, pos)
    implicit none
    integer(kind=kint), intent(in) :: array(:)
    integer(kind=kint), intent(in) :: iS, iE
    integer(kind=kint), intent(in) :: val
    integer(kind=kint), intent(out) :: idx, pos
    integer(kind=kint) :: center, left, right, pivot

    pos = -1
    left = iS
    right = iE

    do
      if (left > right) then
        idx = -1
        pos = right
        exit
      endif
      center = (left + right) / 2
      pivot = array(center)
      if (val < pivot) then
        right = center - 1
        cycle
      elseif (pivot < val) then
        left = center + 1
        cycle
      else
        idx = center
        pos = -1
        exit
      endif
    enddo
  end subroutine monolis_bsearch_int_with_position

  subroutine monolis_uniq_int(array, len, newlen)
    implicit none
    integer(kind=kint), intent(inout) :: array(:)
    integer(kind=kint), intent(in) :: len
    integer(kind=kint), intent(out) :: newlen
    integer(kind=kint) :: i, ndup

    ndup = 0
    do i = 2, len
      if(array(i) == array(i - 1 - ndup))then
        ndup = ndup + 1
      elseif(ndup > 0)then
        array(i - ndup) = array(i)
      endif
    end do
    newlen = len - ndup
  end subroutine monolis_uniq_int

end module mod_monolis_stdlib
