module mod_monolis_stdlib
  use mod_monolis_prm
  implicit none

contains

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

end module mod_monolis_stdlib
