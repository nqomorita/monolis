module mod_monolis_shape_util
  use mod_monolis_prm

contains

  subroutine monolis_get_inverse_matrix_2d(xj, inv, det, is_fail)
    implicit none
    real(kdouble) :: xj(2,2), inv(2,2), det, detinv
    logical, optional :: is_fail

    det = xj(1,1) * xj(2,2) &
        - xj(2,1) * xj(1,2)

    if(det < 0.0d0)then
      if(present(is_fail))then
        is_fail = .true.
      else
        stop "determinant < 0.0"
      endif
    endif

    detinv = 1.0d0/det
    inv(1,1) =  xj(2,2)*detinv
    inv(1,2) = -xj(1,2)*detinv
    inv(2,1) = -xj(2,1)*detinv
    inv(2,2) =  xj(1,1)*detinv
  end subroutine monolis_get_inverse_matrix_2d

  subroutine monolis_get_inverse_matrix_3d(xj, inv, det, is_fail)
    implicit none
    real(kdouble) :: xj(3,3), inv(3,3), det, detinv
    logical, optional :: is_fail

    if(present(is_fail)) is_fail = .false.

    det = xj(1,1) * xj(2,2) * xj(3,3) &
        + xj(2,1) * xj(3,2) * xj(1,3) &
        + xj(3,1) * xj(1,2) * xj(2,3) &
        - xj(3,1) * xj(2,2) * xj(1,3) &
        - xj(2,1) * xj(1,2) * xj(3,3) &
        - xj(1,1) * xj(3,2) * xj(2,3)

    if(det < 0.0d0)then
      if(present(is_fail))then
        is_fail = .true.
      else
        stop "determinant < 0.0"
      endif
    endif

    detinv = 1.0d0/det
    inv(1,1) = detinv * ( xj(2,2)*xj(3,3) - xj(3,2)*xj(2,3))
    inv(1,2) = detinv * (-xj(1,2)*xj(3,3) + xj(3,2)*xj(1,3))
    inv(1,3) = detinv * ( xj(1,2)*xj(2,3) - xj(2,2)*xj(1,3))
    inv(2,1) = detinv * (-xj(2,1)*xj(3,3) + xj(3,1)*xj(2,3))
    inv(2,2) = detinv * ( xj(1,1)*xj(3,3) - xj(3,1)*xj(1,3))
    inv(2,3) = detinv * (-xj(1,1)*xj(2,3) + xj(2,1)*xj(1,3))
    inv(3,1) = detinv * ( xj(2,1)*xj(3,2) - xj(3,1)*xj(2,2))
    inv(3,2) = detinv * (-xj(1,1)*xj(3,2) + xj(3,1)*xj(1,2))
    inv(3,3) = detinv * ( xj(1,1)*xj(2,2) - xj(2,1)*xj(1,2))
  end subroutine monolis_get_inverse_matrix_3d

end module mod_monolis_shape_util
