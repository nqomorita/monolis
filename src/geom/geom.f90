module mod_monolis_geom
  use mod_monolis_prm
  use mod_monolis_c3d8_shape
  use mod_monolis_shape_util

contains

  subroutine monolis_geom_get_local_position(coord, pos, x, is_converge)
    implicit none
    integer(kint) :: i, j
    real(kdouble) :: coord(3,8), pos(3), x(3)
    real(kdouble) :: jacobi(3,3), inJacob(3,3), det
    real(kdouble) :: norm, thresh, threshup, func(8,3), n(8), fr(3), dx(3)
    logical :: is_converge, is_fail

    is_converge = .false.
    thresh = 1.0d-8
    threshup = 4.0d0
    x = 0.0d0
    fr = 0.0d0
    do i = 1, 10
      call monolis_C3D8_shapefunc(x, n)
      fr = matmul(coord, n) - pos
      norm = dsqrt(fr(1)*fr(1) + fr(2)*fr(2) + fr(3)*fr(3))
      if(norm < thresh)then
        is_converge = .true.
        exit
      endif

      call monolis_C3D8_shapefunc_deriv(x, func)
      jacobi = matmul(coord, func)
      call monolis_get_inverse_matrix_3d(jacobi, inJacob, det, is_fail)
      if(is_fail) exit

      dx = - Matmul(inJacob, fr)
      do j = 1, 3
        x(j) = x(j) + dx(j)
      enddo

      norm = dsqrt(x(1)*x(1) + x(2)*x(2) + x(3)*x(3))
      if(threshup < norm)then
        exit
      endif
    enddo
  end subroutine monolis_geom_get_local_position

end module mod_monolis_geom
