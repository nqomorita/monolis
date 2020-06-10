module mod_monolis_c3d10_shape
  use mod_monolis_prm
  use mod_monolis_shape_util
  implicit none

  private

  real(kdouble), parameter :: gsp(3,8) = reshape([ &
    -0.577350269189626d0,-0.577350269189626d0,-0.577350269189626d0, &
     0.577350269189626d0,-0.577350269189626d0,-0.577350269189626d0, &
    -0.577350269189626d0, 0.577350269189626d0,-0.577350269189626d0, &
     0.577350269189626d0, 0.577350269189626d0,-0.577350269189626d0, &
    -0.577350269189626d0,-0.577350269189626d0, 0.577350269189626d0, &
     0.577350269189626d0,-0.577350269189626d0, 0.577350269189626d0, &
    -0.577350269189626d0, 0.577350269189626d0, 0.577350269189626d0, &
     0.577350269189626d0, 0.577350269189626d0, 0.577350269189626d0  &
    ], [3,8])

  real(kdouble), parameter :: np(3,8) = reshape([ &
     -1.0d0, -1.0d0,-1.0d0, &
      1.0d0, -1.0d0,-1.0d0, &
      1.0d0,  1.0d0,-1.0d0, &
     -1.0d0,  1.0d0,-1.0d0, &
     -1.0d0, -1.0d0, 1.0d0, &
      1.0d0, -1.0d0, 1.0d0, &
      1.0d0,  1.0d0, 1.0d0, &
     -1.0d0,  1.0d0, 1.0d0  &
    ], [3,8])

    public :: monolis_C3D10_num_gauss_point
    public :: monolis_C3D10_weight
    public :: monolis_C3D10_integral_point
    !public :: monolis_C3D10_node_point
    public :: monolis_C3D10_get_global_deriv
    public :: monolis_C3D10_shapefunc
    public :: monolis_C3D10_shapefunc_deriv

contains

  function monolis_C3D10_num_gauss_point()
    implicit none
    integer(kint) :: monolis_C3D10_num_gauss_point
    monolis_C3D10_num_gauss_point = 8
  end function monolis_C3D10_num_gauss_point

  function monolis_C3D10_weight(i)
    implicit none
    integer(kint), optional :: i
    real(kdouble) :: monolis_C3D10_weight
    monolis_C3D10_weight = 1.0d0
  end function monolis_C3D10_weight

  subroutine monolis_C3D10_integral_point(i, r)
    implicit none
    integer(kint) :: i
    real(kdouble) :: r(3)

    r(1) = gsp(1,i)
    r(2) = gsp(2,i)
    r(3) = gsp(3,i)
  end subroutine monolis_C3D10_integral_point

  !subroutine monolis_C3D10_node_point(i, r)
  !  implicit none
  !  integer(kint) :: i
  !  real(kdouble) :: r(3)

  !  r(1) = np(1,i)
  !  r(2) = np(2,i)
  !  r(3) = np(3,i)
  !end subroutine monolis_C3D10_node_point

  subroutine monolis_C3D10_get_global_deriv(node, r, dndx, det)
    implicit none
    real(kdouble) :: node(3,10), r(3), dndx(10,3), deriv(10,3)
    real(kdouble) :: xj(3,3), inv(3,3), det

    call monolis_C3D10_shapefunc_deriv(r, deriv)
    xj = matmul(node, deriv)
    call monolis_get_inverse_matrix_3d(xj, inv, det)
    dndx = matmul(deriv, inv)
  end subroutine monolis_C3D10_get_global_deriv

  subroutine monolis_C3D10_shapefunc(local, func)
    implicit none
    real(kdouble) :: local(3), func(10)
    real(kdouble) :: xi, et, st, a

    xi = local(1)
    et = local(2)
    st = local(3)
    a = 1.0d0 - xi - et - st

    func(1) = (2.0d0*a - 1.0d0)*a
    func(2) = xi*(2.0d0*xi - 1.0d0)
    func(3) = et*(2.0d0*et - 1.0d0)
    func(4) = st*(2.0d0*st - 1.0d0)
    func(5) = 4.0d0*xi*a
    func(6) = 4.0d0*xi*et
    func(7) = 4.0d0*et*a
    func(8) = 4.0d0*st*a
    func(9) = 4.0d0*xi*st
    func(10)= 4.0d0*et*st
  end subroutine monolis_C3D10_shapefunc

  subroutine monolis_C3D10_shapefunc_deriv(local, func)
    implicit none
    real(kdouble) :: local(3), func(10,3)
    real(kdouble) :: xi, et, st, a

    xi = local(1)
    et = local(2)
    st = local(3)
    a = 1.0d0 - xi - et - st

    func(1,1) = 1.0d0 - 4.0d0*a
    func(2,1) = 4.0d0*xi - 1.0d0
    func(3,1) = 0.0d0
    func(4,1) = 0.0d0
    func(5,1) = 4.0d0*(1.0d0 - 2.0d0*xi - et - st)
    func(6,1) = 4.0d0*et
    func(7,1) =-4.0d0*et
    func(8,1) =-4.0d0*st
    func(9,1) = 4.0d0*st
    func(10,1)= 0.0d0

    func(1,2) = 1.0d0 - 4.0d0*a
    func(2,2) = 0.0d0
    func(3,2) = 4.0d0*et - 1.0d0
    func(4,2) = 0.0d0
    func(5,2) =-4.0d0*xi
    func(6,2) = 4.0d0*xi
    func(7,2) = 4.0d0*(1.0d0 - xi - 2.0d0*et - st)
    func(8,2) =-4.0d0*st
    func(9,2) = 0.0d0
    func(10,2)= 4.0d0*st

    func(1,3) = 1.0d0 - 4.0d0*a
    func(2,3) = 0.0d0
    func(3,3) = 0.0d0
    func(4,3) = 4.0d0*st - 1.0d0
    func(5,3) =-4.0d0*xi
    func(6,3) = 0.0d0
    func(7,3) =-4.0d0*et
    func(8,3) = 4.0d0*(1.0d0 - xi - et - 2.0d0*st)
    func(9,3) = 4.0d0*xi
    func(10,3)= 4.0d0*et
  end subroutine monolis_C3D10_shapefunc_deriv

end module mod_monolis_c3d10_shape
