module mod_monolis_c3d4_shape
  use mod_monolis_prm

  private

  real(kdouble), parameter :: gsp(3) = [ &
     0.25d0, 0.25d0, 0.25d0  &
    ]

  real(kdouble), parameter :: np(3,4) = reshape([ &
      0.0d0, 0.0d0, 0.0d0, &
      1.0d0, 0.0d0, 0.0d0, &
      0.0d0, 1.0d0, 0.0d0, &
      0.0d0, 0.0d0, 1.0d0  &
    ], [3,4])

    public :: monolis_C3D4_num_gauss_point
    public :: monolis_C3D4_weight
    public :: monolid_C3D4_integral_point
    public :: monolid_C3D4_node_point
    public :: monolid_C3D4_shapefunc
    public :: monolid_C3D4_shapefunc_deriv

contains

  function monolis_C3D4_num_gauss_point()
    implicit none
    integer(kint) :: monolis_C3D4_num_gauss_point
    monolis_C3D4_num_gauss_point = 1
  end function monolis_C3D4_num_gauss_point

  function monolis_C3D4_weight(i)
    implicit none
    integer(kint), optional :: i
    real(kdouble) :: monolis_C3D4_weight
    monolis_C3D4_weight = 0.166666666666666d0
  end function monolis_C3D4_weight

  subroutine monolis_C3D4_integral_point(i, r)
    implicit none
    integer(kint) :: i
    real(kdouble) :: r(3)

    r(1) = gsp(1)
    r(2) = gsp(2)
    r(3) = gsp(3)
  end subroutine monolis_C3D4_integral_point

  subroutine monolis_C3D4_node_point(i, r)
    implicit none
    integer(kint) :: i
    real(kdouble) :: r(3)

    r(1) = np(1,i)
    r(2) = np(2,i)
    r(3) = np(3,i)
  end subroutine monolis_C3D4_node_point

  subroutine monolis_C3D4_shapefunc(local, func)
    implicit none
    real(kdouble) :: local(3), func(4)

    func(1) = 1.0d0 - local(1) - local(2) - local(3)
    func(2) = local(1)
    func(3) = local(2)
    func(4) = local(3)
  end subroutine monolis_C3D4_shapefunc

  subroutine monolis_C3D4_shapefunc_deriv(local, func)
    implicit none
    real(kdouble) :: local(3), func(4,3)

    func(1,1) = -1.d0
    func(2,1) = 1.d0
    func(3,1) = 0.d0
    func(4,1) = 0.d0

    func(1,2) = -1.d0
    func(2,2) = 0.d0
    func(3,2) = 1.d0
    func(4,2) = 0.d0

    func(1,3) = -1.d0
    func(2,3) = 0.d0
    func(3,3) = 0.d0
    func(4,3) = 1.d0
  end subroutine monolis_C3D4_shapefunc_deriv

end module mod_monolis_c3d4_shape
