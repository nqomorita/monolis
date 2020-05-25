module mod_monolis_c2d4_shape
  use mod_monolis_prm

  private

  real(kdouble), parameter :: gsp(2,4) = reshape([ &
    -0.577350269189626d0,-0.577350269189626d0, &
     0.577350269189626d0,-0.577350269189626d0, &
    -0.577350269189626d0, 0.577350269189626d0, &
     0.577350269189626d0, 0.577350269189626d0  &
    ], [2,4])

  real(kdouble), parameter :: np(2,4) = reshape([ &
     -1.0d0, -1.0d0, &
      1.0d0, -1.0d0, &
      1.0d0,  1.0d0, &
     -1.0d0,  1.0d0  &
    ], [2,4])

    public :: monolis_C2D4_num_gauss_point
    public :: monolis_C2D4_weight
    public :: monolid_C2D4_integral_point
    public :: monolid_C2D4_node_point
    public :: monolid_C2D4_shapefunc
    public :: monolid_C2D4_shapefunc_deriv

contains

  function monolis_C2D4_num_gauss_point()
    implicit none
    integer(kint) :: monolis_C2D4_num_gauss_point
    monolis_C2D4_num_gauss_point = 4
  end function monolis_C2D4_num_gauss_point

  function monolis_C2D4_weight(i)
    implicit none
    integer(kint), optional :: i
    real(kdouble) :: monolis_C2D4_weight
    monolis_C2D4_weight = 1.0d0
  end function monolis_C2D4_weight

  subroutine monolis_C2D4_integral_point(i, r)
    implicit none
    integer(kint) :: i
    real(kdouble) :: r(2)

    r(1) = gsp(1,i)
    r(2) = gsp(2,i)
  end subroutine monolis_C2D4_integral_point

  subroutine monolis_C2D4_node_point(i, r)
    implicit none
    integer(kint) :: i
    real(kdouble) :: r(2)

    r(1) = np(1,i)
    r(2) = np(2,i)
  end subroutine monolis_C2D4_node_point

  subroutine monolis_C2D4_shapefunc(local, func)
    implicit none
    real(kdouble) :: local(2), func(4)

    func(1) = 0.25d0*(1.0d0-local(1))*(1.0d0-local(2))
    func(2) = 0.25d0*(1.0d0+local(1))*(1.0d0-local(2))
    func(3) = 0.25d0*(1.0d0+local(1))*(1.0d0+local(2))
    func(4) = 0.25d0*(1.0d0-local(1))*(1.0d0+local(2))
  end subroutine monolis_C2D4_shapefunc

  subroutine monolis_C2D4_shapefunc_deriv(local, func)
    implicit none
    real(kdouble) :: local(2), func(4,2)

    func(1,1) = -0.25d0*(1.0d0-local(2))
    func(2,1) =  0.25d0*(1.0d0-local(2))
    func(3,1) =  0.25d0*(1.0d0+local(2))
    func(4,1) = -0.25d0*(1.0d0+local(2))

    func(1,2) = -0.25d0*(1.0d0-local(1))
    func(2,2) = -0.25d0*(1.0d0+local(1))
    func(3,2) =  0.25d0*(1.0d0+local(1))
    func(4,2) =  0.25d0*(1.0d0-local(1))
  end subroutine monolis_C2D4_shapefunc_deriv

end module mod_monolis_c2d4_shape
