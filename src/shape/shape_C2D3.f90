module mod_monolis_c2d3_shape
  use mod_monolis_prm

  private

  real(kdouble), parameter :: gsp(2,3) = reshape([ &
    0.166666666666667d0, 0.166666666666667d0, &
    0.666666666666667d0, 0.166666666666667d0, &
    0.166666666666667d0, 0.666666666666667d0  &
    ], [2,3])

  real(kdouble), parameter :: np(2,3) = reshape([ &
     0.0d0, 0.0d0, &
     1.0d0, 0.0d0, &
     0.0d0, 1.0d0  &
    ], [2,3])

  real(kdouble), parameter :: weight(3) = [ &
     0.166666666666666d0,0.166666666666666d0,0.166666666666666d0 &
    ]

    public :: monolis_C2D3_num_gauss_point
    public :: monolis_C2D3_weight
    public :: monolid_C2D3_integral_point
    public :: monolid_C2D3_node_point
    public :: monolid_C2D3_shapefunc
    public :: monolid_C2D3_shapefunc_deriv

contains

  function monolis_C2D3_num_gauss_point()
    implicit none
    integer(kint) :: monolis_C2D3_num_gauss_point
    monolis_C2D3_num_gauss_point = 3
  end function monolis_C2D3_num_gauss_point

  function monolis_C2D3_weight(i)
    implicit none
    integer(kint) :: i
    real(kdouble) :: monolis_C2D3_weight
    monolis_C2D3_weight = weight(i)
  end function monolis_C2D3_weight

  subroutine monolis_C2D3_integral_point(i, r)
    implicit none
    integer(kint) :: i
    real(kdouble) :: r(2)

    r(1) = gsp(1,i)
    r(2) = gsp(2,i)
  end subroutine monolis_C2D3_integral_point

  subroutine monolis_C2D3_node_point(i, r)
    implicit none
    integer(kint) :: i
    real(kdouble) :: r(2)

    r(1) = np(1,i)
    r(2) = np(2,i)
  end subroutine monolis_C2D3_node_point

  subroutine monolis_C2D3_shapefunc(local, func)
    implicit none
    real(kdouble) ::  local(2), func(3)

    func(1) = 1.0d0 - local(1) - local(2)
    func(2) = local(1)
    func(3) = local(2)
  end subroutine monolis_C2D3_shapefunc

  subroutine monolis_C2D3_shapefunc_deriv(local, func)
    implicit none
    real(kdouble) :: local(2), func(3,2)

    func(1,1) = -1.0d0
    func(2,1) =  1.0d0
    func(3,1) =  0.0d0

    func(1,2) = -1.0d0
    func(2,2) =  0.0d0
    func(3,2) =  1.0d0
  end subroutine monolis_C2D3_shapefunc_deriv

end module mod_monolis_c2d3_shape
