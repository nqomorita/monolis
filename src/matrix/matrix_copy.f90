module mod_monolis_matrix_copy
  use mod_monolis_mat
  use mod_monolis_com
  use mod_monolis_util
  use mod_monolis_stdlib
  implicit none

contains

  subroutine monolis_copy_mat_by_pointer(min, mout)
    implicit none
    type(monolis_mat) :: min
    type(monolis_mat) :: mout

    mout%N = min%N
    mout%NP = min%NP
    mout%NZ = min%NZ
    mout%NDOF = min%NDOF
    mout%index => min%index
    mout%item => min%item
    mout%indexR => min%indexR
    mout%itemR => min%itemR
    mout%permR => min%permR
    mout%A => min%A
    mout%X => min%X
    mout%B => min%B
  end subroutine monolis_copy_mat_by_pointer

  subroutine monolis_copy_mat_profile(min, mout)
    implicit none
    type(monolis_structure) :: min
    type(monolis_structure) :: mout
    integer(kint) :: ndof2

    call monolis_com_copy(min%COM, mout%COM)

    mout%MAT%index => min%MAT%index
    mout%MAT%item => min%MAT%item
    mout%MAT%indexR => min%MAT%indexR
    mout%MAT%itemR => min%MAT%itemR
    mout%MAT%permR => min%MAT%permR

    mout%MAT%n = min%MAT%n
    mout%MAT%np = min%MAT%np
    mout%MAT%ndof = min%MAT%ndof

    ndof2 = min%MAT%ndof*min%MAT%ndof
    allocate(mout%MAT%A(ndof2*min%MAT%index(min%MAT%np)), source = 0.0d0)
    allocate(mout%MAT%X(min%MAT%np*min%MAT%ndof), source = 0.0d0)
    allocate(mout%MAT%B(min%MAT%np*min%MAT%ndof), source = 0.0d0)
  end subroutine monolis_copy_mat_profile

  subroutine monolis_clear_mat_value(mat)
    implicit none
    type(monolis_structure) :: mat

    call monolis_clear_A_value(mat)
    call monolis_clear_X_value(mat)
    call monolis_clear_B_value(mat)
  end subroutine monolis_clear_mat_value

  subroutine monolis_clear_A_value(mat)
    implicit none
    type(monolis_structure) :: mat
    mat%MAT%A = 0.0d0
  end subroutine monolis_clear_A_value

  subroutine monolis_clear_X_value(mat)
    implicit none
    type(monolis_structure) :: mat
    mat%MAT%X = 0.0d0
  end subroutine monolis_clear_X_value

  subroutine monolis_clear_B_value(mat)
    implicit none
    type(monolis_structure) :: mat
    mat%MAT%B = 0.0d0
  end subroutine monolis_clear_B_value

  subroutine monolis_copy_mat_all(min, mout)
    implicit none
    type(monolis_mat) :: min
    type(monolis_mat) :: mout
    integer(kint) :: i, NZ

    mout%N = min%N
    mout%NP = min%NP
    mout%NZ = min%NZ
    mout%NDOF = min%NDOF

    NZ = min%index(min%NP)
    allocate(mout%index(0:min%NP))
    allocate(mout%item(NZ))
    allocate(mout%indexR(0:min%NP))
    allocate(mout%itemR(NZ))
    allocate(mout%permR(NZ))
    allocate(mout%A(min%NDOF*min%NDOF*NZ))
    allocate(mout%X(min%NDOF*min%NP))
    allocate(mout%B(min%NDOF*min%NP))

    mout%index(0:min%NP) = min%index(0:min%NP)
    mout%item = min%item
    mout%indexR(0:min%NP) = min%indexR(0:min%NP)
    mout%itemR = min%itemR
    mout%permR = min%permR
    mout%A = min%A
    mout%X = min%X
    mout%B = min%B
  end subroutine monolis_copy_mat_all

end module mod_monolis_matrix_copy
