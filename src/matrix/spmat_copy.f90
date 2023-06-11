module mod_monolis_spmat_copy
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc

  implicit none

contains

  subroutine monolis_copy_mat_R(mat_in, mat_out)
    implicit none
    type(monolis_structure) :: mat_in
    type(monolis_structure) :: mat_out

    call monolis_finalize(mat_out)

    call monolis_copy_mat_nonzero_pattern(mat_in, mat_out)
    call monolis_copy_mat_value_R(mat_in, mat_out)
  end subroutine monolis_copy_mat_R

  subroutine monolis_copy_mat_C(mat_in, mat_out)
    implicit none
    type(monolis_structure) :: mat_in
    type(monolis_structure) :: mat_out

    call monolis_finalize(mat_out)

    call monolis_copy_mat_nonzero_pattern(mat_in, mat_out)
    call monolis_copy_mat_value_C(mat_in, mat_out)
  end subroutine monolis_copy_mat_C

  subroutine monolis_copy_mat_nonzero_pattern(mat_in, mat_out)
    implicit none
    type(monolis_structure) :: mat_in
    type(monolis_structure) :: mat_out
  end subroutine monolis_copy_mat_nonzero_pattern

  !> 行列値のコピー（実数型）
  subroutine monolis_copy_mat_value_R(mat_in, mat_out)
    implicit none
    !> monolis 構造体（入力）
    type(monolis_structure) :: mat_in
    !> monolis 構造体（出力）
    type(monolis_structure) :: mat_out
    call monolis_copy_mat_value_matrix_R(mat_in, mat_out)
    call monolis_copy_mat_value_rhs_R(mat_in, mat_out)
    call monolis_copy_mat_value_solution_R(mat_in, mat_out)
  end subroutine monolis_copy_mat_value_R

  !> 行列値のコピー（行列、実数型）
  subroutine monolis_copy_mat_value_matrix_R(mat_in, mat_out)
    implicit none
    !> monolis 構造体（入力）
    type(monolis_structure) :: mat_in
    !> monolis 構造体（出力）
    type(monolis_structure) :: mat_out
    mat_out%MAT%R%A = mat_in%MAT%R%A
  end subroutine monolis_copy_mat_value_matrix_R

  !> 行列値のコピー（右辺ベクトル、実数型）
  subroutine monolis_copy_mat_value_rhs_R(mat_in, mat_out)
    implicit none
    !> monolis 構造体（入力）
    type(monolis_structure) :: mat_in
    !> monolis 構造体（出力）
    type(monolis_structure) :: mat_out
    mat_out%MAT%R%B = mat_in%MAT%R%B
  end subroutine monolis_copy_mat_value_rhs_R

  !> 行列値のコピー（解ベクトル、実数型）
  subroutine monolis_copy_mat_value_solution_R(mat_in, mat_out)
    implicit none
    !> monolis 構造体（入力）
    type(monolis_structure) :: mat_in
    !> monolis 構造体（出力）
    type(monolis_structure) :: mat_out
    mat_out%MAT%R%X = mat_in%MAT%R%X
  end subroutine monolis_copy_mat_value_solution_R

  !> 行列値のコピー（複素数型）
  subroutine monolis_copy_mat_value_C(mat_in, mat_out)
    implicit none
    !> monolis 構造体（入力）
    type(monolis_structure) :: mat_in
    !> monolis 構造体（出力）
    type(monolis_structure) :: mat_out
    call monolis_copy_mat_value_matrix_C(mat_in, mat_out)
    call monolis_copy_mat_value_rhs_C(mat_in, mat_out)
    call monolis_copy_mat_value_solution_C(mat_in, mat_out)
  end subroutine monolis_copy_mat_value_C

  !> 行列値のコピー（行列、複素数型）
  subroutine monolis_copy_mat_value_matrix_C(mat_in, mat_out)
    implicit none
    !> monolis 構造体（入力）
    type(monolis_structure) :: mat_in
    !> monolis 構造体（出力）
    type(monolis_structure) :: mat_out
    mat_out%MAT%C%A = mat_in%MAT%C%A
  end subroutine monolis_copy_mat_value_matrix_C

  !> 行列値のコピー（右辺ベクトル、複素数型）
  subroutine monolis_copy_mat_value_rhs_C(mat_in, mat_out)
    implicit none
    !> monolis 構造体（入力）
    type(monolis_structure) :: mat_in
    !> monolis 構造体（出力）
    type(monolis_structure) :: mat_out
    mat_out%MAT%C%B = mat_in%MAT%C%B
  end subroutine monolis_copy_mat_value_rhs_C

  !> 行列値のコピー（解ベクトル、複素数型）
  subroutine monolis_copy_mat_value_solution_C(mat_in, mat_out)
    implicit none
    !> monolis 構造体（入力）
    type(monolis_structure) :: mat_in
    !> monolis 構造体（出力）
    type(monolis_structure) :: mat_out
    mat_out%MAT%C%X = mat_in%MAT%C%X
  end subroutine monolis_copy_mat_value_solution_C

  !> 行列値のゼロ初期化（実数型）
  subroutine monolis_clear_mat_value_R(mat)
    implicit none
    !> monolis 構造体（入力）
    type(monolis_structure) :: mat
    call monolis_clear_mat_value_matrix_R(mat)
    call monolis_clear_mat_value_rhs_R(mat)
    call monolis_clear_mat_value_solution_R(mat)
  end subroutine monolis_clear_mat_value_R

  !> 行列値のゼロ初期化（行列、実数型）
  subroutine monolis_clear_mat_value_matrix_R(mat)
    implicit none
    !> monolis 構造体（入力）
    type(monolis_structure) :: mat
    mat%MAT%R%A = 0.0d0
  end subroutine monolis_clear_mat_value_matrix_R

  !> 行列値のゼロ初期化（右辺ベクトル、実数型）
  subroutine monolis_clear_mat_value_rhs_R(mat)
    implicit none
    !> monolis 構造体（入力）
    type(monolis_structure) :: mat
    mat%MAT%R%B = 0.0d0
  end subroutine monolis_clear_mat_value_rhs_R

  !> 行列値のゼロ初期化（解ベクトル、実数型）
  subroutine monolis_clear_mat_value_solution_R(mat)
    implicit none
    !> monolis 構造体（入力）
    type(monolis_structure) :: mat
    mat%MAT%R%X = 0.0d0
  end subroutine monolis_clear_mat_value_solution_R

  !> 行列値のゼロ初期化（複素数型）
  subroutine monolis_clear_mat_value_C(mat)
    implicit none
    !> monolis 構造体（入力）
    type(monolis_structure) :: mat
    call monolis_clear_mat_value_matrix_C(mat)
    call monolis_clear_mat_value_rhs_C(mat)
    call monolis_clear_mat_value_solution_C(mat)
  end subroutine monolis_clear_mat_value_C

  !> 行列値のゼロ初期化（行列、複素数型）
  subroutine monolis_clear_mat_value_matrix_C(mat)
    implicit none
    !> monolis 構造体（入力）
    type(monolis_structure) :: mat
    mat%MAT%C%A = 0.0d0
  end subroutine monolis_clear_mat_value_matrix_C

  !> 行列値のゼロ初期化（右辺ベクトル、複素数型）
  subroutine monolis_clear_mat_value_rhs_C(mat)
    implicit none
    !> monolis 構造体（入力）
    type(monolis_structure) :: mat
    mat%MAT%C%B = 0.0d0
  end subroutine monolis_clear_mat_value_rhs_C

  !> 行列値のゼロ初期化（解ベクトル、複素数型）
  subroutine monolis_clear_mat_value_solution_C(mat)
    implicit none
    !> monolis 構造体（入力）
    type(monolis_structure) :: mat
    mat%MAT%C%X = 0.0d0
  end subroutine monolis_clear_mat_value_solution_C
end module mod_monolis_spmat_copy