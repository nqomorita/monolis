module mod_monolis_spmat_copy
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc

  implicit none

contains

  !> 行列構造体のコピー（実数型）
  subroutine monolis_copy_mat_R(mat_in, mat_out)
    implicit none
    !> monolis 構造体（入力）
    type(monolis_structure) :: mat_in
    !> monolis 構造体（出力）
    type(monolis_structure) :: mat_out

    call monolis_finalize(mat_out)

    call monolis_copy_mat_nonzero_pattern_R(mat_in, mat_out)
    call monolis_copy_mat_value_R(mat_in, mat_out)
  end subroutine monolis_copy_mat_R

  !> 行列構造体のコピー（複素数型）
  subroutine monolis_copy_mat_C(mat_in, mat_out)
    implicit none
    !> monolis 構造体（入力）
    type(monolis_structure) :: mat_in
    !> monolis 構造体（出力）
    type(monolis_structure) :: mat_out

    call monolis_finalize(mat_out)

    call monolis_copy_mat_nonzero_pattern_C(mat_in, mat_out)
    call monolis_copy_mat_value_C(mat_in, mat_out)
  end subroutine monolis_copy_mat_C

  !> 行列の非零パターンのコピー（実数型）
  subroutine monolis_copy_mat_nonzero_pattern_R(mat_in, mat_out)
    implicit none
    !> monolis 構造体（入力）
    type(monolis_structure) :: mat_in
    !> monolis 構造体（出力）
    type(monolis_structure) :: mat_out

    mat_out%MAT%N = mat_in%MAT%N
    mat_out%MAT%NP = mat_in%MAT%NP
    mat_out%MAT%NDOF = mat_in%MAT%NDOF

    call monolis_copy_mat_nonzero_pattern_CSR (mat_in%MAT%CSR,  mat_out%MAT%CSR)
    call monolis_copy_mat_nonzero_pattern_CSC (mat_in%MAT%CSC,  mat_out%MAT%CSC)
    call monolis_copy_mat_nonzero_pattern_SCSR(mat_in%MAT%SCSR, mat_out%MAT%SCSR)
    call monolis_copy_mat_nonzero_pattern_val_R(mat_in%MAT%R, mat_out%MAT%R)
  end subroutine monolis_copy_mat_nonzero_pattern_R

  !> 行列の非零パターンのコピー（複素数型）
  subroutine monolis_copy_mat_nonzero_pattern_C(mat_in, mat_out)
    implicit none
    !> monolis 構造体（入力）
    type(monolis_structure) :: mat_in
    !> monolis 構造体（出力）
    type(monolis_structure) :: mat_out

    mat_out%MAT%N = mat_in%MAT%N
    mat_out%MAT%NP = mat_in%MAT%NP
    mat_out%MAT%NDOF = mat_in%MAT%NDOF

    call monolis_copy_mat_nonzero_pattern_CSR (mat_in%MAT%CSR,  mat_out%MAT%CSR)
    call monolis_copy_mat_nonzero_pattern_CSC (mat_in%MAT%CSC,  mat_out%MAT%CSC)
    call monolis_copy_mat_nonzero_pattern_SCSR(mat_in%MAT%SCSR, mat_out%MAT%SCSR)
    call monolis_copy_mat_nonzero_pattern_val_C(mat_in%MAT%C, mat_out%MAT%C)
  end subroutine monolis_copy_mat_nonzero_pattern_C

  !> 行列の非零パターンのコピー（実数型）
  subroutine monolis_copy_mat_nonzero_pattern_val_R(mat_in, mat_out)
    implicit none
    !> monolis 構造体（入力）
    type(monolis_mat_val_R) :: mat_in
    !> monolis 構造体（出力）
    type(monolis_mat_val_R) :: mat_out

  end subroutine monolis_copy_mat_nonzero_pattern_val_R

  !> 行列の非零パターンのコピー（複素数型）
  subroutine monolis_copy_mat_nonzero_pattern_val_C(mat_in, mat_out)
    implicit none
    !> monolis 構造体（入力）
    type(monolis_mat_val_C) :: mat_in
    !> monolis 構造体（出力）
    type(monolis_mat_val_C) :: mat_out

  end subroutine monolis_copy_mat_nonzero_pattern_val_C

  !> 行列の非零パターンのコピー（CSR）
  subroutine monolis_copy_mat_nonzero_pattern_CSR(mat_in, mat_out)
    implicit none
    !> monolis 構造体（入力）
    type(monolis_mat_CSR) :: mat_in
    !> monolis 構造体（出力）
    type(monolis_mat_CSR) :: mat_out

  end subroutine monolis_copy_mat_nonzero_pattern_CSR

  !> 行列の非零パターンのコピー（CSC）
  subroutine monolis_copy_mat_nonzero_pattern_CSC(mat_in, mat_out)
    implicit none
    !> monolis 構造体（入力）
    type(monolis_mat_CSC) :: mat_in
    !> monolis 構造体（出力）
    type(monolis_mat_CSC) :: mat_out

  end subroutine monolis_copy_mat_nonzero_pattern_CSC

  !> 行列の非零パターンのコピー（SCSR）
  subroutine monolis_copy_mat_nonzero_pattern_SCSR(mat_in, mat_out)
    implicit none
    !> monolis 構造体（入力）
    type(monolis_mat_separated_CSR) :: mat_in
    !> monolis 構造体（出力）
    type(monolis_mat_separated_CSR) :: mat_out

  end subroutine monolis_copy_mat_nonzero_pattern_SCSR

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