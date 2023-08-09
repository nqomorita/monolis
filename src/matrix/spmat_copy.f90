module mod_monolis_spmat_copy
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc

  implicit none

contains

  !> @ingroup matrix_copy
  !> 行列構造体のコピー（実数型）
  subroutine monolis_copy_mat_R(mat_in, mat_out)
    implicit none
    !> [in] monolis 構造体（入力）
    type(monolis_structure), intent(in) :: mat_in
    !> [in,out] monolis 構造体（出力）
    type(monolis_structure), intent(inout) :: mat_out

    call monolis_finalize(mat_out)

    call monolis_copy_mat_nonzero_pattern_R(mat_in, mat_out)
    call monolis_copy_mat_value_R(mat_in, mat_out)
  end subroutine monolis_copy_mat_R

  !> @ingroup matrix_copy
  !> 行列構造体のコピー（複素数型）
  subroutine monolis_copy_mat_C(mat_in, mat_out)
    implicit none
    !> [in] monolis 構造体（入力）
    type(monolis_structure), intent(in) :: mat_in
    !> [in,out] monolis 構造体（出力）
    type(monolis_structure), intent(inout) :: mat_out

    call monolis_finalize(mat_out)

    call monolis_copy_mat_nonzero_pattern_C(mat_in, mat_out)
    call monolis_copy_mat_value_C(mat_in, mat_out)
  end subroutine monolis_copy_mat_C

  !> @ingroup matrix_copy
  !> 行列の非零パターンのコピー（実数型）
  subroutine monolis_copy_mat_nonzero_pattern_R(mat_in, mat_out)
    implicit none
    !> [in] monolis 構造体（入力）
    type(monolis_structure), intent(in) :: mat_in
    !> [in,out] monolis 構造体（出力）
    type(monolis_structure), intent(inout) :: mat_out
    integer(kint) :: NP, NZ, NZU, NZL

    mat_out%MAT%N = mat_in%MAT%N
    mat_out%MAT%NP = mat_in%MAT%NP
    mat_out%MAT%NDOF = mat_in%MAT%NDOF

    call monolis_copy_mat_nonzero_pattern_CSR (mat_in%MAT%NP, mat_in%MAT%CSR,  mat_out%MAT%CSR)
    call monolis_copy_mat_nonzero_pattern_CSC (mat_in%MAT%NP, mat_in%MAT%CSC,  mat_out%MAT%CSC)
    call monolis_copy_mat_nonzero_pattern_SCSR(mat_in%MAT%NP, mat_in%MAT%SCSR, mat_out%MAT%SCSR)

    NP = mat_in%MAT%NP
    NZ = mat_in%MAT%CSR%index(NP + 1)

    NZU = 0
    if(associated(mat_in%MAT%SCSR%indexU))then
      NZU = mat_in%MAT%SCSR%indexU(NP + 1)
    endif

    NZL = 0
    if(associated(mat_in%MAT%SCSR%indexL))then
      NZL = mat_in%MAT%SCSR%indexL(NP + 1)
    endif

    call monolis_copy_mat_nonzero_pattern_val_R(mat_in%MAT%NP, mat_in%MAT%NDOF, NZ, NZU, NZL, &
      & mat_in%MAT%R, mat_out%MAT%R)
  end subroutine monolis_copy_mat_nonzero_pattern_R

  !> @ingroup matrix_copy
  !> 行列の非零パターンのコピー（複素数型）
  subroutine monolis_copy_mat_nonzero_pattern_C(mat_in, mat_out)
    implicit none
    !> [in] monolis 構造体（入力）
    type(monolis_structure), intent(in) :: mat_in
    !> [in,out] monolis 構造体（出力）
    type(monolis_structure), intent(inout) :: mat_out
    integer(kint) :: NP, NZ, NZU, NZL

    mat_out%MAT%N = mat_in%MAT%N
    mat_out%MAT%NP = mat_in%MAT%NP
    mat_out%MAT%NDOF = mat_in%MAT%NDOF

    call monolis_copy_mat_nonzero_pattern_CSR (mat_in%MAT%NP, mat_in%MAT%CSR,  mat_out%MAT%CSR)
    call monolis_copy_mat_nonzero_pattern_CSC (mat_in%MAT%NP, mat_in%MAT%CSC,  mat_out%MAT%CSC)
    call monolis_copy_mat_nonzero_pattern_SCSR(mat_in%MAT%NP, mat_in%MAT%SCSR, mat_out%MAT%SCSR)

    NP = mat_in%MAT%NP
    NZ = mat_in%MAT%CSR%index(NP + 1)

    NZU = 0
    if(associated(mat_in%MAT%SCSR%indexU))then
      NZU = mat_in%MAT%SCSR%indexU(NP + 1)
    endif

    NZL = 0
    if(associated(mat_in%MAT%SCSR%indexL))then
      NZL = mat_in%MAT%SCSR%indexL(NP + 1)
    endif

    call monolis_copy_mat_nonzero_pattern_val_C(mat_in%MAT%NP, mat_in%MAT%NDOF, NZ, NZU, NZL, &
      & mat_in%MAT%C, mat_out%MAT%C)
  end subroutine monolis_copy_mat_nonzero_pattern_C

  !> @ingroup matrix_copy
  !> 行列の非零値のコピー（実数型）
  subroutine monolis_copy_mat_nonzero_pattern_val_R(NP, NDOF, NZ, NZU, NZL, mat_in, mat_out)
    implicit none
    !> [in] 全計算点数
    integer(kint), intent(in) :: NP
    !> [in] 計算点が持つ自由度
    integer(kint), intent(in) :: NDOF
    !> [in] 非零要素数
    integer(kint), intent(in) :: NZ
    !> [in] 非零要素数（上三角）
    integer(kint), intent(in) :: NZU
    !> [in] 非零要素数（下三角）
    integer(kint), intent(in) :: NZL
    !> [in] monolis 構造体（入力）
    type(monolis_mat_val_R), intent(in) :: mat_in
    !> [in,out] monolis 構造体（出力）
    type(monolis_mat_val_R), intent(inout) :: mat_out

    call monolis_pdealloc_R_1d(mat_out%X)
    call monolis_palloc_R_1d(mat_out%X, NP*NDOF)
    mat_out%X = mat_in%X

    call monolis_pdealloc_R_1d(mat_out%B)
    call monolis_palloc_R_1d(mat_out%B, NP*NDOF)
    mat_out%B = mat_in%B

    call monolis_pdealloc_R_1d(mat_out%A)
    call monolis_palloc_R_1d(mat_out%A, NZ*NDOF*NDOF)
    mat_out%A = mat_in%A

    if(associated(mat_in%U))then
      call monolis_pdealloc_R_1d(mat_out%U)
      call monolis_palloc_R_1d(mat_out%U, NZU*NDOF*NDOF)
      mat_out%U = mat_in%U
    endif

    if(associated(mat_in%L))then
      call monolis_pdealloc_R_1d(mat_out%L)
      call monolis_palloc_R_1d(mat_out%L, NZL*NDOF*NDOF)
      mat_out%L = mat_in%L
    endif

    if(associated(mat_in%D))then
      call monolis_pdealloc_R_1d(mat_out%D)
      call monolis_palloc_R_1d(mat_out%D, NP*NDOF*NDOF)
      mat_out%D = mat_in%D
    endif
  end subroutine monolis_copy_mat_nonzero_pattern_val_R

  !> @ingroup matrix_copy
  !> 行列の非零パターンのコピー（複素数型）
  subroutine monolis_copy_mat_nonzero_pattern_val_C(NP, NDOF, NZ, NZU, NZL, mat_in, mat_out)
    implicit none
    !> [in] 全計算点数
    integer(kint), intent(in) :: NP
    !> [in] 計算点が持つ自由度
    integer(kint), intent(in) :: NDOF
    !> [in] 非零要素数
    integer(kint), intent(in) :: NZ
    !> [in] 非零要素数（上三角）
    integer(kint), intent(in) :: NZU
    !> [in] 非零要素数（下三角）
    integer(kint), intent(in) :: NZL
    !> [in] monolis 構造体（入力）
    type(monolis_mat_val_C), intent(in) :: mat_in
    !> [in,out] monolis 構造体（出力）
    type(monolis_mat_val_C), intent(inout) :: mat_out

    call monolis_pdealloc_C_1d(mat_out%X)
    call monolis_palloc_C_1d(mat_out%X, NP*NDOF)
    mat_out%X = mat_in%X

    call monolis_pdealloc_C_1d(mat_out%B)
    call monolis_palloc_C_1d(mat_out%B, NP*NDOF)
    mat_out%B = mat_in%B

    call monolis_pdealloc_C_1d(mat_out%A)
    call monolis_palloc_C_1d(mat_out%A, NZ*NDOF*NDOF)
    mat_out%A = mat_in%A

    if(associated(mat_in%U))then
      call monolis_pdealloc_C_1d(mat_out%U)
      call monolis_palloc_C_1d(mat_out%U, NZU*NDOF*NDOF)
      mat_out%U = mat_in%U
    endif

    if(associated(mat_in%L))then
      call monolis_pdealloc_C_1d(mat_out%L)
      call monolis_palloc_C_1d(mat_out%L, NZL*NDOF*NDOF)
      mat_out%L = mat_in%L
    endif

    if(associated(mat_in%D))then
      call monolis_pdealloc_C_1d(mat_out%D)
      call monolis_palloc_C_1d(mat_out%D, NP*NDOF*NDOF)
      mat_out%D = mat_in%D
    endif
  end subroutine monolis_copy_mat_nonzero_pattern_val_C

  !> @ingroup matrix_copy
  !> 行列の非零パターンのコピー（CSR）
  subroutine monolis_copy_mat_nonzero_pattern_CSR(NP, mat_in, mat_out)
    implicit none
    !> [in] 全計算点数
    integer(kint), intent(in) :: NP
    !> [in] monolis 構造体（入力）
    type(monolis_mat_CSR), intent(in) :: mat_in
    !> [in,out] monolis 構造体（出力）
    type(monolis_mat_CSR), intent(inout) :: mat_out
    integer(kint) :: NZ

    NZ = mat_in%index(NP + 1)

    call monolis_pdealloc_I_1d(mat_out%index)
    call monolis_palloc_I_1d(mat_out%index, NP + 1)
    mat_out%index = mat_in%index

    call monolis_pdealloc_I_1d(mat_out%item)
    call monolis_palloc_I_1d(mat_out%item, NZ)
    mat_out%item = mat_in%item
  end subroutine monolis_copy_mat_nonzero_pattern_CSR

  !> @ingroup matrix_copy
  !> 行列の非零パターンのコピー（CSC）
  subroutine monolis_copy_mat_nonzero_pattern_CSC(NP, mat_in, mat_out)
    implicit none
    !> [in] 全計算点数
    integer(kint), intent(in) :: NP
    !> [in] monolis 構造体（入力）
    type(monolis_mat_CSC), intent(in) :: mat_in
    !> [in,out] monolis 構造体（出力）
    type(monolis_mat_CSC), intent(inout) :: mat_out
    integer(kint) :: NZ

    NZ = mat_in%index(NP + 1)

    call monolis_pdealloc_I_1d(mat_out%index)
    call monolis_palloc_I_1d(mat_out%index, NP + 1)
    mat_out%index = mat_in%index

    call monolis_pdealloc_I_1d(mat_out%item)
    call monolis_palloc_I_1d(mat_out%item, NZ)
    mat_out%item = mat_in%item

    call monolis_pdealloc_I_1d(mat_out%perm)
    call monolis_palloc_I_1d(mat_out%perm, NZ)
    mat_out%perm = mat_in%perm
  end subroutine monolis_copy_mat_nonzero_pattern_CSC

  !> @ingroup matrix_copy
  !> 行列の非零パターンのコピー（SCSR）
  subroutine monolis_copy_mat_nonzero_pattern_SCSR(NP, mat_in, mat_out)
    implicit none
    !> [in] 全計算点数
    integer(kint), intent(in) :: NP
    !> [in] monolis 構造体（入力）
    type(monolis_mat_separated_CSR), intent(in) :: mat_in
    !> [in,out] monolis 構造体（出力）
    type(monolis_mat_separated_CSR), intent(inout) :: mat_out
    integer(kint) :: NZ

    if(associated(mat_in%indexU))then
      call monolis_pdealloc_I_1d(mat_out%indexU)
      call monolis_palloc_I_1d(mat_out%indexU, NP + 1)
      mat_out%indexU = mat_in%indexU
    endif

    if(associated(mat_in%itemU))then
      NZ = mat_in%indexU(NP + 1)
      call monolis_pdealloc_I_1d(mat_out%itemU)
      call monolis_palloc_I_1d(mat_out%itemU, NZ)
      mat_out%itemU = mat_in%itemU
    endif

    if(associated(mat_in%indexL))then
      call monolis_pdealloc_I_1d(mat_out%indexL)
      call monolis_palloc_I_1d(mat_out%indexL, NP)
      mat_out%indexL = mat_in%indexL
    endif

    if(associated(mat_in%itemL))then
      NZ = mat_in%indexL(NP + 1)
      call monolis_pdealloc_I_1d(mat_out%itemL)
      call monolis_palloc_I_1d(mat_out%itemL, NZ)
      mat_out%itemL = mat_in%itemL
    endif
  end subroutine monolis_copy_mat_nonzero_pattern_SCSR

  !> @ingroup matrix_copy
  !> 行列値のコピー（実数型）
  subroutine monolis_copy_mat_value_R(mat_in, mat_out)
    implicit none
    !> [in] monolis 構造体（入力）
    type(monolis_structure), intent(in) :: mat_in
    !> [in,out] monolis 構造体（出力）
    type(monolis_structure), intent(inout) :: mat_out
    call monolis_copy_mat_value_matrix_R(mat_in, mat_out)
    call monolis_copy_mat_value_rhs_R(mat_in, mat_out)
    call monolis_copy_mat_value_solution_R(mat_in, mat_out)
  end subroutine monolis_copy_mat_value_R

  !> @ingroup matrix_copy
  !> 行列値のコピー（行列、実数型）
  subroutine monolis_copy_mat_value_matrix_R(mat_in, mat_out)
    implicit none
    !> [in] monolis 構造体（入力）
    type(monolis_structure), intent(in) :: mat_in
    !> [in,out] monolis 構造体（出力）
    type(monolis_structure), intent(inout) :: mat_out
    mat_out%MAT%R%A = mat_in%MAT%R%A
  end subroutine monolis_copy_mat_value_matrix_R

  !> @ingroup matrix_copy
  !> 行列値のコピー（右辺ベクトル、実数型）
  subroutine monolis_copy_mat_value_rhs_R(mat_in, mat_out)
    implicit none
    !> [in] monolis 構造体（入力）
    type(monolis_structure), intent(in) :: mat_in
    !> [in,out] monolis 構造体（出力）
    type(monolis_structure), intent(inout) :: mat_out
    mat_out%MAT%R%B = mat_in%MAT%R%B
  end subroutine monolis_copy_mat_value_rhs_R

  !> @ingroup matrix_copy
  !> 行列値のコピー（解ベクトル、実数型）
  subroutine monolis_copy_mat_value_solution_R(mat_in, mat_out)
    implicit none
    !> [in] monolis 構造体（入力）
    type(monolis_structure), intent(in) :: mat_in
    !> [in,out] monolis 構造体（出力）
    type(monolis_structure), intent(inout) :: mat_out
    mat_out%MAT%R%X = mat_in%MAT%R%X
  end subroutine monolis_copy_mat_value_solution_R

  !> @ingroup matrix_copy
  !> 行列値のコピー（複素数型）
  subroutine monolis_copy_mat_value_C(mat_in, mat_out)
    implicit none
    !> [in] monolis 構造体（入力）
    type(monolis_structure), intent(in) :: mat_in
    !> [in,out] monolis 構造体（出力）
    type(monolis_structure), intent(inout) :: mat_out
    call monolis_copy_mat_value_matrix_C(mat_in, mat_out)
    call monolis_copy_mat_value_rhs_C(mat_in, mat_out)
    call monolis_copy_mat_value_solution_C(mat_in, mat_out)
  end subroutine monolis_copy_mat_value_C

  !> @ingroup matrix_copy
  !> 行列値のコピー（行列、複素数型）
  subroutine monolis_copy_mat_value_matrix_C(mat_in, mat_out)
    implicit none
    !> [in] monolis 構造体（入力）
    type(monolis_structure), intent(in) :: mat_in
    !> [in,out] monolis 構造体（出力）
    type(monolis_structure), intent(inout) :: mat_out
    mat_out%MAT%C%A = mat_in%MAT%C%A
  end subroutine monolis_copy_mat_value_matrix_C

  !> @ingroup matrix_copy
  !> 行列値のコピー（右辺ベクトル、複素数型）
  subroutine monolis_copy_mat_value_rhs_C(mat_in, mat_out)
    implicit none
    !> [in] monolis 構造体（入力）
    type(monolis_structure), intent(in) :: mat_in
    !> [in,out] monolis 構造体（出力）
    type(monolis_structure), intent(inout) :: mat_out
    mat_out%MAT%C%B = mat_in%MAT%C%B
  end subroutine monolis_copy_mat_value_rhs_C

  !> @ingroup matrix_copy
  !> 行列値のコピー（解ベクトル、複素数型）
  subroutine monolis_copy_mat_value_solution_C(mat_in, mat_out)
    implicit none
    !> [in] monolis 構造体（入力）
    type(monolis_structure), intent(in) :: mat_in
    !> [in,out] monolis 構造体（出力）
    type(monolis_structure), intent(inout) :: mat_out
    mat_out%MAT%C%X = mat_in%MAT%C%X
  end subroutine monolis_copy_mat_value_solution_C

  !> @ingroup matrix_copy
  !> 行列値のゼロ初期化（実数型）
  subroutine monolis_clear_mat_value_R(mat)
    implicit none
    !> [in,out] monolis 構造体（入力）
    type(monolis_structure), intent(inout) :: mat
    call monolis_clear_mat_value_matrix_R(mat)
    call monolis_clear_mat_value_rhs_R(mat)
    call monolis_clear_mat_value_solution_R(mat)
  end subroutine monolis_clear_mat_value_R

  !> @ingroup matrix_copy
  !> 行列値のゼロ初期化（行列、実数型）
  subroutine monolis_clear_mat_value_matrix_R(mat)
    implicit none
    !> [in,out] monolis 構造体（入力）
    type(monolis_structure), intent(inout) :: mat
    mat%MAT%R%A = 0.0d0
  end subroutine monolis_clear_mat_value_matrix_R

  !> @ingroup matrix_copy
  !> 行列値のゼロ初期化（右辺ベクトル、実数型）
  subroutine monolis_clear_mat_value_rhs_R(mat)
    implicit none
    !> [in,out] monolis 構造体（入力）
    type(monolis_structure), intent(inout) :: mat
    mat%MAT%R%B = 0.0d0
  end subroutine monolis_clear_mat_value_rhs_R

  !> @ingroup matrix_copy
  !> 行列値のゼロ初期化（解ベクトル、実数型）
  subroutine monolis_clear_mat_value_solution_R(mat)
    implicit none
    !> [in,out] monolis 構造体（入力）
    type(monolis_structure), intent(inout) :: mat
    mat%MAT%R%X = 0.0d0
  end subroutine monolis_clear_mat_value_solution_R

  !> @ingroup matrix_copy
  !> 行列値のゼロ初期化（複素数型）
  subroutine monolis_clear_mat_value_C(mat)
    implicit none
    !> [in,out] monolis 構造体（入力）
    type(monolis_structure), intent(inout) :: mat
    call monolis_clear_mat_value_matrix_C(mat)
    call monolis_clear_mat_value_rhs_C(mat)
    call monolis_clear_mat_value_solution_C(mat)
  end subroutine monolis_clear_mat_value_C

  !> @ingroup matrix_copy
  !> 行列値のゼロ初期化（行列、複素数型）
  subroutine monolis_clear_mat_value_matrix_C(mat)
    implicit none
    !> [in,out] monolis 構造体（入力）
    type(monolis_structure), intent(inout) :: mat
    mat%MAT%C%A = 0.0d0
  end subroutine monolis_clear_mat_value_matrix_C

  !> @ingroup matrix_copy
  !> 行列値のゼロ初期化（右辺ベクトル、複素数型）
  subroutine monolis_clear_mat_value_rhs_C(mat)
    implicit none
    !> [in,out] monolis 構造体（入力）
    type(monolis_structure), intent(inout) :: mat
    mat%MAT%C%B = 0.0d0
  end subroutine monolis_clear_mat_value_rhs_C

  !> @ingroup matrix_copy
  !> 行列値のゼロ初期化（解ベクトル、複素数型）
  subroutine monolis_clear_mat_value_solution_C(mat)
    implicit none
    !> [in,out] monolis 構造体（入力）
    type(monolis_structure), intent(inout) :: mat
    mat%MAT%C%X = 0.0d0
  end subroutine monolis_clear_mat_value_solution_C
end module mod_monolis_spmat_copy