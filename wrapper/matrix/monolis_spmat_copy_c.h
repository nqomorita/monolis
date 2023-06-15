/* monolis_spmat_copy.h */
#ifndef MONOLIS_SPMAT_COPY_H
#define MONOLIS_SPMAT_COPY_H

#ifdef __cplusplus
extern "C" {
#endif

void monolis_copy_mat_R(
  MONOLIS* mat_in,
  MONOLIS* mat_out);

void monolis_copy_mat_C(
  MONOLIS* mat_in,
  MONOLIS* mat_out);

void monolis_copy_mat_nonzero_pattern_R(
  MONOLIS* mat_in,
  MONOLIS* mat_out);

void monolis_copy_mat_nonzero_pattern_C(
  MONOLIS* mat_in,
  MONOLIS* mat_out);

void monolis_copy_mat_nonzero_pattern_val_R(
  int      NP,
  int      NDOF,
  int      NZ,
  int      NZU,
  int      NZL,
  MONOLIS_MAT_VAL_R* mat_in,
  MONOLIS_MAT_VAL_R* mat_out);

void monolis_copy_mat_nonzero_pattern_val_C(
  int      NP,
  int      NDOF,
  int      NZ,
  int      NZU,
  int      NZL,
  MONOLIS_MAT_VAL_C* mat_in,
  MONOLIS_MAT_VAL_C* mat_out);

void monolis_copy_mat_nonzero_pattern_CSR(
  int      NP,
  MONOLIS_MAT_CSR* mat_in,
  MONOLIS_MAT_CSR* mat_out);

void monolis_copy_mat_nonzero_pattern_CSC(
  int      NP,
  MONOLIS_MAT_CSC* mat_in,
  MONOLIS_MAT_CSC* mat_out);

void monolis_copy_mat_nonzero_pattern_SCSR(
  int      NP,
  MONOLIS_MAT_SEPARATED_CSR* mat_in,
  MONOLIS_MAT_SEPARATED_CSR* mat_out);

void monolis_copy_mat_value_R(
  MONOLIS* mat_in,
  MONOLIS* mat_out);

void monolis_copy_mat_value_matrix_R(
  MONOLIS* mat_in,
  MONOLIS* mat_out);

void monolis_copy_mat_value_solution_R(
  MONOLIS* mat_in,
  MONOLIS* mat_out);

void monolis_copy_mat_value_rhs_R(
  MONOLIS* mat_in,
  MONOLIS* mat_out);

void monolis_copy_mat_value_C(
  MONOLIS* mat_in,
  MONOLIS* mat_out);

void monolis_copy_mat_value_matrix_C(
  MONOLIS* mat_in,
  MONOLIS* mat_out);

void monolis_copy_mat_value_solution_C(
  MONOLIS* mat_in,
  MONOLIS* mat_out);

void monolis_copy_mat_value_rhs_C(
  MONOLIS* mat_in,
  MONOLIS* mat_out);

void monolis_clear_mat_value_R(
  MONOLIS* mat);

void monolis_clear_mat_value_matrix_R(
  MONOLIS* mat);

void monolis_clear_mat_value_solution_R(
  MONOLIS* mat);

void monolis_clear_mat_value_rhs_R(
  MONOLIS* mat);

void monolis_clear_mat_value_C(
  MONOLIS* mat);

void monolis_clear_mat_value_matrix_C(
  MONOLIS* mat);

void monolis_clear_mat_value_solution_C(
  MONOLIS* mat);

void monolis_clear_mat_value_rhs_C(
  MONOLIS* mat);

#ifdef __cplusplus
}
#endif

#endif
