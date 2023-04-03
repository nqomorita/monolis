/* monolis_spmat_copy.h */
#ifndef MONOLIS_SPMAT_COPY_H
#define MONOLIS_SPMAT_COPY_H

#ifdef __cplusplus
extern "C" {
#endif

#include <complex.h>

void monolis_copy_mat_value_R(
  MONOLIS* mat_in,
  MONOLIS* mat_out);

void monolis_copy_mat_value_A_R(
  MONOLIS* mat_in,
  MONOLIS* mat_out);

void monolis_copy_mat_value_X_R(
  MONOLIS* mat_in,
  MONOLIS* mat_out);

void monolis_copy_mat_value_B_R(
  MONOLIS* mat_in,
  MONOLIS* mat_out);

void monolis_copy_mat_value_C(
  MONOLIS* mat_in,
  MONOLIS* mat_out);

void monolis_copy_mat_value_A_C(
  MONOLIS* mat_in,
  MONOLIS* mat_out);

void monolis_copy_mat_value_X_C(
  MONOLIS* mat_in,
  MONOLIS* mat_out);

void monolis_copy_mat_value_B_C(
  MONOLIS* mat_in,
  MONOLIS* mat_out);

void monolis_clear_mat_value_R(
  MONOLIS* mat);

void monolis_clear_mat_value_A_R(
  MONOLIS* mat);

void monolis_clear_mat_value_X_R(
  MONOLIS* mat);

void monolis_clear_mat_value_B_R(
  MONOLIS* mat);

void monolis_clear_mat_value_C(
  MONOLIS* mat);

void monolis_clear_mat_value_A_C(
  MONOLIS* mat);

void monolis_clear_mat_value_X_C(
  MONOLIS* mat);

void monolis_clear_mat_value_B_C(
  MONOLIS* mat);

#ifdef __cplusplus
}
#endif

#endif
