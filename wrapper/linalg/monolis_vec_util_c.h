/* monolis_vec_util_c.h */
#ifndef MONOLIS_VEC_UTIL_C_H
#define MONOLIS_VEC_UTIL_C_H

#ifdef __cplusplus
extern "C" {
#endif

#include <complex.h>
#include "monolis_def_struc_c.h"

void monolis_update_R(
  MONOLIS* mat,
  int      n_dof,
  double*  x);

void monolis_update_I(
  MONOLIS* mat,
  int      n_dof,
  int*     x);

void monolis_update_C(
  MONOLIS*        mat,
  int             n_dof,
  double complex* x);

#ifdef __cplusplus
}
#endif

#endif
