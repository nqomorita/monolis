/* monolis_def_struc.h */
#ifndef MONOLIS_DEF_STRUC_H
#define MONOLIS_DEF_STRUC_H

#ifdef __cplusplus
extern "C" {
#endif

#include "monolis_utils.h"
#include "monolis_def_solver_prm_c.h"
#include "monolis_def_mat_c.h"

typedef struct {
  MONOLIS_PRM prm;
  MONOLIS_COM com;
  MONOLIS_MAT mat;
  MONOLIS_MAT prec;
} MONOLIS;

void monolis_global_initialize();

void monolis_global_finalize();

void monolis_initialize(
  MONOLIS*    mat,
  const char* input_file_dir);

void monolis_finalize(
  MONOLIS* mat);

#ifdef __cplusplus
}
#endif

#endif
