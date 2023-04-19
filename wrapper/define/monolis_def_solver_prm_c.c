#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "monolis_def_solver_prm_c.h"
#include "monolis_utils.h"

void monolis_prm_initialize(
  MONOLIS_PRM* prm)
{
  int i;

  strcpy(prm->com_top_dir_name, ".");
  strcpy(prm->com_part_dir_name, "parted.0");
  strcpy(prm->com_file_name, "node.dat");

  for (i = 0; i < MONOLIS_PRM_IARRAY_SIZE; ++i) {
    prm->Iarray[i] = 0;
  }

  for (i = 0; i < MONOLIS_PRM_RARRAY_SIZE; ++i) {
    prm->Rarray[i] = 0.0;
  }

  prm->Iarray[MONOLIS_PRM_I_METHOD] = 1;
  prm->Iarray[MONOLIS_PRM_I_PRECOND] = 1;
  prm->Iarray[MONOLIS_PRM_I_MAX_ITER] = 1000;
  prm->Iarray[MONOLIS_PRM_I_CUR_ITER] = 0;
  prm->Iarray[MONOLIS_PRM_I_IERR] = -1;
  prm->Iarray[MONOLIS_PRM_I_IS_INIT_X] = MONOLIS_I_TRUE;
  prm->Iarray[MONOLIS_PRM_I_IS_SYM_MATRIX] = MONOLIS_I_FALSE;
  prm->Iarray[MONOLIS_PRM_I_IS_DEBUG] = MONOLIS_I_FALSE;
  prm->Iarray[MONOLIS_PRM_I_IS_MEASUREMENT] = MONOLIS_I_FALSE;
  prm->Iarray[MONOLIS_PRM_I_IS_CHECK_DIAG] = MONOLIS_I_FALSE;
  prm->Iarray[MONOLIS_PRM_I_IS_PREC_STORED] = MONOLIS_I_FALSE;
  prm->Iarray[MONOLIS_PRM_I_SHOW_ITERLOG] = MONOLIS_I_TRUE;
  prm->Iarray[MONOLIS_PRM_I_SHOW_TIME] = MONOLIS_I_TRUE;
  prm->Iarray[MONOLIS_PRM_I_SHOW_SUMMARY] = MONOLIS_I_TRUE;
  prm->Iarray[MONOLIS_PRM_I_SHOW_TIME_STATISTICS] = MONOLIS_I_FALSE;

  prm->Rarray[MONOLIS_PRM_R_TOL] = 1.0e-8;
}

void monolis_prm_finalize(
  MONOLIS_PRM* prm)
{
  int i;

  for (i = 0; i < MONOLIS_PRM_IARRAY_SIZE; ++i) {
    prm->Iarray[i] = 0;
  }

  for (i = 0; i < MONOLIS_PRM_RARRAY_SIZE; ++i) {
    prm->Rarray[i] = 0.0;
  }
}
