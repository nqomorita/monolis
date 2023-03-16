#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void monolis_prm_initialize(
  MONOLIS_PRM* prm)
{
  int i;

  for (int i = 0; i < monolis_prm_Iarray_size; ++i) {
    prm->Iarray[i] = 0;
  }

  for (int i = 0; i < monolis_prm_Rarray_size; ++i) {
    prm->Rarray[i] = 0.0;
  }

  prm->Iarray[monolis_prm_I_method] = 1;
  prm->Iarray[monolis_prm_I_precond] = 1;
  prm->Iarray[monolis_prm_I_max_iter] = 1000;
  prm->Iarray[monolis_prm_I_cur_iter] = 0;
  prm->Iarray[monolis_prm_I_ierr] = -1;
  prm->Iarray[monolis_prm_I_is_init_x] = monolis_I_true;
  prm->Iarray[monolis_prm_I_is_sym_matrix] = monolis_I_false;
  prm->Iarray[monolis_prm_I_is_debug] = monolis_I_false;
  prm->Iarray[monolis_prm_I_is_measurement] = monolis_I_false;
  prm->Iarray[monolis_prm_I_is_check_diag] = monolis_I_false;
  prm->Iarray[monolis_prm_I_is_prec_stored] = monolis_I_false;
  prm->Iarray[monolis_prm_I_show_iterlog] = monolis_I_true;
  prm->Iarray[monolis_prm_I_show_time] = monolis_I_true;
  prm->Iarray[monolis_prm_I_show_summary] = monolis_I_true;
  prm->Iarray[monolis_prm_I_show_time_statistics] = monolis_I_false;

  prm->Rarray[monolis_prm_R_tol] = 1.0e-8;
}

void monolis_prm_finalize(
  MONOLIS_PRM* prm)
{
  int i;

  for (int i = 0; i < monolis_prm_Iarray_size; ++i) {
    prm->Iarray[i] = 0;
  }

  for (int i = 0; i < monolis_prm_Rarray_size; ++i) {
    prm->Rarray[i] = 0.0;
  }
}
