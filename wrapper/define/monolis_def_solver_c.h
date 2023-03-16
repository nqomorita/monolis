/* monolis_def_struc.h */
#ifndef MONOLIS_DEF_SOLVER_H
#define MONOLIS_DEF_SOLVER_H

#ifdef __cplusplus
extern "C" {
#endif

const int monolis_prm_Iarray_size = 100;

const int monolis_prm_Rarray_size = 100;

const int monolis_iter_CG       = 1;

const int monolis_iter_GropCG   = 2;

const int monolis_iter_PipeCG   = 3;

const int monolis_iter_PipeCR   = 4;

const int monolis_iter_BiCGSTAB = 5;

const int monolis_iter_PipeBiCGSTAB = 6;

const int monolis_iter_BiCGSTAB_noprec = 7;

//const int monolis_iter_CABiCGSTAB_noprec = 8;

const int monolis_iter_PipeBiCGSTAB_noprec = 9;

//const int monolis_iter_GMRES = 10;

const int monolis_iter_COCG = 11;

const int monolis_prec_NONE   = 0;

const int monolis_prec_DIAG   = 1;

const int monolis_prec_ILU    = 2;

const int monolis_prec_JACOBI = 3;

const int monolis_prec_SOR    = 4;

const int monolis_prec_SPIKE  = 5;

const int monolis_prec_LU     = 6;

const int monolis_prec_MUMPS  = 7;

const int monolis_prec_ROM    = 8;

const int monolis_prec_MF     = 9;

const int monolis_prec_MUMPS_LOCAL = 10;

const int monolis_prm_I_method = 1;

const int monolis_prm_I_precond = 2;

const int monolis_prm_I_max_iter = 3;

const int monolis_prm_I_cur_iter = 4;

const int monolis_prm_I_ierr = 5;

//const int monolis_prm_I_is_scaling = 6;

//const int monolis_prm_I_is_reordering = 7;

const int monolis_prm_I_is_init_x = 8;

const int monolis_prm_I_is_sym_matrix = 9;

const int monolis_prm_I_is_debug = 10;

const int monolis_prm_I_is_measurement = 11;

const int monolis_prm_I_is_check_diag = 12;

const int monolis_prm_I_is_prec_stored = 13;

const int monolis_prm_I_show_iterlog = 14;

const int monolis_prm_I_show_time = 15;

const int monolis_prm_I_show_summary = 16;

const int monolis_prm_I_show_time_statistics = 17;

const int monolis_prm_R_tol = 1;

const int monolis_prm_R_cur_resid = 2;

const int monolis_R_time_sol = 3;

const int monolis_R_time_prep = 4;

const int monolis_R_time_spmv = 5;

const int monolis_R_time_dotp = 6;

const int monolis_R_time_prec = 7;

const int monolis_R_time_comm_dotp = 8;

const int monolis_R_time_comm_spmv = 9;

typedef struct {

  int    Iarray[100];

  double Rarray[100];
} MONOLIS_PRM;

void monolis_prm_initialize(
  MONOLIS_PRM* prm);

void monolis_prm_finalize(
  MONOLIS_PRM* prm);

#ifdef __cplusplus
}
#endif

#endif
