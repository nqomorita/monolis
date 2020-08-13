/* monolis_c.h */
#ifndef MONOLIS_H
#define MONOLIS_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>

const int monolis_iter_CG       = 1;
const int monolis_iter_GropCG   = 2;
const int monolis_iter_PipeCG   = 3;
const int monolis_iter_PipeCR   = 4;
const int monolis_iter_BiCGSTAB = 5;
const int monolis_iter_PipeBiCGSTAB = 6;
const int monolis_iter_BiCGSTAB_noprec = 7;
const int monolis_iter_CABiCGSTAB_noprec = 8;
const int monolis_iter_PipeBiCGSTAB_noprec = 9;
const int monolis_iter_SOR      = 10;
const int monolis_iter_IR       = 11;

const int monolis_prec_NONE   = 0;
const int monolis_prec_DIAG   = 1;
const int monolis_prec_ILU    = 2;
const int monolis_prec_JACOBI = 3;
const int monolis_prec_SOR    = 4;
const int monolis_prec_SAINV  = 5;
const int monolis_prec_RIF    = 6;
const int monolis_prec_SPIKE  = 7;
const int monolis_prec_DIRECT = 8;
const int monolis_prec_MUMPS  = 9;

typedef struct {
  int method;
  int precond;
  int maxiter;
  int curiter;
  int ierr;
  double tol;
  double curresid;

  /* flags */
  bool is_scaling;
  bool is_reordering;
  bool is_init_x;
  bool is_sym_matrix;
  bool is_debug;
  bool is_check_diag;
  bool show_iterlog;
  bool show_timelog;
  bool show_summary;

  /* time */
  double tsol;
  double tprep;
  double tspmv;
  double tdotp;
  double tprec;
  double tcomm;
} MONOLIS_PRM;

typedef struct {
  int myrank;
  int comm;
  int commsize;

  int recv_n_neib;
  int* recv_neib_pe;
  int* recv_index;
  int* recv_item;

  int send_n_neib;
  int* send_neib_pe;
  int* send_index;
  int* send_item;
} MONOLIS_COM;

typedef struct {
  int N, NP, NZ, NDOF;
  int* index;
  int* item;

  /* for CRR format */
  int* indexR;
  int* itemR;
  int* permR;

  double* A;
  double* X;
  double* B;
  //double* diag;
} MONOLIS_MAT;

typedef struct {
  MONOLIS_PRM prm;
  MONOLIS_COM com;
  MONOLIS_MAT mat;
} MONOLIS;

/* getter */
int monolis_get_global_comm();

int monolis_get_global_commsize();

int monolis_get_global_myrank();

double monolis_get_time();

void monolis_get_input_filename();

/* setter */
void monolis_set_method   (MONOLIS* mat, int    flag);

void monolis_set_precond  (MONOLIS* mat, int    flag);

void monolis_set_maxiter  (MONOLIS* mat, int    flag);

void monolis_set_tolerance(MONOLIS* mat, double flag);

void monolis_show_iterlog (MONOLIS* mat, bool   flag);

void monolis_show_timelog (MONOLIS* mat, bool   flag);

void monolis_show_summary (MONOLIS* mat, bool   flag);

/* mat initializer */
void monolis_global_initialize();

void monolis_global_finalize();

void monolis_initialize(
  MONOLIS* mat);

void monolis_finalize(
  MONOLIS* mat);

/* mat clear */
void monolis_clear(
  MONOLIS* mat);

void monolis_clear_mat(
  MONOLIS* mat);

void monolis_clear_rhs(
  MONOLIS* mat);

void monolis_clear_solution(
  MONOLIS* mat);

/* mat copy */
void monolis_copy_all(
  MONOLIS* in,
  MONOLIS* out);

void monolis_copy_nonzero_pattern(
  MONOLIS* in,
  MONOLIS* out);

void monolis_copy_param(
  MONOLIS* in,
  MONOLIS* out);

void monolis_get_CRR_format(
  int      nnode,
  int      nz,
  int*     index,
  int*     item,
  int*     indexR,
  int*     itemR,
  int*     permR);

void monolis_get_nonzero_pattern(
  MONOLIS* mat,
  int      nnode,
  int      nbase_func,
  int      ndof,
  int      nelem,
  int**    elem);

void monolis_add_sparse_matrix(
  MONOLIS* mat,
  int      nbase_func,
  int      *connectivity,
  double** local_mat);

void monolis_add_scalar_to_sparse_matrix(
  MONOLIS* mat,
  double   val,
  int      i,
  int      j,
  int      submat_i,
  int      submat_j);

void monolis_add_sparse_matrix_c_main(
  int      nnode,
  int      nz,
  int      ndof,
  int      nbf,
  int*     index,
  int*     item,
  double*  A,
  int*     con,
  double*  mat);

void monolis_add_scalar_to_sparse_matrix_c_main(
  int      nnode,
  int      nz,
  int      ndof,
  int*     index,
  int*     item,
  double*  A,
  int      i,
  int      j,
  int      submat_i,
  int      submat_j,
  double   val);

void monolis_set_Dirichlet_bc(
  MONOLIS* mat,
  double*  b,
  int      node_id,
  int      ndof_bc,
  double   val);

void monolis_set_Dirichlet_bc_c_main(
  int      nnode,
  int      nz,
  int      ndof,
  int*     index,
  int*     item,
  int*     indexR,
  int*     itemR,
  int*     permR,
  double*  A,
  double*  b,
  int      node_id,
  int      ndof_bc,
  double   val);

void monolis_solve(
  MONOLIS* mat,
  double*  b,
  double*  x);

void monolis_solve_c_main(
  int      n,
  int      np,
  int      nz,
  int      ndof,
  double*  A,
  double*  x,
  double*  b,
  int*     index,
  int*     item,
  int      myrank,
  int      comm,
  int      commsize,
  int      recv_n_neib,
  int      send_n_neib,
  int      method,
  int      precond,
  int      maxiter,
  double   tol,
  int      iterlog,
  int      timelog,
  int      summary);

void monolis_qsort_int(
  int* array,
  int  iS,
  int  iE);

#ifdef __cplusplus
}
#endif

#endif
