/* monolis_struct.h */
#ifndef MONOLIS_STRUCT_H
#define MONOLIS_STRUCT_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>

const int monolis_success = 0;
const int monolis_fail    = 1;

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

const int monolis_sum = 1;
const int monolis_max = 2;
const int monolis_min = 3;

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
  int intrenal_nnode;

  int recv_n_neib;
  int* recv_neib_pe;
  int* recv_index;
  int* recv_item;

  int send_n_neib;
  int* send_neib_pe;
  int* send_index;
  int* send_item;

  int* global_node_id;
  int* global_elem_id;
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

#ifdef __cplusplus
}
#endif

#endif
