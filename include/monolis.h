/* monolis_c.h */
#ifndef MONOLIS_H
#define MONOLIS_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>

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
  bool show_time;
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
  double* diag;
} MONOLIS_MAT;

typedef struct {
  MONOLIS_PRM prm;
  MONOLIS_COM com;
  MONOLIS_MAT mat;
} MONOLIS;

void monolis_global_initialize();

void monolis_global_finalize();

void monolis_initialize(
  MONOLIS* mat);

void monolis_finalize(
  MONOLIS* mat);

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

void monolis_set_Dirichlet_bc(
  MONOLIS* mat,
  double*  b,
  int      node_id,
  int      ndof_bc,
  double   val);

void monolis_solve(
  MONOLIS* mat,
  double*  b,
  double*  x);

#ifdef __cplusplus
}
#endif

#endif
