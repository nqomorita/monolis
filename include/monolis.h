/* monolis_c.h */
#ifndef MONOLIS_H
#define MONOLIS_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  int myrank, comm, commsize, n_neib;
  int *neib_pe;
  int *recv_index;
  int *recv_item;
  int *send_index;
  int *send_item;
} monolis_com;

typedef struct {
  int N, NP, NPU, NPL, NDOF;
  int *indexU;
  int *indexL;
  int *itemU;
  int *itemL;
  double *D;
  double *AU;
  double *AL;
  double *X;
  double *B;
} monolis_mat;

typedef struct {
  int method;
  int precond;
  int maxiter;
  double tol;
  //bool is_scaling;
  //bool is_reordering;
} monolis_prm;

extern void monolis_com_initialize(monolis_com *monoCOM);

extern void monolis(monolis_com *monoCOM, int N, int NP, int NZ, int NDOF,
  double *A, double *X, double *B, int *index, int *item,
  int method, int precond, int maxiter, double tol,
  int is_scaling, int is_reordering, int is_init_x, int show_iterlog, int show_time, int show_summary);

extern void monolis_com_finalize(monolis_com *monoCOM);

extern void monolis_matvec_wrapper(monolis_com *monoCOM, int N, int NP, int NZ, int NDOF, double *A,
  int *index, int *item, double *X, double *B);

#ifdef __cplusplus
}
#endif

#endif
