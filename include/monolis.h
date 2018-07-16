/* monolis_c.h */
#ifndef MONOLIS_H
#define MONOLIS_H

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

extern void monolis(int N, int NP, int NDOF, int NPU, int NPL,
  double *D, double *AU, double *AL, double *X, double *B,
  int *indexU, int *itemU, int *indexL, int *itemL,
  int myrank, int comm, int commsize, int n_neib,
  int *neib_pe, int *recv_index, int *send_index, int *recv_item, int *send_item,
  int method, int precond, int maxiter, double tol, int is_scaling);

extern void monolis_serial(int *N, int *NDOF, int *NPU, int *NPL,
  double *D, double *AU, double *AL, double *X, double *B,
  int *indexU, int *itemU, int *indexL, int *itemL,
  int *method, int *precond, int *maxiter, double *tol, int *is_scaling);

extern void monolis_convert_full_matrix(int Nf, int NDOFf, double *Af, double thresh,
  int N, int NDOF, int NPU, int NPL, double *D, double *AU, double *AL,
  int *indexU, int *itemU, int *indexL, int *itemL);

extern void monolis_convert_coo_get_size(int *N, int *NZ, int *indexI, int *indexJ, int *NPU, int *NPL);

extern void monolis_convert_coo_get_matrix(int *N, int *NZ, int *NDOF, double *A, int *indexI, int *indexJ,
  int *NPU, int *NPL, double *D, double *U, double *L,
  int *indexU, int *itemU, int *indexL, int *itemL);

extern void monolis_convert_csr_get_size(int *N, int *NZ, int *index, int *item, int *NPU, int *NPL);

extern void monolis_convert_csr_get_matrix(int *N, int *NZ, int *NDOF, double *A, int *index, int *item,
  int *NPU, int *NPL, double *D, double *U, double *L,
  int *indexU, int *itemU, int *indexL, int *itemL);

#endif
