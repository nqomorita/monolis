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
  //logical is_scaling
  //logical is_reordering
} monolis_prm;

extern void monolis(int N, int NP, int NDOF, int NPU, int NPL,
  double *D, double *AU, double *AL, double *X, double *B,
  int *indexU, int *itemU, int *indexL, int *itemL,
  int myrank, int comm, int commsize, int n_neib,
  int *neib_pe, int *recv_index, int *send_index, int *recv_item, int *send_item,
  int method, int precond, int maxiter, double tol, int is_scaling);

extern void monolis_serial(int N, int NDOF, int NPU, int NPL,
  double *D, double *AU, double *AL, double *X, double *B,
  int *indexU, int *itemU, int *indexL, int *itemL,
  int method, int precond, int maxiter, double tol, int is_scaling);

extern void monolis_convert_full_matrix(int Nf, int NDOFf, double *Af, double thresh,
  int N, int NDOF, int NPU, int NPL, double *D, double *AU, double *AL,
  int *indexU, int *itemU, int *indexL, int *itemL);

extern void monolis_convert_coo_matrix(int Nf, int NZf, int NDOFf, double *Af, int *indexI, int *indexJ,
  int N, int NDOF, int NPU, int NPL, double *D, double *AU, double *AL,
  int *indexU, int *itemU, int *indexL, int *itemL);

extern void monolis_convert_csr_matrix(int Nf, int NDOFf, double *Af, int *index, int *item,
  int N, int NDOF, int NPU, int NPL, double *D, double *AU, double *AL,
  int *indexU, int *itemU, int *indexL, int *itemL);

extern void monolis_convert_test(int N);

#endif
