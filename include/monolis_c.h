/* monolis_c.h */
#ifndef MONOLIS_H
#define MONOLIS_H

extern void monolis_c_(int N, int NP, int NDOF, int NPU, int NPL,
  double *D, double *AU, double *AL, double *X, double *B,
  int *indexU, int *itemU, int *indexL, int *itemL,
  int myrank, int comm, int commsize, int n_neib,
  int *neib_pe, int *recv_index, int *send_index, int *recv_item, int *send_item,
  int method, int precond, int maxiter, double tol, bool is_scaling);

extern void monolis_serial_c_(int N, int NDOF, int NPU, int NPL,
  double *D, double *AU, double *AL, double *X, double *B,
  int *indexU, int *itemU, int *indexL, int *itemL,
  int method, int precond, int maxiter, double tol, bool is_scaling);

#endif
