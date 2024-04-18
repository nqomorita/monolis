#ifndef MONOLIS_WRAPPER_ML_UTIL_H
#define MONOLIS_WRAPPER_ML_UTIL_H

#include "ml_include.h"

extern void monolis_ml_getrow_nn_(int *id, int *n_requested_rows,
                                int *requested_rows, int *allocated_space,
                                int *cols, double *values, int *row_lengths,
                                int *ierr);
extern void monolis_ml_matvec_nn_(int *id, int *in_length, double *p,
                                int *out_length, double *ap, int *ierr);
extern void monolis_ml_comm_nn_(int *id, double *x, int *ierr);

extern int monolis_ML_getrow_nn(ML_Operator *mat_in, int N_requested_rows,
                              int requested_rows[], int allocated_space,
                              int cols[], double values[], int row_lengths[]);
extern int monolis_ML_matvec_nn(ML_Operator *mat_in, int in_length, double p[],
                              int out_length, double ap[]);
extern int monolis_ML_comm_nn(double x[], void *A_data);
extern int monolis_ML_smoother_diag_apply_nn(ML_Smoother *data, int x_length, double x[],
                                           int rhs_length, double rhs[]);
extern int monolis_ML_smoother_ssor_apply_nn(ML_Smoother *data, int x_length, double x[],
                                           int rhs_length, double rhs[]);

#endif 
