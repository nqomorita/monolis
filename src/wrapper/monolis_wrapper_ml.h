#ifndef MONOLIS_WRAPPER_ML_H
#define MONOLIS_WRAPPER_ML_H

extern void monolis_ml_get_nlocal_(
    int *nlocal, 
    int *nlocal_allcolumns, 
    int *ierr);

extern void monolis_ml_getrow_nn_(
    int *n_requested_rows,
    int *requested_rows, 
    int *allocated_space,
    int *cols, 
    double *values, 
    int *row_lengths,
    int *ierr);

extern void monolis_ml_matvec_nn_(
    int *in_length, 
    double *p,
    int *out_length, 
    double *ap, 
    int *ierr);

extern void monolis_ml_comm_nn_(
    double *x, 
    int *ierr);

#ifdef WITH_ML
extern int monolis_ML_getrow_nn(
    ML_Operator *mat_in, 
    int N_requested_rows,
    int requested_rows[], 
    int allocated_space,
    int cols[], 
    double values[], 
    int row_lengths[]);

extern int monolis_ML_matvec_nn(
    ML_Operator *mat_in, 
    int in_length, 
    double p[],
    int out_length, 
    double ap[]);
#endif

extern int monolis_ML_comm_nn(
    double x[], 
    void *A_data);

#endif
