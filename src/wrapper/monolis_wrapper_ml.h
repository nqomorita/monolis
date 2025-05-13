#ifndef MONOLIS_WRAPPER_ML_H
#define MONOLIS_WRAPPER_ML_H

void monolis_ml_get_nlocal_c(
    int* nlocal, 
    int* nlocal_allcolumns, 
    int* ierr);

void monolis_ml_getrow_nn_c(
    int*    n_requested_rows,
    int*    requested_rows, 
    int*    allocated_space,
    int*    cols, 
    double* values, 
    int*    row_lengths,
    int*    ierr);

void monolis_ml_matvec_nn_c(
    int*    in_length, 
    double* p,
    int*    out_length, 
    double* ap, 
    int*    ierr);

void monolis_ml_comm_nn_c(
    double* x, 
    int*    ierr);

#endif
