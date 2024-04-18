#include "monolis_wrapper_ml_util.h"

int monolis_ML_getrow_nn(ML_Operator *mat_in, int N_requested_rows,
                       int requested_rows[], int allocated_space,
                       int cols[], double values[], int row_lengths[]) {
  int *id, ierr;
  id = (int *)ML_Get_MyGetrowData(mat_in);
  monolis_ml_getrow_nn_(id, &N_requested_rows, requested_rows, &allocated_space,
                      cols, values, row_lengths, &ierr);
  return ierr;
}

int monolis_ML_matvec_nn(ML_Operator *mat_in, int in_length, double p[],
                       int out_length, double ap[]) {
  int *id, ierr;
  id = (int *)ML_Get_MyGetrowData(mat_in);
  monolis_ml_matvec_nn_(id, &in_length, p, &out_length, ap, &ierr);
  return ierr;
}

int monolis_ML_comm_nn(double x[], void *A_data) {
  int *id, ierr;
  id = (int *)A_data;
  monolis_ml_comm_nn_(id, x, &ierr);
  return ierr;
}
