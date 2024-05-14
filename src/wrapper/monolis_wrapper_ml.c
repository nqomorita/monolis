#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include "monolis_utils.h"

#include "Trilinos_version.h"
#include "ml_include.h"
#include "ml_config.h"

#include "monolis_wrapper_ml.h"

int monolis_ML_getrow(ML_Operator *mat_in, int N_requested_rows,
                       int requested_rows[], int allocated_space,
                       int cols[], double values[], int row_lengths[]) {
  int ierr;
  monolis_ml_getrow_nn_(&N_requested_rows, requested_rows, &allocated_space,
                        cols, values, row_lengths, &ierr);
  return ierr;
}

int monolis_ML_matvec(ML_Operator *mat_in, int in_length, double p[],
                       int out_length, double ap[]) {
  int ierr;
  monolis_ml_matvec_nn_(&in_length, p, &out_length, ap, &ierr);
  return ierr;
}

int monolis_ML_comm(double x[], void *A_data) {
  int ierr;
  monolis_ml_comm_nn_(x, &ierr);
  return ierr;
}

struct ml_info {
  ML *ml_object;
  ML_Aggregate *agg_object;
};

static struct ml_info MLInfo;

void monolis_ML_wrapper_setup(int *sym, int *Ndof, int *ierr) {
  int loglevel, myrank, nglobal;
  int N_grids, N_levels;
  int nlocal, nlocal_allcolumns;
  int *id;
  ML *ml_object;
  ML_Aggregate *agg_object;

  N_grids = 100;
  ML_Create(&ml_object, N_grids);

  myrank = monolis_mpi_get_global_my_rank();
  monolis_ml_get_nlocal_(&nlocal, &nlocal_allcolumns, ierr);

  ML_Init_Amatrix(ml_object, 0, nlocal, nlocal, &myrank);
  ML_Set_Amatrix_Getrow(ml_object, 0, monolis_ML_getrow, monolis_ML_comm, nlocal_allcolumns);
  ML_Set_Amatrix_Matvec(ml_object, 0, monolis_ML_matvec);

  ML_Aggregate_Create(&agg_object);

  nglobal = nlocal;
  monolis_allreduce_I(1, &nglobal, MONOLIS_MPI_SUM, monolis_mpi_get_global_comm());
  
  //ML_Aggregate_Set_MaxCoarseSize(agg_object, nglobal - 1);
  //ML_Aggregate_Set_CoarsenScheme_UncoupledMIS(agg_object);
  //ML_Aggregate_Set_Dimensions(agg_object, *Ndof);

  //N_levels = ML_Gen_MGHierarchy_UsingAggregation(ml_object, 0, ML_INCREASING, agg_object); 
  N_levels = ML_Gen_MultiLevelHierarchy_UsingAggregation(ml_object, 0, ML_INCREASING, agg_object); 
  if(myrank == 0){ printf("N_levels %d \n", N_levels); }

  ML_Gen_Smoother_Jacobi(ml_object, ML_ALL_LEVELS, ML_BOTH, 1, ML_DEFAULT);
  //ML_Gen_Smoother_GaussSeidel(ml_object, ML_ALL_LEVELS, ML_BOTH, 1, ML_DEFAULT);
  //ML_Gen_Smoother_Cheby(ml_object, ML_ALL_LEVELS, ML_BOTH, 10.0, 2);

  //ML_Gen_Solver(ml_object, ML_MGV, 0, N_levels - 1);
  //ML_Gen_Solver(ml_object, ML_MGW, 0, N_levels - 1);
  ML_Gen_Solver(ml_object, ML_MGFULLV, 0, N_levels - 1);

  /* Save objects */
  MLInfo.ml_object  = ml_object;
  MLInfo.agg_object = agg_object;
}

void monolis_ML_wrapper_apply(double rhs[], int *ierr) {
  int nlocal, nlocal_allcolumns;
  double *sol;
  int i;
  ML *ml_object;

  ml_object = MLInfo.ml_object;
  monolis_ml_get_nlocal_(&nlocal, &nlocal_allcolumns, ierr);

  sol = monolis_alloc_R_1d(sol, nlocal_allcolumns);

  ML_Solve_MGV(ml_object, rhs, sol);

  for (i = 0; i < nlocal; i++) {
    rhs[i] = sol[i];
  }

  monolis_dealloc_R_1d(&sol);
}

void monolis_ML_wrapper_clear(int *ierr) {
  ML_Aggregate_Destroy(&(MLInfo.agg_object));
  ML_Destroy(&(MLInfo.ml_object));
}

/* Fortran interface */
void monolis_precond_ml_setup_(int *sym, int *ndof, int *ierr) {
  monolis_ML_wrapper_setup(sym, ndof, ierr);
}
void monolis_precond_ml_setup__(int *sym, int *ndof, int *ierr) {
  monolis_ML_wrapper_setup(sym, ndof, ierr);
}
void MONOLIS_PRECOND_ML_SETUP(int *sym, int *ndof, int *ierr) {
  monolis_ML_wrapper_setup(sym, ndof, ierr);
}

void monolis_precond_ml_apply_(double rhs[], int *ierr) {
  monolis_ML_wrapper_apply(rhs, ierr);
}
void monolis_precond_ml_apply__(double rhs[], int *ierr) {
  monolis_ML_wrapper_apply(rhs, ierr);
}
void MONOLIS_PRECOND_ML_APPLY(double rhs[], int *ierr) {
  monolis_ML_wrapper_apply(rhs, ierr);
}

void monolis_precond_ml_clear_(int *ierr) {
  monolis_ML_wrapper_clear(ierr);
}
void monolis_precond_ml_clear__(int *ierr) {
  monolis_ML_wrapper_clear(ierr);
}
void MONOLIS_ML_WRAPPER_CLEAR(int *ierr) {
  monolis_ML_wrapper_clear(ierr);
}
