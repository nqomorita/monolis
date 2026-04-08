#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include "monolis_utils.h"

#ifdef WITH_ML
#include "Trilinos_version.h"
#include "ml_include.h"
#include "ml_config.h"
#endif

#include "monolis_wrapper_ml.h"

#ifdef WITH_ML
int monolis_ML_getrow(
  ML_Operator *mat_in, 
  int    N_requested_rows,
  int    requested_rows[], 
  int    allocated_space,
  int    cols[], 
  double values[], 
  int    row_lengths[])
{
  int ierr;
  monolis_ml_getrow_nn_c(&N_requested_rows, requested_rows, &allocated_space,
                        cols, values, row_lengths, &ierr);
  return ierr;
}

int monolis_ML_matvec(
  ML_Operator *mat_in, 
  int    in_length, 
  double p[],
  int    out_length, 
  double ap[])
{
  int ierr;
  monolis_ml_matvec_nn_c(&in_length, p, &out_length, ap, &ierr);
  return ierr;
}

int monolis_ML_comm(
  double x[], 
  void   *A_data)
{
  int ierr;
  monolis_ml_comm_nn_c(x, &ierr);
  return ierr;
}

struct ml_info {
  ML*           ml_object;
  ML_Aggregate* agg_object;
};

static struct ml_info MLInfo;

void monolis_ML_compute_rigid_body_modes(
  int    ndof,
  int    nnode,
  double coord[],
  int    null_dim,
  int    leng,
  double null_vect[])
{
  int i, base;
  double x, y, z;

  for (i = 0; i < leng * null_dim; i++) {
    null_vect[i] = 0.0;
  }

  if (ndof == 1) {
    for (i = 0; i < nnode; i++) {
      null_vect[i] = 1.0;
    }
  } else if (ndof == 2) {
    for (i = 0; i < nnode; i++) {
      base = i * 2;
      x = coord[i * 2 + 0];
      y = coord[i * 2 + 1];
      /* translation x */
      null_vect[0 * leng + base + 0] = 1.0;
      /* translation y */
      null_vect[1 * leng + base + 1] = 1.0;
      /* rotation z */
      null_vect[2 * leng + base + 0] = -y;
      null_vect[2 * leng + base + 1] =  x;
    }
  } else if (ndof == 3) {
    for (i = 0; i < nnode; i++) {
      base = i * 3;
      x = coord[i * 3 + 0];
      y = coord[i * 3 + 1];
      z = coord[i * 3 + 2];
      /* translation x */
      null_vect[0 * leng + base + 0] = 1.0;
      /* translation y */
      null_vect[1 * leng + base + 1] = 1.0;
      /* translation z */
      null_vect[2 * leng + base + 2] = 1.0;
      /* rotation around x */
      null_vect[3 * leng + base + 1] =  z;
      null_vect[3 * leng + base + 2] = -y;
      /* rotation around y */
      null_vect[4 * leng + base + 0] = -z;
      null_vect[4 * leng + base + 2] =  x;
      /* rotation around z */
      null_vect[5 * leng + base + 0] =  y;
      null_vect[5 * leng + base + 1] = -x;
    }
  }
}
#endif

void monolis_ML_wrapper_setup(
  int*    sym, 
  int*    ndof, 
  int*    nnode,
  double  coord[],
  int*    ierr)
{
  int loglevel, myrank, nglobal;
  int N_grids, N_levels;
  int nlocal, nlocal_allcolumns;
  int null_dim;
  double *null_vect;

#ifdef WITH_ML
  ML *ml_object;
  ML_Aggregate *agg_object;

  loglevel = 1;
  ML_Set_PrintLevel(loglevel);

  /** AMG definition **/
  N_grids = 100;
  ML_Create(&ml_object, N_grids);

  /** matrix definition **/
  myrank = monolis_mpi_get_global_my_rank();
  monolis_ml_get_nlocal_c(&nlocal, &nlocal_allcolumns, ierr);

  ML_Init_Amatrix(ml_object, 0, nlocal, nlocal, &myrank);
  ML_Set_Amatrix_Getrow(ml_object, 0, monolis_ML_getrow, monolis_ML_comm, nlocal_allcolumns);
  ML_Set_Amatrix_Matvec(ml_object, 0, monolis_ML_matvec);

  /** aggregate definition **/
  ML_Aggregate_Create(&agg_object);

  nglobal = nlocal;
  monolis_allreduce_I(1, &nglobal, MONOLIS_MPI_SUM, monolis_mpi_get_global_comm());
  
  //ML_Aggregate_Set_MaxCoarseSize(agg_object, 1);
  ML_Aggregate_Set_MaxCoarseSize(agg_object, nglobal/10);
  
  ML_Aggregate_Set_CoarsenScheme_Uncoupled(agg_object);
  //ML_Aggregate_Set_CoarsenScheme_UncoupledMIS(agg_object);
  //ML_Aggregate_Set_CoarsenScheme_Zoltan(agg_object);

  /** null space definition (rigid body modes) **/
  if (*ndof == 1) {
    null_dim = 1;
  } else if (*ndof == 2) {
    null_dim = 3;
  } else {
    null_dim = 6;
  }

  ML_Aggregate_Set_Dimensions(agg_object, *ndof);

  null_vect = (double*)calloc(nlocal * null_dim, sizeof(double));
  monolis_ML_compute_rigid_body_modes(*ndof, *nnode, coord, null_dim, nlocal, null_vect);
  ML_Aggregate_Set_NullSpace(agg_object, *ndof, null_dim, null_vect, nlocal);
  free(null_vect);

  //N_levels = ML_Gen_MGHierarchy_UsingAggregation(ml_object, 0, ML_INCREASING, agg_object); 
  N_levels = ML_Gen_MultiLevelHierarchy_UsingAggregation(ml_object, 0, ML_INCREASING, agg_object); 
  if(myrank == 0){ printf("N_levels %d \n", N_levels); }

  /** smoother definition **/
  //ML_Gen_Smoother_Jacobi(ml_object, ML_ALL_LEVELS, ML_BOTH, 1, ML_DEFAULT);
  //ML_Gen_Smoother_GaussSeidel(ml_object, ML_ALL_LEVELS, ML_BOTH, 2, ML_DEFAULT);
  ML_Gen_Smoother_Cheby(ml_object, ML_ALL_LEVELS, ML_BOTH, 20.0, 4);

  /** solver definition **/
  ML_Gen_Solver(ml_object, ML_MGV, 0, N_levels - 1);
  //ML_Gen_Solver(ml_object, ML_MGW, 0, N_levels - 1);
  //ML_Gen_Solver(ml_object, ML_MGFULLV, 0, N_levels - 1);

  /* Save objects */
  MLInfo.ml_object  = ml_object;
  MLInfo.agg_object = agg_object;
#else
    printf("%s\n", "* monolis_wrapper_ml.c: ML is NOT enabled");
    exit(1);
#endif
}

void monolis_ML_wrapper_apply(
  double rhs[], 
  int*   ierr)
{
  int nlocal, nlocal_allcolumns;
  double *sol;
  int i;
#ifdef WITH_ML  
  ML *ml_object;

  ml_object = MLInfo.ml_object;

  monolis_ml_get_nlocal_c(&nlocal, &nlocal_allcolumns, ierr);

  sol = monolis_alloc_R_1d(sol, nlocal_allcolumns);

  ML_Solve_MGV(ml_object, rhs, sol);

  for (i = 0; i < nlocal; i++) {
    rhs[i] = sol[i];
  }

  monolis_dealloc_R_1d(&sol);
#endif
}

void monolis_ML_wrapper_clear(
  int* ierr)
{
#ifdef WITH_ML  
  ML_Aggregate_Destroy(&(MLInfo.agg_object));
  ML_Destroy(&(MLInfo.ml_object));
#endif
}

