#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
//#include "monolis_util.h"

#include "Trilinos_version.h"
#include "ml_include.h"
#include "ml_config.h"

#ifdef HAVE_ML_AMESOS
# include "Amesos_config.h"
#endif

#include "monolis_wrapper_ml.h"

/*
 * Options
 */

enum coarse_solver {Smoother, KLU, MUMPS};
enum smoother_type {Cheby};
enum coarsen_scheme {UncoupledMIS, METIS, ParMETIS, Zoltan, DD};

struct ml_options {
  /* Coarse solver
   *  available solvers: Smoother, KLU, MUMPS
   * Note:
   *  - Trilinos must be built with Amesos enabled to use KLU
   *  - Trilinos must be built with Amesos and MPI enabled to use MUMPS
   */
  enum coarse_solver CoarseSolver;

  /* Smoother type
   *  available types: Cheby, SymBlockGaussSeidel, Jacobi
   */
  enum smoother_type SmootherType;

  /* Whether HEC-MW smoother is used at finest level when SmootherType is SymBlockGaussSeidel
   */
  int FlgUsemonolisSmoother;

  /* Solver cycle
   *  available types: ML_MGV (V-cycle), ML_MGW (W-cycle), ML_MGFULLV (Full V-Cycle)
   */
  int MGType;

  /* Max num of levels
   */
  int MaxLevels;

  /* Coarsening scheme
   *  available types: UncoupledMIS, METIS, ParMETIS, Zoltan, DD
   */
  enum coarsen_scheme CoarsenScheme;

  /* Num of smoother sweeps (order of polynomial for Cheby)
   */
  int NumSweeps;

  /* Max coarse size
   */
  int MaxCoarseSize;
};

/* default values */
#ifdef HAVE_ML_AMESOS
# ifdef HAVE_AMESOS_MUMPS
#  define DEFAULT_COARSE_SOLVER MUMPS
# else
#  define DEFAULT_COARSE_SOLVER KLU
# endif
#else
# define DEFAULT_COARSE_SOLVER  Smoother
#endif
#define DEFAULT_SMOOTHER_TYPE   Cheby
#define DEFAULT_MG_TYPE         ML_MGW
#define DEFAULT_MAX_LEVELS      4
#define DEFAULT_COARSEN_SCHEME  UncoupledMIS
#define DEFAULT_NUM_SWEEPS      2
#define MAX_COARSE_SIZE_MUMPS   50000
#define MAX_COARSE_SIZE_KLU     10000

static void ml_options_set(struct ml_options *mlopt, int *ierr) {
  mlopt->CoarseSolver = DEFAULT_COARSE_SOLVER;
  mlopt->SmootherType = DEFAULT_SMOOTHER_TYPE;
  mlopt->MGType = DEFAULT_MG_TYPE;
  mlopt->CoarsenScheme = DEFAULT_COARSEN_SCHEME;
  mlopt->MaxLevels = 2;
  mlopt->NumSweeps = 2;
  mlopt->MaxCoarseSize = 10;
}

int monolis_ML_getrow_nn(ML_Operator *mat_in, int N_requested_rows,
                       int requested_rows[], int allocated_space,
                       int cols[], double values[], int row_lengths[]) {
  int *id, ierr;
  id = (int *)ML_Get_MyGetrowData(mat_in);
  //monolis_ml_getrow_nn_(id, &N_requested_rows, requested_rows, &allocated_space,
  //                    cols, values, row_lengths, &ierr);
  return ierr;
}

int monolis_ML_matvec_nn(ML_Operator *mat_in, int in_length, double p[],
                       int out_length, double ap[]) {
  int *id, ierr;
  id = (int *)ML_Get_MyGetrowData(mat_in);
  //monolis_ml_matvec_nn_(id, &in_length, p, &out_length, ap, &ierr);
  return ierr;
}

int monolis_ML_comm_nn(double x[], void *A_data) {
  int *id, ierr;
  id = (int *)A_data;
  //monolis_ml_comm_nn_(id, x, &ierr);
  return ierr;
}

struct ml_info {
  struct ml_options opt;
  ML *ml_object;
  ML_Aggregate *agg_object;
  int ndof;
};

static struct ml_info MLInfo;

/*
 * public functions
 */

void monolis_ML_wrapper_setup(int *sym, int *Ndof, int *ierr) {
  int loglevel, myrank;
  int N_grids, N_levels;
  int nlocal, nlocal_allcolumns;
  struct ml_options *mlopt;
  ML *ml_object;
  ML_Aggregate *agg_object;

  //monolis_Comm_rank(monolis_comm_get_comm(), &myrank);

  /* Get options */
  mlopt = &(MLInfo.opt);
  ml_options_set(mlopt, ierr);

  /* ML object */
  N_grids = mlopt->MaxLevels;
  ML_Create(&ml_object, N_grids);

  //monolis_ml_get_nlocal_(id, &nlocal, &nlocal_allcolumns, ierr);

  ML_Init_Amatrix(ml_object, 0, nlocal, nlocal, NULL);
  ML_Set_Amatrix_Getrow(ml_object, 0, monolis_ML_getrow_nn, monolis_ML_comm_nn, nlocal_allcolumns);
  ML_Set_Amatrix_Matvec(ml_object, 0, monolis_ML_matvec_nn);

  /* if (!(*sym)) ML_Set_Symmetrize(ml_object, ML_YES); */

  /* Aggregate */
  ML_Aggregate_Create(&agg_object);

  /* Max coarse size */
  {
    int nglobal;
    //monolis_Allreduce(&nlocal, &nglobal, 1, monolis_INT, monolis_SUM, monolis_comm_get_comm());
    if (nglobal <= mlopt->MaxCoarseSize) mlopt->MaxCoarseSize = nglobal - 1; /* coarsen at least once */
    if (mlopt->MaxCoarseSize > 0) ML_Aggregate_Set_MaxCoarseSize(agg_object, mlopt->MaxCoarseSize);
  }

  /* options */
  /* CoarsenScheme */
  {
    if (mlopt->CoarsenScheme == UncoupledMIS) {
      ML_Aggregate_Set_CoarsenScheme_UncoupledMIS(agg_object);
    } else if (mlopt->CoarsenScheme == METIS) {
      ML_Aggregate_Set_CoarsenScheme_METIS(agg_object);
    } else if (mlopt->CoarsenScheme == ParMETIS) {
      ML_Aggregate_Set_CoarsenScheme_ParMETIS(agg_object);
    } else if (mlopt->CoarsenScheme == Zoltan) {
      ML_Aggregate_Set_CoarsenScheme_Zoltan(agg_object);
    } else if (mlopt->CoarsenScheme == DD) {
      ML_Aggregate_Set_CoarsenScheme_DD(agg_object);
    }
    /*
    if (mlopt->MaxLevels == 2) {
      ML_Aggregate_Set_LocalNumber(ml_object, agg_object, ML_ALL_LEVELS, 1);
    } else if (mlopt->MaxLevels == 3) {
      ML_Aggregate_Set_NodesPerAggr(ml_object, agg_object, ML_ALL_LEVELS, 512);
      ML_Aggregate_Set_ReqLocalCoarseSize(ml_object->ML_num_levels, agg_object, ML_ALL_LEVELS, 128);
    }
    */
  }
  /* ML_Aggregate_Set_Threshold(agg_object, threshold); */
  /* ML_Aggregate_Set_DampingFactor(agg_object, dampingfactor); */

  /* eigen-analysis */
  /* ML_Set_SpectralNormScheme_PowerMethod(ml_object); */ /* power-method (default) */
  if (*sym) {
    ML_Set_SpectralNormScheme_Calc(ml_object); /* cg */
  }
  /* ML_Set_SpectralNorm_Iterations(ml_object, 10); */ /* default: 10 */

  ML_Aggregate_Set_Dimensions(agg_object, *Ndof);

  /* Generate MultiGrid */
  /* N_levels = ML_Gen_MGHierarchy_UsingAggregation(ml_object, 0, ML_INCREASING, agg_object); */
  N_levels = ML_Gen_MultiLevelHierarchy_UsingAggregation(ml_object, 0, ML_INCREASING, agg_object);
  if (loglevel >= 1 && myrank == 0) fprintf(stderr, "INFO: ML generated num of levels is %d\n", N_levels);

  /* Smoother */
  /*
   * Set pre- and post-smoother for each level
   *  level      : num in (0, N_levels-1) or ML_ALL_LEVELS
   *  pre-or-post: ML_PRESMOOTHER, ML_POSTSMOOTHER or ML_BOTH
   *  omega      : damping factor for Jacobi, GaussSeidel, etc. (ML_DEFAULT=1.0)
   */
  {
    int level;
    int coarsest_level = N_levels - 1;
    /*
     * levels other than the coarsest level
     */
    for (level = 0; level < coarsest_level; level++) {
      ML_Gen_Smoother_Cheby(ml_object, level, ML_BOTH, 20.0, mlopt->NumSweeps);
    }
    /*
     * coarsest level
     */

    if (mlopt->CoarseSolver == MUMPS) {
      ML_Gen_Smoother_Amesos(ml_object, coarsest_level, ML_AMESOS_MUMPS, 1, 0.0, 1);
    } else if (mlopt->CoarseSolver == KLU) {
      ML_Gen_Smoother_Amesos(ml_object, coarsest_level, ML_AMESOS_KLU, 1, 0.0, 1);
    } else {
      ML_Gen_Smoother_Cheby(ml_object, coarsest_level, ML_BOTH, 20.0, 2);
    }
  }

  /* Solver */
  ML_Gen_Solver(ml_object, mlopt->MGType, 0, N_levels - 1);

  /* Save objects */
  MLInfo.ml_object  = ml_object;
  MLInfo.agg_object = agg_object;
  MLInfo.ndof = *Ndof;
}

void monolis_ML_wrapper_apply(double rhs[], int *ierr) {
  int nlocal, nlocal_allcolumns;
  double *sol;
  int i;
  ML *ml_object;

  ml_object = MLInfo.ml_object;
  //monolis_ml_get_nlocal_(id, &nlocal, &nlocal_allcolumns, ierr);
  //if (*ierr != monolis_SUCCESS) return;
  //sol = (double *)monolis_malloc(sizeof(double) * nlocal_allcolumns);
  
  /* MultiGrid V-cycle */
  ML_Solve_MGV(ml_object, rhs, sol);

  for (i = 0; i < nlocal; i++) {
    rhs[i] = sol[i];
  }

  //monolis_free(sol);
}

void monolis_ML_wrapper_clear(int *ierr) {
  struct ml_options *mlopt = &(MLInfo.opt);

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
