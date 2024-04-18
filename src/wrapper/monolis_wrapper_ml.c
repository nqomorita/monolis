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

#include "monolis_wrapper_ml_util.h"
#include "monolis_wrapper_ml.h"

/*
 * Options
 */

enum coarse_solver {Smoother, KLU, MUMPS};
enum smoother_type {Cheby, SymBlockGaussSeidel, Jacobi};
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


static void ml_options_set(struct ml_options *mlopt, int *id, int myrank, int *ierr) {
  int opt[10];

  //monolis_ml_get_opt_(id, opt, ierr);
  //if (*ierr != monolis_SUCCESS) return;

  switch (opt[0]) {
  case 0:
    mlopt->CoarseSolver = DEFAULT_COARSE_SOLVER;
    break;
  case 1:
    mlopt->CoarseSolver = Smoother;
    break;
#ifdef HAVE_ML_AMESOS
  case 2:
    mlopt->CoarseSolver = KLU;
    break;
  case 3:
# ifdef HAVE_AMESOS_MUMPS
    mlopt->CoarseSolver = MUMPS;
    break;
# else
    if (myrank == 0) fprintf(stderr, "WARNING: MUMPS not available as coarse solver (rebuild Trilinos with MUMPS and MPI enabled)\n");
    break;
# endif
#else
  case 2:
    if (myrank == 0) fprintf(stderr, "WARNING: KLU not available as coarse solver (rebuild Trilinos with Amesos enabled)\n");
    break;
  case 3:
    if (myrank == 0) fprintf(stderr, "WARNING: MUMPS not available as coarse solver (rebuild Trilinos with Amesos, MUMPS and MPI enabled)\n");
    break;
#endif
  default:
    if (myrank == 0) fprintf(stderr, "WARNING: invalid ML_CoarseSolver=%d (ignored)\n", opt[0]);
  }

  switch (opt[1]) {
  case 0:
    mlopt->SmootherType = DEFAULT_SMOOTHER_TYPE;
    break;
  case 1:
    mlopt->SmootherType = Cheby;
    break;
  default:
    if (myrank == 0) fprintf(stderr, "WARNING: invalid ML_Smoother=%d (ignored)\n", opt[1]);
  }

  switch (opt[2]) {
  case 0:
    mlopt->MGType = DEFAULT_MG_TYPE;
    break;
  case 1:
    mlopt->MGType = ML_MGV;
    break;
  case 2:
    mlopt->MGType = ML_MGW;
    break;
  case 3:
    mlopt->MGType = ML_MGFULLV;
    break;
  default:
    if (myrank == 0) fprintf(stderr, "WARNING: invalid ML_MGCycle=%d (ignored)\n", opt[2]);
  }

  if (opt[3] > 0) {
    mlopt->MaxLevels = opt[3];
  } else {
    if (opt[3] < 0) {
      if (myrank == 0) fprintf(stderr, "WARNING: invalid ML_MaxLevels=%d (ignored)\n", opt[3]);
    }
    mlopt->MaxLevels = DEFAULT_MAX_LEVELS;
  }

  switch (opt[4]) {
  case 0:
    mlopt->CoarsenScheme = DEFAULT_COARSEN_SCHEME;
    break;
  case 1:
    mlopt->CoarsenScheme = UncoupledMIS;
    break;
  case 2:
    mlopt->CoarsenScheme = METIS;
    break;
  case 3:
    mlopt->CoarsenScheme = ParMETIS;
    break;
  case 4:
    mlopt->CoarsenScheme = Zoltan;
    break;
  case 5:
    mlopt->CoarsenScheme = DD;
    break;
  default:
    if (myrank == 0) fprintf(stderr, "WARNING: invalid ML_CoarseningScheme=%d (ignored)\n", opt[4]);
  }

  if (opt[5] > 0) {
    mlopt->NumSweeps = opt[5];
  } else {
    if (opt[5] < 0) {
      if (myrank == 0) fprintf(stderr, "WARNING: invalid ML_NumSweep=%d (ignored)\n", opt[5]);
    }
    mlopt->NumSweeps = DEFAULT_NUM_SWEEPS;
    opt[5] = mlopt->NumSweeps;
    //monolis_ml_set_opt_(id, opt, ierr);
    //if (*ierr != monolis_SUCCESS) return;
  }

  if (opt[6] > 0) {
    mlopt->MaxCoarseSize = opt[6];
  } else {
    if (mlopt->CoarseSolver == MUMPS) {
      mlopt->MaxCoarseSize = MAX_COARSE_SIZE_MUMPS;
    } else if (mlopt->CoarseSolver == KLU) {
      mlopt->MaxCoarseSize = MAX_COARSE_SIZE_KLU;
    } else {
      mlopt->MaxCoarseSize = -1; /* use default (128? 32?) */
    }
  }
}

void ml_options_print(struct ml_options *mlopt, FILE *fp, int myrank, int loglevel) {
  char optstr[7][32];
  switch (mlopt->CoarseSolver) {
  case Smoother:
    if (loglevel >= 2 && myrank == 0) fprintf(fp, "INFO: ML coarse solver is smoother\n");
    strcpy(optstr[0], "Smoother");
    break;
  case KLU:
    if (loglevel >= 2 && myrank == 0) fprintf(fp, "INFO: ML coarse solver is KLU\n");
    strcpy(optstr[0], "KLU");
    break;
  case MUMPS:
    if (loglevel >= 2 && myrank == 0) fprintf(fp, "INFO: ML coarse solver is MUMPS\n");
    strcpy(optstr[0], "MUMPS");
    break;
  }
  switch (mlopt->SmootherType) {
  case Cheby:
    if (loglevel >= 2 && myrank == 0) fprintf(fp, "INFO: ML smoother is Cheby\n");
    strcpy(optstr[1], "Cheby");
    break;
  case SymBlockGaussSeidel:
    if (loglevel >= 2 && myrank == 0) fprintf(fp, "INFO: ML smoother is SymBlockGaussSeidel\n");
    strcpy(optstr[1], "SymBlockGaussSeidel");
    break;
  case Jacobi:
    if (loglevel >= 2 && myrank == 0) fprintf(fp, "INFO: ML smoother is Jacobi\n");
    strcpy(optstr[1], "Jacobi");
    break;
  }
  switch (mlopt->MGType) {
  case ML_MGV:
    if (loglevel >= 2 && myrank == 0) fprintf(fp, "INFO: ML multigrid type is V-cycle\n");
    strcpy(optstr[2], "V-cycle");
    break;
  case ML_MGW:
    if (loglevel >= 2 && myrank == 0) fprintf(fp, "INFO: ML multigrid type is W-cycle\n");
    strcpy(optstr[2], "W-cycle");
    break;
  case ML_MGFULLV:
    if (loglevel >= 2 && myrank == 0) fprintf(fp, "INFO: ML multigrid type is Full-V-cycle\n");
    strcpy(optstr[2], "Full-V-cycle");
    break;
  }
  if (loglevel >= 2 && myrank == 0) fprintf(fp, "INFO: ML num of max levels is %d\n", mlopt->MaxLevels);
  sprintf(optstr[3], "MaxLevels=%d", mlopt->MaxLevels);
  switch (mlopt->CoarsenScheme) {
  case UncoupledMIS:
    if (loglevel >= 2 && myrank == 0) fprintf(fp, "INFO: ML coarsening scheme is UncoupledMIS\n");
    strcpy(optstr[4], "UncoupledMIS");
    break;
  case METIS:
    if (loglevel >= 2 && myrank == 0) fprintf(fp, "INFO: ML coarsening scheme is METIS\n");
    strcpy(optstr[4], "METIS");
    break;
  case ParMETIS:
    if (loglevel >= 2 && myrank == 0) fprintf(fp, "INFO: ML coarsening scheme is ParMETIS\n");
    strcpy(optstr[4], "ParMETIS");
    break;
  case Zoltan:
    if (loglevel >= 2 && myrank == 0) fprintf(fp, "INFO: ML coarsening scheme is Zoltan\n");
    strcpy(optstr[4], "Zoltan");
    break;
  case DD:
    if (loglevel >= 2 && myrank == 0) fprintf(fp, "INFO: ML coarsening scheme is DD\n");
    strcpy(optstr[4], "DD");
    break;
  }
  if (loglevel >= 2 && myrank == 0) fprintf(fp, "INFO: ML num of smoother sweeps is %d\n", mlopt->NumSweeps);
  sprintf(optstr[5], "NumSweeps=%d", mlopt->NumSweeps);
  if (loglevel >= 2 && myrank == 0) fprintf(fp, "INFO: ML max coarse size is %d\n", mlopt->MaxCoarseSize);
  sprintf(optstr[6], "MaxCoarseSize=%d", mlopt->MaxCoarseSize);
  if (loglevel >= 1 && myrank == 0) {
    fprintf(fp, "INFO: ML options: %s %s %s %s %s %s %s\n",
            optstr[0], optstr[1], optstr[2], optstr[3], optstr[4], optstr[5], optstr[6]);
  }
}

/*
 * static variable
 */

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

void monolis_ML_wrapper_setup(int *id, int *sym, int *Ndof, int *ierr) {
  int loglevel, myrank;
  int N_grids, N_levels;
  int nlocal, nlocal_allcolumns;
  struct ml_options *mlopt;
  ML *ml_object;
  ML_Aggregate *agg_object;

  //monolis_ml_get_loglevel_(id, &loglevel);
  ML_Set_PrintLevel(loglevel);

  //monolis_Comm_rank(monolis_comm_get_comm(), &myrank);

  /* Get options */
  mlopt = &(MLInfo.opt);
  ml_options_set(mlopt, id, myrank, ierr);
  ml_options_print(mlopt, stderr, myrank, loglevel);

  /* ML object */
  N_grids = mlopt->MaxLevels;
  ML_Create(&ml_object, N_grids);
  //monolis_ml_get_nlocal_(id, &nlocal, &nlocal_allcolumns, ierr);
  //if (*ierr != monolis_SUCCESS) return;
  ML_Init_Amatrix(ml_object, 0, nlocal, nlocal, id);
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

void monolis_ML_wrapper_apply(int *id, double rhs[], int *ierr) {
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

void monolis_ML_wrapper_clear(int *id, int *ierr) {
  struct ml_options *mlopt = &(MLInfo.opt);

  ML_Aggregate_Destroy(&(MLInfo.agg_object));
  ML_Destroy(&(MLInfo.ml_object));
}

/* Fortran interface */
void monolis_precond_ml_setup_(int *id, int *sym, int *ndof, int *ierr) {
  monolis_ML_wrapper_setup(id, sym, ndof, ierr);
}
void monolis_precond_ml_setup__(int *id, int *sym, int *ndof, int *ierr) {
  monolis_ML_wrapper_setup(id, sym, ndof, ierr);
}
void MONOLIS_PRECOND_ML_SETUP(int *id, int *sym, int *ndof, int *ierr) {
  monolis_ML_wrapper_setup(id, sym, ndof, ierr);
}

void monolis_precond_ml_apply_(int *id, double rhs[], int *ierr) {
  monolis_ML_wrapper_apply(id, rhs, ierr);
}
void monolis_precond_ml_apply__(int *id, double rhs[], int *ierr) {
  monolis_ML_wrapper_apply(id, rhs, ierr);
}
void MONOLIS_PRECOND_ML_APPLY(int *id, double rhs[], int *ierr) {
  monolis_ML_wrapper_apply(id, rhs, ierr);
}

void monolis_precond_ml_clear_(int *id, int *ierr) {
  monolis_ML_wrapper_clear(id, ierr);
}
void monolis_precond_ml_clear__(int *id, int *ierr) {
  monolis_ML_wrapper_clear(id, ierr);
}
void MONOLIS_ML_WRAPPER_CLEAR(int *id, int *ierr) {
  monolis_ML_wrapper_clear(id, ierr);
}
