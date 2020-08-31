/* monolis.h */
#ifndef MONOLIS_H
#define MONOLIS_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>
#include "monolis_struct.h"
#include "monolis_c_main.h"

/* initializer */
/* call only once at the beginning of a program */
void monolis_global_initialize();

/* call only once at the end of a program */
void monolis_global_finalize();



/* global getter */
/* get MPI world communicator */
int monolis_get_global_comm();

/* get MPI communicator size */
int monolis_get_global_commsize();

/* get MPI rank */
int monolis_get_global_myrank();

/* get time */
double monolis_get_time();

/* get input file name for parallel computation */
const char* monolis_get_input_filename(
  const char* filename_body);



/* matrix initializer */
/* initialize a monolis structure */
void monolis_initialize(
  MONOLIS* mat);

/* finalize a monolis structure */
void monolis_finalize(
  MONOLIS* mat);



/* parameter setter */
/* set iterative method
  default: monolis_iter_CG
  flags:
    - monolis_iter_CG:           for symmetric positive define matrix
    - monolis_iter_GropCG:       for symmetric positive define matrix
    - monolis_iter_PipeCG:       for symmetric positive define matrix
    - monolis_iter_PipeCR:       for symmetric matrix
    - monolis_iter_BiCGSTAB:     for non-symmetric matrix
    - monolis_iter_PipeBiCGSTAB: for non-symmetric matrix */
void monolis_set_method   (MONOLIS* mat, int    flag);

/* set predonditioner
  default: monolis_prec_DIAG
  flags:
    - monolis_prec_NONE
    - monolis_prec_DIAG
    - monolis_prec_ILU
    - monolis_prec_JACOBI
    - monolis_prec_SOR */
void monolis_set_precond  (MONOLIS* mat, int    flag);

/* set maximum iteration
  default: 10000 */
void monolis_set_maxiter  (MONOLIS* mat, int    flag);

/* set convergence threshold
  default: 1.0e-8 */
void monolis_set_tolerance(MONOLIS* mat, double value);

/* show convergence history (rerative residual per iteration)
  default: false */
void monolis_show_iterlog (MONOLIS* mat, bool   flag);

/* show execution time in details
  default: false */
void monolis_show_timelog (MONOLIS* mat, bool   flag);

/* show summary log
  default: true */
void monolis_show_summary (MONOLIS* mat, bool   flag);



/* matrix clearing handler */
/* set all components of A, x, and b in Ax = b to zero */
void monolis_clear(
  MONOLIS* mat);

/* set all components of A in Ax = b to zero */
void monolis_clear_mat(
  MONOLIS* mat);

/* set all components of b in Ax = b to zero */
void monolis_clear_rhs(
  MONOLIS* mat);

/* set all components of x in Ax = b to zero */
void monolis_clear_solution(
  MONOLIS* mat);



/* matrix copy handler */
/* copy all components of monolis structure */
void monolis_copy_all(
  MONOLIS* in,
  MONOLIS* out);

/* copy all components of monolis structure, and
all components of A, x, and b in Ax = b is set to zero */
void monolis_copy_nonzero_pattern(
  MONOLIS* in,
  MONOLIS* out);



/* sparse matrix handler */
/* get nonzero pattern from connectivity */
void monolis_get_nonzero_pattern(
  MONOLIS* mat,
  int      nnode,
  int      nbase_func,
  int      ndof,
  int      nelem,
  int**    elem);

/* add scalar value to sparse matrix */
void monolis_add_scalar_to_sparse_matrix(
  MONOLIS* mat,
  double   val,
  int      i,
  int      j,
  int      submat_i,
  int      submat_j);

/* set Dirichlet boundary condition to sparse matrix */
void monolis_set_Dirichlet_bc(
  MONOLIS* mat,
  double*  b,
  int      node_id,
  int      ndof_bc,
  double   val);



/* solver */
void monolis_solve(
  MONOLIS* mat,
  double*  b,
  double*  x);

#ifdef __cplusplus
}
#endif

#endif
