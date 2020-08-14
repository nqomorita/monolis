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
void monolis_global_initialize();

void monolis_global_finalize();

/* mat initializer */
void monolis_initialize(
  MONOLIS* mat);

void monolis_finalize(
  MONOLIS* mat);

/* getter */
int monolis_get_global_comm();

int monolis_get_global_commsize();

int monolis_get_global_myrank();

double monolis_get_time();

void monolis_get_input_filename();

/* setter */
void monolis_set_method   (MONOLIS* mat, int    flag);

void monolis_set_precond  (MONOLIS* mat, int    flag);

void monolis_set_maxiter  (MONOLIS* mat, int    flag);

void monolis_set_tolerance(MONOLIS* mat, double flag);

void monolis_show_iterlog (MONOLIS* mat, bool   flag);

void monolis_show_timelog (MONOLIS* mat, bool   flag);

void monolis_show_summary (MONOLIS* mat, bool   flag);

/* mat clear */
void monolis_clear(
  MONOLIS* mat);

void monolis_clear_mat(
  MONOLIS* mat);

void monolis_clear_rhs(
  MONOLIS* mat);

void monolis_clear_solution(
  MONOLIS* mat);

/* mat copy */
void monolis_copy_all(
  MONOLIS* in,
  MONOLIS* out);

void monolis_copy_nonzero_pattern(
  MONOLIS* in,
  MONOLIS* out);

/* sparse matrix handler */
void monolis_get_nonzero_pattern(
  MONOLIS* mat,
  int      nnode,
  int      nbase_func,
  int      ndof,
  int      nelem,
  int**    elem);

void monolis_add_scalar_to_sparse_matrix(
  MONOLIS* mat,
  double   val,
  int      i,
  int      j,
  int      submat_i,
  int      submat_j);

/* boundary condition handler */
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
