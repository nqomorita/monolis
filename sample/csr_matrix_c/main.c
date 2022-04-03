#include <stdio.h>
#include <stdlib.h>
#include "monolis.h"

int main(int argc, char *args[]) {
  MONOLIS monolis;
  const char* dir_name;
  double time[7];

  dir_name = ".";
  printf("* test: monolis in C\n");

  printf("* monolis_global_initialize\n");
  monolis_global_initialize();

  printf("\n");
  printf("* test: first touch \n");
  printf("* monolis_initialize\n");
  monolis_initialize(&monolis, dir_name);

  printf("* monolis_finalize\n");
  monolis_finalize(&monolis);

  printf("\n");
  printf("* test: second touch \n");
  printf("* monolis_initialize\n");
  monolis_initialize(&monolis, dir_name);

  //monolis_set_method   (monolis, solver_type);
  //monolis_set_precond  (monolis, precond_type);
  //monolis_set_maxiter  (monolis, num_max_iters);
  //monolis_set_tolerance(monolis, epsilon);
  //monolis_set_init_x   (monolis, false);

  monolis_solve(
      &monolis,
      monolis.mat.B,
      monolis.mat.X);

  monolis_get_time_solver           (&monolis, &time[0]);
  monolis_get_time_preparing        (&monolis, &time[1]);
  monolis_get_time_spmv             (&monolis, &time[2]);
  monolis_get_time_dot_product      (&monolis, &time[3]);
  monolis_get_time_precondition     (&monolis, &time[4]);
  monolis_get_time_comm_dot_product (&monolis, &time[5]);
  monolis_get_time_comm_spmv        (&monolis, &time[6]);

  printf("* monolis_get_time_solver           %e\n", time[0]);
  printf("* monolis_get_time_preparing        %e\n", time[1]);
  printf("* monolis_get_time_spmv             %e\n", time[2]);
  printf("* monolis_get_time_dot_product      %e\n", time[3]);
  printf("* monolis_get_time_precondition     %e\n", time[4]);
  printf("* monolis_get_time_comm_dot_product %e\n", time[5]);
  printf("* monolis_get_time_comm_spmv        %e\n", time[6]);

  printf("* monolis_finalize\n");
  monolis_finalize(&monolis);

  printf("* monolis_global_finalize\n");
  monolis_global_finalize();
}
