#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "gedatsu.h"
#include "monolis_utils.h"
#include "./define/monolis_def_solver_c_test.h"
#include "./define/monolis_def_solver_util_c_test.h"

int main()
{
  monolis_std_log_string("monolis_c_test");

  monolis_mpi_initialize();

  monolis_def_solver_c_test();

  monolis_def_solver_util_c_test();

  monolis_mpi_finalize();
}
