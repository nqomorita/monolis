program monolis_test
  use mod_monolis
  use mod_monolis_def_mat_test
  use mod_monolis_def_solver_test
  use mod_monolis_def_solver_util_test
  use mod_monolis_def_struc_test
  use mod_monolis_spmat_copy_test
  use mod_monolis_spmat_nonzero_pattern_util_test
  use mod_monolis_spmat_nonzero_pattern_test
  use mod_monolis_spmat_handler_util_test
  use mod_monolis_spmat_handler_test
  use mod_monolis_vec_util_test
  use mod_monolis_linalg_test
  use mod_monolis_converge_test
  use mod_monolis_matvec_test
  use mod_monolis_solver_CG_test
  use mod_monolis_solver_GropCG_test
  use mod_monolis_solver_PipeCG_test
  use mod_monolis_solver_PipeCR_test
  use mod_monolis_solver_BiCGSTAB_test
  use mod_monolis_solver_BiCGSTAB_noprec_test
  use mod_monolis_solver_PipeBiCGSTAB_test
  use mod_monolis_solver_PipeBiCGSTAB_noprec_test
  use mod_monolis_solver_COCG_test
  use mod_monolis_lapack_test
  use mod_monolis_scalapack_test
  use mod_monolis_eigen_lanczos_util_test
  use mod_monolis_eigen_solver_test
  use mod_monolis_precond_test
  use mod_monolis_solve_test
  implicit none

  call monolis_global_initialize()

  call monolis_def_mat_test()
  call monolis_def_solver_test()
  call monolis_def_solver_util_test()
  call monolis_def_struc_test()

  call monolis_spmat_copy_test()
  call monolis_spmat_nonzero_pattern_util_test()
  call monolis_spmat_nonzero_pattern_test()
  call monolis_spmat_handler_util_test()
  call monolis_spmat_handler_test()

  call monolis_vec_util_test()
  call monolis_linalg_test()
  call monolis_converge_test()
  call monolis_matvec_test()

  call monolis_precond_test()
  call monolis_solve_test()

  call monolis_lapack_test()
  call monolis_scalapack_test()

  if(monolis_mpi_get_global_comm_size() == 1)then
    call monolis_solver_CG_test()
    call monolis_solver_GropCG_test()
    call monolis_solver_PipeCG_test()
    call monolis_solver_PipeCR_test()
    call monolis_solver_BiCGSTAB_test()
    call monolis_solver_BiCGSTAB_noprec_test()
    call monolis_solver_PipeBiCGSTAB_test()
    call monolis_solver_PipeBiCGSTAB_noprec_test()
    call monolis_solver_COCG_test()

    call monolis_eigen_lanczos_util_test()
    call monolis_eigen_solver_test()
  endif

  call monolis_global_finalize()

end program monolis_test
