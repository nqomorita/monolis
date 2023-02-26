program monolis_test
  use mod_monolis
  use mod_monolis_def_mat_test
  use mod_monolis_def_solver_test
  use mod_monolis_def_solver_util_test
  use mod_monolis_def_struc_test
  implicit none

  call monolis_def_mat_test()
  call monolis_def_solver_test()
  call monolis_def_solver_util_test()
  call monolis_def_struc_test()

end program monolis_test
