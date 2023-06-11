!> 行列構造体の定義テスト
module mod_monolis_def_mat_test
  use mod_monolis_utils
  implicit none

contains

  subroutine monolis_def_mat_test()
    implicit none

    call monolis_std_global_log_string("monolis_mat_initialize")
    call monolis_std_global_log_string("monolis_mat_initialize_val_R")
    call monolis_std_global_log_string("monolis_mat_initialize_val_C")
    call monolis_std_global_log_string("monolis_mat_initialize_SCSR")
    call monolis_std_global_log_string("monolis_mat_initialize_CSR")
    call monolis_std_global_log_string("monolis_mat_initialize_CSC")
    call monolis_std_global_log_string("monolis_mat_finalize")
    call monolis_std_global_log_string("monolis_mat_finalize_val_R")
    call monolis_std_global_log_string("monolis_mat_finalize_val_C")
    call monolis_std_global_log_string("monolis_mat_finalize_SCSR")
    call monolis_std_global_log_string("monolis_mat_finalize_CSR")
    call monolis_std_global_log_string("monolis_mat_finalize_CSC")

    call monolis_std_global_log_string("monolis_set_RHS_R")
    call monolis_std_global_log_string("monolis_set_initial_solution_R")
    call monolis_std_global_log_string("monolis_set_initial_solution_C")
    call monolis_std_global_log_string("monolis_get_solution_R")
    call monolis_std_global_log_string("monolis_get_solution_C")

    call monolis_std_global_log_string("monolis_check_input_param")
  end subroutine monolis_def_mat_test

end module mod_monolis_def_mat_test