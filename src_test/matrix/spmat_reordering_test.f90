!> 疎行列操作関数群
module mod_monolis_spmat_reordering_test
  use mod_monolis

  implicit none

contains

  subroutine monolis_spmat_reordering_test()
    implicit none

    call monolis_std_global_log_string("monolis_matrix_reordering_fw_R")
    call monolis_std_global_log_string("monolis_matrix_reordering_bk_R")
    call monolis_std_global_log_string("monolis_reorder_vector_fw")
    call monolis_std_global_log_string("monolis_reorder_back_vector_bk")
    call monolis_std_global_log_string("monolis_restruct_matrix")
    call monolis_std_global_log_string("monolis_restruct_matrix_profile")
    call monolis_std_global_log_string("monolis_restruct_matrix_values")
  end subroutine monolis_spmat_reordering_test

end module mod_monolis_spmat_reordering_test
