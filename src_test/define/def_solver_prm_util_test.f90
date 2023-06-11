!> ソルバパラメータの設定関数群テスト
module mod_monolis_def_solver_util_test
  use mod_monolis

  implicit none

contains

  subroutine monolis_def_solver_util_test()
    implicit none
    type(monolis_structure) :: monolis
    integer(kint) :: i_param
    real(kdouble) :: r_param

    call monolis_std_global_log_string("monolis_set_method")
    call monolis_std_global_log_string("monolis_set_precond")
    call monolis_std_global_log_string("monolis_set_maxiter")
    call monolis_std_global_log_string("monolis_set_init_x")
    call monolis_std_global_log_string("monolis_set_sym_matrix")
    call monolis_std_global_log_string("monolis_set_debug")
    call monolis_std_global_log_string("monolis_set_performance_measurement")
    call monolis_std_global_log_string("monolis_set_check_diag")
    call monolis_std_global_log_string("monolis_set_prec_stored")
    call monolis_std_global_log_string("monolis_show_iterlog")
    call monolis_std_global_log_string("monolis_show_timelog")
    call monolis_std_global_log_string("monolis_show_summary")
    call monolis_std_global_log_string("monolis_show_timelog_statistics")
    call monolis_std_global_log_string("monolis_get_converge_iter")
    call monolis_std_global_log_string("monolis_get_error_tag")
    call monolis_std_global_log_string("monolis_set_tolerance")
    call monolis_std_global_log_string("monolis_get_converge_residual")
    call monolis_std_global_log_string("monolis_get_time_solver")
    call monolis_std_global_log_string("monolis_get_time_preparing")
    call monolis_std_global_log_string("monolis_get_time_spmv")
    call monolis_std_global_log_string("monolis_get_time_inner_product")
    call monolis_std_global_log_string("monolis_get_time_precondition")
    call monolis_std_global_log_string("monolis_get_time_comm_inner_product")
    call monolis_std_global_log_string("monolis_get_time_comm_spmv")

    monolis%PRM%Iarray = 0

    call monolis_set_method(monolis, 1)
    call monolis_set_precond(monolis, 2)
    call monolis_set_maxiter(monolis, 3)
    !call monolis_set_scaling(monolis, 4)
    !call monolis_set_reordering(monolis, 5)

    call monolis_test_check_eq_I1("monolis_def_solver_util_test 1", monolis%PRM%Iarray(monolis_prm_I_method), 1)
    call monolis_test_check_eq_I1("monolis_def_solver_util_test 2", monolis%PRM%Iarray(monolis_prm_I_precond), 2)
    call monolis_test_check_eq_I1("monolis_def_solver_util_test 3", monolis%PRM%Iarray(monolis_prm_I_max_iter), 3)
    !call monolis_test_check_eq_I1("monolis_def_solver_util_test 4", monolis%PRM%Iarray(monolis_prm_I_is_scaling), 4)
    !call monolis_test_check_eq_I1("monolis_def_solver_util_test 5", monolis%PRM%Iarray(monolis_prm_I_is_reordering), 5)

    monolis%PRM%Iarray = 0
    call monolis_set_init_x(monolis, .true.)
    call monolis_test_check_eq_I1("monolis_def_solver_util_test 6", monolis%PRM%Iarray(monolis_prm_I_is_init_x), 1)

    monolis%PRM%Iarray = 0
    call monolis_set_sym_matrix(monolis, .true.)
    call monolis_test_check_eq_I1("monolis_def_solver_util_test 7", monolis%PRM%Iarray(monolis_prm_I_is_sym_matrix), 1)

    monolis%PRM%Iarray = 0
    call monolis_set_debug(monolis, .true.)
    call monolis_test_check_eq_I1("monolis_def_solver_util_test 8", monolis%PRM%Iarray(monolis_prm_I_is_debug), 1)

    monolis%PRM%Iarray = 0
    call monolis_set_performance_measurement(monolis, .true.)
    call monolis_test_check_eq_I1("monolis_def_solver_util_test 9", monolis%PRM%Iarray(monolis_prm_I_is_measurement), 1)

    monolis%PRM%Iarray = 0
    call monolis_set_check_diag(monolis, .true.)
    call monolis_test_check_eq_I1("monolis_def_solver_util_test 10", monolis%PRM%Iarray(monolis_prm_I_is_check_diag), 1)

    monolis%PRM%Iarray = 0
    call monolis_set_prec_stored(monolis, .true.)
    call monolis_test_check_eq_I1("monolis_def_solver_util_test 11", monolis%PRM%Iarray(monolis_prm_I_is_prec_stored), 1)

    monolis%PRM%Iarray = 0
    call monolis_show_iterlog(monolis, .true.)
    call monolis_test_check_eq_I1("monolis_def_solver_util_test 12", monolis%PRM%Iarray(monolis_prm_I_show_iterlog), 1)

    monolis%PRM%Iarray = 0
    call monolis_show_timelog(monolis, .true.)
    call monolis_test_check_eq_I1("monolis_def_solver_util_test 13", monolis%PRM%Iarray(monolis_prm_I_show_time), 1)

    monolis%PRM%Iarray = 0
    call monolis_show_summary(monolis, .true.)
    call monolis_test_check_eq_I1("monolis_def_solver_util_test 14", monolis%PRM%Iarray(monolis_prm_I_show_summary), 1)

    monolis%PRM%Iarray = 0
    call monolis_show_timelog_statistics(monolis, .true.)
    call monolis_test_check_eq_I1("monolis_def_solver_util_test 15", monolis%PRM%Iarray(monolis_prm_I_show_time_statistics), 1)

    monolis%PRM%Iarray = 0
    monolis%PRM%Iarray(monolis_prm_I_cur_iter) = 1
    monolis%PRM%Iarray(monolis_prm_I_ierr) = 2

    call monolis_get_converge_iter(monolis, i_param)

    call monolis_test_check_eq_I1("monolis_def_solver_util_test 16", i_param, 1)

    call monolis_get_error_tag(monolis, i_param)

    call monolis_test_check_eq_I1("monolis_def_solver_util_test 17", i_param, 2)

    monolis%PRM%Rarray = 0.0d0

    call monolis_set_tolerance(monolis, 1.0d0)
    call monolis_test_check_eq_R1("monolis_def_solver_util_test 18", monolis%PRM%Rarray(monolis_prm_R_tol), 1.0d0)

    monolis%PRM%Rarray = 0.0d0
    monolis%PRM%Rarray(monolis_prm_R_cur_resid) = 1.0d0
    monolis%PRM%Rarray(monolis_R_time_sol) = 2.0d0
    monolis%PRM%Rarray(monolis_R_time_prep) = 3.0d0
    monolis%PRM%Rarray(monolis_R_time_spmv) = 4.0d0
    monolis%PRM%Rarray(monolis_R_time_dotp) = 5.0d0
    monolis%PRM%Rarray(monolis_R_time_prec) = 6.0d0
    monolis%PRM%Rarray(monolis_R_time_comm_dotp) = 7.0d0
    monolis%PRM%Rarray(monolis_R_time_comm_spmv) = 8.0d0

    call monolis_get_converge_residual(monolis, r_param)
    call monolis_test_check_eq_R1("monolis_def_solver_util_test 19", r_param, 1.0d0)

    call monolis_get_time_solver(monolis, r_param)
    call monolis_test_check_eq_R1("monolis_def_solver_util_test 20", r_param, 2.0d0)

    call monolis_get_time_preparing(monolis, r_param)
    call monolis_test_check_eq_R1("monolis_def_solver_util_test 21", r_param, 3.0d0)

    call monolis_get_time_spmv(monolis, r_param)
    call monolis_test_check_eq_R1("monolis_def_solver_util_test 22", r_param, 4.0d0)

    call monolis_get_time_inner_product(monolis, r_param)
    call monolis_test_check_eq_R1("monolis_def_solver_util_test 23", r_param, 5.0d0)

    call monolis_get_time_precondition(monolis, r_param)
    call monolis_test_check_eq_R1("monolis_def_solver_util_test 24", r_param, 6.0d0)

    call monolis_get_time_comm_inner_product(monolis, r_param)
    call monolis_test_check_eq_R1("monolis_def_solver_util_test 25", r_param, 7.0d0)

    call monolis_get_time_comm_spmv(monolis, r_param)
    call monolis_test_check_eq_R1("monolis_def_solver_util_test 26", r_param, 8.0d0)
  end subroutine monolis_def_solver_util_test

end module mod_monolis_def_solver_util_test
