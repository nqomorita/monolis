!> 収束判定テストモジュール
module mod_monolis_converge_test
  use mod_monolis
  use mod_monolis_converge

  implicit none

contains

  subroutine monolis_converge_test()
    implicit none
    type(monolis_structure) :: monolis
    type(monolis_com) :: com
    real(kdouble) :: R(4)
    real(kdouble) :: B2
    integer(kint) :: iter
    logical :: is_converge
    real(kdouble) :: tdotp
    real(kdouble) :: tcomm

    call monolis_std_global_log_string("monolis_check_converge_R")

    call monolis_initialize(monolis)
    call monolis_com_initialize_by_self(com)

    call monolis_com_set_communicator(com, monolis_mpi_get_global_comm())
    call monolis_com_set_my_rank(com, monolis_mpi_get_global_my_rank())
    call monolis_com_set_comm_size(com, monolis_mpi_get_global_comm_size())

    monolis%MAT%N = 2
    monolis%MAT%NDOF = 2

    B2 = 4.0d0

    R(1) = 2.0d0
    R(2) = 2.0d0
    R(3) = 2.0d0
    R(4) = 2.0d0

    iter = 10

    call monolis_check_converge_R(monolis%PRM, com, monolis%MAT, R, B2, iter, is_converge, tdotp, tcomm)

    if(monolis_mpi_get_global_comm_size() == 2)then
      call monolis_test_check_eq_I1("monolis_converge_test 1", monolis%PRM%Iarray(monolis_prm_I_cur_iter), 10)
      call monolis_test_check_eq_R1("monolis_converge_test 2", monolis%PRM%Rarray(monolis_prm_R_cur_resid), dsqrt(8.0d0))
      call monolis_test_check_eq_L1("monolis_converge_test 3", is_converge, .false.)
    else
      call monolis_test_check_eq_I1("monolis_converge_test 1", monolis%PRM%Iarray(monolis_prm_I_cur_iter), 10)
      call monolis_test_check_eq_R1("monolis_converge_test 2", monolis%PRM%Rarray(monolis_prm_R_cur_resid), 2.0d0)
      call monolis_test_check_eq_L1("monolis_converge_test 3", is_converge, .false.)
      endif
  end subroutine monolis_converge_test

end module mod_monolis_converge_test
