!> ベクトル内積テストモジュール
module mod_monolis_linalg_test
  use mod_monolis
  use mod_monolis_inner_product
  implicit none

contains

  subroutine monolis_linalg_test()
    implicit none
    type(monolis_structure) :: monolis
    type(monolis_com) :: com
    integer(kint) :: n
    integer(kint) :: ndof
    integer(kint) :: iX(4), iY(4), isum
    real(kdouble) :: rX(4), rY(4), rsum
    complex(kdouble) :: cX(4), cY(4), csum

    call monolis_std_global_log_string("monolis_inner_product_main_C")
    call monolis_std_global_log_string("monolis_inner_product_main_I")
    call monolis_std_global_log_string("monolis_inner_product_main_R")
    call monolis_std_global_log_string("monolis_inner_product_main_R_no_comm")

    call monolis_std_global_log_string("monolis_inner_product_I")
    call monolis_std_global_log_string("monolis_inner_product_R")
    call monolis_std_global_log_string("monolis_inner_product_C")
    call monolis_std_global_log_string("monolis_inner_productV_I")
    call monolis_std_global_log_string("monolis_inner_productV_R")
    call monolis_std_global_log_string("monolis_inner_productV_C")

    call monolis_initialize(monolis)
    call monolis_com_initialize_by_self(com)

    call monolis_com_set_communicator(com, monolis_mpi_get_global_comm())
    call monolis_com_set_my_rank(com, monolis_mpi_get_global_my_rank())
    call monolis_com_set_comm_size(com, monolis_mpi_get_global_comm_size())
    call monolis_com_set_n_internal_vertex(com, 2)

    monolis%MAT%N = 2

    ndof = 2

    iX(1) = 1; iY(1) = 1
    iX(2) = 1; iY(2) = 2
    iX(3) = 1; iY(3) = 3
    iX(4) = 1; iY(4) = 4

    call monolis_inner_product_I(monolis, com, ndof, iX, iY, isum)

    if(monolis_mpi_get_global_comm_size() == 2)then
      call monolis_test_check_eq_I1("monolis_linalg_test 1", isum, 20)
    else
      call monolis_test_check_eq_I1("monolis_linalg_test 1", isum, 10)
    endif

    rX(1) = 1.0d0; rY(1) = 1.0d0
    rX(2) = 1.0d0; rY(2) = 2.0d0
    rX(3) = 1.0d0; rY(3) = 3.0d0
    rX(4) = 1.0d0; rY(4) = 4.0d0

    call monolis_inner_product_R(monolis, com, ndof, rX, rY, rsum)

    if(monolis_mpi_get_global_comm_size() == 2)then
      call monolis_test_check_eq_R1("monolis_linalg_test 2", rsum, 20.0d0)
    else
      call monolis_test_check_eq_R1("monolis_linalg_test 2", rsum, 10.0d0)
    endif

    cX(1) = (1.0d0, 0.0d0); cY(1) = (1.0d0, 1.0d0)
    cX(2) = (1.0d0, 0.0d0); cY(2) = (2.0d0, 2.0d0)
    cX(3) = (1.0d0, 0.0d0); cY(3) = (3.0d0, 3.0d0)
    cX(4) = (1.0d0, 0.0d0); cY(4) = (4.0d0, 4.0d0)

    call monolis_inner_product_C(monolis, com, ndof, cX, cY, csum)

    if(monolis_mpi_get_global_comm_size() == 2)then
      call monolis_test_check_eq_C1("monolis_linalg_test 3", csum, (20.0d0, 20.0d0))
    else
      call monolis_test_check_eq_C1("monolis_linalg_test 3", csum, (10.0d0, 10.0d0))
    endif

    n = 2

    ndof = 2

    iX(1) = 1; iY(1) = 1
    iX(2) = 1; iY(2) = 2
    iX(3) = 1; iY(3) = 3
    iX(4) = 1; iY(4) = 4

    call monolis_inner_productV_I(monolis, com, n, ndof, iX, iY, isum)

    if(monolis_mpi_get_global_comm_size() == 2)then
      call monolis_test_check_eq_I1("monolis_linalg_test 4", isum, 20)
    else
      call monolis_test_check_eq_I1("monolis_linalg_test 4", isum, 10)
    endif

    rX(1) = 1.0d0; rY(1) = 1.0d0
    rX(2) = 1.0d0; rY(2) = 2.0d0
    rX(3) = 1.0d0; rY(3) = 3.0d0
    rX(4) = 1.0d0; rY(4) = 4.0d0

    call monolis_inner_productV_R(monolis, com, n, ndof, rX, rY, rsum)

    if(monolis_mpi_get_global_comm_size() == 2)then
      call monolis_test_check_eq_R1("monolis_linalg_test 5", rsum, 20.0d0)
    else
      call monolis_test_check_eq_R1("monolis_linalg_test 5", rsum, 10.0d0)
    endif

    cX(1) = (1.0d0, 0.0d0); cY(1) = (1.0d0, 1.0d0)
    cX(2) = (1.0d0, 0.0d0); cY(2) = (2.0d0, 2.0d0)
    cX(3) = (1.0d0, 0.0d0); cY(3) = (3.0d0, 3.0d0)
    cX(4) = (1.0d0, 0.0d0); cY(4) = (4.0d0, 4.0d0)

    call monolis_inner_productV_C(monolis, com, n, ndof, cX, cY, csum)

    if(monolis_mpi_get_global_comm_size() == 2)then
      call monolis_test_check_eq_C1("monolis_linalg_test 6", csum, (20.0d0, 20.0d0))
    else
      call monolis_test_check_eq_C1("monolis_linalg_test 6", csum, (10.0d0, 10.0d0))
    endif

    rX(1) = 1.0d0; rY(1) = 1.0d0
    rX(2) = 1.0d0; rY(2) = 2.0d0
    rX(3) = 1.0d0; rY(3) = 3.0d0
    rX(4) = 1.0d0; rY(4) = 4.0d0

    call monolis_inner_product_main_R_no_comm(n, ndof, rX, rY, rsum)

    call monolis_test_check_eq_R1("monolis_linalg_test 7", rsum, 10.0d0)
  end subroutine monolis_linalg_test
end module mod_monolis_linalg_test
