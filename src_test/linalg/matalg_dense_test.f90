!> 疎行列ベクトル積テストモジュール
module mod_monolis_matalg_dense_test
  use mod_monolis
  use mod_monolis_utils
  use mod_monolis_matalg_dense
  implicit none

contains

  subroutine monolis_matalg_dense_test()
    implicit none

    if(monolis_mpi_get_global_comm_size() == 2) return

    call monolis_dense_matvec_local_R_test()
    call monolis_dense_matmul_local_R_test()

    call monolis_std_global_log_string("monolis_dense_matvec_local_R")
    call monolis_std_global_log_string("monolis_dense_matmul_local_R")
  end subroutine monolis_matalg_dense_test

  subroutine monolis_dense_matvec_local_R_test()
    implicit none
    real(kdouble) :: a(5), b(5), b_th(5), mat_dense(5,5), time

    call monolis_std_global_log_string("monolis_dense_matvec_local_R")

    mat_dense(1, 1) = 2.0d0
    mat_dense(1, 2) = 1.0d0
    mat_dense(2, 1) = 1.0d0
    mat_dense(2, 2) = 2.0d0
    mat_dense(2, 3) = 3.0d0
    mat_dense(3, 2) = 1.0d0
    mat_dense(3, 3) = 2.0d0
    mat_dense(3, 4) = 4.0d0
    mat_dense(4, 3) = 1.0d0
    mat_dense(4, 4) = 2.0d0
    mat_dense(4, 5) = 5.0d0
    mat_dense(5, 4) = 1.0d0
    mat_dense(5, 5) = 2.0d0

    a(1) = 1.0d0
    a(2) = 1.0d0
    a(3) = 1.0d0
    a(4) = 1.0d0
    a(5) = 1.0d0

    call monolis_dense_matvec_local_R(5, 5, mat_dense, a, b, time)

    b_th = matmul(mat_dense, a)

    call monolis_test_check_eq_R("monolis_dense_matvec_local_R", b, b_th)
  end subroutine monolis_dense_matvec_local_R_test

  subroutine monolis_dense_matmul_local_R_test()
    implicit none
    real(kdouble) :: a(5,3), b(5,3), b_th(5,3), mat_dense(5,5), time

    call monolis_std_global_log_string("monolis_dense_matmul_local_R")

    mat_dense(1, 1) = 2.0d0
    mat_dense(1, 2) = 1.0d0
    mat_dense(2, 1) = 1.0d0
    mat_dense(2, 2) = 2.0d0
    mat_dense(2, 3) = 3.0d0
    mat_dense(3, 2) = 1.0d0
    mat_dense(3, 3) = 2.0d0
    mat_dense(3, 4) = 4.0d0
    mat_dense(4, 3) = 1.0d0
    mat_dense(4, 4) = 2.0d0
    mat_dense(4, 5) = 5.0d0
    mat_dense(5, 4) = 1.0d0
    mat_dense(5, 5) = 2.0d0

    a(1,1) = 1.0d0; a(1,2) = 2.0d0; a(1,3) = 3.0d0
    a(2,1) = 1.0d0; a(2,2) = 2.0d0; a(2,3) = 3.0d0
    a(3,1) = 1.0d0; a(3,2) = 2.0d0; a(3,3) = 3.0d0
    a(4,1) = 1.0d0; a(4,2) = 2.0d0; a(4,3) = 3.0d0
    a(5,1) = 1.0d0; a(5,2) = 2.0d0; a(5,3) = 3.0d0

    call monolis_dense_matmul_local_R(5, 5, 3, mat_dense, a, b, time)

    b_th = matmul(mat_dense, a)

    call monolis_test_check_eq_R("monolis_dense_matmul_local_R 1", b(:,1), b_th(:,1))
    call monolis_test_check_eq_R("monolis_dense_matmul_local_R 2", b(:,2), b_th(:,2))
    call monolis_test_check_eq_R("monolis_dense_matmul_local_R 3", b(:,3), b_th(:,3))
  end subroutine monolis_dense_matmul_local_R_test
end module mod_monolis_matalg_dense_test
