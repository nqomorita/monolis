!> 疎行列ベクトル積テストモジュール
module mod_monolis_matalg_dense_test
  use mod_monolis
  use mod_monolis_utils
  implicit none

contains

  subroutine monolis_matalg_dense_test()
    implicit none

    if(monolis_mpi_get_global_comm_size() == 2) return

    call monolis_dense_matvec_local_R_test()

    call monolis_std_global_log_string("monolis_dense_matvec_local_R")
    call monolis_std_global_log_string("monolis_dense_matmul_local_R")
  end subroutine monolis_matalg_dense_test

  subroutine monolis_dense_matvec_local_R_test()
    implicit none
    type(monolis_structure) :: mat
    type(monolis_com) :: com
    integer(kint) :: n_node, n_elem, elem(2,4)
    real(kdouble) :: a(5), b(5), b_th(5), mat_dense(5,5)

    call monolis_std_global_log_string("monolis_dense_matvec_local_R")

    !call monolis_dense_matvec_local_R(N, M, MAT, X, Y, tdemv)

    !call monolis_test_check_eq_R("monolis_dense_matvec_local_R", b, b_th)

    call monolis_finalize(mat)
  end subroutine monolis_dense_matvec_local_R_test

  subroutine monolis_dense_matmul_local_R_test()
    implicit none
    type(monolis_structure) :: mat
    type(monolis_com) :: com
    integer(kint) :: n_node, n_elem, elem(2,4)
    real(kdouble) :: a(5), b(5), b_th(5), mat_dense(5,5)

    call monolis_std_global_log_string("monolis_dense_matmul_local_R")

    !call monolis_dense_matmul_local_R(N1, N2, N3, MAT1, MAT2, Y, tdemv)

    !call monolis_test_check_eq_R("monolis_dense_matmul_local_R", b, b_th)

    call monolis_finalize(mat)
  end subroutine monolis_dense_matmul_local_R_test

end module mod_monolis_matalg_dense_test
