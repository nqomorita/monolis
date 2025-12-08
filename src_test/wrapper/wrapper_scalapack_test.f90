!> ScaLAPACK ラッパーモジュール
module mod_monolis_scalapack_test
  use mod_monolis_utils
  use mod_monolis_scalapack

  implicit none

contains

  subroutine monolis_scalapack_test()
    implicit none
    call monolis_scalapack_test_1()
    call monolis_scalapack_test_2()
    call monolis_scalapack_test_3()
    call monolis_scalapack_test_4()
    call monolis_scalapack_test_5()
    call monolis_scalapack_test_6()
  end subroutine monolis_scalapack_test

  subroutine monolis_scalapack_test_1()
    implicit none
    !> 行列の大きさ（行数 N_loc）
    integer(kint) :: N_loc = 4
    !> 行列の大きさ（列数 M）
    integer(kint) :: M = 3
    !> 入力行列（N_loc x M）
    real(kdouble) :: A(4,3)
    !> 左特異行列（N_loc x P, P = min(N, M)）
    real(kdouble) :: S(4,3)
    !> 特異値（P）
    real(kdouble) :: V(3)
    !> 右特異行列（P x M）
    real(kdouble) :: D(3,3)
    integer(kint) :: comm
    integer(kint) :: scalapack_comm

    real(kdouble) :: Vt(3,3)
    real(kdouble) :: VD(3,3)
    real(kdouble) :: A_res(4,3)

    call monolis_std_global_log_string("monolis_scalapack_gesvd_R")

    comm = monolis_mpi_get_global_comm()

    A = 0.0d0
    S = 0.0d0
    V = 0.0d0
    D = 0.0d0

    if(monolis_mpi_get_global_my_rank() == 0)then
      A(1,1) = 1.0d0
      A(2,1) = 2.0d0
      A(3,1) = 3.0d0
      A(4,1) = 4.0d0

      A(1,2) = 11.0d0
      A(2,2) = 12.0d0
      A(3,2) = 13.0d0
      A(4,2) = 14.0d0

      A(1,3) = 21.0d0
      A(2,3) = 22.0d0
      A(3,3) = 23.0d0
      A(4,3) = 24.0d0
    else
      A(1,1) = 5.0d0
      A(2,1) = 6.0d0
      A(3,1) = 7.0d0
      A(4,1) = 8.0d0

      A(1,2) = 15.0d0
      A(2,2) = 16.0d0
      A(3,2) = 17.0d0
      A(4,2) = 18.0d0

      A(1,3) = 25.0d0
      A(2,3) = 26.0d0
      A(3,3) = 27.0d0
      A(4,3) = 28.0d0
    endif

    call monolis_scalapack_comm_initialize(comm, scalapack_comm)
    call monolis_scalapack_gesvd_R(N_loc, M, A, S, V, D, comm, scalapack_comm)
    call monolis_scalapack_comm_finalize(scalapack_comm)

    Vt = 0.0d0
    Vt(1,1) = V(1)
    Vt(2,2) = V(2)
    Vt(3,3) = V(3)

    VD = matmul(Vt, D)
    A_res = matmul(S,VD)

    call monolis_test_check_eq_R("monolis_scalapack_gesvd_R 1-1", A(:,1), A_res(:,1))
    call monolis_test_check_eq_R("monolis_scalapack_gesvd_R 1-2", A(:,2), A_res(:,2))
    call monolis_test_check_eq_R("monolis_scalapack_gesvd_R 1-3", A(:,3), A_res(:,3))
  end subroutine monolis_scalapack_test_1

  subroutine monolis_scalapack_test_2()
    implicit none
    !> 行列の大きさ（行数 N）
    integer(kint) :: N_loc = 2
    !> 行列の大きさ（列数 M）
    integer(kint) :: M = 6
    !> 入力行列（N_loc x M）
    real(kdouble) :: A(2,6)
    !> 左特異行列（N_loc x P, P = min(N, M)）
    real(kdouble) :: S(2,4)
    !> 特異値（P）
    real(kdouble) :: V(4)
    !> 右特異行列（P x M）
    real(kdouble) :: D(4,6)
    integer(kint) :: comm
    integer(kint) :: scalapack_comm

    real(kdouble) :: Vt(4,4)
    real(kdouble) :: VD(4,6)
    real(kdouble) :: A_res(2,6)

    call monolis_std_global_log_string("monolis_scalapack_gesvd_R")

    comm = monolis_mpi_get_global_comm()

    A = 0.0d0
    S = 0.0d0
    V = 0.0d0
    D = 0.0d0

    if(monolis_mpi_get_global_my_rank() == 0)then
      A(1,1) = 1.0d0
      A(2,1) = 2.0d0

      A(1,2) = 13.0d0
      A(2,2) = 14.0d0

      A(1,3) = 25.0d0
      A(2,3) = 26.0d0

      A(1,4) = 37.0d0
      A(2,4) = 38.0d0

      A(1,5) = 49.0d0
      A(2,5) = 50.0d0

      A(1,6) = 51.0d0
      A(2,6) = 52.0d0
    else
      A(1,1) = 3.0d0
      A(2,1) = 4.0d0

      A(1,2) = 15.0d0
      A(2,2) = 16.0d0

      A(1,3) = 27.0d0
      A(2,3) = 28.0d0

      A(1,4) = 39.0d0
      A(2,4) = 40.0d0

      A(1,5) = 41.0d0
      A(2,5) = 42.0d0

      A(1,6) = 53.0d0
      A(2,6) = 54.0d0
    endif

    call monolis_scalapack_comm_initialize(comm, scalapack_comm)
    call monolis_scalapack_gesvd_R(N_loc, M, A, S, V, D, comm, scalapack_comm)
    call monolis_scalapack_comm_finalize(scalapack_comm)

    if(monolis_mpi_get_global_comm_size() == 1)then
      Vt = 0.0d0
      Vt(1,1) = V(1)
      Vt(2,2) = V(2)

      VD(1:2,1:6) = matmul(Vt(1:2,1:2), D(1:2,1:6))
      A_res = matmul(S(1:2,1:2),VD(1:2,1:6))
    else
      Vt = 0.0d0
      Vt(1,1) = V(1)
      Vt(2,2) = V(2)
      Vt(3,3) = V(3)
      Vt(4,4) = V(4)

      VD = matmul(Vt, D)
      A_res = matmul(S,VD)
    endif

    call monolis_test_check_eq_R("monolis_scalapack_gesvd_R 2-1", A(:,1), A_res(:,1))
    call monolis_test_check_eq_R("monolis_scalapack_gesvd_R 2-2", A(:,2), A_res(:,2))
    call monolis_test_check_eq_R("monolis_scalapack_gesvd_R 2-3", A(:,3), A_res(:,3))
    call monolis_test_check_eq_R("monolis_scalapack_gesvd_R 2-4", A(:,4), A_res(:,4))
    call monolis_test_check_eq_R("monolis_scalapack_gesvd_R 2-5", A(:,5), A_res(:,5))
    call monolis_test_check_eq_R("monolis_scalapack_gesvd_R 2-6", A(:,6), A_res(:,6))
  end subroutine monolis_scalapack_test_2

 subroutine monolis_scalapack_test_3()
    implicit none
    !> 行列の大きさ（行数 N）
    integer(kint) :: N_loc = 2
    !> 行列の大きさ（列数 M）
    integer(kint) :: M = 6
    !> 入力行列（N_loc x M）
    real(kdouble) :: A(2,6)
    !> 左特異行列（N_loc x P, P = min(N, M)）
    real(kdouble) :: S(2,2)
    !> 特異値（P）
    real(kdouble) :: V(2)
    !> 右特異行列（P x M）
    real(kdouble) :: D(2,6)
    integer(kint) :: comm
    integer(kint) :: scalapack_comm

    real(kdouble) :: Vt(2,2)
    real(kdouble) :: VD(2,6)
    real(kdouble) :: A_res(2,6)

    call monolis_std_global_log_string("monolis_scalapack_gesvd_R")

    comm = monolis_mpi_get_self_comm()

    A = 0.0d0
    S = 0.0d0
    V = 0.0d0
    D = 0.0d0

    if(monolis_mpi_get_global_my_rank() == 0)then
      A(1,1) = 1.0d0
      A(2,1) = 2.0d0

      A(1,2) = 13.0d0
      A(2,2) = 14.0d0

      A(1,3) = 25.0d0
      A(2,3) = 26.0d0

      A(1,4) = 37.0d0
      A(2,4) = 38.0d0

      A(1,5) = 49.0d0
      A(2,5) = 50.0d0

      A(1,6) = 51.0d0
      A(2,6) = 52.0d0
    else
      A(1,1) = 3.0d0
      A(2,1) = 4.0d0

      A(1,2) = 15.0d0
      A(2,2) = 16.0d0

      A(1,3) = 27.0d0
      A(2,3) = 28.0d0

      A(1,4) = 39.0d0
      A(2,4) = 40.0d0

      A(1,5) = 41.0d0
      A(2,5) = 42.0d0

      A(1,6) = 53.0d0
      A(2,6) = 54.0d0
    endif

    call monolis_scalapack_comm_initialize(comm, scalapack_comm)
    call monolis_scalapack_gesvd_R(N_loc, M, A, S, V, D, comm, scalapack_comm)
    call monolis_scalapack_comm_finalize(scalapack_comm)

    Vt = 0.0d0
    Vt(1,1) = V(1)
    Vt(2,2) = V(2)

    VD = matmul(Vt, D)
    A_res = matmul(S,VD)

    call monolis_test_check_eq_R("monolis_scalapack_gesvd_R 3-1", A(:,1), A_res(:,1))
    call monolis_test_check_eq_R("monolis_scalapack_gesvd_R 3-2", A(:,2), A_res(:,2))
    call monolis_test_check_eq_R("monolis_scalapack_gesvd_R 3-3", A(:,3), A_res(:,3))
    call monolis_test_check_eq_R("monolis_scalapack_gesvd_R 3-4", A(:,4), A_res(:,4))
    call monolis_test_check_eq_R("monolis_scalapack_gesvd_R 3-5", A(:,5), A_res(:,5))
    call monolis_test_check_eq_R("monolis_scalapack_gesvd_R 3-6", A(:,6), A_res(:,6))
  end subroutine monolis_scalapack_test_3

 subroutine monolis_scalapack_test_4()
    implicit none
    !> 行列の大きさ（行数 N）
    integer(kint) :: N_loc = 2
    !> 行列の大きさ（列数 M）
    integer(kint) :: M = 6
    !> 入力行列（N_loc x M）
    real(kdouble) :: A(2,6)
    !> 左特異行列（N_loc x P, P = min(N, M)）
    real(kdouble) :: S(2,2)
    !> 特異値（P）
    real(kdouble) :: V(2)
    !> 右特異行列（P x M）
    real(kdouble) :: D(2,6)
    integer(kint) :: comm
    integer(kint) :: scalapack_comm

    real(kdouble) :: Vt(2,2)
    real(kdouble) :: VD(2,6)
    real(kdouble) :: A_res(2,6)

    call monolis_std_global_log_string("monolis_scalapack_gesvd_R section 4")

    comm = monolis_mpi_get_self_comm()

    A = 0.0d0
    S = 0.0d0
    V = 0.0d0
    D = 0.0d0

    call monolis_scalapack_comm_initialize(comm, scalapack_comm)

    if(monolis_mpi_get_global_my_rank() == 0)then
      A(1,1) = 1.0d0
      A(2,1) = 2.0d0

      A(1,2) = 13.0d0
      A(2,2) = 14.0d0

      A(1,3) = 25.0d0
      A(2,3) = 26.0d0

      A(1,4) = 37.0d0
      A(2,4) = 38.0d0

      A(1,5) = 49.0d0
      A(2,5) = 50.0d0

      A(1,6) = 51.0d0
      A(2,6) = 52.0d0

      call monolis_scalapack_gesvd_R(N_loc, M, A, S, V, D, comm, scalapack_comm)

      Vt = 0.0d0
      Vt(1,1) = V(1)
      Vt(2,2) = V(2)

      VD = matmul(Vt, D)
      A_res = matmul(S,VD)

      call monolis_test_check_eq_R("monolis_scalapack_gesvd_R 4-1", A(:,1), A_res(:,1))
      call monolis_test_check_eq_R("monolis_scalapack_gesvd_R 4-2", A(:,2), A_res(:,2))
      call monolis_test_check_eq_R("monolis_scalapack_gesvd_R 4-3", A(:,3), A_res(:,3))
      call monolis_test_check_eq_R("monolis_scalapack_gesvd_R 4-4", A(:,4), A_res(:,4))
      call monolis_test_check_eq_R("monolis_scalapack_gesvd_R 4-5", A(:,5), A_res(:,5))
      call monolis_test_check_eq_R("monolis_scalapack_gesvd_R 4-6", A(:,6), A_res(:,6))
    endif

    call monolis_scalapack_comm_finalize(scalapack_comm)
  end subroutine monolis_scalapack_test_4

  subroutine monolis_scalapack_test_5()
    implicit none
    !> 行列の大きさ（ローカル行数）
    integer(kint) :: N_loc = 2
    !> 行列の大きさ（全体のサイズ N x N）
    integer(kint) :: N = 4
    !> 入力行列（N_loc x N）
    real(kdouble) :: A(2,4)
    !> ピボット情報（N_loc）
    integer(kint) :: ipiv(2)
    !> 右辺ベクトル（N_loc x 1）
    real(kdouble) :: B(2,1)
    !> 解ベクトル（N_loc x 1）
    real(kdouble) :: X(2,1), X_ref(2,1)
    !> 検算用
    real(kdouble) :: AX(2,1)
    integer(kint) :: comm
    integer(kint) :: scalapack_comm

    call monolis_std_global_log_string("monolis_scalapack_getrf_R/getrs_R")

    comm = monolis_mpi_get_global_comm()

    A = 0.0d0
    B = 0.0d0
    ipiv = 0

    if(monolis_mpi_get_global_my_rank() == 0)then
      ! 行列 A の設定 (正定値対称行列)
      A(1,1) = 4.0d0
      A(1,2) = 1.0d0
      A(1,3) = 2.0d0
      A(1,4) = 1.0d0

      A(2,1) = 1.0d0
      A(2,2) = 3.0d0
      A(2,3) = 1.0d0
      A(2,4) = 2.0d0

      ! 右辺ベクトル B の設定
      B(1,1) = 1.0d0
      B(2,1) = 2.0d0
    else
      A(1,1) = 2.0d0
      A(1,2) = 1.0d0
      A(1,3) = 5.0d0
      A(1,4) = 1.0d0

      A(2,1) = 1.0d0
      A(2,2) = 2.0d0
      A(2,3) = 1.0d0
      A(2,4) = 4.0d0

      B(1,1) = 3.0d0
      B(2,1) = 4.0d0
    endif

    call monolis_scalapack_comm_initialize(comm, scalapack_comm)

    ! LU分解
    call monolis_scalapack_getrf_R(N_loc, N, A, ipiv, scalapack_comm)

    ! 線形方程式の求解
    X = B
    call monolis_scalapack_getrs_R(N_loc, N, 1, A, ipiv, X, comm, scalapack_comm)

    call monolis_scalapack_comm_finalize(scalapack_comm)

    if(monolis_mpi_get_global_my_rank() == 0)then
      X_ref(1,1) = -0.23008849557522
      X_ref(2,1) = -0.070796460176991
    else
      X_ref(1,1) = 0.51327433628319
      X_ref(2,1) = 0.9646017699115
    endif

    call monolis_test_check_eq_R("monolis_scalapack_getrf_R/getrs_R 5", X(:,1), X_ref(:,1))
  end subroutine monolis_scalapack_test_5

  subroutine monolis_scalapack_test_6()
    implicit none
    !> 行列の大きさ（ローカル行数）
    integer(kint) :: N_loc = 2
    !> 行列の大きさ（全体のサイズ N x N）
    integer(kint) :: N = 4
    !> 入力行列（N_loc x N）
    real(kdouble) :: A(2,4)
    !> ピボット情報（N_loc）
    integer(kint) :: ipiv(2)
    !> 右辺ベクトル（N_loc x 2）、2つの右辺
    real(kdouble) :: B(2,2)
    !> 解ベクトル（N_loc x 2）
    real(kdouble) :: X(2,2), X_ref(2,2)
    !> 検算用
    real(kdouble) :: AX(2,2)
    integer(kint) :: comm
    integer(kint) :: scalapack_comm

    call monolis_std_global_log_string("monolis_scalapack_getrf_R/getrs_R (self_comm)")

    comm = monolis_mpi_get_self_comm()

    A = 0.0d0
    B = 0.0d0
    ipiv = 0

    if(monolis_mpi_get_global_my_rank() == 0)then
      ! 行列 A の設定
      A(1,1) = 4.0d0
      A(1,2) = 1.0d0
      A(1,3) = 2.0d0
      A(1,4) = 1.0d0

      A(2,1) = 1.0d0
      A(2,2) = 3.0d0
      A(2,3) = 1.0d0
      A(2,4) = 2.0d0

      B(1,1) = 1.0d0
      B(2,1) = 2.0d0
      B(1,2) = 3.0d0
      B(2,2) = 4.0d0
    else
      A(1,1) = 2.0d0
      A(1,2) = 1.0d0
      A(1,3) = 5.0d0
      A(1,4) = 1.0d0

      A(2,1) = 1.0d0
      A(2,2) = 2.0d0
      A(2,3) = 1.0d0
      A(2,4) = 4.0d0

      B(1,1) = 5.0d0
      B(2,1) = 6.0d0
      B(1,2) = 7.0d0
      B(2,2) = 8.0d0
    endif

    call monolis_scalapack_comm_initialize(comm, scalapack_comm)

    ! LU分解
    call monolis_scalapack_getrf_R(N_loc, N, A, ipiv, scalapack_comm)

    ! 線形方程式の求解（2つの右辺）
    X = B
    call monolis_scalapack_getrs_R(N_loc, N, 2, A, ipiv, X, comm, scalapack_comm)

    call monolis_scalapack_comm_finalize(scalapack_comm)

    if(monolis_mpi_get_global_my_rank() == 0)then
      X_ref(1,1) = -0.23008849557522
      X_ref(2,1) = -0.070796460176991

      X_ref(1,2) = -0.23008849557522
      X_ref(2,2) = -0.070796460176991
    else
      X_ref(1,1) = 0.51327433628319
      X_ref(2,1) = 0.9646017699115

      X_ref(1,2) = 0.51327433628319
      X_ref(2,2) = 0.9646017699115
    endif

    call monolis_test_check_eq_R("monolis_scalapack_getrf_R/getrs_R 6-1", X(:,1), X_ref(:,1))
    call monolis_test_check_eq_R("monolis_scalapack_getrf_R/getrs_R 6-2", X(:,2), X_ref(:,2))
  end subroutine monolis_scalapack_test_6
end module mod_monolis_scalapack_test
