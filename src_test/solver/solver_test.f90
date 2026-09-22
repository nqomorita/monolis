!> 線形ソルバテストモジュール
module mod_monolis_solve_test
  use mod_monolis

  implicit none

  !> GPU 動作確認テストで走査する反復解法数
  integer(kint), parameter :: monolis_solve_test_n_method = 9
  !> GPU 動作確認テストにおける 1 ランクあたりの内部計算点数
  integer(kint), parameter :: monolis_solve_test_n_internal = 4

contains

  subroutine monolis_solve_test()
    implicit none

    call monolis_std_global_log_string("monolis_solve_R")
    call monolis_std_global_log_string("monolis_solve_C")
    call monolis_std_global_log_string("monolis_solve_main_R")
    call monolis_std_global_log_string("monolis_solve_main_C")
    call monolis_std_global_log_string("monolis_solver_select_R")
    call monolis_std_global_log_string("monolis_solver_select_C")

    call monolis_solve_gpu_test()
  end subroutine monolis_solve_test

  !> GPU 動作確認テスト
  !> 構成 (b) シミュレータ CPU + monolis GPU（ホスト同期 ON）と
  !> 構成 (c) シミュレータ GPU + monolis GPU（ホスト同期 OFF）を検証する。
  !> CPU ビルドでは OpenACC 指示文がコメントとなり、袖領域の整合性テストとして機能する
  subroutine monolis_solve_gpu_test()
    implicit none
    integer(kint) :: i, n_dof
    integer(kint) :: method(monolis_solve_test_n_method)

    method(1) = monolis_iter_CG
    method(2) = monolis_iter_GropCG
    method(3) = monolis_iter_PipeCG
    method(4) = monolis_iter_PipeCR
    method(5) = monolis_iter_BiCGSTAB
    method(6) = monolis_iter_PipeBiCGSTAB
    method(7) = monolis_iter_BiCGSTAB_noprec
    method(8) = monolis_iter_PipeBiCGSTAB_noprec
    method(9) = monolis_iter_BiCGSAFE

    do n_dof = 1, 3
      do i = 1, monolis_solve_test_n_method
        call monolis_solve_host_sync_test_main(n_dof, method(i), monolis_prec_NONE)
        call monolis_solve_host_sync_test_main(n_dof, method(i), monolis_prec_DIAG)
        call monolis_solve_host_sync_test_main(n_dof, method(i), monolis_prec_SOR)
      enddo
    enddo

    do i = 1, monolis_solve_test_n_method
      call monolis_solve_device_resident_test_main(1, method(i), monolis_prec_DIAG)
      call monolis_solve_device_resident_test_main(3, method(i), monolis_prec_DIAG)
    enddo
  end subroutine monolis_solve_gpu_test

  !> 構成 (b)：ホスト同期 ON（既定）の検証
  !> 求解後のホスト側解ベクトルが袖（共有計算点）を含めて正しいことを確認する
  subroutine monolis_solve_host_sync_test_main(n_dof, method, prec)
    implicit none
    !> [in] 計算点が持つ自由度
    integer(kint), intent(in) :: n_dof
    !> [in] 反復解法
    integer(kint), intent(in) :: method
    !> [in] 前処理
    integer(kint), intent(in) :: prec
    type(monolis_structure) :: mat
    type(monolis_com) :: com
    integer(kint) :: NPNDOF
    real(kdouble), allocatable :: b(:), x(:), x_comm(:), x_ans(:)

    call monolis_std_log_I1("DOF", n_dof)
    call monolis_std_log_I1("METHOD", method)
    call monolis_std_log_I1("PRECOND", prec)

    call monolis_solve_test_get_mat(mat, com, n_dof, NPNDOF, b)

    call monolis_alloc_R_1d(x, NPNDOF)
    call monolis_alloc_R_1d(x_comm, NPNDOF)
    call monolis_alloc_R_1d(x_ans, NPNDOF)
    x_ans = 1.0d0

    call monolis_set_method(mat, method)
    call monolis_set_precond(mat, prec)
    call monolis_set_maxiter(mat, 1000)
    call monolis_set_tolerance(mat, 1.0d-10)
    call monolis_prm_enable_host_sync(mat, .true.)

    call monolis_solve_R(mat, com, b, x)

    !# 袖成分を含む全成分が厳密解に一致すること
    call monolis_test_check_eq_R("monolis_solve_gpu_test host sync solution", x, x_ans)

    !# 袖成分が所有ランクの値と一致すること（再通信で値が変化しない）
    x_comm = x
    call monolis_mpi_update_R(com, n_dof, x_comm)
    call monolis_test_check_eq_R("monolis_solve_gpu_test host sync ghost", x_comm, x)

    call monolis_dealloc_R_1d(b)
    call monolis_dealloc_R_1d(x)
    call monolis_dealloc_R_1d(x_comm)
    call monolis_dealloc_R_1d(x_ans)
    call monolis_finalize(mat)
  end subroutine monolis_solve_host_sync_test_main

  !> 構成 (c)：ホスト同期 OFF（デバイス常駐モード）の検証
  !> 呼出し側が解・右辺ベクトルのデバイス常駐を所有し、求解後もデバイス上の値が正しいことを確認する
  subroutine monolis_solve_device_resident_test_main(n_dof, method, prec)
    implicit none
    !> [in] 計算点が持つ自由度
    integer(kint), intent(in) :: n_dof
    !> [in] 反復解法
    integer(kint), intent(in) :: method
    !> [in] 前処理
    integer(kint), intent(in) :: prec
    type(monolis_structure) :: mat
    type(monolis_com) :: com
    integer(kint) :: NPNDOF, i
    real(kdouble), allocatable :: b(:), x(:), x_comm(:), x_ans(:)

    call monolis_std_log_I1("DOF", n_dof)
    call monolis_std_log_I1("METHOD", method)
    call monolis_std_log_I1("PRECOND", prec)

    call monolis_solve_test_get_mat(mat, com, n_dof, NPNDOF, b)

    call monolis_alloc_R_1d(x, NPNDOF)
    call monolis_alloc_R_1d(x_comm, NPNDOF)
    call monolis_alloc_R_1d(x_ans, NPNDOF)
    x_ans = 1.0d0

    call monolis_set_method(mat, method)
    call monolis_set_precond(mat, prec)
    call monolis_set_maxiter(mat, 1000)
    call monolis_set_tolerance(mat, 1.0d-10)
    call monolis_prm_enable_host_sync(mat, .false.)

    do i = 1, NPNDOF
      mat%MAT%R%B(i) = b(i)
      mat%MAT%R%X(i) = 0.0d0
    enddo

    !# シミュレータ側が右辺・解ベクトルのデバイス常駐を所有する
    !$acc enter data copyin(mat%MAT%R%X, mat%MAT%R%B)

    call monolis_solve_main_R(mat%PRM, com, mat%MAT, mat%PREC)

    !# ホスト同期 OFF ではデバイス上の値が正。ホストへの取り出しは呼出し側の責任で行う
    !$acc update self(mat%MAT%R%X)
    !$acc exit data delete(mat%MAT%R%X, mat%MAT%R%B)

    do i = 1, NPNDOF
      x(i) = mat%MAT%R%X(i)
    enddo

    call monolis_test_check_eq_R("monolis_solve_gpu_test device resident solution", x, x_ans)

    x_comm = x
    call monolis_mpi_update_R(com, n_dof, x_comm)
    call monolis_test_check_eq_R("monolis_solve_gpu_test device resident ghost", x_comm, x)

    call monolis_dealloc_R_1d(b)
    call monolis_dealloc_R_1d(x)
    call monolis_dealloc_R_1d(x_comm)
    call monolis_dealloc_R_1d(x_ans)
    call monolis_finalize(mat)
  end subroutine monolis_solve_device_resident_test_main

  !> GPU 動作確認テスト用の分割行列・右辺ベクトルの生成
  !> 各ランクが内部計算点 4 個を持つ 1 次元チェーンを構成し、隣接ランクの端点を袖として保持する。
  !> 厳密解 x = 1 に対応する右辺を返す
  subroutine monolis_solve_test_get_mat(mat, com, n_dof, NPNDOF, b)
    implicit none
    !> [out] monolis 構造体
    type(monolis_structure), intent(inout) :: mat
    !> [out] 通信テーブル構造体
    type(monolis_com), intent(inout) :: com
    !> [in] 計算点が持つ自由度
    integer(kint), intent(in) :: n_dof
    !> [out] 袖を含むベクトルサイズ
    integer(kint), intent(out) :: NPNDOF
    !> [out] 右辺ベクトル
    real(kdouble), allocatable, intent(out) :: b(:)
    integer(kint) :: n_internal, n_node, n_elem, comm_size, my_rank
    integer(kint) :: i, j, k, id_l, id_r
    integer(kint), allocatable :: global_id(:), elem(:,:)
    real(kdouble), allocatable :: a(:)

    comm_size = monolis_mpi_get_global_comm_size()
    my_rank = monolis_mpi_get_global_my_rank()
    n_internal = monolis_solve_test_n_internal

    id_l = 0
    id_r = 0
    n_node = n_internal
    if(my_rank > 0)then
      n_node = n_node + 1
      id_l = n_node
    endif
    if(my_rank < comm_size - 1)then
      n_node = n_node + 1
      id_r = n_node
    endif

    call monolis_alloc_I_1d(global_id, n_node)
    do i = 1, n_internal
      global_id(i) = n_internal*my_rank + i
    enddo
    if(id_l > 0) global_id(id_l) = n_internal*my_rank
    if(id_r > 0) global_id(id_r) = n_internal*(my_rank + 1) + 1

    n_elem = n_internal - 1
    if(id_l > 0) n_elem = n_elem + 1
    if(id_r > 0) n_elem = n_elem + 1

    call monolis_alloc_I_2d(elem, 2, n_elem)
    do i = 1, n_internal - 1
      elem(1,i) = i
      elem(2,i) = i + 1
    enddo
    k = n_internal - 1
    if(id_l > 0)then
      k = k + 1
      elem(1,k) = id_l
      elem(2,k) = 1
    endif
    if(id_r > 0)then
      k = k + 1
      elem(1,k) = n_internal
      elem(2,k) = id_r
    endif

    call monolis_com_initialize_by_global_id(com, monolis_mpi_get_global_comm(), &
      & n_internal, n_node, global_id)

    call monolis_initialize(mat)
    call monolis_get_nonzero_pattern_by_simple_mesh_R(mat, n_node, 2, n_dof, n_elem, elem)

    !# 対角優位な対称行列（対角 10、隣接 -1）
    do i = 1, n_node
      do j = 1, n_dof
        call monolis_add_scalar_to_sparse_matrix_R(mat, i, i, j, j, 1.0d1)
      enddo
    enddo
    do i = 1, n_elem
      do j = 1, n_dof
        call monolis_add_scalar_to_sparse_matrix_R(mat, elem(1,i), elem(2,i), j, j, -1.0d0)
        call monolis_add_scalar_to_sparse_matrix_R(mat, elem(2,i), elem(1,i), j, j, -1.0d0)
      enddo
    enddo

    if(comm_size > 1) mat%MAT%N = com%n_internal_vertex

    NPNDOF = n_dof*n_node
    call monolis_alloc_R_1d(a, NPNDOF)
    call monolis_alloc_R_1d(b, NPNDOF)
    a = 1.0d0

    call monolis_matvec_product_R(mat, com, a, b)

    call monolis_dealloc_R_1d(a)
    call monolis_dealloc_I_1d(global_id)
    call monolis_dealloc_I_2d(elem)
  end subroutine monolis_solve_test_get_mat

end module mod_monolis_solve_test
