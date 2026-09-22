!> IDR(s) 法テストモジュール
module mod_monolis_solver_IDRS_test
  use mod_monolis
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

  implicit none

contains

  subroutine monolis_solver_IDRS_test()
    implicit none
    integer(kint) :: n_dof

    do n_dof = 1, 3
      call monolis_solver_IDRS_test_main(n_dof, monolis_prec_NONE)
      call monolis_solver_IDRS_test_main(n_dof, monolis_prec_DIAG)
      call monolis_solver_IDRS_test_main(n_dof, monolis_prec_SOR)
      call monolis_solver_IDRS_test_main(n_dof, monolis_prec_LU)
    enddo
  end subroutine monolis_solver_IDRS_test

  subroutine monolis_solver_IDRS_test_main(n_dof, prec)
    implicit none
    type(monolis_structure) :: mat
    type(monolis_com) :: com
    integer(kint) :: nnode, nelem, elem(2,9)
    integer(kint) :: i1, i2, j2
    integer(kint) :: n_dof, prec
    real(kdouble) :: val
    real(kdouble) :: a(n_dof*10), b(n_dof*10)

    call monolis_std_global_log_string("monolis_solver_IDRS")
    call monolis_std_log_I1("DOF", n_dof)
    call monolis_std_log_I1("PRECOND", prec)

    call monolis_initialize(mat)
    call monolis_com_initialize_by_self(com)

    nnode = 10

    nelem = 9

    elem(1,1) = 1; elem(2,1) = 2;
    elem(1,2) = 2; elem(2,2) = 3;
    elem(1,3) = 3; elem(2,3) = 4;
    elem(1,4) = 4; elem(2,4) = 5;
    elem(1,5) = 5; elem(2,5) = 6;
    elem(1,6) = 6; elem(2,6) = 7;
    elem(1,7) = 7; elem(2,7) = 8;
    elem(1,8) = 8; elem(2,8) = 9;
    elem(1,9) = 9; elem(2,9) =10;

    call monolis_get_nonzero_pattern_by_simple_mesh_R(mat, nnode, 2, n_dof, nelem, elem)

    do i1 = 1, 10
      do i2 = 1, n_dof
        call random_number(val)
        val = val + 2.0d0*n_dof
        call monolis_add_scalar_to_sparse_matrix_R(mat, i1, i1, i2, i2, val)
      enddo
    enddo

    do i1 = 1, 9
      do i2 = 1, n_dof
      do j2 = 1, n_dof
        call random_number(val)
        call monolis_add_scalar_to_sparse_matrix_R(mat, elem(1,i1), elem(2,i1), i2, j2, val)
        call random_number(val)
        call monolis_add_scalar_to_sparse_matrix_R(mat, elem(2,i1), elem(1,i1), j2, i2, val)
      enddo
      enddo
    enddo

    a = 1.0d0

    call monolis_matvec_product_R(mat, com, a, b)

    call monolis_set_method(mat, monolis_iter_IDRS)
    call monolis_set_precond(mat, prec)
    call monolis_set_tolerance(mat, 1.0d-12)
    call monolis_show_timelog_statistics(mat, .true.)

    a = 0.0d0

    call monolis_solve_R(mat, com, b, a)

    b = 1.0d0

    call monolis_test_check_eq_R("monolis_solver_IDRS_test_main", a, b)

    call monolis_finalize(mat)
  end subroutine monolis_solver_IDRS_test_main

  !> 非一定解に対する IDR(s) 法の初期値・乱数依存性の回帰テスト
  !> 1・2 ランクでは C 並列テストと同じ 10 計算点の行列を使用する
  subroutine monolis_solver_IDRS_regression_test()
    implicit none
    integer(kint) :: n_dof, n_seed
    integer(kint), allocatable :: saved_seed(:)

    call monolis_std_global_log_string("monolis_solver_IDRS")
    call random_seed(size = n_seed)
    call monolis_alloc_I_1d(saved_seed, n_seed)
    call random_seed(get = saved_seed)

    do n_dof = 1, 3
      call monolis_solver_IDRS_regression_main(n_dof, monolis_prec_NONE, .false.)
      call monolis_solver_IDRS_regression_main(n_dof, monolis_prec_DIAG, .false.)
      call monolis_solver_IDRS_regression_main(n_dof, monolis_prec_SOR, .false.)
      !# 対角前処理後に単位行列となる系でも breakdown せずに収束すること
      call monolis_solver_IDRS_regression_main(n_dof, monolis_prec_DIAG, .true.)
    enddo
    call monolis_solver_IDRS_regression_small_test()

    call random_seed(put = saved_seed)
    call monolis_dealloc_I_1d(saved_seed)
  end subroutine monolis_solver_IDRS_regression_test

  !> 固定した行列に対して乱数 seed と初期解を変え、真の残差も確認する
  subroutine monolis_solver_IDRS_regression_main(n_dof, prec, is_diagonal)
    implicit none
    !> [in] 計算点が持つ自由度
    integer(kint), intent(in) :: n_dof
    !> [in] 前処理
    integer(kint), intent(in) :: prec
    !> [in] 単位行列相当の前処理を検証する場合は真
    logical, intent(in) :: is_diagonal
    type(monolis_structure) :: mat
    type(monolis_com) :: com
    integer(kint), parameter :: seed_id(8) = [1, 2, 3, 7, 17, 53, 101, 509]
    integer(kint), parameter :: n_basis(8) = [1, 2, 4, 1, 2, 4, 2, 4]
    integer(kint) :: i, i_seed, i_init, n_seed, nn_dof, iter_conv
    integer(kint), allocatable :: seed(:)
    real(kdouble) :: norms(2), residual
    real(kdouble), allocatable :: x(:), x_ans(:), b(:), ax(:)
    character(monolis_charlen) :: header

    call monolis_solver_IDRS_regression_get_mat(mat, com, n_dof, is_diagonal, x_ans)
    nn_dof = n_dof*mat%MAT%N
    call monolis_alloc_R_1d(x, size(x_ans))
    call monolis_alloc_R_1d(b, size(x_ans))
    call monolis_alloc_R_1d(ax, size(x_ans))
    call monolis_matvec_product_R(mat, com, x_ans, b)

    call monolis_set_method(mat, monolis_iter_IDRS)
    call monolis_set_precond(mat, prec)
    call monolis_set_maxiter(mat, 1000)
    call monolis_set_tolerance(mat, 1.0d-12)
    call monolis_show_timelog(mat, .false.)
    call monolis_show_iterlog(mat, .false.)
    call monolis_show_summary(mat, .false.)

    call random_seed(size = n_seed)
    call monolis_alloc_I_1d(seed, n_seed)
    do i_seed = 1, size(seed_id)
      call monolis_set_solver_IDRS_num_basis(mat, n_basis(i_seed))
      !# 残差置換後の再始動も、短い周期で明示的に検証する
      call monolis_set_iter_RR(mat, 50)
      if(i_seed == 7) call monolis_set_iter_RR(mat, 1)
      do i = 1, n_seed
        seed(i) = 104729*seed_id(i_seed) + 7919*i + 1049*com%my_rank
      enddo
      do i_init = 0, 2
        write(header, '(a,4(a,i0),a,l1,a,i0)') "monolis_solver_IDRS_regression", &
          " DOF=", n_dof, " PREC=", prec, " SEED=", seed_id(i_seed), " S=", n_basis(i_seed), &
          " DIAGONAL=", is_diagonal, " INIT=", i_init
        call monolis_std_global_log_string(trim(header))
        call random_seed(put = seed)

        !# 既定のゼロ初期化、非ゼロ初期解、既に収束した初期解を区別する
        call monolis_set_init_x(mat, i_init == 0)
        select case(i_init)
        case(0)
          x = 0.0d0
        case(1)
          x = 0.25d0 + 0.1d0*x_ans
        case(2)
          x = x_ans
        end select

        call monolis_solve_R(mat, com, b, x)

        !# 等値比較だけでは NaN を検出できないため、有限性を先に検証する
        call monolis_test_check_eq_L1(trim(header)//" finite solution", all(ieee_is_finite(x)), .true.)
        call monolis_test_check_eq_R(trim(header)//" solution including ghosts", x, x_ans)

        !# ソルバ内部の再帰残差とは独立に b-Ax を再計算する
        call monolis_matvec_product_R(mat, com, x, ax)
        ax = b - ax
        norms(1) = dot_product(ax(1:nn_dof), ax(1:nn_dof))
        norms(2) = dot_product(b(1:nn_dof), b(1:nn_dof))
        call monolis_allreduce_R(2, norms, monolis_mpi_sum, com%comm)
        residual = sqrt(norms(1)/norms(2))
        call monolis_test_check_eq_L1(trim(header)//" finite true residual", ieee_is_finite(residual), .true.)
        call monolis_test_check_eq_L1(trim(header)//" true residual", residual <= 1.0d-11, .true.)

        call monolis_get_converge_residual(mat, residual)
        call monolis_test_check_eq_L1(trim(header)//" finite reported residual", ieee_is_finite(residual), .true.)
        call monolis_test_check_eq_L1(trim(header)//" reported residual", residual <= 1.0d-12, .true.)
        call monolis_get_converge_iter(mat, iter_conv)
        if(i_init == 2)then
          call monolis_test_check_eq_I1(trim(header)//" initially converged", iter_conv, 0)
        endif
      enddo
    enddo

    call monolis_dealloc_I_1d(seed)
    call monolis_dealloc_R_1d(x)
    call monolis_dealloc_R_1d(x_ans)
    call monolis_dealloc_R_1d(b)
    call monolis_dealloc_R_1d(ax)
    call monolis_finalize(mat)
    call monolis_com_finalize(com)
  end subroutine monolis_solver_IDRS_regression_main

  !> 真の残差と再帰残差が乖離する小規模系、および自由度数を超える基底数の検証
  subroutine monolis_solver_IDRS_regression_small_test()
    implicit none
    type(monolis_structure) :: mat
    type(monolis_com) :: com
    integer(kint), parameter :: n_nodes(3) = [4, 2, 1]
    integer(kint), parameter :: preconds(2) = [monolis_prec_NONE, monolis_prec_DIAG]
    integer(kint) :: i_case, n_node, n_elem, elem(2,3), i, j, i_prec, n_basis, n_seed, seed_id, iter_conv
    integer(kint), allocatable :: seed(:)
    real(kdouble) :: x(4), x_ans(4), b(4), ax(4), diagonal, residual
    real(kdouble) :: p(2), q(2), value
    character(monolis_charlen) :: header

    call random_seed(size = n_seed)
    call monolis_alloc_I_1d(seed, n_seed)
    do i_case = 1, size(n_nodes)
      n_node = n_nodes(i_case)
      n_elem = max(1, n_node - 1)
      do i = 1, n_elem
        elem(:,i) = [i, min(i + 1, n_node)]
      enddo
      call monolis_initialize(mat)
      call monolis_com_initialize_by_self(com)
      call monolis_get_nonzero_pattern_by_simple_mesh_R(mat, n_node, 2, 1, n_elem, elem(:,1:n_elem))

      !# n=s=4 では旧実装が真の相対残差約 0.04 に対して収束を報告した
      seed_id = 24
      x_ans = 1.0d0
      if(n_node <= 2)then
        !# 2x2 では旧実装が seed=784、s=4 で偽収束した
        seed_id = 784
        x_ans(2) = 1.05d0
      endif
      do i = 1, n_node
        diagonal = 2.0d0
        if(n_node == 4)then
          diagonal = 4.0d0
          if(i == 1 .or. i == n_node) diagonal = 5.0d0
        endif
        call monolis_add_scalar_to_sparse_matrix_R(mat, i, i, 1, 1, diagonal)
      enddo
      do i = 1, n_node - 1
        call monolis_add_scalar_to_sparse_matrix_R(mat, i, i + 1, 1, 1, 1.0d0)
        call monolis_add_scalar_to_sparse_matrix_R(mat, i + 1, i, 1, 1, 1.0d0)
      enddo
      call monolis_matvec_product_R(mat, com, x_ans(1:n_node), b(1:n_node))
      do i = 1, n_seed
        seed(i) = seed_id*104729 + i*7919
      enddo
      call monolis_set_method(mat, monolis_iter_IDRS)
      call monolis_set_tolerance(mat, 1.0d-12)
      call monolis_show_timelog(mat, .false.)
      call monolis_show_iterlog(mat, .false.)
      call monolis_show_summary(mat, .false.)

      do i_prec = 1, size(preconds)
        call monolis_set_precond(mat, preconds(i_prec))
        do n_basis = 2, 4, 2
          write(header, '(a,4(a,i0))') "monolis_solver_IDRS_regression_small", &
            " N=", n_node, " PREC=", preconds(i_prec), " S=", n_basis, " SEED=", seed_id
          call monolis_std_global_log_string(trim(header))
          call monolis_set_solver_IDRS_num_basis(mat, n_basis)
          call random_seed(put = seed)
          x = 0.0d0
          call monolis_solve_R(mat, com, b(1:n_node), x(1:n_node))

          call monolis_test_check_eq_L1(trim(header)//" finite solution", &
            all(ieee_is_finite(x(1:n_node))), .true.)
          call monolis_test_check_eq_R(trim(header)//" solution", x(1:n_node), x_ans(1:n_node))
          call monolis_matvec_product_R(mat, com, x(1:n_node), ax(1:n_node))
          ax(1:n_node) = b(1:n_node) - ax(1:n_node)
          residual = sqrt(dot_product(ax(1:n_node), ax(1:n_node))/dot_product(b(1:n_node), b(1:n_node)))
          call monolis_test_check_eq_L1(trim(header)//" finite true residual", ieee_is_finite(residual), .true.)
          call monolis_test_check_eq_L1(trim(header)//" true residual", residual <= 1.0d-11, .true.)
          call monolis_get_converge_residual(mat, residual)
          call monolis_test_check_eq_L1(trim(header)//" finite reported residual", ieee_is_finite(residual), .true.)
          call monolis_test_check_eq_L1(trim(header)//" reported residual", residual <= 1.0d-12, .true.)
        enddo
      enddo
      call monolis_finalize(mat)
      call monolis_com_finalize(com)
    enddo

    !# s=1 の shadow 方向を p とし、q をその直交方向とする。
    !# A=I+p*p^T、厳密解 p+q なら最初の更新後の残差は q/2 となり、
    !# omega=1 の更新で収束する。最大反復数 1 でも途中で打ち切らないこと。
    do i = 1, n_seed
      seed(i) = 104729 + i*7919
    enddo
    call random_seed(put = seed)
    call random_number(p)
    p = p - 0.5d0
    p = p/sqrt(dot_product(p, p))
    q = [-p(2), p(1)]
    x_ans(1:2) = p + q
    elem(:,1) = [1, 2]
    call monolis_initialize(mat)
    call monolis_com_initialize_by_self(com)
    call monolis_get_nonzero_pattern_by_simple_mesh_R(mat, 2, 2, 1, 1, elem(:,1:1))
    do i = 1, 2
      do j = 1, 2
        value = p(i)*p(j)
        if(i == j) value = value + 1.0d0
        call monolis_add_scalar_to_sparse_matrix_R(mat, i, j, 1, 1, value)
      enddo
    enddo
    call monolis_matvec_product_R(mat, com, x_ans(1:2), b(1:2))
    call monolis_set_method(mat, monolis_iter_IDRS)
    call monolis_set_precond(mat, monolis_prec_NONE)
    call monolis_set_solver_IDRS_num_basis(mat, 1)
    call monolis_set_maxiter(mat, 1)
    call monolis_set_tolerance(mat, 1.0d-12)
    call monolis_show_timelog(mat, .false.)
    call monolis_show_iterlog(mat, .false.)
    call monolis_show_summary(mat, .false.)
    call random_seed(put = seed)
    x = 0.0d0
    call monolis_solve_R(mat, com, b(1:2), x(1:2))

    header = "monolis_solver_IDRS_regression_small omega at maxiter"
    call monolis_test_check_eq_L1(trim(header)//" finite solution", all(ieee_is_finite(x(1:2))), .true.)
    call monolis_test_check_eq_R(trim(header)//" solution", x(1:2), x_ans(1:2))
    call monolis_matvec_product_R(mat, com, x(1:2), ax(1:2))
    ax(1:2) = b(1:2) - ax(1:2)
    residual = sqrt(dot_product(ax(1:2), ax(1:2))/dot_product(b(1:2), b(1:2)))
    call monolis_test_check_eq_L1(trim(header)//" finite true residual", ieee_is_finite(residual), .true.)
    call monolis_test_check_eq_L1(trim(header)//" true residual", residual <= 1.0d-12, .true.)
    call monolis_get_converge_iter(mat, iter_conv)
    call monolis_test_check_eq_I1(trim(header)//" iterations", iter_conv, 1)
    call monolis_get_converge_residual(mat, residual)
    call monolis_test_check_eq_L1(trim(header)//" finite reported residual", ieee_is_finite(residual), .true.)
    call monolis_test_check_eq_L1(trim(header)//" reported residual", residual <= 1.0d-12, .true.)

    call monolis_finalize(mat)
    call monolis_com_finalize(com)
    call monolis_dealloc_I_1d(seed)
  end subroutine monolis_solver_IDRS_regression_small_test

  !> C 並列テストの対角 4・隣接 1 の鎖を、ファイルに依存せず分割構築する
  subroutine monolis_solver_IDRS_regression_get_mat(mat, com, n_dof, is_diagonal, x_ans)
    implicit none
    !> [out] monolis 構造体
    type(monolis_structure), intent(inout) :: mat
    !> [out] 通信テーブル構造体
    type(monolis_com), intent(inout) :: com
    !> [in] 計算点が持つ自由度
    integer(kint), intent(in) :: n_dof
    !> [in] 対角行列のみを構築する場合は真
    logical, intent(in) :: is_diagonal
    !> [out] 袖を含む非一定の厳密解
    real(kdouble), allocatable, intent(out) :: x_ans(:)
    integer(kint) :: n_global, n_internal, n_node, n_elem, comm_size, my_rank
    integer(kint) :: first, last, id_l, id_r, i, j, k, global_dof
    integer(kint), allocatable :: global_id(:), elem(:,:)
    real(kdouble) :: diagonal

    comm_size = monolis_mpi_get_global_comm_size()
    my_rank = monolis_mpi_get_global_my_rank()
    n_global = max(10, 5*comm_size)
    first = n_global*my_rank/comm_size + 1
    last = n_global*(my_rank + 1)/comm_size
    n_internal = last - first + 1
    n_node = n_internal
    id_l = 0
    id_r = 0
    if(first > 1)then
      n_node = n_node + 1
      id_l = n_node
    endif
    if(last < n_global)then
      n_node = n_node + 1
      id_r = n_node
    endif

    call monolis_alloc_I_1d(global_id, n_node)
    do i = 1, n_internal
      global_id(i) = first + i - 1
    enddo
    if(id_l > 0) global_id(id_l) = first - 1
    if(id_r > 0) global_id(id_r) = last + 1

    n_elem = n_node - 1
    call monolis_alloc_I_2d(elem, 2, n_elem)
    do i = 1, n_internal - 1
      elem(:,i) = [i, i + 1]
    enddo
    k = n_internal - 1
    if(id_l > 0)then
      k = k + 1
      elem(:,k) = [id_l, 1]
    endif
    if(id_r > 0)then
      k = k + 1
      elem(:,k) = [n_internal, id_r]
    endif

    call monolis_com_initialize_by_global_id(com, monolis_mpi_get_global_comm(), &
      n_internal, n_node, global_id)
    call monolis_initialize(mat)
    call monolis_get_nonzero_pattern_by_simple_mesh_R(mat, n_node, 2, n_dof, n_elem, elem)
    mat%MAT%N = n_internal
    call monolis_alloc_R_1d(x_ans, n_node*n_dof)
    do i = 1, n_node
      do j = 1, n_dof
        global_dof = (global_id(i) - 1)*n_dof + j
        x_ans((i - 1)*n_dof + j) = 1.0d0 + 0.05d0*(global_dof - 1)
        diagonal = 4.0d0
        !# 2 の整数乗なら前処理による単位行列化に丸め誤差が入らない
        if(is_diagonal) diagonal = 2.0d0**mod(global_dof - 1, 4)
        call monolis_add_scalar_to_sparse_matrix_R(mat, i, i, j, j, diagonal)
      enddo
    enddo
    if(.not. is_diagonal)then
      do i = 1, n_elem
        do j = 1, n_dof
          call monolis_add_scalar_to_sparse_matrix_R(mat, elem(1,i), elem(2,i), j, j, 1.0d0)
          call monolis_add_scalar_to_sparse_matrix_R(mat, elem(2,i), elem(1,i), j, j, 1.0d0)
        enddo
      enddo
    endif

    call monolis_dealloc_I_1d(global_id)
    call monolis_dealloc_I_2d(elem)
  end subroutine monolis_solver_IDRS_regression_get_mat
end module mod_monolis_solver_IDRS_test
