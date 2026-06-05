!> 3 次元構造格子 SPD 行列を CG + DIAG で「MPI 並列 + GPU」で解くサンプル
!>
!> 分割済みメッシュ (parted.0/) を読み込み、各プロセスが局所行列を組み立て、
!> ハロ通信を伴う CG 法を GPU 上で実行する。
!>
!> 使い方:
!>   mpirun -np P ./solver -ndof NDOF
!>     P    : MPI プロセス数 (= メッシュ分割数)
!>     NDOF : 1 節点あたりのブロックサイズ
program main
  use mod_monolis
  use mod_monolis_utils
  implicit none

  type(monolis_structure) :: mat
  type(monolis_com) :: com
  integer(kint) :: NDOF, NLOOP
  integer(kint) :: n_node, n_elem, n_base
  integer(kint) :: i, idof, n_internal
  integer(kint), allocatable :: elem(:,:)
  real(kdouble), allocatable :: a(:), b(:), x(:)
  real(kdouble) :: max_err
  character(monolis_charlen) :: fname
  real(kdouble), allocatable :: node(:,:)
  logical :: is_get

  call monolis_global_initialize()

  NDOF = 1
  NLOOP = 1
  call monolis_get_arg_input_I("-ndof", NDOF, is_get)
  call monolis_get_arg_input_I("-nloop", NLOOP, is_get)

  if(NDOF < 1)then
    if(monolis_mpi_get_global_my_rank() == 0) write(*,*)"NDOF must be >= 1"
    call monolis_global_finalize()
    stop
  endif

  !> 分割済み (または単一の) メッシュを読み込み
  fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "node.dat")
  call monolis_input_node(fname, n_node, node)

  fname = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "elem.dat")
  call monolis_input_elem(fname, n_elem, n_base, elem)

  if(monolis_mpi_get_global_my_rank() == 0)then
    write(*,"(a,i0)")" ** MPI size = ", monolis_mpi_get_global_comm_size()
    write(*,"(a,i0)")" ** NDOF     = ", NDOF
    write(*,"(a,i0)")" ** NLOOP    = ", NLOOP
  endif

  call monolis_initialize(mat)

  !> 通信テーブルの初期化 (並列時は分割ファイルから、単一時は自プロセスのみ)
  if(monolis_mpi_get_global_comm_size() > 1)then
    call monolis_com_initialize_by_parted_files(com, &
      monolis_mpi_get_global_comm(), &
      MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "node.dat")
  else
    call monolis_com_initialize_by_self(com)
  endif

  call monolis_get_nonzero_pattern_by_simple_mesh_R(mat, n_node, 2, NDOF, n_elem, elem)

  !> 対角ブロック: 5.95*I (3 次元 7 点ステンシル + 微小シフトで SPD)
  do i = 1, n_node
    do idof = 1, NDOF
      call monolis_add_scalar_to_sparse_matrix_R(mat, i, i, idof, idof, 5.95d0)
    enddo
  enddo

  !> 非対角ブロック: -I (各エッジについて両方向に追加)
  do i = 1, n_elem
    do idof = 1, NDOF
      call monolis_add_scalar_to_sparse_matrix_R(mat, elem(1,i), elem(2,i), idof, idof, -1.0d0)
      call monolis_add_scalar_to_sparse_matrix_R(mat, elem(2,i), elem(1,i), idof, idof, -1.0d0)
    enddo
  enddo

  !> 厳密解 x = 1 を仮定し、b = A*x を生成 (matvec はハロ通信を伴う)
  call monolis_alloc_R_1d(x, n_node*NDOF)
  call monolis_alloc_R_1d(b, n_node*NDOF)
  call monolis_alloc_R_1d(a, n_node*NDOF)

  x = 1.0d0
  call monolis_matvec_product_R(mat, com, x, b)

  !> CG + DIAG 設定
  call monolis_set_method   (mat, monolis_iter_CG)
  call monolis_set_precond  (mat, monolis_prec_DIAG)
  call monolis_set_maxiter  (mat, 10000)
  call monolis_set_tolerance(mat, 1.0d-10)
  call monolis_show_iterlog (mat, .false.)
  call monolis_show_timelog (mat, .true.)
  call monolis_show_summary (mat, .true.)

  do i = 1, NLOOP
    a = 0.0d0
    call monolis_solve_R(mat, com, b, a)
  enddo

  !> 誤差確認 (内部自由度のみ)
  if(monolis_mpi_get_global_comm_size() > 1)then
    n_internal = com%n_internal_vertex
  else
    n_internal = n_node
  endif

  max_err = 0.0d0
  do i = 1, n_internal*NDOF
    max_err = max(max_err, dabs(a(i) - 1.0d0))
  enddo
  call monolis_allreduce_R1(max_err, monolis_mpi_max, monolis_mpi_get_global_comm())

  if(monolis_mpi_get_global_my_rank() == 0)then
    write(*,"(a,1pe12.4)")" ** max error |x - 1| = ", max_err
  endif

  call monolis_finalize(mat)
  call monolis_global_finalize()

end program main
