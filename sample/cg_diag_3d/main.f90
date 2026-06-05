!> 3 次元構造格子から SPD 行列を生成し、CG 法 + 対角前処理で解くサンプル
!>
!> 使い方:
!>   ./solver -n N -ndof NDOF
!>   N    : 1 軸方向の節点数 (総節点数は N*N*N)
!>   NDOF : 1 節点あたりのブロックサイズ
program main
  use mod_monolis
  use mod_monolis_utils
  implicit none

  type(monolis_structure) :: mat
  type(monolis_com) :: com
  integer(kint) :: N, NDOF, NLOOP
  integer(kint) :: n_node, n_elem
  integer(kint) :: ix, iy, iz, i, idof, eid
  integer(kint), allocatable :: elem(:,:)
  real(kdouble), allocatable :: a(:), b(:), x(:)
  character(monolis_charlen) :: argv
  logical :: is_get

  call monolis_global_initialize()

  !> 引数取得
  N = 10
  NDOF = 1
  NLOOP = 1

  call monolis_get_arg_input_I("-n", N, is_get)
  call monolis_get_arg_input_I("-ndof", NDOF, is_get)
  call monolis_get_arg_input_I("-nloop", NLOOP, is_get)

  if(N < 2)then
    if(monolis_mpi_get_global_my_rank() == 0) write(*,*)"N must be >= 2"
    call monolis_global_finalize()
    stop
  endif

  if(NDOF < 1)then
    if(monolis_mpi_get_global_my_rank() == 0) write(*,*)"NDOF must be >= 1"
    call monolis_global_finalize()
    stop
  endif

  if(NLOOP < 1)then
    if(monolis_mpi_get_global_my_rank() == 0) write(*,*)"NLOOP must be >= 1"
    call monolis_global_finalize()
    stop
  endif

  !> 3 次元構造格子の節点数および辺数 (= 2 節点要素数)
  n_node = N*N*N
  n_elem = 3*(N-1)*N*N

  if(monolis_mpi_get_global_my_rank() == 0)then
    write(*,"(a,i0)")" ** N     = ", n_node
    write(*,"(a,i0)")" ** NDOF  = ", NDOF
    write(*,"(a,i0)")" ** NLOOP = ", NLOOP
  endif

  call monolis_alloc_I_2d(elem, 2, n_elem)

  !> 構造格子のエッジ生成 (x, y, z 各方向)
  eid = 0
  do iz = 1, N
    do iy = 1, N
      do ix = 1, N
        if(ix < N)then
          eid = eid + 1
          elem(1,eid) = node_id(ix,   iy,   iz,   N)
          elem(2,eid) = node_id(ix+1, iy,   iz,   N)
        endif
        if(iy < N)then
          eid = eid + 1
          elem(1,eid) = node_id(ix,   iy,   iz,   N)
          elem(2,eid) = node_id(ix,   iy+1, iz,   N)
        endif
        if(iz < N)then
          eid = eid + 1
          elem(1,eid) = node_id(ix,   iy,   iz,   N)
          elem(2,eid) = node_id(ix,   iy,   iz+1, N)
        endif
      enddo
    enddo
  enddo

  !> 行列構造体の初期化と非ゼロパターンの構築
  call monolis_initialize(mat)
  call monolis_com_initialize_by_self(com)

  call monolis_get_nonzero_pattern_by_simple_mesh_R(mat, n_node, 2, NDOF, n_elem, elem)

  !> 対角ブロック: 7*I
  !>   3 次元 7 点ステンシルラプラシアン (中心 6) に微小シフトを加え、
  !>   全節点で対角優位 → SPD を保証
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

  !> 厳密解 x = 1 を仮定し、b = A*x を生成
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

  !> 解
  do i = 1, NLOOP
    a = 0.0d0
    call monolis_solve_R(mat, com, b, a)
  enddo

  !> 誤差確認
  if(monolis_mpi_get_global_my_rank() == 0)then
    write(*,"(a,1pe12.4)")" ** max error |x - 1| = ", maxval(dabs(a - 1.0d0))
  endif

  call monolis_finalize(mat)
  call monolis_global_finalize()

contains

  !> 構造格子節点 (ix, iy, iz) の 1-based 通し番号
  integer(kint) function node_id(ix, iy, iz, N)
    implicit none
    integer(kint), intent(in) :: ix, iy, iz, N
    node_id = (iz - 1)*N*N + (iy - 1)*N + ix
  end function node_id

end program main
