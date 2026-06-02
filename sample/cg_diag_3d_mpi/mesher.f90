!> 3 次元構造格子のメッシュ (node.dat / elem.dat) を出力するメッシャ
!>
!> 使い方:
!>   ./mesher -n N
!>   N : 1 軸方向の節点数 (総節点数は N*N*N, 要素は 2 節点エッジ)
program mesher
  use mod_monolis
  use mod_monolis_utils
  implicit none

  integer(kint) :: N
  integer(kint) :: n_node, n_elem
  integer(kint) :: ix, iy, iz, eid
  integer(kint), allocatable :: elem(:,:)
  real(kdouble), allocatable :: node(:,:)
  logical :: is_get

  call monolis_global_initialize()

  N = 10
  call monolis_get_arg_input_I("-n", N, is_get)

  if(N < 2)then
    write(*,*)"N must be >= 2"
    call monolis_global_finalize()
    stop
  endif

  n_node = N*N*N
  n_elem = 3*(N-1)*N*N

  !> 節点座標 (格子点)
  call monolis_alloc_R_2d(node, 3, n_node)
  do iz = 1, N
    do iy = 1, N
      do ix = 1, N
        node(1, node_id(ix, iy, iz, N)) = dble(ix - 1)
        node(2, node_id(ix, iy, iz, N)) = dble(iy - 1)
        node(3, node_id(ix, iy, iz, N)) = dble(iz - 1)
      enddo
    enddo
  enddo

  !> エッジ (2 節点要素) の生成
  call monolis_alloc_I_2d(elem, 2, n_elem)
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

  call monolis_output_node("node.dat", n_node, node)
  call monolis_output_elem("elem.dat", n_elem, 2, elem)

  write(*,"(a,i0,a,i0)")" ** mesher: n_node = ", n_node, ", n_elem = ", n_elem

  call monolis_global_finalize()

contains

  integer(kint) function node_id(ix, iy, iz, N)
    implicit none
    integer(kint), intent(in) :: ix, iy, iz, N
    node_id = (iz - 1)*N*N + (iy - 1)*N + ix
  end function node_id

end program mesher
