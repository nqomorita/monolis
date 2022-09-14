module mod_monolis_comm_overlap
  use mod_monolis_util
  use mod_monolis_util_debug
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mesh
  use mod_monolis_comm_util
  implicit none

  private

  public :: monolis_get_overlap_commtable_graph
  public :: monolis_get_overlap_commtable
  public :: get_commnication_table

contains

  subroutine monolis_get_overlap_commtable_graph(graph, comm, node_list, n_domain)
    implicit none
    type(monolis_graph) :: graph
    type(monolis_com), allocatable :: comm(:)
    type(monolis_node_list), allocatable :: node_list(:)
    integer(kint) :: n_domain

    call get_commnication_boundary_graph(graph, node_list, n_domain)
    call get_commnication_table(graph, node_list, comm, n_domain)
  end subroutine monolis_get_overlap_commtable_graph

  subroutine get_commnication_boundary_graph(graph, node_list, n_domain)
    implicit none
    type(monolis_graph) :: graph
    type(monolis_node_list), allocatable :: node_list(:)
    integer(kint) :: n_domain
    integer(kint) :: nid, in, j, k, avg

    call monolis_debug_header("get_commnication_boundary_graph")

    allocate(node_list(n_domain))

    !> get node and nid
    avg = 0
    call monolis_debug_header("get_nnode_and_nid_at_subdomain")
    write(*,"(a)")"**     nid,    total, internal,     comm"
    do nid =  1, n_domain
      call get_nnode_and_nid_at_subdomain_graph(graph, node_list(nid), nid, avg)
    enddo

    write(*,"(4i10)") avg/n_domain
  end subroutine get_commnication_boundary_graph

  subroutine monolis_get_overlap_commtable(graph, comm, node_list, n_domain)
    implicit none
    type(monolis_graph) :: graph
    type(monolis_com), allocatable :: comm(:)
    type(monolis_node_list), allocatable :: node_list(:)
    integer(kint) :: n_domain

    call get_overlap_domain(graph, n_domain)
    call get_commnication_boundary(graph, node_list, n_domain)
    call get_commnication_table(graph, node_list, comm, n_domain)
  end subroutine monolis_get_overlap_commtable

  subroutine get_commnication_boundary(graph, node_list, n_domain)
    implicit none
    type(monolis_graph) :: graph
    type(monolis_node_list), allocatable :: node_list(:)
    integer(kint) :: n_domain
    integer(kint) :: nid, in, j, k, avg

    call monolis_debug_header("get_commnication_boundary")

    allocate(node_list(n_domain))

    !> get local elem and eid
    call monolis_debug_header("get_nelem_and_eid_at_subdomain")
    do nid =  1, n_domain
      call get_nelem_and_eid_at_subdomain(graph, node_list(nid), nid)
    enddo

    !> get node and nid
    avg = 0
    call monolis_debug_header("get_nnode_and_nid_at_subdomain")
    write(*,"(a)")"**     nid,    total, internal,     comm"
    do nid =  1, n_domain
      call get_nnode_and_nid_at_subdomain(graph, node_list(nid), nid, avg)
    enddo

    write(*,"(4i10)") avg/n_domain
  end subroutine get_commnication_boundary

  subroutine get_commnication_table(graph, node_list, comm, n_domain)
    implicit none
    type(monolis_graph) :: graph
    type(monolis_node_list) :: node_list(:)
    type(monolis_com), allocatable :: comm(:)
    integer(kint) :: nnode, nelem, n_domain

    call monolis_debug_header("get_commnication_table")

    allocate(comm(n_domain))

    call get_neib_PE(node_list, graph, comm, n_domain)
    call get_recv_table(node_list, graph, comm, n_domain)
    call get_send_table(node_list, graph, comm, n_domain)
  end subroutine get_commnication_table
end module mod_monolis_comm_overlap
