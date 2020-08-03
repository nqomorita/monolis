program monolis_partitioner
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_graph
  use mod_monolis_mesh
  use mod_monolis_comm_overlap
  use mod_monolis_io
  implicit none
  type(monolis_mesh) :: mesh
  type(monolis_graph) :: graph
  type(monolis_com), allocatable :: comm(:)
  type(monolis_node_list), allocatable :: node_list(:)
  integer(kint) :: n_domain
  logical :: is_format_id, is_overlap

  call monolis_get_part_arg(n_domain, is_format_id, is_overlap)

  call monolis_input_mesh(mesh, is_format_id)

  call monolis_part_graph(mesh, graph, n_domain)

  if(is_overlap)then
    call monolis_get_overlap_commtable(mesh, graph, comm, node_list, n_domain)
  else
    stop "monolis_partitioner: nonoverlapping partition is not supported"
  endif

  call monolis_visual_parted_mesh(mesh%nnode, mesh%node, mesh%nelem, mesh%nbase_func, mesh%elem, &
    graph%node_domid_raw, graph%elem_domid)

  call monolis_output_mesh(mesh, graph, comm, node_list, n_domain)

  call monolis_part_finalize()

contains

  subroutine monolis_part_finalize()
    implicit none

    if(allocated(mesh%node)) deallocate(mesh%node)
    if(allocated(mesh%nid)) deallocate(mesh%nid)
    if(allocated(mesh%elem)) deallocate(mesh%elem)
    if(allocated(mesh%eid)) deallocate(mesh%eid)
    if(allocated(comm)) deallocate(comm)
    if(allocated(node_list)) deallocate(node_list)
    !if(associated(mesh%elem_domid_raw)) deallocate(mesh%elem_domid_raw)
    !if(associated(mesh%elem_domid)) deallocate(mesh%elem_domid)
  end subroutine monolis_part_finalize

end program monolis_partitioner
