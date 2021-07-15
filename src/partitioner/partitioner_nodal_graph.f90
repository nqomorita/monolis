program monolis_partitioner_nodal_graph
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_graph
  use mod_monolis_comm_overlap
  use mod_monolis_io
  implicit none
  type(monolis_graph) :: graph
  type(monolis_graph_format) :: graph_format
  !type(monolis_com), allocatable :: comm(:)
  !type(monolis_node_list), allocatable :: node_list(:)
  integer(kint) :: n_domain
  character :: fname*100

  call monolis_get_nodal_graph_part_arg(fname, n_domain)

  call monolis_input_graph(fname, graph_format)

  call monolis_get_graph_from_graph_format(graph_format, graph)

  call monolis_part_nodal_graph(graph, n_domain)

  !call monolis_get_overlap_commtable(mesh, graph, comm, node_list, n_domain)

  !call monolis_output_parted_graph(graph, comm, node_list, n_domain)

  call monolis_part_finalize()

contains

  subroutine monolis_part_finalize()
    implicit none
    !if(allocated(comm)) deallocate(comm)
    !if(allocated(node_list)) deallocate(node_list)
  end subroutine monolis_part_finalize

end program monolis_partitioner_nodal_graph
