program monolis_partitioner
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_graph
  use mod_monolis_mesh
  use mod_monolis_io
  implicit none
  type(monolis_mesh) :: mesh
  type(monolis_graph) :: graph
  type(monolis_com), allocatable :: comm(:)
  type(monolis_node_list), allocatable :: node_list(:)
  integer(kint) :: n_domain
  logical :: is_format_id

  call monolis_get_arg(n_domain, is_format_id)

end program monolis_partitioner
