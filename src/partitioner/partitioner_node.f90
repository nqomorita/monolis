program monolis_partitioner_node
  use mod_monolis_mesh
  use mod_monolis_io
  implicit none
  integer(kint) :: n_domain, nnode, ndof
  type(monolis_mesh), allocatable :: mesh(:)
  real(kdouble), allocatable :: node(:,:)
  character :: fname*100

  call monolis_get_part_bc_arg(n_domain, fname)

  call monolis_input_mesh_node(fname, nnode, node)

  call monolis_par_input_node_id(n_domain, mesh)

  call monolis_par_output_node(n_domain, mesh, fname, node)

end program monolis_partitioner_node
