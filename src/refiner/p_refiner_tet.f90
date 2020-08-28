program monolis_p_refiner_tet
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mesh
  use mod_monolis_io
  use mod_monolis_p_refiner
  implicit none
  type(monolis_mesh) :: mesh, mesh_ref
  character :: output_dir*100, fname*100
  logical :: is_format_id

  call monolis_set_debug(.true.)
  call monolis_debug_header("monolis_p_refiner for TET element")

  !call monolis_get_refiner_arg(is_format_id)
  is_format_id = .false.

  call monolis_input_mesh(mesh, is_format_id)

  if(mesh%nbase_func /= 4)then
    call monolis_debug_header("ERROR: please input 1st order tetra mesh")
    stop 1
  endif

  call monolis_p_refine(mesh, mesh_ref)

  output_dir = "refined/"
  call system('if [ ! -d refined ]; then (echo "** create refined"; mkdir -p refined); fi')

  fname = trim(output_dir)//"node.dat"
  call monolis_output_mesh_node(fname, mesh_ref%nnode, mesh_ref%nnode, mesh%node)

  fname = trim(output_dir)//"elem.dat"
  call monolis_output_mesh_elem(fname, mesh_ref%nelem, mesh%nbase_func, mesh_ref%elem)
end program monolis_p_refiner_tet

