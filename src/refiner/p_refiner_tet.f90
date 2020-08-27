program monolis_refiner
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mesh
  use mod_monolis_io
  use mod_monolis_p_refiner
  implicit none
  type(monolis_mesh) :: mesh, mesh_ref
  logical :: is_format_id

  call monolis_set_debug(.true.)
  call monolis_debug_header("monolis_p_refiner for TET element")

  !call monolis_get_refiner_arg(is_format_id)
  is_format_id = .false.

  call monolis_input_mesh(mesh, is_format_id)
  call monolis_p_refine(mesh, mesh_ref)
  call monolis_output_refined_mesh(mesh_ref, is_format_id)
end program monolis_refiner
