program monolis_dbc_all_tet
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mesh
  use mod_monolis_io
  use mod_monolis_dbc_all_util
  implicit none
  type(monolis_mesh) :: mesh
  character :: output_dir*100, fname*100
  logical :: is_format_id
  integer(kint), allocatable :: is_surf_node(:)

  call monolis_set_debug(.true.)
  call monolis_debug_header("monolis_dbc_all_tet")
  is_format_id = .false.

  call monolis_input_mesh(mesh, is_format_id)

  if(mesh%nbase_func /= 4)then
    call monolis_debug_header("ERROR: please input 1st order hex mesh")
    stop 1
  endif

  call monolis_get_surf_node(mesh, 4, 4, 4, is_surf_node)

  !call output_dbc(mesh, is_surf_node)
end program monolis_dbc_all_tet

