program monolis_dbc_all_tet
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mesh
  use mod_monolis_io
  implicit none
  type(monolis_mesh) :: mesh
  character :: output_dir*100, fname*100
  logical :: is_format_id

  call monolis_set_debug(.true.)
  call monolis_debug_header("monolis_dbc_all_tet")
  is_format_id = .false.

  call monolis_input_mesh(mesh, is_format_id)

  if(mesh%nbase_func /= 4)then
    call monolis_debug_header("ERROR: please input 1st order hex mesh")
    stop 1
  endif
end program monolis_dbc_all_tet

