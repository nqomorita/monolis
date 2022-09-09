program gedatsu_conveter_mesh2graph
  use mod_monolis_mesh
  use mod_monolis_io
  use mod_monolis_io_arg
  use mod_monolis_graph
  use iso_c_binding
  implicit none
  type(monolis_mesh) :: mesh
  type(monolis_graph) :: graph
  integer(c_int), pointer :: index(:), item(:)
  integer(kint) :: shift
  logical :: is_format_id
  character :: fname*100

  is_format_id = .false.
  call monolis_input_mesh(mesh, is_format_id)

  call monolis_convert_mesh_to_connectivity &
   & (mesh%nelem, mesh%nbase_func, mesh%elem, graph%ebase_func, graph%connectivity)

  call monolis_convert_connectivity_to_nodal_graph &
    & (mesh%nnode, mesh%nelem, graph%ebase_func, graph%connectivity, index, item)

  shift = 0
  if(minval(mesh%nid) == 0) shift = -1 !> for C binding

  fname = "graph.dat"
  call output_graph_format(fname, mesh%nnode, mesh%nid, index, item, shift)

contains

  subroutine output_graph_format(fname, nnode, nid, ebase_func, connectivity, shift)
    implicit none
    integer(kint) :: nnode
    integer(kint) :: i, in, j, jn, k, jS, shift, id
    integer(kint) :: nid(:), ebase_func(:), connectivity(:)
    character :: fname*100

    open(20, file = fname, status = "replace")
      write(20,"(i0)")nnode
      do i = 1, nnode
        jS = ebase_func(i)
        in = ebase_func(i+1) - ebase_func(i)
        write(20,"(i0,x,i0,$)") i+shift, in
        do j = 1, in
          k = connectivity(jS+j)
          write(20,"(x,i0,$)") k+shift
        enddo
        write(20,*)""
      enddo
    close(20)
  end subroutine output_graph_format

end program gedatsu_conveter_mesh2graph
