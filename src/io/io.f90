module mod_monolis_solver_io
  use mod_monolis_prm
  use mod_monolis_util
  use mod_monolis_mesh
  !use mod_monolis_io_arg
  implicit none

contains

!  subroutine monolis_input_mesh(mesh, is_format_id)
!    implicit none
!    type(monolis_mesh) :: mesh
!    integer(kint) :: i, shift
!    logical :: is_format_id
!    character :: fname*100
!
!    if(is_format_id)then
!      fname = "node.dat"
!      call monolis_input_mesh_node(fname, mesh%nnode, mesh%node, mesh%nid)
!
!      fname = "elem.dat"
!      call monolis_input_mesh_elem(fname, mesh%nelem, mesh%nbase_func, mesh%elem, mesh%eid)
!    else
!      fname = "node.dat"
!      call monolis_input_mesh_node(fname, mesh%nnode, mesh%node)
!
!      fname = "elem.dat"
!      call monolis_input_mesh_elem(fname, mesh%nelem, mesh%nbase_func, mesh%elem)
!
!      shift = 0
!      if(minval(mesh%elem) == 0) shift = -1 !> for C binding
!
!      allocate(mesh%nid(mesh%nnode), source = 0)
!      do i = 1, mesh%nnode
!        mesh%nid(i) = i + shift
!      enddo
!      allocate(mesh%eid(mesh%nelem), source = 0)
!      do i = 1, mesh%nelem
!        mesh%eid(i) = i + shift
!      enddo
!    endif
!
!    call monolis_global_to_local_elem(mesh%nnode, mesh%nid, mesh%nelem, mesh%elem, mesh%nbase_func)
!  end subroutine monolis_input_mesh
!
!  subroutine monolis_output_mesh(mesh, graph, comm, node_list, n_domain)
!    implicit none
!    type(monolis_mesh) :: mesh
!    type(monolis_graph) :: graph
!    type(monolis_com) :: comm(:)
!    type(monolis_node_list) :: node_list(:)
!    integer(kint) :: i, n_domain, shift
!    character :: cnum*5, output_dir*100, fname*100
!
!    call monolis_debug_header("monolis_output_mesh")
!
!    output_dir = "parted.0/"
!    call system('if [ ! -d parted.0 ]; then (echo "** create parted.0"; mkdir -p parted.0); fi')
!
!    shift = 0
!    if(minval(mesh%nid) == 0) shift = -1 !> for C binding
!
!    do i = 1, n_domain
!      write(cnum,"(i0)") i-1
!      fname = trim(output_dir)//"node.dat."//trim(cnum)
!      call monolis_output_mesh_node(fname, node_list(i)%nnode, mesh%node, node_list(i)%nid)
!
!      fname = trim(output_dir)//"elem.dat."//trim(cnum)
!      call monolis_output_mesh_elem_ref(fname, node_list(i)%nelem, mesh%nbase_func, node_list(i)%eid, &
!        graph%ebase_func, graph%connectivity, node_list(i), shift)
!
!      fname = trim(output_dir)//"node.id."//trim(cnum)
!      call monolis_output_mesh_global_nid(fname, node_list(i)%nnode, mesh%nid, node_list(i)%nid)
!
!      fname = trim(output_dir)//"connectivity.id."//trim(cnum)
!      call monolis_output_mesh_global_eid(fname, mesh%nelem, node_list(i)%nelem, node_list(i)%nelem_in, &
!      &  mesh%eid, node_list(i)%eid, mesh%nbase_func, graph%ebase_func)
!
!      fname = trim(output_dir)//"monolis.send."//trim(cnum)
!      call monolis_output_mesh_comm(fname, comm(i)%send_n_neib, comm(i)%send_neib_pe, &
!        comm(i)%send_index, comm(i)%send_item)
!
!      fname = trim(output_dir)//"monolis.recv."//trim(cnum)
!      call monolis_output_mesh_comm(fname, comm(i)%recv_n_neib, comm(i)%recv_neib_pe, &
!        comm(i)%recv_index, comm(i)%recv_item)
!
!      fname = trim(output_dir)//"node.n_internal."//trim(cnum)
!      call monolis_output_mesh_n_internal(fname, node_list(i)%nnode_in)
!
!      fname = trim(output_dir)//"connectivity.n_internal."//trim(cnum)
!      call monolis_output_mesh_n_internal(fname, node_list(i)%nelem_in)
!    enddo
!  end subroutine monolis_output_mesh
!
!  subroutine monolis_output_parted_nodal_graph(fmain, graph, graph_format, comm, node_list, n_domain)
!    implicit none
!    type(monolis_graph) :: graph
!    type(monolis_graph_format) :: graph_format
!    type(monolis_com) :: comm(:)
!    type(monolis_node_list) :: node_list(:)
!    integer(kint) :: i, n_domain, shift
!    character :: cnum*5, output_dir*100, fname*100, fmain*100
!
!    call monolis_debug_header("monolis_output_parted_nodal_graph")
!
!    output_dir = "parted.0/"
!    call system('if [ ! -d parted.0 ]; then (echo "** create parted.0"; mkdir -p parted.0); fi')
!
!    shift = 0
!    if(minval(graph_format%point_id) == 0) shift = -1 !> for C binding
!
!    do i = 1, n_domain
!      write(cnum,"(i0)") i-1
!
!      fname = trim(output_dir)//trim(fmain)//"."//trim(cnum)
!      call monolis_output_graph_format(fname, node_list(i)%nnode, node_list(i)%nid, &
!        node_list(i)%index, node_list(i)%item, shift)
!
!      fname = trim(output_dir)//"node.id."//trim(cnum)
!      call monolis_output_mesh_global_nid(fname, node_list(i)%nnode, graph_format%point_id, node_list(i)%nid)
!
!      fname = trim(output_dir)//"connectivity.id."//trim(cnum)
!      call monolis_output_mesh_global_eid_null(fname)
!
!      fname = trim(output_dir)//"monolis.send."//trim(cnum)
!      call monolis_output_mesh_comm(fname, comm(i)%send_n_neib, comm(i)%send_neib_pe, &
!        comm(i)%send_index, comm(i)%send_item)
!
!      fname = trim(output_dir)//"monolis.recv."//trim(cnum)
!      call monolis_output_mesh_comm(fname, comm(i)%recv_n_neib, comm(i)%recv_neib_pe, &
!        comm(i)%recv_index, comm(i)%recv_item)
!
!      fname = trim(output_dir)//"node.n_internal."//trim(cnum)
!      call monolis_output_mesh_n_internal(fname, node_list(i)%nnode_in)
!
!      fname = trim(output_dir)//"connectivity.n_internal."//trim(cnum)
!      call monolis_output_mesh_n_internal(fname, 0)
!    enddo
!  end subroutine monolis_output_parted_nodal_graph
!
!  subroutine monolis_output_parted_connectivity_graph(fmain, mesh, graph, graph_format, node_list, n_domain)
!    implicit none
!    type(monolis_mesh) :: mesh(:)
!    type(monolis_graph) :: graph
!    type(monolis_graph_format) :: graph_format
!    type(monolis_node_list) :: node_list(:)
!    integer(kint) :: i, n_domain, shift, in, j, jn, kS, kE, k, kn, idx, nnode, nmin
!    integer(kint), allocatable :: perm(:)
!    character :: cnum*5, output_dir*100, fname*100, fmain*100
!
!    call monolis_debug_header("monolis_output_parted_connectivity_graph")
!
!    output_dir = "parted.0/"
!    call system('if [ ! -d parted.0 ]; then (echo "** create parted.0"; mkdir -p parted.0); fi')
!
!    shift = 0
!    nmin = 1
!    do i = 1, n_domain
!      in = minval(mesh(i)%nid)
!      if(in < nmin) nmin = in
!    enddo
!    if(nmin == 0) shift = -1 !> for C binding
!
!    do i = 1, n_domain
!      write(cnum,"(i0)") i-1
!
!      nnode = mesh(i)%nnode
!      allocate(perm(nnode), source = 0)
!      do j = 1, nnode
!        perm(j) = j
!      enddo
!      call monolis_qsort_int_with_perm(mesh(i)%nid, 1, nnode, perm)
!
!      fname = trim(output_dir)//"connectivity.dat."//trim(cnum)
!      !call monolis_output_connectivity_graph(fname, node_list(i)%nnode, node_list(i)%nid, node_list(i)%nid_perm, &
!      !  node_list(i)%index, node_list(i)%item, shift)
!      open(20, file = fname, status = "replace")
!        in = node_list(i)%nelem_in + node_list(i)%nelem_out
!        write(20,"(i0)")in
!        do j = 1, in
!          jn = node_list(i)%eid(j)
!          kS = graph%index(jn) + 1
!          kE = graph%index(jn+1)
!          write(20,"(i0,x,i0,$)")j, kE-kS+1
!          do k = kS, kE
!            kn = graph%item(k) + shift
!            call monolis_bsearch_int(mesh(i)%nid, 1, nnode, kn, idx)
!            write(20,"(x,i0,$)")perm(idx) + shift
!          enddo
!          write(20,"(a)")""
!        enddo
!      close(20)
!
!      fname = trim(output_dir)//"connectivity.id."//trim(cnum)
!      !call monolis_output_connectivity_global_eid(fname, mesh%nelem, node_list(i)%nelem, node_list(i)%nelem_in, &
!      !&  mesh%eid, node_list(i)%eid, mesh%nbase_func, graph%ebase_func)
!      open(20, file = fname, status = "replace")
!        in = node_list(i)%nelem_in + node_list(i)%nelem_out
!        write(20,"(i0)")in
!        do j = 1, in
!          write(20,"(i0,x,i0,x,i0)") j, 1, node_list(i)%eid(j) + shift
!        enddo
!      close(20)
!
!      fname = trim(output_dir)//"connectivity.n_internal."//trim(cnum)
!      !call monolis_output_connectivity_global_eid(fname, mesh%nelem, node_list(i)%nelem, node_list(i)%nelem_in, &
!      !&  mesh%eid, node_list(i)%eid, mesh%nbase_func, graph%ebase_func)
!      open(20, file = fname, status = "replace")
!        write(20,"(i0)")node_list(i)%nelem_in
!      close(20)
!
!      deallocate(perm)
!    enddo
!  end subroutine monolis_output_parted_connectivity_graph
!
!  subroutine monolis_visual_parted_mesh(nnode, node, nelem, nbase, elem, nodeid, elemid)
!    implicit none
!    integer(kint) :: i, j
!    integer(kint) :: nnode, nelem, nbase, elem(:,:), nodeid(:), elemid(:)
!    real(kdouble) :: node(:,:)
!    character :: etype*6, output_dir*100
!
!    call monolis_debug_header("monolis_visual_parted_mesh")
!
!    output_dir = "visual/"
!    call system('if [ ! -d visual ]; then (echo "** create visual"; mkdir -p visual); fi')
!
!    open(20, file = trim(output_dir)//"mesh.parted.inp", status = "replace")
!      write(20,"(5i12)") nnode, nelem, 1, 1, 0
!      do i = 1, nnode
!        write(20,"(i0,1p3e12.5)") i, node(1,i), node(2,i), node(3,i)
!      enddo
!
!      if(nbase == 3) etype = " tri  "
!      if(nbase == 4) etype = " tet  "
!      if(nbase == 8) etype = " hex  "
!      if(nbase ==10) etype = " tet2 "
!
!      do i = 1, nelem
!        write(20,"(i0,i4,a,$)") i, 0, etype
!        do j = 1, nbase
!          write(20,"(i12,$)") elem(j,i)
!        enddo
!        write(20,*)""
!      enddo
!
!      write(20,"(a)")"1 1"
!      write(20,"(a)")"node_domid, unknown"
!      do i = 1, nnode
!        write(20,"(i0,x,i0,x,i0,x,i0)") i, nodeid(i)!, bp_graph%elem_domid(i), in
!      enddo
!
!      write(20,"(a)")"1 1"
!      write(20,"(a)")"elem_domid, unknown"
!      do i = 1, nelem
!        write(20,"(i0,x,i0,x,i0,x,i0)") i, elemid(i)!, bp_graph%elem_domid(i), in
!      enddo
!    close(20)
!  end subroutine monolis_visual_parted_mesh
!
!  subroutine monolis_output_graph_format(fname, nnode, nid, ebase_func, connectivity, shift)
!    implicit none
!    integer(kint) :: nnode
!    integer(kint) :: i, in, j, jn, k, jS, shift, id
!    integer(kint) :: nid(:), ebase_func(:), connectivity(:)
!    integer(kint), allocatable :: t1(:), p1(:)
!    character :: fname*100
!
!    allocate(t1(nnode), source = 0)
!    allocate(p1(nnode), source = 0)
!    t1 = nid
!    do i = 1, nnode
!      p1(i) = i
!    enddo
!    call monolis_qsort_int_with_perm(t1, 1, nnode, p1)
!
!    open(20, file = fname, status = "replace")
!      write(20,"(i0)")nnode
!      do i = 1, nnode
!        jS = ebase_func(i)
!        in = ebase_func(i+1) - ebase_func(i)
!        write(20,"(i0,x,i0,$)") i+shift, in
!        do j = 1, in
!          k = connectivity(jS+j)
!          call monolis_bsearch_int(t1, 1, nnode, k, id)
!          write(20,"(x,i0,$)") p1(id)+shift
!        enddo
!        write(20,*)""
!      enddo
!    close(20)
!  end subroutine monolis_output_graph_format

end module mod_monolis_solver_io
