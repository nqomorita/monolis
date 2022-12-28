module mod_monolis_io
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_util
  use mod_monolis_util_debug
  use mod_monolis_mesh
  use mod_monolis_stdlib
  use mod_monolis_io_arg
  implicit none

contains

  subroutine monolis_input_mesh(mesh, is_format_id)
    implicit none
    type(monolis_mesh) :: mesh
    integer(kint) :: i, shift
    logical :: is_format_id
    character :: fname*100

    if(is_format_id)then
      fname = "node.dat"
      call monolis_input_mesh_node(fname, mesh%nnode, mesh%node, mesh%nid)

      fname = "elem.dat"
      call monolis_input_mesh_elem(fname, mesh%nelem, mesh%nbase_func, mesh%elem, mesh%eid)
    else
      fname = "node.dat"
      call monolis_input_mesh_node(fname, mesh%nnode, mesh%node)

      fname = "elem.dat"
      call monolis_input_mesh_elem(fname, mesh%nelem, mesh%nbase_func, mesh%elem)

      shift = 0
      if(minval(mesh%elem) == 0) shift = -1 !> for C binding

      allocate(mesh%nid(mesh%nnode), source = 0)
      do i = 1, mesh%nnode
        mesh%nid(i) = i + shift
      enddo
      allocate(mesh%eid(mesh%nelem), source = 0)
      do i = 1, mesh%nelem
        mesh%eid(i) = i + shift
      enddo
    endif

    call monolis_global_to_local_elem(mesh%nnode, mesh%nid, mesh%nelem, mesh%elem, mesh%nbase_func)
  end subroutine monolis_input_mesh

  subroutine monolis_output_mesh(mesh, graph, comm, node_list, n_domain)
    implicit none
    type(monolis_mesh) :: mesh
    type(monolis_graph) :: graph
    type(monolis_com) :: comm(:)
    type(monolis_node_list) :: node_list(:)
    integer(kint) :: i, n_domain, shift
    character :: cnum*5, output_dir*100, fname*100

    call monolis_debug_header("monolis_output_mesh")

    output_dir = "parted.0/"
    call system('if [ ! -d parted.0 ]; then (echo "** create parted.0"; mkdir -p parted.0); fi')

    shift = 0
    if(minval(mesh%nid) == 0) shift = -1 !> for C binding

    do i = 1, n_domain
      write(cnum,"(i0)") i-1
      fname = trim(output_dir)//"node.dat."//trim(cnum)
      call monolis_output_mesh_node(fname, node_list(i)%nnode, mesh%node, node_list(i)%nid)

      fname = trim(output_dir)//"elem.dat."//trim(cnum)
      call monolis_output_mesh_elem_ref(fname, node_list(i)%nelem, mesh%nbase_func, node_list(i)%eid, &
        graph%ebase_func, graph%connectivity, node_list(i), shift)

      fname = trim(output_dir)//"node.id."//trim(cnum)
      call monolis_output_mesh_global_nid(fname, node_list(i)%nnode, mesh%nid, node_list(i)%nid)

      fname = trim(output_dir)//"connectivity.id."//trim(cnum)
      call monolis_output_mesh_global_eid(fname, mesh%nelem, node_list(i)%nelem, node_list(i)%nelem_in, &
      &  mesh%eid, node_list(i)%eid, mesh%nbase_func, graph%ebase_func)

      fname = trim(output_dir)//"monolis.send."//trim(cnum)
      call monolis_output_mesh_comm(fname, comm(i)%send_n_neib, comm(i)%send_neib_pe, &
        comm(i)%send_index, comm(i)%send_item)

      fname = trim(output_dir)//"monolis.recv."//trim(cnum)
      call monolis_output_mesh_comm(fname, comm(i)%recv_n_neib, comm(i)%recv_neib_pe, &
        comm(i)%recv_index, comm(i)%recv_item)

      fname = trim(output_dir)//"node.n_internal."//trim(cnum)
      call monolis_output_mesh_n_internal(fname, node_list(i)%nnode_in)

      fname = trim(output_dir)//"connectivity.n_internal."//trim(cnum)
      call monolis_output_mesh_n_internal(fname, node_list(i)%nelem_in)
    enddo
  end subroutine monolis_output_mesh

  subroutine monolis_output_parted_nodal_graph(fmain, graph, graph_format, comm, node_list, n_domain)
    implicit none
    type(monolis_graph) :: graph
    type(monolis_graph_format) :: graph_format
    type(monolis_com) :: comm(:)
    type(monolis_node_list) :: node_list(:)
    integer(kint) :: i, n_domain, shift
    character :: cnum*5, output_dir*100, fname*100, fmain*100

    call monolis_debug_header("monolis_output_parted_nodal_graph")

    output_dir = "parted.0/"
    call system('if [ ! -d parted.0 ]; then (echo "** create parted.0"; mkdir -p parted.0); fi')

    shift = 0
    if(minval(graph_format%point_id) == 0) shift = -1 !> for C binding

    do i = 1, n_domain
      write(cnum,"(i0)") i-1

      fname = trim(output_dir)//trim(fmain)//"."//trim(cnum)
      call monolis_output_graph_format(fname, node_list(i)%nnode, node_list(i)%nid, &
        node_list(i)%index, node_list(i)%item, shift)

      fname = trim(output_dir)//"node.id."//trim(cnum)
      call monolis_output_mesh_global_nid(fname, node_list(i)%nnode, graph_format%point_id, node_list(i)%nid)

      fname = trim(output_dir)//"connectivity.id."//trim(cnum)
      call monolis_output_mesh_global_eid_null(fname)

      fname = trim(output_dir)//"monolis.send."//trim(cnum)
      call monolis_output_mesh_comm(fname, comm(i)%send_n_neib, comm(i)%send_neib_pe, &
        comm(i)%send_index, comm(i)%send_item)

      fname = trim(output_dir)//"monolis.recv."//trim(cnum)
      call monolis_output_mesh_comm(fname, comm(i)%recv_n_neib, comm(i)%recv_neib_pe, &
        comm(i)%recv_index, comm(i)%recv_item)

      fname = trim(output_dir)//"node.n_internal."//trim(cnum)
      call monolis_output_mesh_n_internal(fname, node_list(i)%nnode_in)

      fname = trim(output_dir)//"connectivity.n_internal."//trim(cnum)
      call monolis_output_mesh_n_internal(fname, 0)
    enddo
  end subroutine monolis_output_parted_nodal_graph

  subroutine monolis_output_parted_connectivity_graph(fmain, mesh, graph, graph_format, node_list, n_domain)
    implicit none
    type(monolis_mesh) :: mesh(:)
    type(monolis_graph) :: graph
    type(monolis_graph_format) :: graph_format
    type(monolis_node_list) :: node_list(:)
    integer(kint) :: i, n_domain, shift, in, j, jn, kS, kE, k, kn, idx, nnode, nmin
    integer(kint), allocatable :: perm(:)
    character :: cnum*5, output_dir*100, fname*100, fmain*100

    call monolis_debug_header("monolis_output_parted_connectivity_graph")

    output_dir = "parted.0/"
    call system('if [ ! -d parted.0 ]; then (echo "** create parted.0"; mkdir -p parted.0); fi')

    shift = 0
    nmin = 1
    do i = 1, n_domain
      in = minval(mesh(i)%nid)
      if(in < nmin) nmin = in
    enddo
    if(nmin == 0) shift = -1 !> for C binding

    do i = 1, n_domain
      write(cnum,"(i0)") i-1

      nnode = mesh(i)%nnode
      allocate(perm(nnode), source = 0)
      do j = 1, nnode
        perm(j) = j
      enddo
      call monolis_qsort_int_with_perm(mesh(i)%nid, 1, nnode, perm)

      fname = trim(output_dir)//"connectivity.dat."//trim(cnum)
      !call monolis_output_connectivity_graph(fname, node_list(i)%nnode, node_list(i)%nid, node_list(i)%nid_perm, &
      !  node_list(i)%index, node_list(i)%item, shift)
      open(20, file = fname, status = "replace")
        in = node_list(i)%nelem_in + node_list(i)%nelem_out
        write(20,"(i0)")in
        do j = 1, in
          jn = node_list(i)%eid(j)
          kS = graph%index(jn) + 1
          kE = graph%index(jn+1)
          write(20,"(i0,x,i0,$)")j, kE-kS+1
          do k = kS, kE
            kn = graph%item(k) + shift
            call monolis_bsearch_int(mesh(i)%nid, 1, nnode, kn, idx)
            write(20,"(x,i0,$)")perm(idx) + shift
          enddo
          write(20,"(a)")""
        enddo
      close(20)

      fname = trim(output_dir)//"connectivity.id."//trim(cnum)
      !call monolis_output_connectivity_global_eid(fname, mesh%nelem, node_list(i)%nelem, node_list(i)%nelem_in, &
      !&  mesh%eid, node_list(i)%eid, mesh%nbase_func, graph%ebase_func)
      open(20, file = fname, status = "replace")
        in = node_list(i)%nelem_in + node_list(i)%nelem_out
        write(20,"(i0)")in
        do j = 1, in
          write(20,"(i0,x,i0,x,i0)") j, 1, node_list(i)%eid(j) + shift
        enddo
      close(20)

      fname = trim(output_dir)//"connectivity.n_internal."//trim(cnum)
      !call monolis_output_connectivity_global_eid(fname, mesh%nelem, node_list(i)%nelem, node_list(i)%nelem_in, &
      !&  mesh%eid, node_list(i)%eid, mesh%nbase_func, graph%ebase_func)
      open(20, file = fname, status = "replace")
        write(20,"(i0)")node_list(i)%nelem_in
      close(20)

      deallocate(perm)
    enddo
  end subroutine monolis_output_parted_connectivity_graph

  subroutine monolis_input_mesh_node(fname, nnode, node, nid)
    implicit none
    integer(kint) :: nnode, i, t(4)
    real(kdouble), allocatable :: node(:,:)
    integer(kint), optional, allocatable :: nid(:)
    character :: fname*100

    open(20, file = fname, status = "old")
      if(present(nid))then
        read(20,*) nnode
        allocate(node(3,nnode), source = 0.0d0)
        allocate(nid(nnode), source = 0)
        do i = 1, nnode
          read(20,*) nid(i), node(1,i), node(2,i), node(3,i)
        enddo
      else
        read(20,*) nnode
        allocate(node(3,nnode), source = 0.0d0)
        do i = 1, nnode
          read(20,*) node(1,i), node(2,i), node(3,i)
        enddo
      endif
    close(20)

    call monolis_debug_int("nnode", nnode)
  end subroutine monolis_input_mesh_node

  subroutine monolis_input_mesh_distval_r(fname, nnode, ndof, val, label)
    implicit none
    integer(kint) :: nnode, ndof, i, j
    real(kdouble), allocatable :: val(:,:)
    character :: fname*100, label*100

    open(20, file = fname, status = "old")
      read(20,"(a100)") label
      read(20,*) nnode, ndof
      allocate(val(ndof,nnode), source = 0.0d0)
      do i = 1, nnode
        read(20,*) (val(j,i), j = 1, ndof)
      enddo
    close(20)

    call monolis_debug_int("nnode", nnode)
  end subroutine monolis_input_mesh_distval_r

  subroutine monolis_input_mesh_distval_i(fname, nnode, ndof, val, label)
    implicit none
    integer(kint) :: nnode, ndof, i, j
    integer(kint), allocatable :: val(:,:)
    character :: fname*100, label*100

    open(20, file = fname, status = "old")
      read(20,"(a100)") label
      read(20,*) nnode, ndof
      allocate(val(ndof,nnode), source = 0)
      do i = 1, nnode
        read(20,*) (val(j,i), j = 1, ndof)
      enddo
    close(20)

    call monolis_debug_int("nnode", nnode)
  end subroutine monolis_input_mesh_distval_i

  subroutine monolis_input_mesh_elem(fname, nelem, nbase, elem, eid)
    implicit none
    integer(kint) :: nelem, nbase, i, j
    integer(kint), allocatable :: elem(:,:)
    integer(kint), optional, allocatable :: eid(:)
    character :: fname*100

    open(20, file = fname, status = "old")
      read(20,*) nelem, nbase
      call monolis_debug_int("nelem", nelem)
      call monolis_debug_int("nbase_func", nbase)

      allocate(elem(nbase,nelem), source = 0)

      if(present(eid))then
        allocate(eid(nelem), source = 0)
        do i = 1, nelem
          read(20,*) eid(i), (elem(j,i), j = 1, nbase)
        enddo
      else
        do i = 1, nelem
          read(20,*) (elem(j,i), j = 1, nbase)
        enddo
      endif
    close(20)
  end subroutine monolis_input_mesh_elem

  subroutine monolis_input_id(fname, id, nid_out)
    implicit none
    integer(kint) :: nid, i, j
    integer(kint), optional :: nid_out
    integer(kint), allocatable :: id(:)
    character :: fname*100

    open(20, file = fname, status = "old")
      read(20,*) nid
      call monolis_debug_int("nid", nid)

      allocate(id(nid), source = 0)

      do i = 1, nid
        read(20,*) j, j, id(i)
      enddo
    close(20)

    if(present(nid_out)) nid_out = nid
  end subroutine monolis_input_id

  subroutine monolis_par_input_node_id(n_domain, mesh)
    implicit none
    type(monolis_mesh), allocatable :: mesh(:)
    integer(kint) :: n_domain, i
    character :: fname*100, cnum*5

    do i = 1, n_domain
      write(cnum,"(i0)") i-1
      fname = "parted.0/node.id."//trim(cnum)
      call monolis_input_id(fname, mesh(i)%nid, mesh(i)%nnode)
    enddo
  end subroutine monolis_par_input_node_id

  subroutine monolis_par_input_connectivity_id(n_domain, mesh)
    implicit none
    type(monolis_mesh), allocatable :: mesh(:)
    integer(kint) :: n_domain, i
    character :: fname*100, cnum*5

    do i = 1, n_domain
      write(cnum,"(i0)") i-1
      fname = "parted.0/connectivity.id."//trim(cnum)
      call monolis_input_id(fname, mesh(i)%eid, mesh(i)%nelem)
    enddo
  end subroutine monolis_par_input_connectivity_id

  subroutine monolis_par_input_node_n_internal(n_domain, mesh)
    implicit none
    type(monolis_mesh), allocatable :: mesh(:)
    integer(kint) :: n_domain, i
    character :: fname*100, cnum*5

    do i = 1, n_domain
      write(cnum,"(i0)") i-1
      fname = "parted.0/node.n_internal."//trim(cnum)
      call monolis_input_n_internal(fname, mesh(i)%nnode_in)
    enddo
  end subroutine monolis_par_input_node_n_internal

  subroutine monolis_input_condition(fname, ncond, ndof, icond, cond)
    implicit none
    integer(kint) :: ncond, ndof, i, j
    integer(kint), allocatable :: icond(:,:)
    real(kdouble), allocatable :: cond(:)
    character :: fname*100

    open(20, file = fname, status = "old")
      read(20,*) ncond, ndof
      call monolis_debug_int("ncond", ncond)
      call monolis_debug_int("ndof", ndof)

      allocate(icond(2,ncond), source = 0)
      allocate(cond(ncond), source = 0.0d0)

      do i = 1, ncond
        read(20,*) icond(1,i), icond(2,i), cond(i)
      enddo
    close(20)
  end subroutine monolis_input_condition

  subroutine monolis_input_mesh_restart_data(fname, n, nbase, var, gid)
    implicit none
    integer(kint) :: n, nbase, i, j, in, nid, id
    integer(kint) :: gid(:)
    real(kdouble) :: r(nbase)
    character :: fname*100
    integer(kint), allocatable :: perm(:), temp(:)
    real(kdouble), allocatable :: var(:)

    nid = size(gid)
    allocate(temp(nid), source = 0)
    allocate(perm(nid), source = 0)
    do i = 1, nid
      temp(i) = gid(i)
      perm(i) = i
    enddo
    call monolis_qsort_int_with_perm(temp, 1, nid, perm)

    open(20, file = fname, status = "old")
      read(20,*) n
      call monolis_debug_int("nvar", n)
      call monolis_debug_int("nvar_dof", nbase)

      allocate(var(nbase*n), source = 0.0d0)

      do i = 1, n
        read(20,*) in, (r(j), j = 1, nbase)
        call monolis_bsearch_int(temp, 1, nid, in, id)
        do j = 1, nbase
          var(nbase*id - nbase + j) = r(j)
        enddo
      enddo
    close(20)
  end subroutine monolis_input_mesh_restart_data

  subroutine monolis_input_graph(fname, graph_format)
    implicit none
    type(monolis_graph_format) :: graph_format
    integer(kint) :: i, in, j, tmp, NZ
    character :: fname*100

    NZ = 0
    open(20, file = fname, status = "old")
      read(20,*) graph_format%n_point
      do i = 1, graph_format%n_point
        read(20,*) tmp, in, (tmp, j = 1, in)
        NZ = NZ + in
      enddo
    close(20)

    allocate(graph_format%point_id(graph_format%n_point), source = 0)
    allocate(graph_format%n_adjacent(graph_format%n_point), source = 0)
    allocate(graph_format%adjacent_id(NZ), source = 0)

    NZ = 0
    open(20, file = fname, status = "old")
      read(20,*) graph_format%n_point
      do i = 1, graph_format%n_point
        read(20,*) graph_format%point_id(i), in, (graph_format%adjacent_id(NZ+j), j = 1, in)
        graph_format%n_adjacent(i) = in
        NZ = NZ + in
      enddo
    close(20)
  end subroutine monolis_input_graph

  subroutine monolis_output_mesh_n_internal(fname, n_in)
    implicit none
    integer(kint) :: n_in
    character :: fname*100

    open(20, file = fname, status = "replace")
      write(20,"(i0)")n_in
    close(20)
  end subroutine monolis_output_mesh_n_internal

  subroutine monolis_input_n_internal(fname, n_in)
    implicit none
    integer(kint) :: n_in
    character :: fname*100

    open(20, file = fname, status = "old")
      read(20,*)n_in
    close(20)
  end subroutine monolis_input_n_internal

  subroutine monolis_output_mesh_node(fname, nnode, node, nid)
    implicit none
    integer(kint) :: nnode
    integer(kint) :: i, in
    integer(kint), optional :: nid(:)
    real(kdouble) :: node(:,:)
    character :: fname*100

    open(20, file = fname, status = "replace")
    if(present(nid))then
      write(20,"(i0,x,i0)")nnode
      do i = 1, nnode
        in = nid(i)
        write(20,"(1p3e22.14)") node(1,in), node(2,in), node(3,in)
      enddo
    else
      write(20,"(i0)")nnode
      do i = 1, nnode
        write(20,"(1p3e22.14)") node(1,i), node(2,i), node(3,i)
      enddo
    endif
    close(20)
  end subroutine monolis_output_mesh_node

  subroutine monolis_output_mesh_elem(fname, nelem, nbase, elem)
    implicit none
    integer(kint) :: i, j, nelem, nbase, elem(:,:)
    character :: fname*100

    open(20, file = fname, status = "replace")
      write(20,"(i0,x,i0)")nelem, nbase
      do i = 1, nelem
        do j = 1, nbase
          write(20,"(x,i0,$)") elem(j,i)
        enddo
        write(20,*)""
      enddo
    close(20)
  end subroutine monolis_output_mesh_elem

  subroutine monolis_output_mesh_elem_ref(fname, nelem, nbase, eid, ebase_func, conn, node_list, shift)
    implicit none
    type(monolis_node_list) :: node_list
    integer(kint) :: i, j, in, jn, nelem, nbase, eid(:), ebase_func(:), conn(:), shift, jS, jE
    character :: fname*100

    in = 0
    do i = 1, nelem
      jn = eid(i)
      jS = ebase_func(jn) + 1
      jE = ebase_func(jn+1)
      if(jE-jS+1 == nbase) in = in + 1
    enddo

    open(20, file = fname, status = "replace")
      write(20,"(i0,x,i0)")in, nbase
      do i = 1, nelem
        jn = eid(i)
        !do j = 1, nbase
        jS = ebase_func(jn) + 1
        jE = ebase_func(jn+1)
        if(jE-jS+1 /= nbase) cycle
        do j = jS, jE
          in = conn(j)
          write(20,"(x,i0,$)") node_list%nid_perm(in) + shift
        enddo
        write(20,*)""
      enddo
    close(20)
  end subroutine monolis_output_mesh_elem_ref

  subroutine monolis_output_mesh_global_nid(fname, nnode, global_nid, nid)
    implicit none
    integer(kint) :: i, in, nnode, global_nid(:), nid(:)
    character :: fname*100

    open(20, file = fname, status = "replace")
      write(20,"(i0)")nnode
      do i = 1, nnode
        in = nid(i)
        write(20,"(i0,x,i0,x,i0)") i, 1, global_nid(in)
      enddo
    close(20)
  end subroutine monolis_output_mesh_global_nid

  subroutine monolis_output_mesh_global_eid(fname, nelem_all, nelem, nelem_in, global_eid, eid, &
  &  nbase, ebase_func)
    implicit none
    integer(kint) :: i, in, jn, kn, nelem, nelem_in, global_eid(:), eid(:), nelem_all, jE, jS
    integer(kint) :: ebase_func(:), nbase
    character :: fname*100

    in = 0
    kn = 0
    do i = 1, nelem
      jn = eid(i)
      jS = ebase_func(jn) + 1
      jE = ebase_func(jn+1)
      if(jE-jS+1 == nbase) in = in + 1
      if(jE-jS+1 == nbase .and. i <= nelem_in) kn = kn + 1
    enddo
    nelem_in = kn

    open(20, file = fname, status = "replace")
      write(20,"(i0,x,i0)")in
      do i = 1, nelem
        in = eid(i)
        if(in > nelem_all) cycle
        write(20,"(i0,x,i0,x,i0)") i, 1, global_eid(in)
      enddo
    close(20)
  end subroutine monolis_output_mesh_global_eid

  subroutine monolis_output_mesh_global_eid_null(fname)
    implicit none
    character :: fname*100

    open(20, file = fname, status = "replace")
      write(20,"(i0,x,i0)")0
    close(20)
  end subroutine monolis_output_mesh_global_eid_null

  subroutine monolis_par_output_node(n_domain, mesh, fname_body, node)
    use mod_monolis_mesh
    implicit none
    type(monolis_mesh) :: mesh(:)
    integer(kint) :: n_domain, i, j, in, nnode, shift, nmin
    real(kdouble) :: node(:,:)
    character :: fname*100, fname_body*100, cnum*5

    shift = 0
    nmin = 1
    do i = 1, n_domain
      in = minval(mesh(i)%nid)
      if(in < nmin) nmin = in
    enddo
    if(nmin == 0) shift = 1

    do i = 1, n_domain
      write(cnum,"(i0)") i-1
      fname = "parted.0/"//trim(fname_body)//"."//trim(cnum)

      open(20, file = fname, status = "replace")
        write(20,"(i0)") mesh(i)%nnode
        do j = 1, mesh(i)%nnode
          in = mesh(i)%nid(j) + shift
          write(20,"(1p3e22.12)") node(1,in), node(2,in), node(3,in)
        enddo
      close(20)
    enddo
  end subroutine monolis_par_output_node

  subroutine monolis_par_output_nodal_val_r(n_domain, mesh, fname_body, ndof, val, label)
    use mod_monolis_mesh
    implicit none
    type(monolis_mesh) :: mesh(:)
    integer(kint) :: n_domain, i, j, k, in, nnode, shift, nmin, ndof
    real(kdouble) :: val(:,:)
    character :: fname*100, fname_body*100, cnum*5, label*100

    shift = 0
    nmin = 1
    do i = 1, n_domain
      in = minval(mesh(i)%nid)
      if(in < nmin) nmin = in
    enddo
    if(nmin == 0) shift = 1

    do i = 1, n_domain
      write(cnum,"(i0)") i-1
      fname = "parted.0/"//trim(fname_body)//"."//trim(cnum)

      open(20, file = fname, status = "replace")
        write(20,"(a)") trim(label)
        write(20,"(i0,x,i0)") mesh(i)%nnode, ndof
        do j = 1, mesh(i)%nnode
          in = mesh(i)%nid(j) + shift
          write(20,"(1p20e22.12)") (val(k,in), k = 1, ndof)
        enddo
      close(20)
    enddo
  end subroutine monolis_par_output_nodal_val_r

  subroutine monolis_par_output_nodal_val_i(n_domain, mesh, fname_body, ndof, val, label)
    use mod_monolis_mesh
    implicit none
    type(monolis_mesh) :: mesh(:)
    integer(kint) :: n_domain, i, j, k, in, nnode, shift, nmin, ndof
    integer(kint) :: val(:,:)
    character :: fname*100, fname_body*100, cnum*5, label*100

    shift = 0
    nmin = 1
    do i = 1, n_domain
      in = minval(mesh(i)%eid)
      if(in < nmin) nmin = in
    enddo
    if(nmin == 0) shift = 1

    do i = 1, n_domain
      write(cnum,"(i0)") i-1
      fname = "parted.0/"//trim(fname_body)//"."//trim(cnum)

      open(20, file = fname, status = "replace")
        write(20,"(a)") trim(label)
        write(20,"(i0,x,i0)") mesh(i)%nnode, ndof
        do j = 1, mesh(i)%nnode
          in = mesh(i)%eid(j) + shift
          write(20,"(20(i0,x))") (val(k,in), k = 1, ndof)
        enddo
      close(20)
    enddo
  end subroutine monolis_par_output_nodal_val_i

  subroutine monolis_par_output_connectivity_val_r(n_domain, mesh, fname_body, ndof, val, label)
    use mod_monolis_mesh
    implicit none
    type(monolis_mesh) :: mesh(:)
    integer(kint) :: n_domain, i, j, k, in, nelem, shift, nmin, ndof
    real(kdouble) :: val(:,:)
    character :: fname*100, fname_body*100, cnum*5, label*100

    shift = 0
    nmin = 1
    do i = 1, n_domain
      in = minval(mesh(i)%eid)
      if(in < nmin) nmin = in
    enddo
    if(nmin == 0) shift = 1

    do i = 1, n_domain
      write(cnum,"(i0)") i-1
      fname = "parted.0/"//trim(fname_body)//"."//trim(cnum)

      open(20, file = fname, status = "replace")
        write(20,"(a)") trim(label)
        write(20,"(i0,x,i0)") mesh(i)%nelem, ndof
        do j = 1, mesh(i)%nelem
          in = mesh(i)%eid(j) + shift
          write(20,"(1p20e22.12)") (val(k,in), k = 1, ndof)
        enddo
      close(20)
    enddo
  end subroutine monolis_par_output_connectivity_val_r

  subroutine monolis_par_output_connectivity_val_i(n_domain, mesh, fname_body, ndof, val, label)
    use mod_monolis_mesh
    implicit none
    type(monolis_mesh) :: mesh(:)
    integer(kint) :: n_domain, i, j, k, in, nelem, shift, nmin, ndof
    integer(kint) :: val(:,:)
    character :: fname*100, fname_body*100, cnum*5, label*100

    shift = 0
    nmin = 1
    do i = 1, n_domain
      in = minval(mesh(i)%nid)
      if(in < nmin) nmin = in
    enddo
    if(nmin == 0) shift = 1

    do i = 1, n_domain
      write(cnum,"(i0)") i-1
      fname = "parted.0/"//trim(fname_body)//"."//trim(cnum)

      open(20, file = fname, status = "replace")
        write(20,"(a)") trim(label)
        write(20,"(i0,x,i0)") mesh(i)%nelem, ndof
        do j = 1, mesh(i)%nelem
          in = mesh(i)%nid(j) + shift
          write(20,"(20(i0,x))") (val(k,in), k = 1, ndof)
        enddo
      close(20)
    enddo
  end subroutine monolis_par_output_connectivity_val_i

  subroutine monolis_par_output_condition(n_domain, mesh, fname_body, ncond_all, ndof, icond, cond)
    use mod_monolis_mesh
    implicit none
    type(monolis_mesh) :: mesh(:)
    integer(kint) :: n_domain, ncond_all, ndof, i, j, k, in, nid, ncond, shift, nmin
    integer(kint) :: icond(:,:)
    real(kdouble) :: cond(:)
    integer(kint), allocatable :: temp(:), perm(:)
    character :: fname*100, fname_body*100, cnum*5

    nmin = 1
    shift = 0
    do i = 1, n_domain
      in = minval(mesh(i)%nid)
      if(in < nmin) nmin = in
    enddo
    if(nmin == 0) shift = -1

    do i = 1, n_domain
      allocate(temp(mesh(i)%nnode), source = 0)
      allocate(perm(mesh(i)%nnode), source = 0)
      do j = 1, mesh(i)%nnode
        perm(j) = j
      enddo
      temp = mesh(i)%nid(:)
      call monolis_qsort_int_with_perm(temp, 1, mesh(i)%nnode, perm)

      write(cnum,"(i0)") i-1
      fname = "parted.0/"//trim(fname_body)//"."//trim(cnum)

      ncond = 0
      do j = 1, ncond_all
        nid = icond(1,j)
        call monolis_bsearch_int(temp, 1, mesh(i)%nnode, nid, in)
        if(in /= -1) ncond = ncond + 1
      enddo

      open(20, file = fname, status = "replace")
        write(20,"(i0,x,i0)") ncond, ndof
        do j = 1, ncond_all
          nid = icond(1,j)
          call monolis_bsearch_int(temp, 1, mesh(i)%nnode, nid, in)
          if(in == -1) cycle
          write(20,"(i0,x,i0,x,1pe12.5)") perm(in) + shift, icond(2,j), cond(j)
        enddo
      close(20)
      deallocate(temp)
      deallocate(perm)
    enddo
  end subroutine monolis_par_output_condition

  subroutine monolis_output_mesh_comm(fname, n_neib, neib_pe, index, item)
    implicit none
    integer(kint) :: n_neib
    integer(kint) :: neib_pe(:)
    integer(kint) :: index(0:n_neib)
    integer(kint) :: item(:)
    integer(kint) :: i
    character :: fname*100

    open(20, file = fname, status = "replace")
      !> for overlap
      write(20,"(i0,x,i0)")n_neib, index(n_neib)
      do i = 1, n_neib
        write(20,"(i0)")neib_pe(i) - 1 !> for MPI
      enddo

      do i = 0, n_neib
        write(20,"(i0)")index(i)
      enddo
      do i = 1, index(n_neib)
        write(20,"(i0)")item(i)
      enddo

      !> for non-overlap
      write(20,"(i0,x,i0)")0, 0
    close(20)
  end subroutine monolis_output_mesh_comm

  subroutine monolis_visual_parted_mesh(nnode, node, nelem, nbase, elem, nodeid, elemid)
    implicit none
    integer(kint) :: i, j
    integer(kint) :: nnode, nelem, nbase, elem(:,:), nodeid(:), elemid(:)
    real(kdouble) :: node(:,:)
    character :: etype*6, output_dir*100

    call monolis_debug_header("monolis_visual_parted_mesh")

    output_dir = "visual/"
    call system('if [ ! -d visual ]; then (echo "** create visual"; mkdir -p visual); fi')

    open(20, file = trim(output_dir)//"mesh.parted.inp", status = "replace")
      write(20,"(5i12)") nnode, nelem, 1, 1, 0
      do i = 1, nnode
        write(20,"(i0,1p3e12.5)") i, node(1,i), node(2,i), node(3,i)
      enddo

      if(nbase == 3) etype = " tri  "
      if(nbase == 4) etype = " tet  "
      if(nbase == 8) etype = " hex  "
      if(nbase ==10) etype = " tet2 "

      do i = 1, nelem
        write(20,"(i0,i4,a,$)") i, 0, etype
        do j = 1, nbase
          write(20,"(i12,$)") elem(j,i)
        enddo
        write(20,*)""
      enddo

      write(20,"(a)")"1 1"
      write(20,"(a)")"node_domid, unknown"
      do i = 1, nnode
        write(20,"(i0,x,i0,x,i0,x,i0)") i, nodeid(i)!, bp_graph%elem_domid(i), in
      enddo

      write(20,"(a)")"1 1"
      write(20,"(a)")"elem_domid, unknown"
      do i = 1, nelem
        write(20,"(i0,x,i0,x,i0,x,i0)") i, elemid(i)!, bp_graph%elem_domid(i), in
      enddo
    close(20)
  end subroutine monolis_visual_parted_mesh

  subroutine monolis_output_graph_format(fname, nnode, nid, ebase_func, connectivity, shift)
    implicit none
    integer(kint) :: nnode
    integer(kint) :: i, in, j, jn, k, jS, shift, id
    integer(kint) :: nid(:), ebase_func(:), connectivity(:)
    integer(kint), allocatable :: t1(:), p1(:)
    character :: fname*100

    allocate(t1(nnode), source = 0)
    allocate(p1(nnode), source = 0)
    t1 = nid
    do i = 1, nnode
      p1(i) = i
    enddo
    call monolis_qsort_int_with_perm(t1, 1, nnode, p1)

    open(20, file = fname, status = "replace")
      write(20,"(i0)")nnode
      do i = 1, nnode
        jS = ebase_func(i)
        in = ebase_func(i+1) - ebase_func(i)
        write(20,"(i0,x,i0,$)") i+shift, in
        do j = 1, in
          k = connectivity(jS+j)
          call monolis_bsearch_int(t1, 1, nnode, k, id)
          write(20,"(x,i0,$)") p1(id)+shift
        enddo
        write(20,*)""
      enddo
    close(20)
  end subroutine monolis_output_graph_format

end module mod_monolis_io
