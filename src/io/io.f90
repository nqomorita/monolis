module mod_monolis_io
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_util
  use mod_monolis_mesh
  use mod_monolis_stdlib
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
      call monolis_input_mesh_node(fname, mesh%nnode_in, mesh%nnode, mesh%node, mesh%nid)

      fname = "elem.dat"
      call monolis_input_mesh_elem(fname, mesh%nelem, mesh%nbase_func, mesh%elem, mesh%eid)
    else
      fname = "node.dat"
      call monolis_input_mesh_node(fname, mesh%nnode_in, mesh%nnode, mesh%node)

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
    integer(kint) :: i, n_domain
    character :: cnum*5, output_dir*100, fname*100

    call monolis_debug_header("monolis_output_mesh")

    output_dir = "parted/"
    call system('if [ ! -d parted ]; then (echo "** create parted"; mkdir -p parted); fi')

    do i = 1, n_domain
      write(cnum,"(i0)") i-1
      fname = trim(output_dir)//"node.dat."//trim(cnum)
      call monolis_output_mesh_node(fname, node_list(i)%nnode, node_list(i)%nnode_in, mesh%node, node_list(i)%nid)

      fname = trim(output_dir)//"elem.dat."//trim(cnum)
      call monolis_output_mesh_elem_ref(fname, node_list(i)%nelem, mesh%nbase_func, node_list(i)%eid, &
        mesh%elem, node_list(i))

      fname = trim(output_dir)//"monolis.send."//trim(cnum)
      call monolis_output_mesh_comm(fname, comm(i)%send_n_neib, comm(i)%send_neib_pe, &
        comm(i)%send_index, comm(i)%send_item)

      fname = trim(output_dir)//"monolis.recv."//trim(cnum)
      call monolis_output_mesh_comm(fname, comm(i)%recv_n_neib, comm(i)%recv_neib_pe, &
        comm(i)%recv_index, comm(i)%recv_item)

      fname = trim(output_dir)//"node.id."//trim(cnum)
      call monolis_output_mesh_global_nid(fname, node_list(i)%nnode, mesh%nid, node_list(i)%nid)

      fname = trim(output_dir)//"elem.id."//trim(cnum)
      call monolis_output_mesh_global_eid(fname, node_list(i)%nelem, mesh%nid, node_list(i)%eid)
    enddo
  end subroutine monolis_output_mesh

  subroutine monolis_input_mesh_node(fname, nnode_in, nnode, node, nid)
    implicit none
    integer(kint) :: nnode_in, nnode, i, ierr, t(4)
    real(kdouble), allocatable :: node(:,:)
    integer(kint), optional, allocatable :: nid(:)
    character :: fname*100

    open(20, file = fname, status = "old")
      if(present(nid))then
        read(20,*, iostat = ierr) t(1), t(2), t(3)
        !read(20,*, iostat = ierr) nnode_in, nnode
        if(ierr /= 0)then
          rewind(20)
          read(20,*) nnode
          nnode_in = nnode
        else
          rewind(20)
          read(20,*) nnode_in, nnode
        endif

        allocate(node(3,nnode), source = 0.0d0)
        allocate(nid(nnode), source = 0)
        do i = 1, nnode
          read(20,*) nid(i), node(1,i), node(2,i), node(3,i)
        enddo
      else
        read(20,*, iostat = ierr) t(1), t(2)
        !read(20,*, iostat = ierr) nnode_in, nnode
        if(ierr /= 0)then
          rewind(20)
          read(20,*) nnode
          nnode_in = nnode
        else
          rewind(20)
          read(20,*) nnode_in, nnode
        endif

        allocate(node(3,nnode), source = 0.0d0)
        do i = 1, nnode
          read(20,*) node(1,i), node(2,i), node(3,i)
        enddo
      endif
    close(20)

    call monolis_debug_int("nnode", nnode)
  end subroutine monolis_input_mesh_node

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

  subroutine monolis_input_id(fname, id)
    implicit none
    integer(kint) :: nid, i, j
    integer(kint), allocatable :: id(:)
    character :: fname*100

    open(20, file = fname, status = "old")
      read(20,*) nid
      call monolis_debug_int("nid", nid)

      allocate(id(nid), source = 0)

      do i = 1, nid
        read(20,*) id(i)
      enddo
    close(20)
  end subroutine monolis_input_id

  subroutine monolis_input_condition(fname, ncond, icond, cond)
    implicit none
    integer(kint) :: ncond, i, j
    integer(kint), allocatable :: icond(:,:)
    real(kdouble), allocatable :: cond(:)
    character :: fname*100

    open(20, file = fname, status = "old")
      read(20,*) ncond
      call monolis_debug_int("ncond", ncond)

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

  subroutine monolis_output_mesh_node(fname, nnode, nnode_in, node, nid)
    implicit none
    integer(kint) :: nnode, nnode_in
    integer(kint) :: i, in
    integer(kint), optional :: nid(:)
    real(kdouble) :: node(:,:)
    character :: fname*100

    open(20, file = fname, status = "replace")
    if(present(nid))then
      write(20,"(i0,x,i0)")nnode_in, nnode
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
    integer(kint) :: i, j, in, jn, nelem, nbase, elem(:,:)
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

  subroutine monolis_output_mesh_elem_ref(fname, nelem, nbase, eid, elem, node_list)
    implicit none
    type(monolis_node_list) :: node_list
    integer(kint) :: i, j, in, jn, nelem, nbase, eid(:), elem(:,:)
    character :: fname*100

    open(20, file = fname, status = "replace")
      write(20,"(i0,x,i0)")nelem, nbase
      do i = 1, nelem
        jn = eid(i)
        do j = 1, nbase
          in = elem(j,jn)
          write(20,"(x,i0,$)") node_list%nid_perm(in)
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
        write(20,"(i0)") global_nid(in)
      enddo
    close(20)
  end subroutine monolis_output_mesh_global_nid

  subroutine monolis_output_mesh_global_eid(fname, nelem, global_eid, eid)
    implicit none
    integer(kint) :: i, in, nelem, global_eid(:), eid(:)
    character :: fname*100

    open(20, file = fname, status = "replace")
      write(20,"(i0,x,i0)")nelem
      do i = 1, nelem
        in = eid(i)
        write(20,"(i0,x,i0)") global_eid(in)
      enddo
    close(20)
  end subroutine monolis_output_mesh_global_eid

  subroutine monolis_output_mesh_comm(fname, n_neib, neib_pe, index, item)
    implicit none
    integer(kint) :: n_neib
    integer(kint) :: neib_pe(:)
    integer(kint) :: index(0:n_neib)
    integer(kint) :: item(:)
    integer(kint) :: i
    character :: fname*100

    open(20, file = fname, status = "replace")
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

      write(20,"(i0,x,i0)")0, 0 !> for non-over
    close(20)
  end subroutine monolis_output_mesh_comm

  subroutine monolis_visual_parted_mesh(nnode, node, nelem, nbase, elem, nodeid, elemid)
    implicit none
    integer(kint) :: i, j, in
    integer(kint) :: nnode, nelem, nbase, elem(:,:), nodeid(:), elemid(:)
    real(kdouble) :: node(:,:)
    character :: etype*6, output_dir*100, fname*100

    call monolis_debug_header("monolis_visual_parted_mesh")

    output_dir = "visual/"
    call system('if [ ! -d visual ]; then (echo "** create visual"; mkdir -p visual); fi')

    open(20, file = trim(output_dir)//"mesh.parted.inp", status = "replace")
      write(20,"(5i12)") nnode, nelem, 1, 1, 0
      do i = 1, nnode
        write(20,"(i0,1p3e12.5)") i, node(1,i), node(2,i), node(3,i)
      enddo

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

  subroutine monolis_get_part_arg(n_domain, is_format_id, is_overlap)
    implicit none
    integer(kint) :: i, count, n, n_domain
    character :: argc1*128, argc2*128
    logical :: is_format_id, is_overlap

    call monolis_debug_header("monolis_get_part_arg")

    call monolis_set_debug(.true.)

    count = iargc()
    if(count == 1)then
      call getarg(1, argc1)
      if(trim(argc1) == "-h")then
        write(*,"(a)")"-n {num of subdomain}: the number of subdomain"
        write(*,"(a)")"-t {N/O}: type of domain decomposition (N:non-overlapping, O:overlapping)"
        write(*,"(a)")"--with-id {Y/N}: node or elem id appears at the beginning of each line"
        write(*,"(a)")"-h: help"
        stop
      endif
    endif

    n_domain = 1
    is_overlap = .true.
    is_format_id = .false.

    if(mod(count,2) /= 0) stop "* monolis partitioner input arg error"
    do i = 1, count/2
      call getarg(2*i-1, argc1)
      call getarg(2*i  , argc2)
      if(trim(argc1) == "-n")then
        read(argc2,*) n
        n_domain = n

      elseif(trim(argc1) == "-t")then
        if(trim(argc2) == "O") is_overlap = .true.
        if(trim(argc2) == "N") is_overlap = .false.

      elseif(trim(argc1) == "--with-id")then
        if(trim(argc2) == "Y") is_format_id = .true.
        if(trim(argc2) == "N") is_format_id = .false.

      else
        write(*,"(a)")"* monolis input arg error"
        stop
      endif
    enddo

    call monolis_debug_int("n_domain", n_domain)
  end subroutine monolis_get_part_arg

end module mod_monolis_io
