module mod_monolis_p_refiner
  use mod_monolis_util
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mesh
  use mod_monolis_hash
  use mod_monolis_stdlib
  implicit none

  type(type_monolis_hash_tree) :: hash_tree

contains

  subroutine monolis_p_refine(mesh, mesh_ref)
    use mod_monolis_mesh
    implicit none
    type(monolis_mesh) :: mesh, mesh_ref
    integer(kint) :: itable(3,6), conn(mesh%nbase_func)
    integer(kint) :: i, eid, i1, i2, tmp, newid, nid(6), maxid
    character :: ckey*27
    logical :: is_exist, is_pushed
    integer(kint), allocatable :: tet2conn(:,:), tet2nid(:)
    real(kdouble), allocatable :: pnode(:,:)

    itable(1,1) = 1; itable(2,1) = 2
    itable(1,2) = 2; itable(2,2) = 3
    itable(1,3) = 3; itable(2,3) = 1
    itable(1,4) = 1; itable(2,4) = 4
    itable(1,5) = 2; itable(2,5) = 4
    itable(1,6) = 3; itable(2,6) = 4

    allocate(tet2conn(10,mesh%nelem), source = 0)
    allocate(tet2nid(6*mesh%nelem), source = 0)
    allocate(pnode(3,6*mesh%nelem), source = 0.0d0)

    call monolis_hash_init(hash_tree)

    newid = 0
    maxid = maxval(mesh%nid)
    do eid = 1, mesh%nelem
      call monolis_get_connectivity(mesh, eid, mesh%nbase_func, conn)

      do i = 1, 6
        i1 = conn(itable(1, i))
        i2 = conn(itable(2, i))
        ckey = get_key(i1, i2, 0)

        is_exist = .false.
        call monolis_hash_get(hash_tree, ckey, tmp, is_exist)
        if(is_exist)then
          nid(i) = tmp
        else
          newid = newid + 1
          nid(i) = mesh%nnode + newid
          tet2nid(newid) = maxid + newid
          pnode(:,newid) = get_mid_point_position(mesh, i1, i2)
          call monolis_hash_push(hash_tree, ckey, nid(i), is_pushed, is_exist)
        endif
      enddo

      tet2conn( 1,eid) = conn(1)
      tet2conn( 2,eid) = conn(2)
      tet2conn( 3,eid) = conn(3)
      tet2conn( 4,eid) = conn(4)
      tet2conn( 5,eid) = nid(1)
      tet2conn( 6,eid) = nid(2)
      tet2conn( 7,eid) = nid(3)
      tet2conn( 8,eid) = nid(4)
      tet2conn( 9,eid) = nid(5)
      tet2conn(10,eid) = nid(6)
    enddo

    !> construct p-refined mesh
    call alloc_mesh(mesh, mesh_ref, newid)
    call copy_node(mesh, mesh_ref, pnode)
    call copy_elem(mesh_ref, tet2conn)

    call monolis_hash_finalize(hash_tree)
  end subroutine monolis_p_refine

  subroutine alloc_mesh(mesh, mesh_ref, newid)
    use mod_monolis_mesh
    implicit none
    type(monolis_mesh) :: mesh, mesh_ref
    integer(kint) :: newid

    mesh_ref%nnode = mesh%nnode + newid
    mesh_ref%nelem = mesh%nelem
    mesh_ref%nbase_func = 10
    allocate(mesh_ref%node( 3,mesh_ref%nnode), source = 0.0d0)
    allocate(mesh_ref%elem(10,mesh_ref%nelem), source = 0)
  end subroutine alloc_mesh

  subroutine copy_node(mesh, mesh_ref, pnode)
    use mod_monolis_mesh
    implicit none
    type(monolis_mesh) :: mesh, mesh_ref
    integer(kint) :: i, j, in
    real(kdouble) :: pnode(:,:)

    do i = 1, mesh%nnode
      do j = 1, 3
        mesh_ref%node(j,i) = mesh%node(j,i)
      enddo
    enddo

    do i = mesh%nnode+1, mesh_ref%nnode
      in = i - mesh%nnode
      do j = 1, 3
        mesh_ref%node(j,i) = pnode(j,in)
      enddo
    enddo
  end subroutine copy_node

  subroutine copy_elem(mesh_ref, tet2conn)
    use mod_monolis_mesh
    implicit none
    type(monolis_mesh) :: mesh_ref
    integer(kint) :: i, j, tet2conn(:,:)

    do i = 1, mesh_ref%nelem
      do j = 1, 10
        mesh_ref%elem(j,i) = tet2conn(j,i)
      enddo
    enddo
  end subroutine copy_elem

  function get_mid_point_position(mesh, i1, i2)
    use mod_monolis_mesh
    implicit none
    type(monolis_mesh) :: mesh
    integer(kint) :: i1, i2
    real(kdouble) :: get_mid_point_position(3)

    get_mid_point_position(1) = (mesh%node(1,i1) + mesh%node(1,i2))/2.0d0
    get_mid_point_position(2) = (mesh%node(2,i1) + mesh%node(2,i2))/2.0d0
    get_mid_point_position(3) = (mesh%node(3,i1) + mesh%node(3,i2))/2.0d0
  end function

  function get_key(i1, i2, i3)
    implicit none
    integer(kint) :: i1, i2, i3, array(3)
    character :: c1*9, c2*9, c3*9, get_key*27

    array(1) = i1
    array(2) = i2
    array(3) = i3
    call monolis_qsort_int(array, 1, 3)

    write(c1,"(i9.9)")array(1)
    write(c2,"(i9.9)")array(2)
    write(c3,"(i9.9)")array(3)
    get_key = c1//c2//c3
  end function get_key

end module mod_monolis_p_refiner
