module mod_monolis_sparse_util
  use mod_monolis_util
  use mod_monolis_stdlib
  use mod_monolis_graph
  implicit none

contains

!> to be deleted >>>
  subroutine monolis_get_nonzero_pattern(monolis, nnode, nbase_func, ndof, nelem, elem)
    use iso_c_binding
    implicit none
    type(monolis_structure) :: monolis
    integer(kint) :: nnode, nbase_func, ndof, nelem, elem(:,:)
    call monolis_get_nonzero_pattern_by_simple_mesh(monolis, nnode, nbase_func, ndof, nelem, elem)
  end subroutine monolis_get_nonzero_pattern
!> <<< to be deleted

  subroutine monolis_get_nonzero_pattern_by_simple_mesh(monolis, nnode, nbase_func, ndof, nelem, elem)
    use iso_c_binding
    implicit none
    type(monolis_structure) :: monolis
    integer(kint) :: nnode, nbase_func, ndof, nelem, elem(:,:)
    integer(kint), pointer :: ebase_func(:), connectivity(:)
    integer(c_int), pointer :: index(:), item(:)

    call monolis_convert_mesh_to_connectivity &
      & (nelem, nbase_func, elem, ebase_func, connectivity)

    call monolis_convert_connectivity_to_nodal &
      & (nnode, nelem, ebase_func, connectivity, index, item)

    call monolis_get_nonzero_pattern_by_nodal_graph &
      & (monolis, nnode, ndof, index, item)

    deallocate(ebase_func); nullify(ebase_func)
    deallocate(connectivity); nullify(connectivity)
    deallocate(index); nullify(index)
    deallocate(item); nullify(item)
  end subroutine monolis_get_nonzero_pattern_by_simple_mesh

  subroutine monolis_get_nonzero_pattern_with_arbitrary_dof &
    (monolis, nnode, nbase_func, n_dof_list, nelem, elem)
    use iso_c_binding
    implicit none
    type(monolis_structure) :: monolis
    integer(kint) :: nnode, nbase_func, n_dof_list(:), nelem, elem(:,:)
    integer(kint), pointer :: ebase_func(:), connectivity(:)
    integer(c_int), pointer :: index(:), item(:)

    call monolis_convert_mesh_to_connectivity &
      & (nelem, nbase_func, elem, ebase_func, connectivity)

    call monolis_convert_connectivity_to_nodal &
      & (nnode, nelem, ebase_func, connectivity, index, item)

     call monolis_get_nonzero_pattern_by_nodal_graph_with_arbitrary_dof &
       & (monolis, nnode, n_dof_list, index, item)
  end subroutine monolis_get_nonzero_pattern_with_arbitrary_dof

  subroutine monolis_get_nonzero_pattern_by_connectivity &
      & (monolis, nnode, ndof, nelem, ebase_func, connectivity)
    use iso_c_binding
    implicit none
    type(monolis_structure) :: monolis
    integer(kint) :: nnode, ndof, nelem
    integer(kint), pointer :: ebase_func(:), connectivity(:)
    integer(c_int), pointer :: index(:), item(:)

    call monolis_convert_connectivity_to_nodal &
      & (nnode, nelem, ebase_func, connectivity, index, item)

     call monolis_get_nonzero_pattern_by_nodal_graph &
      & (monolis, nnode, ndof, index, item)
  end subroutine monolis_get_nonzero_pattern_by_connectivity

  subroutine monolis_get_nonzero_pattern_by_nodal_graph(monolis, nnode, ndof, index, item)
    use iso_c_binding
    implicit none
    type(monolis_structure) :: monolis
    integer(kint) :: nnode, ndof
    integer(kint) :: i, j, nz, jS, jE
    integer(c_int), pointer :: index(:), item(:)

    monolis%MAT%N = nnode
    monolis%MAT%NP = nnode
    monolis%MAT%NDOF = ndof
    allocate(monolis%MAT%X(ndof*nnode), source = 0.0d0)
    allocate(monolis%MAT%B(ndof*nnode), source = 0.0d0)
    allocate(monolis%MAT%index(0:nnode), source = 0)
    do i = 1, nnode
      monolis%MAT%index(i) = index(i+1) + i
    enddo

    nz = monolis%MAT%index(nnode)
    monolis%MAT%NZ = nz
    allocate(monolis%MAT%A(ndof*ndof*nz), source = 0.0d0)
    allocate(monolis%MAT%item(nz), source = 0)
    do i = 1, nnode
      jS = monolis%MAT%index(i-1) + 1
      jE = monolis%MAT%index(i)
      monolis%MAT%item(jS) = i
      do j = jS+1, jE
        monolis%MAT%item(j) = item(j-i)
      enddo
      call monolis_qsort_int(monolis%MAT%item(jS:jE), 1, jE - jS + 1)
    enddo

    allocate(monolis%MAT%indexR(0:nnode), source = 0)
    allocate(monolis%MAT%itemR(nz), source = 0)
    allocate(monolis%MAT%permR(nz), source = 0)

    call monolis_get_CRR_format(monolis%MAT%N, monolis%MAT%N, nz, &
      & monolis%MAT%index, monolis%MAT%item, &
      & monolis%MAT%indexR, monolis%MAT%itemR, monolis%MAT%permR)
  end subroutine monolis_get_nonzero_pattern_by_nodal_graph

  subroutine monolis_get_nonzero_pattern_by_nodal_graph_with_arbitrary_dof &
    (monolis, nnode, n_dof_list, index, item)
    use iso_c_binding
    implicit none
    type(monolis_structure) :: monolis
    integer(kint) :: nnode, n_dof_list(:)
    integer(kint) :: i, j, k, nz, jS, jE, kS, kE
    integer(kint) :: total_dof, in, jn, kn, nrow, ncol, l
    integer(kint), allocatable :: n_dof_index(:)
    integer(c_int), pointer :: index(:), item(:)

    total_dof = 0
    do i = 1, nnode
      total_dof = total_dof + n_dof_list(i)
    enddo

    allocate(n_dof_index(nnode), source = 0)
    call monolis_get_n_dof_index(nnode, n_dof_list, n_dof_index)

    monolis%MAT%N = total_dof
    monolis%MAT%NP = total_dof
    monolis%MAT%NDOF = 1
    allocate(monolis%MAT%X(total_dof), source = 0.0d0)
    allocate(monolis%MAT%B(total_dof), source = 0.0d0)
    allocate(monolis%MAT%index(0:total_dof), source = 0)

    !> count nz
    nz = 0
    do i = 1, nnode
      jS = index(i) + 1
      jE = index(i+1)
      nz = nz + n_dof_list(i)*n_dof_list(i)
      do j = jS, jE
        jn = item(j)
        nz = nz +  n_dof_list(i)*n_dof_list(jn)
      enddo
    enddo

    monolis%MAT%NZ = nz
    allocate(monolis%MAT%A(nz), source = 0.0d0)
    allocate(monolis%MAT%item(nz), source = 0)

    !> construct index and item
    in = 0
    ncol = 0
    do i = 1, nnode
      jS = index(i) + 1
      jE = index(i+1)
      do k = 1, n_dof_list(i)
        ncol = ncol + 1
        nrow = n_dof_list(i)
        do j = 1, n_dof_list(i)
          in = in + 1
          monolis%MAT%item(in) = n_dof_index(i) + j
        enddo
        do j = jS, jE
          jn = item(j)
          nrow = nrow + n_dof_list(jn)
          do l = 1, n_dof_list(jn)
            in = in + 1
            monolis%MAT%item(in) = n_dof_index(jn) + l
          enddo
        enddo
        monolis%MAT%index(ncol) = monolis%MAT%index(ncol-1) + nrow

        kS = monolis%MAT%index(ncol-1) + 1
        kE = monolis%MAT%index(ncol)
        call monolis_qsort_int(monolis%MAT%item(kS:kE), 1, kE-kS+1)
      enddo
    enddo

    allocate(monolis%MAT%indexR(0:total_dof), source = 0)
    allocate(monolis%MAT%itemR(nz), source = 0)
    allocate(monolis%MAT%permR(nz), source = 0)

    call monolis_get_CRR_format(monolis%MAT%N, monolis%MAT%N, nz, &
      & monolis%MAT%index, monolis%MAT%item, &
      & monolis%MAT%indexR, monolis%MAT%itemR, monolis%MAT%permR)

    nullify(index)
    nullify(item)
  end subroutine monolis_get_nonzero_pattern_by_nodal_graph_with_arbitrary_dof

  subroutine monolis_get_n_dof_index(n_node, n_dof_list, n_dof_index)
    implicit none
    integer(kint), intent(in) :: n_node, n_dof_list(:)
    integer(kint) :: i, n_dof_index(:)

    do i = 1, n_node - 1
      n_dof_index(i+1) = n_dof_index(i) + n_dof_list(i)
    enddo
  end subroutine monolis_get_n_dof_index

  !> setter
  subroutine monolis_set_scalar_to_sparse_matrix(monolis, i, j, sub_i, sub_j, val)
    implicit none
    type(monolis_structure) :: monolis
    integer(kint), intent(in) :: i, j, sub_i, sub_j
    real(kdouble), intent(in) :: val

    call monolis_sparse_matrix_set_value(monolis%MAT%index, monolis%MAT%item, monolis%MAT%A, &
      & monolis%MAT%ndof, i, j, sub_i, sub_j, val)
  end subroutine monolis_set_scalar_to_sparse_matrix

  subroutine monolis_set_matrix_to_sparse_matrix(monolis, i, j, val)
    implicit none
    type(monolis_structure) :: monolis
    integer(kint), intent(in) :: i, j
    real(kdouble), intent(in) :: val(:,:)

    call monolis_sparse_matrix_set_block(monolis%MAT%index, monolis%MAT%item, monolis%MAT%A, &
      & monolis%MAT%ndof, i, j, val)
  end subroutine monolis_set_matrix_to_sparse_matrix

  !> getter
  subroutine monolis_get_scalar_from_sparse_matrix(monolis, i, j, sub_i, sub_j, val, is_find)
    implicit none
    type(monolis_structure) :: monolis
    integer(kint), intent(in) :: i, j, sub_i, sub_j
    real(kdouble) :: val
    logical :: is_find

    call monolis_sparse_matrix_get_value(monolis%MAT%index, monolis%MAT%item, monolis%MAT%A, &
      & monolis%MAT%ndof, i, j, sub_i, sub_j, val, is_find)
  end subroutine monolis_get_scalar_from_sparse_matrix

  subroutine monolis_sparse_matrix_get_value(index, item, A, ndof, ci, cj, csub_i, csub_j, val, is_find)
    implicit none
    integer(kint), intent(in) :: ndof
    integer(kint), intent(in) :: index(0:), item(:), ci, cj, csub_i, csub_j
    real(kdouble), intent(inout) :: A(:)
    real(kdouble), intent(out) :: val
    integer(kint) :: j, jn, im, jS, jE, NDOF2
    logical :: is_find

    val = 0.0d0
    is_find = .false.

    NDOF2 = ndof*ndof
    if(ndof < csub_i) return
    if(ndof < csub_j) return

    jS = index(ci-1) + 1
    jE = index(ci)
    do j = jS, jE
      jn = item(j)
      if(jn == cj)then
        im = NDOF2*(j-1) + ndof*(csub_i-1) + csub_j
        val = A(im)
        is_find = .true.
        return
      endif
    enddo
  end subroutine monolis_sparse_matrix_get_value

  !> adder
  subroutine monolis_add_scalar_to_sparse_matrix(monolis, i, j, sub_i, sub_j, val)
    implicit none
    type(monolis_structure) :: monolis
    integer(kint), intent(in) :: i, j, sub_i, sub_j
    real(kdouble), intent(in) :: val

    call monolis_sparse_matrix_add_value(monolis%MAT%index, monolis%MAT%item, monolis%MAT%A, &
      & monolis%MAT%ndof, i, j, sub_i, sub_j, val)
  end subroutine monolis_add_scalar_to_sparse_matrix

  subroutine monolis_add_matrix_to_sparse_matrix(monolis, nbase_func, connectivity, stiff)
    implicit none
    type(monolis_structure) :: monolis
    integer(kint), intent(in) :: nbase_func, connectivity(nbase_func)
    real(kdouble), intent(in) :: stiff(:,:)

    call monolis_sparse_matrix_add_matrix(monolis%MAT%index, monolis%MAT%item, monolis%MAT%A, &
      & nbase_func, nbase_func, monolis%MAT%ndof, connectivity, connectivity, stiff)
  end subroutine monolis_add_matrix_to_sparse_matrix

  subroutine monolis_add_matrix_to_sparse_matrix_offdiag(monolis, n1, n2, c1, c2, stiff)
    implicit none
    type(monolis_structure) :: monolis
    integer(kint), intent(in) :: n1, n2, c1(n1), c2(n2)
    real(kdouble), intent(in) :: stiff(:,:)

    call monolis_sparse_matrix_add_matrix(monolis%MAT%index, monolis%MAT%item, monolis%MAT%A, &
      & n1, n2, monolis%MAT%ndof, c1, c2, stiff)
  end subroutine monolis_add_matrix_to_sparse_matrix_offdiag

  subroutine monolis_sparse_matrix_add_matrix(index, item, A, n1, n2, ndof, e1t, e2t, stiff)
    implicit none
    integer(kint), intent(in) :: n1, n2, ndof
    integer(kint), intent(in) :: index(0:), item(:), e1t(n1), e2t(n2)
    real(kdouble), intent(inout) :: A(:)
    real(kdouble), intent(in) :: stiff(n1*ndof,n2*ndof)
    integer(kint) :: e1(n1), e2(n2)
    integer(kint) :: i, j, k, in, jn, im, jS, jE, i2, j2, i1, j1, NDOF2
    integer(kint) :: eperm1(n1), eperm2(n2)
    real(kdouble) :: temp(n1*ndof,n2*ndof)

    NDOF2 = ndof*ndof

    e1 = e1t
    do i = 1, n1
      eperm1(i) = i
    enddo
    call monolis_qsort_int_with_perm(e1, 1, n1, eperm1)

    e2 = e2t
    do i = 1, n2
      eperm2(i) = i
    enddo
    call monolis_qsort_int_with_perm(e2, 1, n2, eperm2)

    temp = 0.0d0
    do i = 1, n2
      i1 = eperm2(i)
      do j = 1, n1
        j1 = eperm1(j)
        do i2 = 1, ndof
          do j2 = 1, ndof
            temp(ndof*(i-1)+i2, ndof*(j-1)+j2) = stiff(ndof*(j1-1)+j2, ndof*(i1-1)+i2)
          enddo
        enddo
      enddo
    enddo

    do i = 1, n1
      in = e1(i)
      jS = index(in-1) + 1
      jE = index(in)
      aa:do j = 1, n2
        do k = jS, jE
          jn = item(k)
          if(jn == e2(j))then
            do i1 = 1, ndof
            do i2 = 1, ndof
              im = NDOF2*(k-1) + ndof*(i1-1) + i2
              A(im) = A(im) + temp(ndof*(j-1)+i2, ndof*(i-1)+i1)
            enddo
            enddo
            jS = k + 1
            cycle aa
          endif
        enddo
      call monolis_stop_by_matrix_assemble(e1(i), e2(j))
      enddo aa
    enddo
  end subroutine monolis_sparse_matrix_add_matrix

  subroutine monolis_sparse_matrix_set_block(index, item, A, ndof, e1, e2, val)
    implicit none
    integer(kint), intent(in) :: e1, e2, ndof
    integer(kint), intent(in) :: index(0:), item(:)
    real(kdouble), intent(inout) :: A(:)
    real(kdouble), intent(in) :: val(ndof,ndof)
    integer(kint) :: i, j, k, in, jn, im, jS, jE, i2, i1, NDOF2

    NDOF2 = ndof*ndof
    in = e1
    jS = index(in-1) + 1
    jE = index(in)
    do k = jS, jE
      jn = item(k)
      if(jn == e2)then
        do i1 = 1, ndof
        do i2 = 1, ndof
          im = NDOF2*(k-1) + ndof*(i1-1) + i2
          A(im) = val(i2, i1)
        enddo
        enddo
        return
      endif
    enddo
    call monolis_stop_by_matrix_assemble(e1, e2)
  end subroutine monolis_sparse_matrix_set_block

  subroutine monolis_sparse_matrix_add_value(index, item, A, ndof, ci, cj, csub_i, csub_j, val)
    implicit none
    integer(kint), intent(in) :: ndof
    integer(kint), intent(in) :: index(0:), item(:), ci, cj, csub_i, csub_j
    real(kdouble), intent(inout) :: A(:)
    real(kdouble), intent(in) :: val
    integer(kint) :: j, jn, im, jS, jE, NDOF2
    character :: cerr*128

    NDOF2 = ndof*ndof
    if(ndof < csub_i) call monolis_stop_by_submatrix_access(ndof, csub_i)
    if(ndof < csub_j) call monolis_stop_by_submatrix_access(ndof, csub_j)

    jS = index(ci-1) + 1
    jE = index(ci)
    do j = jS, jE
      jn = item(j)
      if(jn == cj)then
        im = NDOF2*(j-1) + ndof*(csub_i-1) + csub_j
        A(im) = A(im) + val
        return
      endif
    enddo

    call monolis_stop_by_matrix_assemble(ci, cj)
  end subroutine monolis_sparse_matrix_add_value

  subroutine monolis_sparse_matrix_set_value(index, item, A, ndof, ci, cj, csub_i, csub_j, val)
    implicit none
    integer(kint), intent(in) :: ndof
    integer(kint), intent(in) :: index(0:), item(:), ci, cj, csub_i, csub_j
    real(kdouble), intent(inout) :: A(:)
    real(kdouble), intent(in) :: val
    integer(kint) :: j, jn, im, jS, jE, NDOF2

    NDOF2 = ndof*ndof
    if(ndof < csub_i) call monolis_stop_by_submatrix_access(ndof, csub_i)
    if(ndof < csub_j) call monolis_stop_by_submatrix_access(ndof, csub_j)

    jS = index(ci-1) + 1
    jE = index(ci)
    do j = jS, jE
      jn = item(j)
      if(jn == cj)then
        im = NDOF2*(j-1) + ndof*(csub_i-1) + csub_j
        A(im) = val
        return
      endif
      call monolis_stop_by_matrix_assemble(ci, cj)
    enddo
  end subroutine monolis_sparse_matrix_set_value

  !> CSR data setter
  subroutine monolis_set_matrix_BCSR(monolis, N, NP, NDOF, NZ, A, index, item)
    implicit none
    type(monolis_structure) :: monolis
    integer(kint), intent(in) :: N, NP, NDOF, NZ
    integer(kint), intent(in) :: index(0:), item(:)
    real(kdouble), intent(in) :: A(:)
    integer(kint) :: i

    if(associated(monolis%MAT%index))  deallocate(monolis%MAT%index)
    if(associated(monolis%MAT%item))   deallocate(monolis%MAT%item)
    if(associated(monolis%MAT%A))      deallocate(monolis%MAT%A)
    if(associated(monolis%MAT%X))      deallocate(monolis%MAT%X)
    if(associated(monolis%MAT%B))      deallocate(monolis%MAT%B)
    if(associated(monolis%MAT%indexR)) deallocate(monolis%MAT%indexR)
    if(associated(monolis%MAT%itemR))  deallocate(monolis%MAT%itemR)
    if(associated(monolis%MAT%permR))  deallocate(monolis%MAT%permR)

    monolis%MAT%N = N
    monolis%MAT%NP = NP
    monolis%MAT%NDOF = NDOF

    allocate(monolis%MAT%index(0:NP), source = 0)
    allocate(monolis%MAT%item(NZ), source = 0)
    allocate(monolis%MAT%A(NDOF*NDOF*NZ), source = 0.0d0)
    allocate(monolis%MAT%X(NDOF*NP), source = 0.0d0)
    allocate(monolis%MAT%B(NDOF*NP), source = 0.0d0)

    do i = 1, NP
      monolis%MAT%index(i) = index(i)
    enddo
    do i = 1, NZ
      monolis%MAT%item(i) = item(i)
    enddo
    do i = 1, NDOF*NDOF*NZ
      monolis%MAT%A(i) = A(i)
    enddo

    allocate(monolis%MAT%indexR(0:NP), source = 0)
    allocate(monolis%MAT%itemR(NZ), source = 0)
    allocate(monolis%MAT%permR(NZ), source = 0)

    call monolis_get_CRR_format(N, NP, NZ, &
      & monolis%MAT%index, monolis%MAT%item, &
      & monolis%MAT%indexR, monolis%MAT%itemR, monolis%MAT%permR)
  end subroutine monolis_set_matrix_BCSR

  !> boundary condition
  subroutine monolis_set_Dirichlet_bc(monolis, B, node_id, ndof_bc, val)
    implicit none
    type(monolis_structure) :: monolis
    integer(kint), intent(in) :: node_id, ndof_bc
    real(kdouble), intent(in) :: val
    real(kdouble) :: B(:)

    call monolis_sparse_matrix_add_bc(monolis%MAT%index, monolis%MAT%item, monolis%MAT%A, B, &
      & monolis%MAT%indexR, monolis%MAT%itemR, monolis%MAT%permR, &
      & monolis%MAT%ndof, node_id, ndof_bc, val)
  end subroutine monolis_set_Dirichlet_bc

  subroutine monolis_sparse_matrix_add_bc(index, item, A, B, indexR, itemR, permA, &
    & ndof, nnode, idof, val)
    implicit none
    integer(kint), intent(in) :: nnode, ndof, idof
    integer(kint), intent(in) :: index(0:), item(:), indexR(0:), itemR(:), permA(:)
    real(kdouble), intent(inout) :: A(:), B(:)
    real(kdouble), intent(in) :: val
    integer(kint) :: j, k, jn, kn, jS, jE, NDOF2
    logical :: is_add

    if(ndof < idof) call monolis_stop_by_submatrix_access(ndof, idof)

    is_add = .false.
    NDOF2 = ndof*ndof

    jS = indexR(nnode-1) + 1
    jE = indexR(nnode)
    do j = jS, jE
      jn = itemR(j)
      kn = permA(j)
      do k = 1, ndof
        B(ndof*(jn-1)+k) = B(ndof*(jn-1)+k) - val*A(NDOF2*(kn-1) + ndof*(k-1) + idof)
        A(NDOF2*(kn-1) + ndof*(k-1) + idof) = 0.0d0
      enddo
    enddo

    jS = index(nnode-1) + 1
    jE = index(nnode)
    do j = jS, jE
      do k = 1, ndof
        A(NDOF2*(j-1) + ndof*(idof-1) + k) = 0.0d0
      enddo

      jn = item(j)
      if(jn == nnode)then
        A(NDOF2*(j-1) + (ndof+1)*(idof-1) + 1) = 1.0d0
        is_add = .true.
      endif
    enddo

    if(.not. is_add) stop "error: not find a diagonal element in monolis_sparse_matrix_add_bc"

    B(ndof*nnode-ndof+idof) = val
  end subroutine monolis_sparse_matrix_add_bc

  subroutine monolis_get_CRR_format(NC, NR, NZ, index, item, indexR, itemR, permR)
    implicit none
    integer(kint), intent(in) :: NC, NR, NZ, index(0:), item(:)
    integer(kint), pointer :: indexR(:), itemR(:), permR(:)
    integer(kint), allocatable :: temp(:)
    integer(kint) :: i, j, in, jS, jE, m, p

    allocate(temp(NR), source = 0)
    do i = 1, NC
      jS = index(i-1) + 1
      jE = index(i)
      do j = jS, jE
        in = item(j)
        temp(in) = temp(in) + 1
      enddo
    enddo

    do i = 1, NR
      indexR(i) = indexR(i-1) + temp(i)
    enddo

    temp = 0
    do i = 1, NC
      jS = index(i-1) + 1
      jE = index(i)
      do j = jS, jE
        in = item(j)
        m = indexR(in-1)
        temp(in) = temp(in) + 1
        p = temp(in)
        itemR(m + p) = i
        permR(m + p) = j
      enddo
    enddo

    do i = 1, NR
      jS = indexR(i-1) + 1
      jE = indexR(i)
      call monolis_qsort_int(itemR(jS:jE), 1, jE-jS+1)
    enddo
    deallocate(temp)
  end subroutine monolis_get_CRR_format

  subroutine monolis_stop_by_matrix_assemble(ci, cj)
    integer(kint), intent(in) :: ci, cj
    write(*,"(a,i0,a,i0,a)") "error: The non-zero element at (", ci, ", ", cj, &
      & ") is not allocated. The value is not accessible."
    stop
  end subroutine monolis_stop_by_matrix_assemble

  subroutine monolis_stop_by_submatrix_access(ndof, sub_dof)
    integer(kint), intent(in) :: ndof, sub_dof
    write(*,"(a)")   "error: set value greater than the DoF of submatrix."
    write(*,"(a,i8)")"       the DoF of submatrix: ", ndof
    write(*,"(a,i8)")"       value:                ", sub_dof
    stop
  end subroutine monolis_stop_by_submatrix_access

end module mod_monolis_sparse_util
