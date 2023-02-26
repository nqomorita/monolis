module mod_monolis_spmat_handler
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc

  implicit none

contains

  !> setter
  subroutine monolis_set_scalar_to_sparse_matrix(monolis, i, j, sub_i, sub_j, val)
    implicit none
    type(monolis_structure) :: monolis
    integer(kint), intent(in) :: i, j, sub_i, sub_j
    real(kdouble), intent(in) :: val

!    call monolis_sparse_matrix_set_value(monolis%MAT%index, monolis%MAT%item, monolis%MAT%A, &
!      & monolis%MAT%ndof, i, j, sub_i, sub_j, val)
  end subroutine monolis_set_scalar_to_sparse_matrix

  subroutine monolis_set_matrix_to_sparse_matrix(monolis, i, j, val)
    implicit none
    type(monolis_structure) :: monolis
    integer(kint), intent(in) :: i, j
    real(kdouble), intent(in) :: val(:,:)

!    call monolis_sparse_matrix_set_block(monolis%MAT%index, monolis%MAT%item, monolis%MAT%A, &
!      & monolis%MAT%ndof, i, j, val)
  end subroutine monolis_set_matrix_to_sparse_matrix

  !> getter
  subroutine monolis_get_scalar_from_sparse_matrix(monolis, i, j, sub_i, sub_j, val, is_find)
    implicit none
    type(monolis_structure) :: monolis
    integer(kint), intent(in) :: i, j, sub_i, sub_j
    real(kdouble) :: val
    logical :: is_find

!    call monolis_sparse_matrix_get_value(monolis%MAT%index, monolis%MAT%item, monolis%MAT%A, &
!      & monolis%MAT%ndof, i, j, sub_i, sub_j, val, is_find)
  end subroutine monolis_get_scalar_from_sparse_matrix

  subroutine monolis_sparse_matrix_get_value(index, item, A, ndof, ci, cj, csub_i, csub_j, val, is_find)
    implicit none
    integer(kint), intent(in) :: ndof
    integer(kint), intent(in) :: index(0:), item(:), ci, cj, csub_i, csub_j
    real(kdouble), intent(inout) :: A(:)
    real(kdouble), intent(out) :: val
    integer(kint) :: j, jn, im, jS, jE, NDOF2
    logical :: is_find

!    val = 0.0d0
!    is_find = .false.

!    NDOF2 = ndof*ndof
!    if(ndof < csub_i) return
!    if(ndof < csub_j) return

!    jS = index(ci-1) + 1
!    jE = index(ci)
!    do j = jS, jE
!      jn = item(j)
!      if(jn == cj)then
!        im = NDOF2*(j-1) + ndof*(csub_i-1) + csub_j
!        val = A(im)
!        is_find = .true.
!        return
!      endif
!    enddo
  end subroutine monolis_sparse_matrix_get_value

  !> adder
  subroutine monolis_add_scalar_to_sparse_matrix(monolis, i, j, sub_i, sub_j, val)
    implicit none
    type(monolis_structure) :: monolis
    integer(kint), intent(in) :: i, j, sub_i, sub_j
    real(kdouble), intent(in) :: val

!    call monolis_sparse_matrix_add_value(monolis%MAT%index, monolis%MAT%item, monolis%MAT%A, &
!      & monolis%MAT%ndof, i, j, sub_i, sub_j, val)
  end subroutine monolis_add_scalar_to_sparse_matrix

  subroutine monolis_add_matrix_to_sparse_matrix(monolis, nbase_func, connectivity, stiff)
    implicit none
    type(monolis_structure) :: monolis
    integer(kint), intent(in) :: nbase_func, connectivity(nbase_func)
    real(kdouble), intent(in) :: stiff(:,:)

!    call monolis_sparse_matrix_add_matrix(monolis%MAT%index, monolis%MAT%item, monolis%MAT%A, &
!      & nbase_func, nbase_func, monolis%MAT%ndof, connectivity, connectivity, stiff)
  end subroutine monolis_add_matrix_to_sparse_matrix

  subroutine monolis_add_matrix_to_sparse_matrix_offdiag(monolis, n1, n2, c1, c2, stiff)
    implicit none
    type(monolis_structure) :: monolis
    integer(kint), intent(in) :: n1, n2, c1(n1), c2(n2)
    real(kdouble), intent(in) :: stiff(:,:)

!    call monolis_sparse_matrix_add_matrix(monolis%MAT%index, monolis%MAT%item, monolis%MAT%A, &
!      & n1, n2, monolis%MAT%ndof, c1, c2, stiff)
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

!    NDOF2 = ndof*ndof
!
!    e1 = e1t
!    do i = 1, n1
!      eperm1(i) = i
!    enddo
!    call monolis_qsort_I_2d(e1, eperm1, 1, n1)
!
!    e2 = e2t
!    do i = 1, n2
!      eperm2(i) = i
!    enddo
!    call monolis_qsort_I_2d(e2, eperm2, 1, n2)
!
!    temp = 0.0d0
!    do i = 1, n2
!      i1 = eperm2(i)
!      do j = 1, n1
!        j1 = eperm1(j)
!        do i2 = 1, ndof
!          do j2 = 1, ndof
!            temp(ndof*(i-1)+i2, ndof*(j-1)+j2) = stiff(ndof*(j1-1)+j2, ndof*(i1-1)+i2)
!          enddo
!        enddo
!      enddo
!    enddo
!
!    do i = 1, n1
!      in = e1(i)
!      jS = index(in-1) + 1
!      jE = index(in)
!      aa:do j = 1, n2
!        do k = jS, jE
!          jn = item(k)
!          if(jn == e2(j))then
!            do i1 = 1, ndof
!            do i2 = 1, ndof
!              im = NDOF2*(k-1) + ndof*(i1-1) + i2
!              A(im) = A(im) + temp(ndof*(j-1)+i2, ndof*(i-1)+i1)
!            enddo
!            enddo
!            jS = k + 1
!            cycle aa
!          endif
!        enddo
!      call monolis_stop_by_matrix_assemble(e1(i), e2(j))
!      enddo aa
!    enddo
  end subroutine monolis_sparse_matrix_add_matrix

  subroutine monolis_sparse_matrix_set_block(index, item, A, ndof, e1, e2, val)
    implicit none
    integer(kint), intent(in) :: e1, e2, ndof
    integer(kint), intent(in) :: index(0:), item(:)
    real(kdouble), intent(inout) :: A(:)
    real(kdouble), intent(in) :: val(ndof,ndof)
    integer(kint) :: i, j, k, in, jn, im, jS, jE, i2, i1, NDOF2

!    NDOF2 = ndof*ndof
!    in = e1
!    jS = index(in-1) + 1
!    jE = index(in)
!    do k = jS, jE
!      jn = item(k)
!      if(jn == e2)then
!        do i1 = 1, ndof
!        do i2 = 1, ndof
!          im = NDOF2*(k-1) + ndof*(i1-1) + i2
!          A(im) = val(i2, i1)
!        enddo
!        enddo
!        return
!      endif
!    enddo
!    call monolis_stop_by_matrix_assemble(e1, e2)
  end subroutine monolis_sparse_matrix_set_block

  subroutine monolis_sparse_matrix_add_value(index, item, A, ndof, ci, cj, csub_i, csub_j, val)
    implicit none
    integer(kint), intent(in) :: ndof
    integer(kint), intent(in) :: index(0:), item(:), ci, cj, csub_i, csub_j
    real(kdouble), intent(inout) :: A(:)
    real(kdouble), intent(in) :: val
    integer(kint) :: j, jn, im, jS, jE, NDOF2
    character :: cerr*128

!    NDOF2 = ndof*ndof
!    if(ndof < csub_i) call monolis_stop_by_submatrix_access(ndof, csub_i)
!    if(ndof < csub_j) call monolis_stop_by_submatrix_access(ndof, csub_j)
!
!    jS = index(ci-1) + 1
!    jE = index(ci)
!    do j = jS, jE
!      jn = item(j)
!      if(jn == cj)then
!        im = NDOF2*(j-1) + ndof*(csub_i-1) + csub_j
!        A(im) = A(im) + val
!        return
!      endif
!    enddo
!
!    call monolis_stop_by_matrix_assemble(ci, cj)
  end subroutine monolis_sparse_matrix_add_value

  subroutine monolis_sparse_matrix_set_value(index, item, A, ndof, ci, cj, csub_i, csub_j, val)
    implicit none
    integer(kint), intent(in) :: ndof
    integer(kint), intent(in) :: index(0:), item(:), ci, cj, csub_i, csub_j
    real(kdouble), intent(inout) :: A(:)
    real(kdouble), intent(in) :: val
    integer(kint) :: j, jn, im, jS, jE, NDOF2

!    NDOF2 = ndof*ndof
!    if(ndof < csub_i) call monolis_stop_by_submatrix_access(ndof, csub_i)
!    if(ndof < csub_j) call monolis_stop_by_submatrix_access(ndof, csub_j)
!
!    jS = index(ci-1) + 1
!    jE = index(ci)
!    do j = jS, jE
!      jn = item(j)
!      if(jn == cj)then
!        im = NDOF2*(j-1) + ndof*(csub_i-1) + csub_j
!        A(im) = val
!        return
!      endif
!      call monolis_stop_by_matrix_assemble(ci, cj)
!    enddo
  end subroutine monolis_sparse_matrix_set_value

  !> CSR data setter
  subroutine monolis_set_matrix_BCSR(monolis, N, NP, NDOF, NZ, A, index, item)
    implicit none
    type(monolis_structure) :: monolis
    integer(kint), intent(in) :: N, NP, NDOF, NZ
    integer(kint), intent(in) :: index(0:), item(:)
    real(kdouble), intent(in) :: A(:)
    integer(kint) :: i

!    if(associated(monolis%MAT%index))  deallocate(monolis%MAT%index)
!    if(associated(monolis%MAT%item))   deallocate(monolis%MAT%item)
!    if(associated(monolis%MAT%A))      deallocate(monolis%MAT%A)
!    if(associated(monolis%MAT%X))      deallocate(monolis%MAT%X)
!    if(associated(monolis%MAT%B))      deallocate(monolis%MAT%B)
!    if(associated(monolis%MAT%indexR)) deallocate(monolis%MAT%indexR)
!    if(associated(monolis%MAT%itemR))  deallocate(monolis%MAT%itemR)
!    if(associated(monolis%MAT%permR))  deallocate(monolis%MAT%permR)
!
!    monolis%MAT%N = N
!    monolis%MAT%NP = NP
!    monolis%MAT%NDOF = NDOF
!
!    allocate(monolis%MAT%index(0:NP), source = 0)
!    allocate(monolis%MAT%item(NZ), source = 0)
!    allocate(monolis%MAT%A(NDOF*NDOF*NZ), source = 0.0d0)
!    allocate(monolis%MAT%X(NDOF*NP), source = 0.0d0)
!    allocate(monolis%MAT%B(NDOF*NP), source = 0.0d0)
!
!    do i = 1, NP
!      monolis%MAT%index(i) = index(i)
!    enddo
!    do i = 1, NZ
!      monolis%MAT%item(i) = item(i)
!    enddo
!    do i = 1, NDOF*NDOF*NZ
!      monolis%MAT%A(i) = A(i)
!    enddo
!
!    allocate(monolis%MAT%indexR(0:NP), source = 0)
!    allocate(monolis%MAT%itemR(NZ), source = 0)
!    allocate(monolis%MAT%permR(NZ), source = 0)
!
!    call monolis_get_CRR_format(N, NP, NZ, &
!      & monolis%MAT%index, monolis%MAT%item, &
!      & monolis%MAT%indexR, monolis%MAT%itemR, monolis%MAT%permR)
  end subroutine monolis_set_matrix_BCSR

  !> boundary condition
  subroutine monolis_set_Dirichlet_bc(monolis, B, node_id, ndof_bc, val)
    implicit none
    type(monolis_structure) :: monolis
    integer(kint), intent(in) :: node_id, ndof_bc
    real(kdouble), intent(in) :: val
    real(kdouble) :: B(:)

!    call monolis_sparse_matrix_add_bc(monolis%MAT%index, monolis%MAT%item, monolis%MAT%A, B, &
!      & monolis%MAT%indexR, monolis%MAT%itemR, monolis%MAT%permR, &
!      & monolis%MAT%ndof, node_id, ndof_bc, val)
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

!    if(ndof < idof) call monolis_stop_by_submatrix_access(ndof, idof)
!
!    is_add = .false.
!    NDOF2 = ndof*ndof
!
!    jS = indexR(nnode-1) + 1
!    jE = indexR(nnode)
!    do j = jS, jE
!      jn = itemR(j)
!      kn = permA(j)
!      do k = 1, ndof
!        B(ndof*(jn-1)+k) = B(ndof*(jn-1)+k) - val*A(NDOF2*(kn-1) + ndof*(k-1) + idof)
!        A(NDOF2*(kn-1) + ndof*(k-1) + idof) = 0.0d0
!      enddo
!    enddo
!
!    jS = index(nnode-1) + 1
!    jE = index(nnode)
!    do j = jS, jE
!      do k = 1, ndof
!        A(NDOF2*(j-1) + ndof*(idof-1) + k) = 0.0d0
!      enddo
!
!      jn = item(j)
!      if(jn == nnode)then
!        A(NDOF2*(j-1) + (ndof+1)*(idof-1) + 1) = 1.0d0
!        is_add = .true.
!      endif
!    enddo
!
!    if(.not. is_add) stop "error: not find a diagonal element in monolis_sparse_matrix_add_bc"
!
!    B(ndof*nnode-ndof+idof) = val
  end subroutine monolis_sparse_matrix_add_bc


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

  subroutine monolis_get_penalty_value(monoMAT)
    implicit none
    type(monolis_mat) :: monoMAT
    integer(kint) :: i, j, k, jS, jE, in, kn, NP, NDOF, NDOF2
    real(kdouble) :: max

!    NP =  monoMAT%NP
!    NDOF  = monoMAT%NDOF
!    NDOF2 = NDOF*NDOF
!    max = 0.0d0
!
!    do i = 1, NP
!      jS = monoMAT%index(i-1) + 1
!      jE = monoMAT%index(i)
!      do j = jS, jE
!        in = monoMAT%item(j)
!        if(i == in)then
!          do k = 1, NDOF
!            kn = NDOF2*(j-1) + (NDOF+1)*(k-1) + 1
!            if(max < monoMAT%A(kn)) max = monoMAT%A(kn)
!          enddo
!        endif
!      enddo
!    enddo
!    monolis_get_penalty_value = max
  end subroutine monolis_get_penalty_value

  subroutine monolis_check_diagonal(monoPRM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_mat) :: monoMAT
    integer(kint) :: i, j, k, jS, jE, in, kn, N, NDOF, NDOF2
    real(kdouble) :: t1, t2

!    if(.not. monoPRM%is_check_diag) return
!    t1 = monolis_get_time()
!
!    N =  monoMAT%N
!    NDOF  = monoMAT%NDOF
!    NDOF2 = NDOF*NDOF
!
!    do i = 1, N
!      jS = monoMAT%index(i-1) + 1
!      jE = monoMAT%index(i)
!      do j = jS, jE
!        in = monoMAT%item(j)
!        if(i == in)then
!          do k = 1, NDOF
!            kn = NDOF2*(j-1) + (NDOF+1)*(k-1) + 1
!            if(monoMAT%A(kn) == 0.0d0)then
!              write(*,"(a,i8,a,i8)")" ** monolis error: zero diagonal at node:", i, " , dof: ", k
!              stop
!            endif
!          enddo
!        endif
!      enddo
!    enddo
!
!    t2 = monolis_get_time()
!    monoPRM%tprep = monoPRM%tprep + t2 - t1
  end subroutine monolis_check_diagonal

end module mod_monolis_spmat_handler
