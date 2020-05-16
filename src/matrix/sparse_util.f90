module mod_monolis_sparse_util
  use mod_monolis_util
  use mod_monolis_stdlib
  implicit none

contains

  subroutine monolis_sparse_matrix_assemble(index, item, A, nnode, ndof, e1t, e2t, stiff)
    implicit none
    integer(kint), intent(in) :: nnode, ndof
    integer(kint), intent(in) :: index(0:), item(:), e1t(nnode), e2t(nnode)
    real(kdouble), intent(inout) :: A(:)
    real(kdouble), intent(in) :: stiff(nnode*ndof,nnode*ndof)
    integer(kint) :: e1(nnode), e2(nnode)
    integer(kint) :: i, j, k, in, jn, im, jS, jE, i2, j2, i1, j1, NDOF2
    integer(kint) :: eperm1(nnode), eperm2(nnode)
    real(kdouble) :: temp(nnode*ndof,nnode*ndof)

    NDOF2 = ndof*ndof
    e1 = e1t
    e2 = e2t
    do i = 1, nnode
      eperm1(i) = i
      eperm2(i) = i
    enddo
    call monolis_qsort_int_with_perm(e1, 1, nnode, eperm1)
    call monolis_qsort_int_with_perm(e2, 1, nnode, eperm2)

    temp = 0.0d0
    do i = 1, nnode
      i1 = eperm2(i)
      do j = 1, nnode
        j1 = eperm1(j)
        do i2 = 1, ndof
          do j2 = 1, ndof
            temp(ndof*(i-1)+i2, ndof*(j-1)+j2) = stiff(ndof*(j1-1)+j2, ndof*(i1-1)+i2)
          enddo
        enddo
      enddo
    enddo

    do i = 1, nnode
      in = e1(i)
      jS = index(in-1) + 1
      jE = index(in)
      aa:do j = 1, nnode
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
        stop "error: merge"
      enddo aa
    enddo
  end subroutine monolis_sparse_matrix_assemble

  subroutine monolis_sparse_matrix_add_bc(index, item, A, B, indexR, itemR, permA, &
    & ndof, nnode, idof, val)
    implicit none
    integer(kint), intent(in) :: nnode, ndof, idof
    integer(kint), intent(in) :: index(0:), item(:), indexR(0:), itemR(:), permA(:)
    real(kdouble), intent(inout) :: A(:), B(:)
    real(kdouble), intent(in) :: val
    integer(kint) :: j, k, jn, kn, jS, jE, NDOF2

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
      endif
    enddo

    B(3*nnode-3+idof) = val
  end subroutine monolis_sparse_matrix_add_bc

  subroutine monolis_get_CRR_format(N, index, item, indexR, itemR, permA)
    implicit none
    integer(kint), intent(in) :: N, index(0:), item(:)
    integer(kint), allocatable :: indexR(:), itemR(:), temp(:), permA(:)
    integer(kint) :: i, j, in, jS, jE, nz, m, p

    allocate(temp(N), source = 0)
    do i = 1, N
      jS = index(i-1) + 1
      jE = index(i)
      do j = jS, jE
        in = item(j)
        temp(in) = temp(in) + 1
      enddo
    enddo

    nz = index(N)
    allocate(indexR(0:N), source = 0)
    allocate(itemR(nz), source = 0)
    allocate(permA(nz), source = 0)

    do i = 1, N
      indexR(i) = indexR(i-1) + temp(i)
    enddo

    temp = 0
    do i = 1, N
      jS = index(i-1) + 1
      jE = index(i)
      do j = jS, jE
        in = item(j)
        m = indexR(in-1)
        temp(in) = temp(in) + 1
        p = temp(in)
        itemR(m + p) = i
        permA(m + p) = j
      enddo
    enddo

    do i = 1, N
      jS = indexR(i-1) + 1
      jE = indexR(i)
      call monolis_qsort_int(itemR(jS:jE), 1, jE-jS+1)
    enddo
  end subroutine monolis_get_CRR_format

end module mod_monolis_sparse_util
