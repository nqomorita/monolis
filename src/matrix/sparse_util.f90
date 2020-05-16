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
    integer(kint) :: i, j, k, in, jn, jS, jE, i2, j2, i1, j1
    integer(kint) :: eperm1(nnode), eperm2(nnode)
    real(kdouble) :: temp(nnode*ndof,nnode*ndof)

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
            temp(ndof*(j-1)+j2, ndof*(i-1)+i2) = stiff(ndof*(j1-1)+j2, ndof*(i1-1)+i2)
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
            A(9*k-8) = A(9*k-8) + temp(3*i-2,3*j-2)
            A(9*k-7) = A(9*k-7) + temp(3*i-2,3*j-1)
            A(9*k-6) = A(9*k-6) + temp(3*i-2,3*j  )
            A(9*k-5) = A(9*k-5) + temp(3*i-1,3*j-2)
            A(9*k-4) = A(9*k-4) + temp(3*i-1,3*j-1)
            A(9*k-3) = A(9*k-3) + temp(3*i-1,3*j  )
            A(9*k-2) = A(9*k-2) + temp(3*i  ,3*j-2)
            A(9*k-1) = A(9*k-1) + temp(3*i  ,3*j-1)
            A(9*k  ) = A(9*k  ) + temp(3*i  ,3*j  )
            jS = k + 1
            cycle aa
          endif
        enddo
        stop "error: merge"
      enddo aa
    enddo
  end subroutine monolis_sparse_matrix_assemble

end module mod_monolis_sparse_util
