module mod_monolis_reorder
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_restruct
  implicit none

  type monolis_edge_info
    integer(kind=kint) :: N = 0
    integer(kind=kint), pointer :: node(:) => null()
  endtype monolis_edge_info

  integer(kind=kint), save, pointer ::  perm(:) => null()
  integer(kind=kint), save, pointer :: iperm(:) => null()

contains

  subroutine monolis_reorder_matrix_fw(monoPRM, monoCOM, monoCOM_reorder, monoMAT, monoMAT_reorder)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_com) :: monoCOM_reorder
    type(monolis_mat) :: monoMAT
    type(monolis_mat) :: monoMAT_reorder

    if(monoPRM%is_reordering)then
#ifdef WITH_METIS
      allocate( perm(monoMAT%NP))
      allocate(iperm(monoMAT%NP))
      call monolis_reorder_matrix_metis(monoMAT, monoMAT_reorder)
      call monolis_restruct_matrix(monoMAT, monoMAT_reorder, perm, iperm)
      call monolis_restruct_comm(monoCOM, monoCOM_reorder, iperm)
      call monolis_reorder_vector_fw(monoMAT%NP, monoMAT%NDOF, monoMAT%B, monoMAT_reorder%B)
#else
      call monolis_mat_copy(monoMAT, monoMAT_reorder)
      call monolis_com_copy(monoCOM, monoCOM_reorder)
#endif
    else
      call monolis_mat_copy(monoMAT, monoMAT_reorder)
      call monolis_com_copy(monoCOM, monoCOM_reorder)
    endif
  end subroutine monolis_reorder_matrix_fw

  subroutine monolis_reorder_matrix_bk(monoPRM, monoCOM, monoMAT_reorder, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    type(monolis_mat) :: monoMAT_reorder

    if(monoPRM%is_reordering)then
#ifdef WITH_METIS
      call monolis_reorder_back_vector_bk(monoMAT%NP, monoMAT%NDOF, monoMAT_reorder%X, monoMAT%X)
      deallocate( perm)
      deallocate(iperm)
#endif
    endif
  end subroutine monolis_reorder_matrix_bk

  subroutine monolis_reorder_vector_fw(N, NDOF, A, B)
    implicit none
    integer(kind=kint) :: N, NDOF
    real(kind=kdouble) :: A(:)
    real(kind=kdouble) :: B(:)
    integer(kind=kint) :: i, in, jn, jo, j
    do i = 1, N
      in = iperm(i)
      jn = (i -1)*NDOF
      jo = (in-1)*NDOF
      do j = 1, NDOF
        B(jn + j) = A(jo + j)
      enddo
    enddo
  end subroutine monolis_reorder_vector_fw

  subroutine monolis_reorder_back_vector_bk(N, NDOF, B, A)
    implicit none
    integer(kind=kint) :: N, NDOF
    real(kind=kdouble) :: B(:)
    real(kind=kdouble) :: A(:)
    integer(kind=kint) :: i, in, jn, jo, j
    do i = 1, N
      in = perm(i)
      jn = (i -1)*NDOF
      jo = (in-1)*NDOF
      do j = 1, NDOF
        A(jo + j) = B(jn + j)
      enddo
    enddo
  end subroutine monolis_reorder_back_vector_bk

  subroutine monolis_reorder_matrix_metis(monoMAT, monoMAT_reorder)
    implicit none
    type(monolis_edge_info), allocatable :: edge(:)
    type(monolis_mat) :: monoMAT
    type(monolis_mat) :: monoMAT_reorder
    integer(kind=kint) :: N, NP
    integer(kind=kint) :: i, j, k, iS, jS, jE, nu, nl
    integer(kind=kint) :: icel, nedge, ic_type, in, jn, kn, nn, ne
    integer(kind=kint) :: imax, imin
    integer(kind=kint) :: nlocal(20)
    integer(kind=kint), allocatable :: check(:), nozero(:)
    integer(kind=kint) :: ierr
    integer(kind=kint) :: nvtxs
    integer(kind=kint), pointer :: indexL(:)      => null()
    integer(kind=kint), pointer :: indexU(:)      => null()
    integer(kind=kint), pointer :: itemL(:)       => null()
    integer(kind=kint), pointer :: itemU(:)       => null()
    integer(kind=kint), pointer :: xadj(:)        => null()
    integer(kind=kint), pointer :: adjncy(:)      => null()
    integer(kind=kint), pointer :: vwgt(:)        => null()
    integer(kind=kint), pointer :: options(:)     => null()
    integer(kind=kint), pointer :: metis_perm(:)  => null()
    integer(kind=kint), pointer :: metis_iperm(:) => null()

    N = monoMAT%N
    NP = monoMAT%NP
    indexL => monoMAT%indexL
    indexU => monoMAT%indexU
    itemL  => monoMAT%itemL
    itemU  => monoMAT%itemU

    nvtxs = N
    allocate(edge(N))

    do i = 1, N
      nl = indexL(i) - indexL(i-1)
      nu = indexU(i) - indexU(i-1)
      in = nl + nu

      if(0 < in)then
        allocate(nozero(in))
        nozero = 0

        jn = 0
        jS = indexL(i-1) + 1
        jE = indexL(i  )
        do j = jS, jE
          jn = jn + 1
          nozero(jn) = itemL(j)
        enddo
        jS = indexU(i-1) + 1
        jE = indexU(i  )
        do j = jS, jE
          if(itemU(j) <= N)then
            jn = jn + 1
            nozero(jn) = itemU(j)
          endif
        enddo

        call reallocate_array(edge(i)%N, jn, edge(i)%node)
        edge(i)%N = jn

        do j = 1, jn
          edge(i)%node(j) = nozero(j)
        enddo
        deallocate(nozero)
      endif
    enddo

    allocate(xadj(N+1))
    xadj(1) = 0
    do i = 1, N
      xadj(i+1) = xadj(i) + edge(i)%N
    enddo

    nedge = xadj(N+1)
    allocate(adjncy(nedge))
    in = 1
    do i = 1, N
      do j = 1, edge(i)%N
        adjncy(in) = edge(i)%node(j) - 1
        in = in + 1
      enddo
    enddo

#ifdef WITH_METIS
    call METIS_NodeND(nvtxs, xadj, adjncy, vwgt, options, perm, iperm)
#endif

    do i = 1, N
       perm(i) =  perm(i) + 1
      iperm(i) = iperm(i) + 1
    enddo
    do i = N+1, NP
       perm(i) = i
      iperm(i) = i
    enddo

    do i = 1, N
      if(associated(edge(i)%node)) deallocate(edge(i)%node)
    enddo
    deallocate(edge)
    deallocate(xadj)
    deallocate(adjncy)
  end subroutine monolis_reorder_matrix_metis

  subroutine reallocate_array(in, inew, x)
    implicit none
    integer(kind=kint), intent(in) :: in, inew
    integer(kind=kint), pointer :: x(:), t(:)
    integer(kind=kint) :: i

    if(.not. associated(x))then
      allocate(x(inew))
    else
      t => x
      x => null()
      allocate(x(inew))
      do i=1,in
        x(i) = t(i)
      enddo
      deallocate(t)
    endif
  end subroutine reallocate_array
end module mod_monolis_reorder