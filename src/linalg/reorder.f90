module mod_monolis_reorder
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_restruct
  use mod_monolis_util
  implicit none

  type monolis_edge_info
    integer(kind=kint) :: N = 0
    integer(kind=kint), pointer :: node(:) => null()
  endtype monolis_edge_info

contains

  subroutine monolis_reorder_matrix_fw(monoPRM, monoCOM, monoCOM_reorder, monoMAT, monoMAT_reorder)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_com) :: monoCOM_reorder
    type(monolis_mat) :: monoMAT
    type(monolis_mat) :: monoMAT_reorder
    integer(kind=kint), pointer ::  perm(:), iperm(:)

    if(monoPRM%is_debug) call monolis_debug_header("monolis_reorder_matrix_fw")

    if(monoPRM%is_reordering)then
#ifdef WITH_METIS
      allocate(monoMAT%perm(monoMAT%NP))
      allocate(monoMAT%iperm(monoMAT%NP))
      perm => monoMAT%perm
      iperm => monoMAT%iperm
      call monolis_reorder_matrix_metis(monoMAT, monoMAT_reorder)
      call monolis_restruct_matrix(monoMAT, monoMAT_reorder, perm, iperm)
      call monolis_restruct_comm(monoCOM, monoCOM_reorder, iperm)
      call monolis_reorder_vector_fw(monoMAT, monoMAT%NP, monoMAT%NDOF, monoMAT%B, monoMAT_reorder%B)
      if(.not. monoPRM%is_init_x)then
        call monolis_reorder_vector_fw(monoMAT, monoMAT%NP, monoMAT%NDOF, monoMAT%X, monoMAT_reorder%X)
      endif
#else
      call monolis_mat_copy_by_pointer(monoMAT, monoMAT_reorder)
      call monolis_com_copy(monoCOM, monoCOM_reorder)
#endif
    else
      call monolis_mat_copy_by_pointer(monoMAT, monoMAT_reorder)
      call monolis_com_copy(monoCOM, monoCOM_reorder)
    endif
  end subroutine monolis_reorder_matrix_fw

  subroutine monolis_reorder_matrix_bk(monoPRM, monoCOM, monoMAT_reorder, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    type(monolis_mat) :: monoMAT_reorder

    if(monoPRM%is_debug) call monolis_debug_header("monolis_reorder_matrix_bk")

    if(monoPRM%is_reordering)then
#ifdef WITH_METIS
      call monolis_reorder_back_vector_bk(monoMAT, monoMAT%NP, monoMAT%NDOF, monoMAT_reorder%X, monoMAT%X)
      deallocate(monoMAT%perm)
      deallocate(monoMAT%iperm)
#endif
    endif
  end subroutine monolis_reorder_matrix_bk

  subroutine monolis_reorder_vector_fw(monoMAT, N, NDOF, A, B)
    implicit none
    type(monolis_mat) :: monoMAT
    integer(kind=kint) :: N, NDOF
    real(kind=kdouble) :: A(:)
    real(kind=kdouble) :: B(:)
    integer(kind=kint) :: i, in, jn, jo, j
    do i = 1, N
      in = monoMAT%iperm(i)
      jn = (in-1)*NDOF
      jo = (i -1)*NDOF
      do j = 1, NDOF
        B(jn + j) = A(jo + j)
      enddo
    enddo
  end subroutine monolis_reorder_vector_fw

  subroutine monolis_reorder_back_vector_bk(monoMAT, N, NDOF, B, A)
    implicit none
    type(monolis_mat) :: monoMAT
    integer(kind=kint) :: N, NDOF
    real(kind=kdouble) :: B(:)
    real(kind=kdouble) :: A(:)
    integer(kind=kint) :: i, in, jn, jo, j
    do i = 1, N
      in = monoMAT%perm(i)
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
    integer(kind=kint) :: i, j, jS, jE
    integer(kind=kint) :: nedge, in, jn
    integer(kind=kint), allocatable :: nozero(:)
    integer(kind=kint) :: nvtxs
    integer(kind=kint), pointer :: index(:)   => null()
    integer(kind=kint), pointer :: item(:)    => null()
    integer(kind=kint), pointer :: xadj(:)    => null()
    integer(kind=kint), pointer :: adjncy(:)  => null()
    integer(kind=kint), pointer :: vwgt(:)    => null()
    integer(kind=kint), pointer :: options(:) => null()
    integer(kind=kint), pointer :: perm(:), iperm(:)

    N = monoMAT%N
    NP = monoMAT%NP
    index => monoMAT%index
    item  => monoMAT%item
    perm => monoMAT%perm
    iperm => monoMAT%iperm

    nvtxs = N
    allocate(edge(N))

    do i = 1, N
      in = index(i) - index(i-1)

      if(0 < in)then
        allocate(nozero(in-1))
        nozero = 0

        jn = 0
        jS = index(i-1) + 1
        jE = index(i  )
        do j = jS, jE
          if(item(j) /= i)then
            jn = jn + 1
            nozero(jn) = item(j)
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