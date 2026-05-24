!*****************************************************************************
! SPACE (SPArse Cholesky Elimination) Library: mod_monolis_pord_graph.f90
!
! Fortran implementation of graph.c functions needed by SPACE_ordering.
!*****************************************************************************

module mod_monolis_pord_graph
  use mod_monolis_pord_const
  use mod_monolis_pord_types
  implicit none
  private

  public :: newGraph, freeGraph, setupSubgraph, compressGraph
  public :: setupGridGraph

contains

  !===========================================================================
  subroutine newGraph(G, nvtx, nedges)
    type(graph_t), intent(out) :: G
    integer(kint),   intent(in)  :: nvtx, nedges
    integer(kint) :: i

    G%nvtx     = nvtx
    G%nedges   = nedges
    G%gtype    = MONOLIS_PORD_UNWEIGHTED
    G%totvwght = nvtx
    allocate(G%xadj(0:nvtx))
    allocate(G%adjncy(0:max(nedges-1, 0)))
    allocate(G%vwght(0:nvtx-1))
    do i = 0, nvtx-1
      G%vwght(i) = 1
    end do
  end subroutine newGraph

  !===========================================================================
  subroutine freeGraph(G)
    type(graph_t), intent(inout) :: G
    if (allocated(G%xadj))   deallocate(G%xadj)
    if (allocated(G%adjncy)) deallocate(G%adjncy)
    if (allocated(G%vwght))  deallocate(G%vwght)
    G%nvtx = 0; G%nedges = 0; G%totvwght = 0
  end subroutine freeGraph

  !===========================================================================
  ! Extracts the subgraph induced by intvertex(0:nvint-1) from G.
  ! vtxmap(0:G%nvtx-1) maps G vertices to subgraph vertices.
  !===========================================================================
  subroutine setupSubgraph(Gsub, G, intvertex, nvint, vtxmap)
    type(graph_t), intent(out) :: Gsub
    type(graph_t), intent(in)  :: G
    integer(kint),   intent(in)  :: nvint
    integer(kint),   intent(in)  :: intvertex(0:nvint-1)
    integer(kint),   intent(out) :: vtxmap(0:G%nvtx-1)

    integer(kint) :: nedgesGsub, totvwght, u, v, i, j, jstart, jstop, ptr

    ! --- count edges and set up local indices ---
    nedgesGsub = 0
    do i = 0, nvint-1
      u = intvertex(i)
      jstart = G%xadj(u); jstop = G%xadj(u+1)
      do j = jstart, jstop-1
        vtxmap(G%adjncy(j)) = -1
      end do
      nedgesGsub = nedgesGsub + (jstop - jstart)
    end do
    do i = 0, nvint-1
      vtxmap(intvertex(i)) = i
    end do

    call newGraph(Gsub, nvint, nedgesGsub)

    ! --- build induced subgraph ---
    totvwght = 0; ptr = 0
    do i = 0, nvint-1
      u = intvertex(i)
      Gsub%xadj(i) = ptr
      Gsub%vwght(i) = G%vwght(u)
      totvwght = totvwght + G%vwght(u)
      jstart = G%xadj(u); jstop = G%xadj(u+1)
      do j = jstart, jstop-1
        v = G%adjncy(j)
        if (vtxmap(v) >= 0) then
          Gsub%adjncy(ptr) = vtxmap(v)
          ptr = ptr + 1
        end if
      end do
    end do
    Gsub%xadj(nvint) = ptr
    Gsub%gtype    = G%gtype
    Gsub%totvwght = totvwght
  end subroutine setupSubgraph

  !===========================================================================
  ! compressGraph: returns .true. if compression was found.
  ! vtxmap(0:G%nvtx-1) maps each vertex to its representative.
  ! Gc is the compressed graph (only valid when function returns .true.)
  !===========================================================================
  subroutine compressGraph(Gc, G, vtxmap, compressed)
    type(graph_t), intent(out) :: Gc
    type(graph_t), intent(in)  :: G
    integer(kint),   intent(out) :: vtxmap(0:G%nvtx-1)
    logical,       intent(out) :: compressed

    integer(kint) :: nvtx, nvtxGc, nedgesGc
    integer(kint), allocatable :: perm(:)
    integer(kint) :: u, v, i, istart, istop

    nvtx = G%nvtx
    call indNodes(G, vtxmap, nvtxGc)

    if (real(nvtxGc, kdouble) > MONOLIS_PORD_COMPRESS_FRACTION * real(nvtx, kdouble)) then
      compressed = .false.
      return
    end if

    allocate(perm(0:nvtx-1))
    compressed = .true.

    ! count edges of compressed graph
    nedgesGc = 0
    do u = 0, nvtx-1
      if (vtxmap(u) == u) then
        istart = G%xadj(u); istop = G%xadj(u+1)
        do i = istart, istop-1
          v = G%adjncy(i)
          if (vtxmap(v) == v) nedgesGc = nedgesGc + 1
        end do
      end if
    end do

    call newGraph(Gc, nvtxGc, nedgesGc)

    nvtxGc = 0; nedgesGc = 0
    do u = 0, nvtx-1
      if (vtxmap(u) == u) then
        istart = G%xadj(u); istop = G%xadj(u+1)
        Gc%xadj(nvtxGc) = nedgesGc
        Gc%vwght(nvtxGc) = 0
        perm(u) = nvtxGc
        nvtxGc = nvtxGc + 1
        do i = istart, istop-1
          v = G%adjncy(i)
          if (vtxmap(v) == v) then
            Gc%adjncy(nedgesGc) = v
            nedgesGc = nedgesGc + 1
          end if
        end do
      end if
    end do
    Gc%xadj(nvtxGc) = nedgesGc

    do i = 0, nedgesGc-1
      Gc%adjncy(i) = perm(Gc%adjncy(i))
    end do
    do u = 0, nvtx-1
      vtxmap(u) = perm(vtxmap(u))
      Gc%vwght(vtxmap(u)) = Gc%vwght(vtxmap(u)) + G%vwght(u)
    end do
    Gc%gtype    = MONOLIS_PORD_WEIGHTED
    Gc%totvwght = G%totvwght

    deallocate(perm)
  end subroutine compressGraph

  !===========================================================================
  ! indNodes: finds indistinguishable vertices; fills vtxmap.
  ! Returns cnvtx = number of representative vertices.
  !===========================================================================
  subroutine indNodes(G, vtxmap, cnvtx)
    type(graph_t), intent(in)  :: G
    integer(kint),   intent(out) :: vtxmap(0:G%nvtx-1)
    integer(kint),   intent(out) :: cnvtx

    integer(kint) :: nvtx
    integer(kint), allocatable :: deg(:), checksum(:), tmp(:)
    integer(kint) :: u, v, i, j, istart, istop, jstart, jstop
    logical :: fail

    nvtx = G%nvtx
    allocate(deg(0:nvtx-1), checksum(0:nvtx-1), tmp(0:nvtx-1))

    do u = 0, nvtx-1
      istart = G%xadj(u); istop = G%xadj(u+1)
      deg(u) = istop - istart
      checksum(u) = u
      tmp(u) = -1
      vtxmap(u) = u
      do i = istart, istop-1
        checksum(u) = checksum(u) + G%adjncy(i)
      end do
    end do

    cnvtx = nvtx
    do u = 0, nvtx-1
      if (vtxmap(u) == u) then
        tmp(u) = u
        istart = G%xadj(u); istop = G%xadj(u+1)
        do i = istart, istop-1
          tmp(G%adjncy(i)) = u
        end do
        do i = istart, istop-1
          v = G%adjncy(i)
          if ((v > u) .and. (checksum(v) == checksum(u)) .and. &
              (deg(v) == deg(u)) .and. (vtxmap(v) == v)) then
            jstart = G%xadj(v); jstop = G%xadj(v+1)
            fail = .false.
            do j = jstart, jstop-1
              if (tmp(G%adjncy(j)) /= u) then
                fail = .true.
                exit
              end if
            end do
            if (.not. fail) then
              vtxmap(v) = u
              cnvtx = cnvtx - 1
            end if
          end if
        end do
      end if
    end do

    deallocate(deg, checksum, tmp)
  end subroutine indNodes

  !===========================================================================
  ! setupGridGraph: creates an unweighted grid (MONOLIS_PORD_GRID), mesh (MONOLIS_PORD_MESH), or
  ! torus (MONOLIS_PORD_TORUS) graph of dimensions dimX x dimY.
  !===========================================================================
  subroutine setupGridGraph(G, dimX, dimY, gtype_in)
    type(graph_t), intent(out) :: G
    integer(kint),   intent(in)  :: dimX, dimY, gtype_in

    integer(kint) :: nvtx, nedges, knz, k

    nvtx = dimX * dimY

    select case (gtype_in)

      case (MONOLIS_PORD_GRID)
        nedges = 8 &
               + 6 * (dimX-2 + dimY-2) &
               + 4 * (dimX-2) * (dimY-2)
        call newGraph(G, nvtx, nedges)
        knz = 0
        do k = 0, nvtx-1
          G%xadj(k) = knz
          if (mod(k+1, dimX) /= 0) then
            G%adjncy(knz) = k+1;   knz = knz+1
          end if
          if (mod(k, dimX) /= 0) then
            G%adjncy(knz) = k-1;   knz = knz+1
          end if
          if (k+dimX < nvtx) then
            G%adjncy(knz) = k+dimX; knz = knz+1
          end if
          if (k-dimX >= 0) then
            G%adjncy(knz) = k-dimX; knz = knz+1
          end if
        end do
        G%xadj(nvtx) = knz

      case (MONOLIS_PORD_MESH)
        nedges = 8 &
               + 6 * (dimX-2 + dimY-2) &
               + 4 * (dimX-2) * (dimY-2) &
               + 4 * (dimX-1) * (dimY-1)
        call newGraph(G, nvtx, nedges)
        knz = 0
        do k = 0, nvtx-1
          G%xadj(k) = knz
          if (mod(k+1, dimX) /= 0) then
            G%adjncy(knz) = k+1; knz = knz+1
            if (k+1+dimX < nvtx) then
              G%adjncy(knz) = k+1+dimX; knz = knz+1
            end if
            if (k+1-dimX >= 0) then
              G%adjncy(knz) = k+1-dimX; knz = knz+1
            end if
          end if
          if (mod(k, dimX) /= 0) then
            G%adjncy(knz) = k-1; knz = knz+1
            if (k-1+dimX < nvtx) then
              G%adjncy(knz) = k-1+dimX; knz = knz+1
            end if
            if (k-1-dimX >= 0) then
              G%adjncy(knz) = k-1-dimX; knz = knz+1
            end if
          end if
          if (k+dimX < nvtx) then
            G%adjncy(knz) = k+dimX;  knz = knz+1
          end if
          if (k-dimX >= 0) then
            G%adjncy(knz) = k-dimX;  knz = knz+1
          end if
        end do
        G%xadj(nvtx) = knz

      case (MONOLIS_PORD_TORUS)
        nedges = 4 * nvtx
        call newGraph(G, nvtx, nedges)
        knz = 0
        do k = 0, nvtx-1
          G%xadj(k) = knz
          if (mod(k+1, dimX) == 0) then
            G%adjncy(knz) = k+1-dimX
          else
            G%adjncy(knz) = k+1
          end if
          knz = knz + 1
          if (mod(k, dimX) == 0) then
            G%adjncy(knz) = k-1+dimX
          else
            G%adjncy(knz) = k-1
          end if
          knz = knz + 1
          G%adjncy(knz) = mod(k+dimX, nvtx)
          knz = knz + 1
          G%adjncy(knz) = mod(k + dimX*(dimY-1), nvtx)
          knz = knz + 1
        end do
        G%xadj(nvtx) = knz

      case default
        write(*,*) 'setupGridGraph: unknown gtype_in', gtype_in
        stop

    end select
  end subroutine setupGridGraph

end module mod_monolis_pord_graph
!*****************************************************************************
! SPACE (SPArse Cholesky Elimination) Library: mod_monolis_pord_tree.f90
!
! Fortran implementation of tree.c functions needed by SPACE_ordering.
! (newElimTree, freeElimTree, initFchSilbRoot, firstPostorder,
!  nextPostorder, permFromElimTree, expandElimTree)
!*****************************************************************************

module mod_monolis_pord_tree
  use mod_monolis_pord_const
  use mod_monolis_pord_types
  implicit none
  private

  public :: newElimTree, freeElimTree
  public :: initFchSilbRoot
  public :: firstPostorder, nextPostorder
  public :: permFromElimTree
  public :: expandElimTree

contains

  !===========================================================================
  subroutine newElimTree(T, nvtx, nfronts)
    type(elimtree_t), intent(out) :: T
    integer(kint),      intent(in)  :: nvtx, nfronts

    T%nvtx    = nvtx
    T%nfronts = nfronts
    T%root    = -1
    allocate(T%ncolfactor(0:nfronts-1))
    allocate(T%ncolupdate(0:nfronts-1))
    allocate(T%parent(0:nfronts-1))
    allocate(T%firstchild(0:nfronts-1))
    allocate(T%silbings(0:nfronts-1))
    allocate(T%vtx2front(0:nvtx-1))
    T%ncolfactor = 0
    T%ncolupdate = 0
    T%parent     = -1
    T%firstchild = -1
    T%silbings   = -1
    T%vtx2front  = -1
  end subroutine newElimTree

  !===========================================================================
  subroutine freeElimTree(T)
    type(elimtree_t), intent(inout) :: T
    if (allocated(T%ncolfactor)) deallocate(T%ncolfactor)
    if (allocated(T%ncolupdate)) deallocate(T%ncolupdate)
    if (allocated(T%parent))     deallocate(T%parent)
    if (allocated(T%firstchild)) deallocate(T%firstchild)
    if (allocated(T%silbings))   deallocate(T%silbings)
    if (allocated(T%vtx2front))  deallocate(T%vtx2front)
    T%nvtx = 0; T%nfronts = 0; T%root = -1
  end subroutine freeElimTree

  !===========================================================================
  ! Sets firstchild and silbings from parent, also sets T%root.
  !===========================================================================
  subroutine initFchSilbRoot(T)
    type(elimtree_t), intent(inout) :: T

    integer(kint) :: nfronts, J, pJ

    nfronts = T%nfronts
    T%firstchild = -1
    T%silbings   = -1
    T%root       = -1

    do J = nfronts-1, 0, -1
      pJ = T%parent(J)
      if (pJ /= -1) then
        T%silbings(J) = T%firstchild(pJ)
        T%firstchild(pJ) = J
      else
        T%silbings(J) = T%root
        T%root = J
      end if
    end do
  end subroutine initFchSilbRoot

  !===========================================================================
  ! firstPostorder: returns the leftmost leaf in a post-order traversal.
  !===========================================================================
  integer(kint) function firstPostorder(T)
    type(elimtree_t), intent(in) :: T
    integer(kint) :: J

    J = T%root
    if (J /= -1) then
      do while (T%firstchild(J) /= -1)
        J = T%firstchild(J)
      end do
    end if
    firstPostorder = J
  end function firstPostorder

  !===========================================================================
  ! nextPostorder: returns next node after J in post-order, or -1.
  !===========================================================================
  integer(kint) function nextPostorder(T, J_in)
    type(elimtree_t), intent(in) :: T
    integer(kint),      intent(in) :: J_in
    integer(kint) :: J

    J = J_in
    if (T%silbings(J) /= -1) then
      J = T%silbings(J)
      do while (T%firstchild(J) /= -1)
        J = T%firstchild(J)
      end do
    else
      J = T%parent(J)
    end if
    nextPostorder = J
  end function nextPostorder

  !===========================================================================
  ! permFromElimTree: fills perm(0:nvtx-1) from T via post-order traversal.
  !===========================================================================
  subroutine permFromElimTree(T, perm)
    type(elimtree_t), intent(in)  :: T
    integer(kint),      intent(out) :: perm(0:T%nvtx-1)

    integer(kint) :: nvtx, nfronts, K, u, count
    integer(kint), allocatable :: first(:), link(:)

    nvtx   = T%nvtx
    nfronts = T%nfronts

    allocate(first(0:nfronts-1), link(0:nvtx-1))
    first = -1

    do u = nvtx-1, 0, -1
      K = T%vtx2front(u)
      link(u)  = first(K)
      first(K) = u
    end do

    count = 0
    K = firstPostorder(T)
    do while (K /= -1)
      u = first(K)
      do while (u /= -1)
        perm(u) = count
        count = count + 1
        u = link(u)
      end do
      K = nextPostorder(T, K)
    end do

    deallocate(first, link)
  end subroutine permFromElimTree

  !===========================================================================
  ! expandElimTree: expands T (for compressed graph) back to nvtxorg vertices.
  ! vtxmap(0:nvtxorg-1) maps original vertices to compressed vertices.
  !===========================================================================
  subroutine expandElimTree(T2, T, vtxmap, nvtxorg)
    type(elimtree_t), intent(out) :: T2
    type(elimtree_t), intent(in)  :: T
    integer(kint),      intent(in)  :: nvtxorg
    integer(kint),      intent(in)  :: vtxmap(0:nvtxorg-1)

    integer(kint) :: nfronts, J, u

    nfronts = T%nfronts
    call newElimTree(T2, nvtxorg, nfronts)
    T2%root = T%root

    do J = 0, nfronts-1
      T2%ncolfactor(J) = T%ncolfactor(J)
      T2%ncolupdate(J) = T%ncolupdate(J)
      T2%parent(J)     = T%parent(J)
      T2%firstchild(J) = T%firstchild(J)
      T2%silbings(J)   = T%silbings(J)
    end do

    do u = 0, nvtxorg-1
      T2%vtx2front(u) = T%vtx2front(vtxmap(u))
    end do
  end subroutine expandElimTree

end module mod_monolis_pord_tree
