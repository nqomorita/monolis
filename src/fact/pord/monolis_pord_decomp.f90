!*****************************************************************************
! SPACE (SPArse Cholesky Elimination) Library: mod_monolis_pord_nestdiss.f90
!
! Fortran implementation of nestdiss.c
!*****************************************************************************

module mod_monolis_pord_nestdiss
  use mod_monolis_pord_const
  use mod_monolis_pord_types
  use mod_monolis_pord_graph
  use mod_monolis_pord_gbisect
  implicit none
  private

  public :: newNDnode, freeNDnode, freeNDtree
  public :: setupNDroot
  public :: buildNDtree

contains

  !===========================================================================
  subroutine newNDnode(nd, G, map_arr, nvint)
    type(nestdiss_t), pointer,  intent(out) :: nd
    type(graph_t),    target,   intent(in)  :: G
    integer(kint),                intent(in)  :: map_arr(0:*)
    integer(kint),                intent(in)  :: nvint

    integer(kint) :: dummy_map

    allocate(nd)
    dummy_map = map_arr(0)  ! suppress unused-argument warning
    nd%G       => G
    nd%depth   = 0
    nd%nvint   = nvint
    nd%cwght   = 0
    allocate(nd%intvertex(0:max(1,nvint)-1))
    allocate(nd%intcolor(0:max(1,nvint)-1))
    nd%parent => null(); nd%childB => null(); nd%childW => null()
    ! map_arr pointer is stored as a separate pointer component in nestdiss_t.
    ! However, since map_arr is not owned, we just store a reference via G.
    ! In practice map_arr is the same array shared by the whole ND tree.
    ! We store nothing here; the caller passes map_arr to setupSubgraph.
  end subroutine newNDnode

  !===========================================================================
  subroutine freeNDnode(nd)
    type(nestdiss_t), pointer, intent(inout) :: nd
    if (.not. associated(nd)) return
    if (allocated(nd%intvertex)) deallocate(nd%intvertex)
    if (allocated(nd%intcolor))  deallocate(nd%intcolor)
    deallocate(nd)
  end subroutine freeNDnode

  !===========================================================================
  subroutine freeNDtree(ndroot)
    type(nestdiss_t), pointer, intent(inout) :: ndroot

    type(nestdiss_t), pointer :: nd, parent_nd

    nd => ndroot
    do while (associated(nd%childB))
      nd => nd%childB
    end do

    do while (.not. associated(nd, ndroot))
      parent_nd => nd%parent
      if (associated(parent_nd%childB, nd)) then
        call freeNDnode(nd)
        nd => parent_nd%childW
        do while (associated(nd%childB)); nd => nd%childB; end do
      else
        call freeNDnode(nd)
        nd => parent_nd
      end if
    end do
  end subroutine freeNDtree

  !===========================================================================
  subroutine setupNDroot(ndroot, G, map_arr, nvtx)
    type(nestdiss_t), pointer, intent(out) :: ndroot
    type(graph_t),    target,  intent(in)  :: G
    integer(kint),               intent(out) :: map_arr(0:nvtx-1)
    integer(kint),               intent(in)  :: nvtx

    integer(kint) :: i

    call newNDnode(ndroot, G, map_arr, nvtx)
    do i = 0, nvtx-1
      ndroot%intvertex(i) = i
    end do
  end subroutine setupNDroot

  !===========================================================================
  ! splitNDnode: compute bisection for the subgraph at nd.
  ! map_arr is the global shared mapping array (0:G%nvtx-1).
  !===========================================================================
  subroutine splitNDnode(nd, map_arr, options, cpus)
    type(nestdiss_t), pointer, intent(inout) :: nd
    integer(kint),               intent(inout) :: map_arr(0:nd%G%nvtx-1)
    integer(kint),               intent(in)    :: options(0:*)
    real(kdouble),                  intent(inout) :: cpus(0:*)

    type(graph_t), target  :: Gsub_local
    type(graph_t), pointer :: Gptr
    type(gbisect_t) :: Gbisect
    type(nestdiss_t), pointer :: b_nd, w_nd
    integer(kint), allocatable :: b_intv(:), w_intv(:)
    integer(kint) :: nvint, b_nvint, w_nvint, u, i
    logical :: built_sub

    nvint     = nd%nvint
    built_sub = .false.

    if (nd%G%nvtx == nvint) then
      Gptr => nd%G
      do u = 0, nvint-1
        map_arr(u) = u
      end do
    else
      call setupSubgraph(Gsub_local, nd%G, nd%intvertex, nvint, map_arr)
      Gptr => Gsub_local
      built_sub = .true.
    end if

    call newGbisect(Gbisect, Gptr)

    call pord_starttimer(cpus(TIME_MULTILEVEL))
    call constructSeparator(Gbisect, options, cpus)
    call pord_stoptimer(cpus(TIME_MULTILEVEL))

    call pord_starttimer(cpus(TIME_SMOOTH))
    if (Gbisect%cwght(GRAY) > 0) call smoothSeparator(Gbisect, options)
    call pord_stoptimer(cpus(TIME_SMOOTH))

    ! copy bisection back
    nd%cwght(GRAY)  = Gbisect%cwght(GRAY)
    nd%cwght(BLACK) = Gbisect%cwght(BLACK)
    nd%cwght(WHITE) = Gbisect%cwght(WHITE)

    b_nvint = 0; w_nvint = 0
    do i = 0, nvint-1
      u = nd%intvertex(i)
      nd%intcolor(i) = Gbisect%color(map_arr(u))
      select case (nd%intcolor(i))
        case (BLACK); b_nvint = b_nvint + 1
        case (WHITE); w_nvint = w_nvint + 1
      end select
    end do

    ! build child nodes
    allocate(b_intv(0:max(1,b_nvint)-1), w_intv(0:max(1,w_nvint)-1))
    b_nvint = 0; w_nvint = 0
    do i = 0, nvint-1
      u = nd%intvertex(i)
      if (nd%intcolor(i) == BLACK) then
        b_intv(b_nvint) = u; b_nvint = b_nvint + 1
      else if (nd%intcolor(i) == WHITE) then
        w_intv(w_nvint) = u; w_nvint = w_nvint + 1
      end if
    end do

    call newNDnode(b_nd, nd%G, map_arr, b_nvint)
    call newNDnode(w_nd, nd%G, map_arr, w_nvint)
    if (b_nvint > 0) b_nd%intvertex(0:b_nvint-1) = b_intv(0:b_nvint-1)
    if (w_nvint > 0) w_nd%intvertex(0:w_nvint-1) = w_intv(0:w_nvint-1)
    ! correct nvint after allocation with max(1,...)
    b_nd%nvint = b_nvint; w_nd%nvint = w_nvint
    deallocate(b_intv, w_intv)

    nd%childB => b_nd; b_nd%parent => nd
    nd%childW => w_nd; w_nd%parent => nd
    b_nd%depth = nd%depth + 1
    w_nd%depth = nd%depth + 1

    call freeGbisect(Gbisect)
    if (built_sub) call freeGraph(Gsub_local)
  end subroutine splitNDnode

  !===========================================================================
  ! buildNDtree: build the full nested dissection tree.
  !===========================================================================
  subroutine buildNDtree(ndroot, map_arr, options, cpus)
    type(nestdiss_t), pointer, intent(inout) :: ndroot
    integer(kint),               intent(inout) :: map_arr(0:ndroot%G%nvtx-1)
    integer(kint),               intent(in)    :: options(0:*)
    real(kdouble),                  intent(inout) :: cpus(0:*)

    ! Use an allocatable array of pointers (Fortran has no pointer arrays directly,
    ! so we use an auxiliary derived type for the queue).
    type :: nd_ptr_t
      type(nestdiss_t), pointer :: p => null()
    end type nd_ptr_t

    type(nd_ptr_t), allocatable :: queue(:)
    type(nestdiss_t), pointer   :: nd
    integer(kint) :: maxseps, seps, domainsize, qhead, qtail, qsize

    maxseps    = MAX_SEPS
    domainsize = options(OPTION_DOMAIN_SIZE)
    if (domainsize == 1) maxseps = DEFAULT_SEPS

    qsize = 2*MAX_SEPS + 2
    allocate(queue(0:qsize-1))
    queue(0)%p => ndroot
    qhead = 0; qtail = 1; seps = 0

    do while ((qhead /= qtail) .and. (seps < maxseps))
      seps = seps + 1
      nd => queue(qhead)%p; qhead = qhead + 1

      call splitNDnode(nd, map_arr, options, cpus)

      if (options(OPTION_MSGLVL) > 1) then
        write(*,'(I4,A,I6,A,I6,A,I6)') &
          seps, '. S ', nd%cwght(GRAY), ', B ', nd%cwght(BLACK), &
          ', W ', nd%cwght(WHITE)
      end if

      if ((nd%childB%nvint > MIN_NODES) .and. &
          ((nd%cwght(BLACK) > domainsize) .or. (qtail < DEFAULT_SEPS))) then
        queue(qtail)%p => nd%childB; qtail = qtail + 1
      end if
      if ((nd%childW%nvint > MIN_NODES) .and. &
          ((nd%cwght(WHITE) > domainsize) .or. (qtail < DEFAULT_SEPS))) then
        queue(qtail)%p => nd%childW; qtail = qtail + 1
      end if
    end do

    deallocate(queue)
  end subroutine buildNDtree

end module mod_monolis_pord_nestdiss
!*****************************************************************************
! SPACE (SPArse Cholesky Elimination) Library: mod_monolis_pord_multisector.f90
!
! Fortran implementation of multisector.c
!*****************************************************************************

module mod_monolis_pord_multisector
  use mod_monolis_pord_const
  use mod_monolis_pord_types
  use mod_monolis_pord_nestdiss
  implicit none
  private

  public :: newMultisector, freeMultisector
  public :: trivialMultisector
  public :: constructMultisector

contains

  !===========================================================================
  subroutine newMultisector(ms, G)
    type(multisector_t), intent(out) :: ms
    type(graph_t),       target, intent(inout) :: G

    ms%G => G
    ms%nstages   = 0
    ms%nnodes    = 0
    ms%totmswght = 0
    allocate(ms%stage(0:G%nvtx-1))
    ms%stage = 0
  end subroutine newMultisector

  !===========================================================================
  subroutine freeMultisector(ms)
    type(multisector_t), intent(inout) :: ms
    if (allocated(ms%stage)) deallocate(ms%stage)
    ms%G => null()
  end subroutine freeMultisector

  !===========================================================================
  subroutine trivialMultisector(ms, G)
    type(multisector_t), intent(out) :: ms
    type(graph_t),       target, intent(inout) :: G

    call newMultisector(ms, G)
    ms%stage     = 0
    ms%nstages   = 1
    ms%nnodes    = 0
    ms%totmswght = 0
  end subroutine trivialMultisector

  !===========================================================================
  subroutine constructMultisector(ms, G, options, cpus)
    type(multisector_t), intent(out) :: ms
    type(graph_t),       target, intent(inout) :: G
    integer(kint),         intent(inout) :: options(0:*)
    real(kdouble),            intent(inout) :: cpus(0:*)

    type(nestdiss_t), pointer :: ndroot
    integer(kint), allocatable, target :: map_arr(:)
    integer(kint) :: nvtx, ordtype

    nvtx = G%nvtx

    if ((nvtx <= MIN_NODES) .and. (options(OPTION_ORDTYPE) /= MINIMUM_PRIORITY) &
        .and. (options(OPTION_MSGLVL) > 0)) then
      write(*,'(A,I0,A)') &
        "Warning in constructMultisector: graph has less than ", MIN_NODES, &
        " nodes, using MINIMUM_PRIORITY"
      options(OPTION_ORDTYPE) = MINIMUM_PRIORITY
    end if

    ordtype = options(OPTION_ORDTYPE)

    select case (ordtype)
      case (MINIMUM_PRIORITY)
        call trivialMultisector(ms, G)

      case (INCOMPLETE_ND, MULTISECTION, TRISTAGE_MULTISECTION)
        allocate(map_arr(0:nvtx-1))
        call setupNDroot(ndroot, G, map_arr, nvtx)
        call buildNDtree(ndroot, map_arr, options, cpus)
        if (ordtype == MULTISECTION) then
          call extractMS2stage(ms, ndroot, G)
        else
          call extractMSmultistage(ms, ndroot, G)
        end if
        call freeNDtree(ndroot)
        call freeNDnode(ndroot)
        deallocate(map_arr)

      case default
        write(*,'(A,I0)') "Error in constructMultisector: unknown ordtype ", ordtype
        stop
    end select
  end subroutine constructMultisector

  !===========================================================================
  subroutine extractMS2stage(ms, ndroot, G)
    type(multisector_t), intent(out) :: ms
    type(nestdiss_t),    pointer, intent(in) :: ndroot
    type(graph_t),       target, intent(inout) :: G

    type(nestdiss_t), pointer :: nd, parent_nd
    integer(kint) :: nvint, totmswght, nnodes, i

    call trivialMultisector(ms, G)

    nnodes = 0; totmswght = 0

    ! post-order traversal: start from leftmost leaf
    nd => ndroot
    do while (associated(nd%childB))
      nd => nd%childB
    end do

    do while (.not. associated(nd, ndroot))
      parent_nd => nd%parent
      if (associated(parent_nd%childB, nd)) then
        ! left subtree visited; go to right subtree's leftmost leaf
        nd => parent_nd%childW
        do while (associated(nd%childB)); nd => nd%childB; end do
      else
        ! right subtree visited; move up to parent and process its separator
        nd => parent_nd
        totmswght = totmswght + nd%cwght(GRAY)
        nvint = nd%nvint
        do i = 0, nvint-1
          if (nd%intcolor(i) == GRAY) then
            nnodes = nnodes + 1
            ms%stage(nd%intvertex(i)) = 1
          end if
        end do
        ! loop continues, will check nd /= ndroot at the top
      end if
    end do

    ms%nstages   = 2
    ms%nnodes    = nnodes
    ms%totmswght = totmswght
  end subroutine extractMS2stage

  !===========================================================================
  subroutine extractMSmultistage(ms, ndroot, G)
    type(multisector_t), intent(out) :: ms
    type(nestdiss_t),    pointer, intent(in) :: ndroot
    type(graph_t),       target, intent(inout) :: G

    type(nestdiss_t), pointer :: nd, parent_nd
    integer(kint) :: nvtx, nvint, maxstage, istage, nnodes, totmswght, i, u

    call trivialMultisector(ms, G)
    nvtx = G%nvtx; maxstage = 0; nnodes = 0; totmswght = 0

    nd => ndroot
    do while (associated(nd%childB))
      nd => nd%childB
    end do

    do while (.not. associated(nd, ndroot))
      parent_nd => nd%parent
      if (associated(parent_nd%childB, nd)) then
        nd => parent_nd%childW
        do while (associated(nd%childB)); nd => nd%childB; end do
      else
        nd => parent_nd
        istage = nd%depth + 1
        if (istage > maxstage) maxstage = istage
        totmswght = totmswght + nd%cwght(GRAY)
        nvint = nd%nvint
        do i = 0, nvint-1
          if (nd%intcolor(i) == GRAY) then
            nnodes = nnodes + 1
            ms%stage(nd%intvertex(i)) = istage
          end if
        end do
        ! loop continues
      end if
    end do

    ! reverse stage numbers so that leaf separator = 1, root = maxstage
    do u = 0, nvtx-1
      if (ms%stage(u) > 0) ms%stage(u) = maxstage - ms%stage(u) + 1
    end do

    ms%nstages   = maxstage + 1
    ms%nnodes    = nnodes
    ms%totmswght = totmswght
  end subroutine extractMSmultistage

end module mod_monolis_pord_multisector
