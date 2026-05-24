!*****************************************************************************
!
! SPACE (SPArse Cholesky Elimination) Library: space_ordering.f90
!
! Pure Fortran implementation of SPACE_ordering.
! All algorithms are implemented in Fortran; no C calls are made.
!
! Public subroutines
! ------------------
!   monolis_pord_ordering(G_in, options, use_defaults, T)
!     Computes an ordering for graph G_in and returns an elimination tree T.
!     options(1..6) maps to internal opts(0..5):
!       options(1)=OPTION_ORDTYPE, (2)=OPTION_NODE_SELECTION1,
!       (3)=OPTION_NODE_SELECTION2, (4)=OPTION_NODE_SELECTION3,
!       (5)=OPTION_DOMAIN_SIZE, (6)=OPTION_MSGLVL
!     use_defaults: if .true., options(:) is ignored.
!     T: caller must call freeElimTree when done.
!
!   monolis_pord_perm_from_elimtree(T, nvtx, perm)
!     Extracts 1-based permutation array perm(1:nvtx).
!
!*****************************************************************************

module mod_monolis_pord_ordering
  use mod_monolis_pord_const
  use mod_monolis_pord_types
  use mod_monolis_pord_graph
  use mod_monolis_pord_tree
  use mod_monolis_pord_multisector
  use mod_monolis_pord_minpriority
  implicit none
  private

  public :: monolis_pord_ordering
  public :: monolis_pord_perm_from_elimtree

contains

  !===========================================================================
  subroutine monolis_pord_ordering(G_in, options, use_defaults, T)
    type(graph_t),    intent(in)    :: G_in
    integer(ip),      intent(in)    :: options(:)
    logical,          intent(in)    :: use_defaults
    type(elimtree_t), intent(out)   :: T

    ! -----------------------------------------------------------------------
    ! Local variables
    ! -----------------------------------------------------------------------
    integer(ip) :: opts(0:ORD_OPTION_SLOTS-1)
    real(dp)    :: cpusOrd(0:ORD_TIME_SLOTS-1)

    ! Gwork is the graph used for ordering (either Gc or a copy of G_in).
    ! We keep it heap-allocated so that ms%G and minprior%ms%G remain valid.
    type(graph_t), pointer :: Gwork => null()

    ! vtxmap: maps original vertices to compressed vertices (allocated here)
    integer(ip), allocatable, target :: vtxmap(:)

    ! Multisector and minimum-priority objects. Both declared with TARGET so
    ! that internal pointer components remain valid after the calls return.
    type(multisector_t), target :: ms
    type(minprior_t)            :: minprior
    type(elimtree_t)            :: T_inner

    integer(ip) :: nvtx_G, msglvl
    integer(ip) :: totnstep, totnzf
    real(dp)    :: totops
    logical     :: compressed

    ! -----------------------------------------------------------------------
    ! 1. Set options
    ! -----------------------------------------------------------------------
    if (use_defaults) then
      opts(OPTION_ORDTYPE)         = SPACE_ORDTYPE
      opts(OPTION_NODE_SELECTION1) = SPACE_NODE_SELECTION1
      opts(OPTION_NODE_SELECTION2) = SPACE_NODE_SELECTION2
      opts(OPTION_NODE_SELECTION3) = SPACE_NODE_SELECTION3
      opts(OPTION_DOMAIN_SIZE)     = SPACE_DOMAIN_SIZE
      opts(OPTION_MSGLVL)          = SPACE_MSGLVL
    else
      opts(OPTION_ORDTYPE)         = options(1)
      opts(OPTION_NODE_SELECTION1) = options(2)
      opts(OPTION_NODE_SELECTION2) = options(3)
      opts(OPTION_NODE_SELECTION3) = options(4)
      opts(OPTION_DOMAIN_SIZE)     = options(5)
      opts(OPTION_MSGLVL)          = options(6)
    end if
    msglvl = opts(OPTION_MSGLVL)
    cpusOrd = 0.0_dp

    ! -----------------------------------------------------------------------
    ! 2. Compress graph
    !    compressGraph may allocate a new compressed graph into Gwork, or
    !    signal compressed=.false. in which case we use a heap copy of G_in.
    ! -----------------------------------------------------------------------
    nvtx_G = G_in%nvtx
    allocate(vtxmap(0:nvtx_G-1))
    allocate(Gwork)

    call pord_starttimer(cpusOrd(TIME_COMPRESS))
    call compressGraph(Gwork, G_in, vtxmap, compressed)
    call pord_stoptimer(cpusOrd(TIME_COMPRESS))

    if (.not. compressed) then
      ! compressGraph left Gwork uninitialised; fill it with a copy of G_in
      call freeGraph(Gwork)   ! clean up any partial state
      call copy_graph(Gwork, G_in)
      if (msglvl > 0) write(*,'(A)') "no compressed graph constructed"
    else
      if (msglvl > 0) then
        write(*,'(A,I0,A,I0,A)') &
          "compressed graph constructed (#nodes ", Gwork%nvtx, &
          ", #edges ", Gwork%nedges/2, ")"
      end if
    end if

    ! -----------------------------------------------------------------------
    ! 3. Compute multisector
    !    opts may be modified (ORDTYPE downgrade) inside constructMultisector.
    ! -----------------------------------------------------------------------
    call pord_starttimer(cpusOrd(TIME_MS))
    call constructMultisector(ms, Gwork, opts, cpusOrd)
    call pord_stoptimer(cpusOrd(TIME_MS))

    if (msglvl > 0) then
      write(*,'(A,I0,A,I0,A,I0)') &
        "quality of multisector: #stages ", ms%nstages, &
        ", #nodes ", ms%nnodes, &
        ", weight ", ms%totmswght
    end if

    ! -----------------------------------------------------------------------
    ! 4. Minimum-priority ordering
    !    ms is declared TARGET so minprior%ms remains valid after the call.
    ! -----------------------------------------------------------------------
    call pord_starttimer(cpusOrd(TIME_BOTTOMUP))
    call setupMinPriority(minprior, ms)
    call orderMinPriority(T_inner, minprior, opts, cpusOrd)
    call pord_stoptimer(cpusOrd(TIME_BOTTOMUP))

    if (msglvl > 0) then
      call stageinfo_sum(minprior, ms%nstages, totnstep, totnzf, totops)
      write(*,'(A,I0,A,I0,A,ES12.4)') &
        "quality of ordering: #steps ", totnstep, &
        ", nzl ", totnzf, ", ops ", totops
    end if

    ! -----------------------------------------------------------------------
    ! 5. Expand elimination tree (undo compression mapping)
    ! -----------------------------------------------------------------------
    if (compressed) then
      call expandElimTree(T, T_inner, vtxmap, nvtx_G)
      call freeElimTree(T_inner)
    else
      T = T_inner
      T_inner%nvtx = 0; T_inner%nfronts = 0  ! prevent double free
    end if

    ! -----------------------------------------------------------------------
    ! 6. Free working objects (in reverse order of construction)
    ! -----------------------------------------------------------------------
    call freeMinPriority(minprior)
    call freeMultisector(ms)
    call freeGraph(Gwork)
    deallocate(Gwork)
    deallocate(vtxmap)

  contains

    ! Heap copy: allocates new arrays inside dst
    subroutine copy_graph(dst, src)
      type(graph_t), intent(out) :: dst
      type(graph_t), intent(in)  :: src
      integer(ip) :: n, ne
      n  = src%nvtx
      ne = src%nedges
      call newGraph(dst, n, max(ne, 0))
      dst%gtype    = src%gtype
      dst%totvwght = src%totvwght
      dst%xadj(0:n)    = src%xadj(0:n)
      if (ne > 0) dst%adjncy(0:ne-1) = src%adjncy(0:ne-1)
      dst%vwght(0:n-1) = src%vwght(0:n-1)
    end subroutine copy_graph

  end subroutine monolis_pord_ordering


  !===========================================================================
  ! monolis_pord_perm_from_elimtree: fills perm(1:nvtx) (Fortran 1-based).
  !===========================================================================
  subroutine monolis_pord_perm_from_elimtree(T, nvtx, perm)
    type(elimtree_t), intent(in)  :: T
    integer(ip),      intent(in)  :: nvtx
    integer(ip),      intent(out) :: perm(nvtx)

    integer(ip), allocatable :: perm0(:)
    allocate(perm0(0:nvtx-1))
    call permFromElimTree(T, perm0)
    perm(1:nvtx) = perm0(0:nvtx-1)
    deallocate(perm0)
  end subroutine monolis_pord_perm_from_elimtree

end module mod_monolis_pord_ordering
