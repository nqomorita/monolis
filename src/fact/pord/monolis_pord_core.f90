!*****************************************************************************
! SPACE (SPArse Cholesky Elimination) Library: mod_monolis_pord_const.f90
!
! Fortran module: constants, parameters and separator evaluation function.
! Translated from include/const.h, include/params.h and include/eval.h.
!*****************************************************************************

module mod_monolis_pord_const
  use iso_c_binding, only: c_int, c_double
  implicit none

  integer, parameter :: dp = c_double
  integer, parameter :: ip = c_int

  ! --- matrix/graph topology types (from const.h) ---
  integer(ip), parameter :: GRID  = 0
  integer(ip), parameter :: MESH  = 1
  integer(ip), parameter :: TORUS = 2
  integer(ip), parameter :: HB    = 3

  ! --- graph types ---
  integer(ip), parameter :: UNWEIGHTED      = 0
  integer(ip), parameter :: WEIGHTED        = 1

  ! --- ordering types ---
  integer(ip), parameter :: MINIMUM_PRIORITY      = 0
  integer(ip), parameter :: INCOMPLETE_ND         = 1
  integer(ip), parameter :: MULTISECTION          = 2
  integer(ip), parameter :: TRISTAGE_MULTISECTION = 3

  ! --- node selection strategies (interior ordering) ---
  integer(ip), parameter :: AMD   = 0
  integer(ip), parameter :: AMF   = 1
  integer(ip), parameter :: AMMF  = 2
  integer(ip), parameter :: AMIND = 3

  ! --- node selection strategies (domain decomposition) ---
  integer(ip), parameter :: QMD   = 0
  integer(ip), parameter :: QMRDV = 1
  integer(ip), parameter :: QRAND = 2

  ! --- default options ---
  integer(ip), parameter :: SPACE_ORDTYPE         = MULTISECTION
  integer(ip), parameter :: SPACE_NODE_SELECTION1 = AMMF
  integer(ip), parameter :: SPACE_NODE_SELECTION2 = AMMF
  integer(ip), parameter :: SPACE_NODE_SELECTION3 = QMRDV
  integer(ip), parameter :: SPACE_DOMAIN_SIZE     = 200
  integer(ip), parameter :: SPACE_MSGLVL          = 2
  integer(ip), parameter :: SPACE_ETREE_NONZ      = 256
  integer(ip), parameter :: SPACE_ETREE_BAL       = 5
  integer(ip), parameter :: SPACE_MASK_OFFSET      = 2

  ! --- misc ---
  integer(ip), parameter :: PORD_TRUE  = 1
  integer(ip), parameter :: PORD_FALSE = 0
  integer(ip), parameter :: PORD_ERR   = -1
  integer(ip), parameter :: PORD_NOERR = 0
  integer(ip), parameter :: MAX_INT    = 2**30 - 1
  real(dp),    parameter :: MAX_FLOAT  = 1.0e31_dp
  real(dp),    parameter :: EPS        = 0.001_dp

  ! --- color constants (also used as array indices, 0-based) ---
  integer(ip), parameter :: GRAY  = 0
  integer(ip), parameter :: BLACK = 1
  integer(ip), parameter :: WHITE = 2

  ! --- Dulmage-Mendelsohn flags ---
  integer(ip), parameter :: SI = 0
  integer(ip), parameter :: SX = 1
  integer(ip), parameter :: SR = 2
  integer(ip), parameter :: BI = 3
  integer(ip), parameter :: BX = 4
  integer(ip), parameter :: BR = 5

  ! --- option array indices (0-based) ---
  integer(ip), parameter :: ORD_OPTION_SLOTS       = 7
  integer(ip), parameter :: OPTION_ORDTYPE         = 0
  integer(ip), parameter :: OPTION_NODE_SELECTION1 = 1
  integer(ip), parameter :: OPTION_NODE_SELECTION2 = 2
  integer(ip), parameter :: OPTION_NODE_SELECTION3 = 3
  integer(ip), parameter :: OPTION_DOMAIN_SIZE      = 4
  integer(ip), parameter :: OPTION_MSGLVL           = 5
  integer(ip), parameter :: OPTION_ETREE_NONZ       = 6

  ! --- timing array indices (0-based) ---
  integer(ip), parameter :: ORD_TIME_SLOTS    = 12
  integer(ip), parameter :: TIME_COMPRESS     = 0
  integer(ip), parameter :: TIME_MS           = 1
  integer(ip), parameter :: TIME_MULTILEVEL   = 2
  integer(ip), parameter :: TIME_INITDOMDEC   = 3
  integer(ip), parameter :: TIME_COARSEDOMDEC = 4
  integer(ip), parameter :: TIME_INITSEP      = 5
  integer(ip), parameter :: TIME_REFINESEP    = 6
  integer(ip), parameter :: TIME_SMOOTH       = 7
  integer(ip), parameter :: TIME_BOTTOMUP     = 8
  integer(ip), parameter :: TIME_UPDADJNCY    = 9
  integer(ip), parameter :: TIME_FINDINODES   = 10
  integer(ip), parameter :: TIME_UPDSCORE     = 11

  ! --- params.h ---
  integer(ip), parameter :: MAX_BAD_FLIPS        = 100
  real(dp),    parameter :: COMPRESS_FRACTION     = 0.75_dp
  integer(ip), parameter :: MIN_NODES            = 100
  integer(ip), parameter :: DEFAULT_SEPS         = 31
  integer(ip), parameter :: MAX_SEPS             = 255
  integer(ip), parameter :: MIN_DOMAINS          = 100
  integer(ip), parameter :: MAX_COARSENING_STEPS = 10

  ! --- eval.h constants (eval1 function) ---
  real(dp), parameter :: TOL1 = 0.50_dp
  real(dp), parameter :: PEN1 = 100.0_dp

contains

  ! Separator evaluation function F(S, B, W) = eval1(S, B, W)
  ! (from eval.h, the default separator cost function)
  pure real(dp) function eval_sep(S, B, W)
    integer(ip), intent(in) :: S, B, W
    integer(ip) :: bmax, bmin
    bmax = max(B, W)
    bmin = min(B, W)
    ! Match C eval1 exactly: S + PEN1*max(0, max(W,B)*(1-TOL1) - min(W,B))
    !                          + (max-min)/max   (all done in floating point)
    eval_sep = real(S, dp) &
             + PEN1 * max(0.0_dp, real(bmax, dp) * (1.0_dp - TOL1) - real(bmin, dp)) &
             + real(bmax - bmin, dp) / real(max(1_ip, bmax), dp)
  end function eval_sep

  ! Timer helpers (using cpu_time)
  subroutine pord_resettimer(t)
    real(dp), intent(out) :: t
    t = 0.0_dp
  end subroutine pord_resettimer

  subroutine pord_starttimer(t)
    real(dp), intent(inout) :: t
    real :: tc
    call cpu_time(tc)
    t = t - real(tc, dp)
  end subroutine pord_starttimer

  subroutine pord_stoptimer(t)
    real(dp), intent(inout) :: t
    real :: tc
    call cpu_time(tc)
    t = t + real(tc, dp)
  end subroutine pord_stoptimer

end module mod_monolis_pord_const
!*****************************************************************************
! SPACE (SPArse Cholesky Elimination) Library: mod_monolis_pord_types.f90
!
! All Fortran derived-type definitions corresponding to C structs.
! Uses 0-based array indexing to match the C source directly.
!*****************************************************************************

module mod_monolis_pord_types
  use mod_monolis_pord_const
  implicit none

  !---------------------------------------------------------------------------
  ! graph_t  (owns all data)
  !---------------------------------------------------------------------------
  type :: graph_t
    integer(ip) :: nvtx     = 0
    integer(ip) :: nedges   = 0
    integer(ip) :: gtype    = UNWEIGHTED   ! "type" is reserved in some compilers
    integer(ip) :: totvwght = 0
    integer(ip), allocatable :: xadj(:)    ! 0:nvtx
    integer(ip), allocatable :: adjncy(:)  ! 0:nedges-1  (may grow in gelim)
    integer(ip), allocatable :: vwght(:)   ! 0:nvtx-1
  end type graph_t

  !---------------------------------------------------------------------------
  ! elimtree_t  (owns all data)
  !---------------------------------------------------------------------------
  type :: elimtree_t
    integer(ip) :: nvtx    = 0
    integer(ip) :: nfronts = 0
    integer(ip) :: root    = -1
    integer(ip), allocatable :: ncolfactor(:)  ! 0:nfronts-1
    integer(ip), allocatable :: ncolupdate(:)  ! 0:nfronts-1
    integer(ip), allocatable :: parent(:)      ! 0:nfronts-1
    integer(ip), allocatable :: firstchild(:)  ! 0:nfronts-1
    integer(ip), allocatable :: silbings(:)    ! 0:nfronts-1
    integer(ip), allocatable :: vtx2front(:)   ! 0:nvtx-1
  end type elimtree_t

  !---------------------------------------------------------------------------
  ! bucket_t  (owns all data)
  !---------------------------------------------------------------------------
  type :: bucket_t
    integer(ip) :: maxbin  = 0
    integer(ip) :: maxitem = 0
    integer(ip) :: offset  = 0
    integer(ip) :: nobj    = 0
    integer(ip) :: minbin  = 0    ! = MAX_INT when empty
    integer(ip), allocatable :: bin(:)    ! 0:maxbin
    integer(ip), allocatable :: nextp(:)  ! 0:maxitem  (C: next)
    integer(ip), allocatable :: lastp(:)  ! 0:maxitem  (C: last)
    integer(ip), allocatable :: key(:)    ! 0:maxitem
  end type bucket_t

  !---------------------------------------------------------------------------
  ! gelim_t  (owns graph and arrays)
  !---------------------------------------------------------------------------
  type :: gelim_t
    type(graph_t) :: G
    integer(ip) :: maxedges = 0
    integer(ip), allocatable :: len(:)    ! 0:nvtx-1
    integer(ip), allocatable :: elen(:)   ! 0:nvtx-1
    integer(ip), allocatable :: parent(:) ! 0:nvtx-1
    integer(ip), allocatable :: degree(:) ! 0:nvtx-1
    integer(ip), allocatable :: score(:)  ! 0:nvtx-1
  end type gelim_t

  !---------------------------------------------------------------------------
  ! stageinfo_t
  !---------------------------------------------------------------------------
  type :: stageinfo_t
    integer(ip) :: nstep = 0
    integer(ip) :: welim = 0
    integer(ip) :: nzf   = 0
    real(dp)    :: ops   = 0.0_dp
  end type stageinfo_t

  !---------------------------------------------------------------------------
  ! multisector_t
  ! Does NOT own the graph; holds a pointer to the external graph.
  !---------------------------------------------------------------------------
  type :: multisector_t
    type(graph_t), pointer :: G => null()   ! not owned
    integer(ip), allocatable :: stage(:)    ! 0:nvtx-1
    integer(ip) :: nstages   = 0
    integer(ip) :: nnodes    = 0
    integer(ip) :: totmswght = 0
  end type multisector_t

  !---------------------------------------------------------------------------
  ! minprior_t
  !---------------------------------------------------------------------------
  type :: minprior_t
    type(gelim_t),       pointer :: Gelim  => null()
    type(multisector_t), pointer :: ms     => null()
    type(bucket_t),      pointer :: bucket => null()
    type(stageinfo_t), allocatable :: stageinfo(:)  ! 0:nstages-1
    integer(ip), allocatable :: reachset(:)  ! 0:nvtx-1
    integer(ip), allocatable :: auxaux(:)    ! 0:nvtx-1
    integer(ip), allocatable :: auxbin(:)    ! 0:nvtx-1
    integer(ip), allocatable :: auxtmp(:)    ! 0:nvtx-1
    integer(ip) :: nreach = 0
    integer(ip) :: flag   = 1
  end type minprior_t

  !---------------------------------------------------------------------------
  ! domdec_t  (owns graph; linked list via pointer components)
  !---------------------------------------------------------------------------
  type :: domdec_t
    type(graph_t) :: G
    integer(ip) :: ndom    = 0
    integer(ip) :: domwght = 0
    integer(ip) :: cwght(0:2) = 0
    integer(ip), allocatable :: vtype(:)   ! 0:nvtx-1
    integer(ip), allocatable :: color(:)   ! 0:nvtx-1
    integer(ip), allocatable :: map_(:)    ! 0:nvtx-1  (C: map)
    type(domdec_t), pointer  :: prev => null()
    type(domdec_t), pointer  :: next_ptr => null()
  end type domdec_t

  !---------------------------------------------------------------------------
  ! gbisect_t
  ! Does NOT own the graph; pointer to external graph.
  !---------------------------------------------------------------------------
  type :: gbisect_t
    type(graph_t), pointer :: G => null()  ! not owned
    integer(ip), allocatable :: color(:)   ! 0:nvtx-1
    integer(ip) :: cwght(0:2) = 0
  end type gbisect_t

  !---------------------------------------------------------------------------
  ! gbipart_t  (owns graph)
  !---------------------------------------------------------------------------
  type :: gbipart_t
    type(graph_t) :: G
    integer(ip) :: nX = 0
    integer(ip) :: nY = 0
  end type gbipart_t

  !---------------------------------------------------------------------------
  ! nestdiss_t  (recursive binary tree; G and map_ are NOT owned)
  !---------------------------------------------------------------------------
  type :: nestdiss_t
    type(graph_t), pointer   :: G       => null()   ! not owned
    integer(ip),   pointer   :: map_arr(:) => null()  ! not owned (0:nvtx-1)
    integer(ip) :: depth = 0
    integer(ip) :: nvint = 0
    integer(ip), allocatable :: intvertex(:)  ! 0:nvint-1
    integer(ip), allocatable :: intcolor(:)   ! 0:nvint-1
    integer(ip) :: cwght(0:2) = 0
    type(nestdiss_t), pointer :: parent => null()
    type(nestdiss_t), pointer :: childB => null()
    type(nestdiss_t), pointer :: childW => null()
  end type nestdiss_t

end module mod_monolis_pord_types
!*****************************************************************************
! SPACE (SPArse Cholesky Elimination) Library: mod_monolis_pord_sort.f90
!
! Fortran implementation of sort.c (only distributionCounting is needed
! by constructDomainDecomposition).
!*****************************************************************************

module mod_monolis_pord_sort
  use mod_monolis_pord_const
  use mod_monolis_pord_types
  implicit none
  private

  public :: distributionCounting
  public :: insertUpInts
  public :: insertUpIntsWithStaticIntKeys
  public :: qsortUpInts

  integer(ip), parameter :: THRES = 10

contains

  !===========================================================================
  ! distributionCounting: sorts node(0:n-1) ascending by key(node(i)).
  !===========================================================================
  subroutine distributionCounting(n, node, key)
    integer(ip), intent(in)    :: n
    integer(ip), intent(inout) :: node(0:n-1)
    integer(ip), intent(in)    :: key(0:*)   ! key(u) for any u in node

    integer(ip), allocatable :: tmp(:), count(:)
    integer(ip) :: minkey, maxkey, l, i, u, vk

    if (n <= 0) return

    ! determine key range
    minkey = MAX_INT; maxkey = 0
    do i = 0, n-1
      u = node(i)
      if (key(u) > maxkey) maxkey = key(u)
      if (key(u) < minkey) minkey = key(u)
    end do
    l = maxkey - minkey

    allocate(count(0:l), tmp(0:n-1))
    count = 0

    do i = 0, n-1
      u = node(i)
      vk = key(u) - minkey
      count(vk) = count(vk) + 1
    end do

    ! prefix sum
    do i = 1, l
      count(i) = count(i) + count(i-1)
    end do

    ! fill tmp in reverse to be stable
    do i = n-1, 0, -1
      u = node(i)
      vk = key(u) - minkey
      count(vk) = count(vk) - 1
      tmp(count(vk)) = u
    end do

    do i = 0, n-1
      node(i) = tmp(i)
    end do

    deallocate(count, tmp)
  end subroutine distributionCounting

  !===========================================================================
  ! myrandom: returns a pseudo-random integer in [0, range-1]
  !===========================================================================
  integer(ip) function pord_random(range)
    integer(ip), intent(in) :: range
    real :: r
    call random_number(r)
    pord_random = int(r * real(range))
    if (pord_random >= range) pord_random = range - 1
    if (pord_random < 0) pord_random = 0
  end function pord_random

  !===========================================================================
  ! insertUpInts: insertion sort upwards (without keys).
  !===========================================================================
  subroutine insertUpInts(n, array)
    integer(ip), intent(in)    :: n
    integer(ip), intent(inout) :: array(0:n-1)
    integer(ip) :: i, j, v
    do i = 1, n-1
      v = array(i); j = i
      do while (j > 0)
        if (array(j-1) <= v) exit
        array(j) = array(j-1)
        j = j - 1
      end do
      array(j) = v
    end do
  end subroutine insertUpInts

  !===========================================================================
  ! insertUpIntsWithStaticIntKeys: insertion sort upwards keyed by key(array(i)).
  !===========================================================================
  subroutine insertUpIntsWithStaticIntKeys(n, array, key)
    integer(ip), intent(in)    :: n
    integer(ip), intent(inout) :: array(0:n-1)
    integer(ip), intent(in)    :: key(0:*)
    integer(ip) :: i, j, e, ke
    do i = 1, n-1
      e = array(i); ke = key(e); j = i
      do while (j > 0)
        if (key(array(j-1)) <= ke) exit
        array(j) = array(j-1)
        j = j - 1
      end do
      array(j) = e
    end do
  end subroutine insertUpIntsWithStaticIntKeys

  !===========================================================================
  ! qsortUpInts: median-of-three quicksort upwards (without keys).
  ! stack must be at least 2*log2(n)+4 long; we allocate internally to keep
  ! the C signature transparent.
  !===========================================================================
  subroutine qsortUpInts(n, array, stack)
    integer(ip), intent(in)    :: n
    integer(ip), intent(inout) :: array(0:n-1)
    integer(ip), intent(inout) :: stack(0:*)

    integer(ip) :: i, j, t, l, m, r, p, mv

    l = 0; r = n-1; p = 2
    do while (p > 0)
      if ((r - l) > THRES) then
        m = l + ishft(r - l, -1)
        if (array(l) > array(r)) then; t = array(l); array(l) = array(r); array(r) = t; end if
        if (array(l) > array(m)) then; t = array(l); array(l) = array(m); array(m) = t; end if
        if (array(r) > array(m)) then; t = array(m); array(m) = array(r); array(r) = t; end if
        mv = array(r); i = l - 1; j = r
        do
          do
            i = i + 1
            if (array(i) >= mv) exit
          end do
          do
            j = j - 1
            if (array(j) <= mv) exit
          end do
          if (i >= j) exit
          t = array(i); array(i) = array(j); array(j) = t
        end do
        t = array(i); array(i) = array(r); array(r) = t
        if ((i - l) > (r - i)) then
          stack(p) = l;     p = p + 1
          stack(p) = i - 1; p = p + 1
          l = i + 1
        else
          stack(p) = i + 1; p = p + 1
          stack(p) = r;     p = p + 1
          r = i - 1
        end if
      else
        p = p - 1; r = stack(p)
        p = p - 1; l = stack(p)
      end if
    end do
    if (THRES > 0) call insertUpInts(n, array)
  end subroutine qsortUpInts

end module mod_monolis_pord_sort
!*****************************************************************************
! SPACE (SPArse Cholesky Elimination) Library: mod_monolis_pord_bucket.f90
!
! Fortran implementation of bucket.c
!*****************************************************************************

module mod_monolis_pord_bucket
  use mod_monolis_pord_const
  use mod_monolis_pord_types
  implicit none
  private

  public :: setupBucket, freeBucket, minBucket, insertBucket, removeBucket

contains

  !===========================================================================
  subroutine setupBucket(bucket, maxbin, maxitem, offset)
    type(bucket_t), intent(out) :: bucket
    integer(ip),    intent(in)  :: maxbin, maxitem, offset

    integer(ip) :: i

    bucket%maxbin  = maxbin
    bucket%maxitem = maxitem
    bucket%offset  = offset
    bucket%nobj    = 0
    bucket%minbin  = MAX_INT

    allocate(bucket%bin(0:maxbin))
    allocate(bucket%nextp(0:maxitem))
    allocate(bucket%lastp(0:maxitem))
    allocate(bucket%key(0:maxitem))

    do i = 0, maxbin
      bucket%bin(i) = -1
    end do
    do i = 0, maxitem
      bucket%nextp(i) = -1
      bucket%lastp(i) = -1
      bucket%key(i)   = MAX_INT
    end do
  end subroutine setupBucket

  !===========================================================================
  subroutine freeBucket(bucket)
    type(bucket_t), intent(inout) :: bucket
    if (allocated(bucket%bin))   deallocate(bucket%bin)
    if (allocated(bucket%nextp)) deallocate(bucket%nextp)
    if (allocated(bucket%lastp)) deallocate(bucket%lastp)
    if (allocated(bucket%key))   deallocate(bucket%key)
    bucket%nobj = 0; bucket%minbin = MAX_INT
  end subroutine freeBucket

  !===========================================================================
  ! minBucket: returns the item with minimum key, or -1 if empty.
  !===========================================================================
  integer(ip) function minBucket(bucket)
    type(bucket_t), intent(inout) :: bucket

    integer(ip) :: maxbin, minbin, nobj
    integer(ip) :: item, bestitem, bestkey

    maxbin = bucket%maxbin
    nobj   = bucket%nobj
    minbin = bucket%minbin

    if (nobj <= 0) then
      minBucket = -1
      return
    end if

    ! advance minbin to leftmost non-empty bin
    do while (bucket%bin(minbin) == -1)
      minbin = minbin + 1
    end do
    bucket%minbin = minbin

    bestitem = bucket%bin(minbin)
    bestkey  = minbin

    ! bins 0 and maxbin may hold items with different actual keys
    if ((minbin == 0) .or. (minbin == maxbin)) then
      item = bucket%nextp(bestitem)
      do while (item /= -1)
        if (bucket%key(item) < bestkey) then
          bestitem = item
          bestkey  = bucket%key(item)
        end if
        item = bucket%nextp(item)
      end do
    end if

    minBucket = bestitem
  end function minBucket

  !===========================================================================
  subroutine insertBucket(bucket, k, item)
    type(bucket_t), intent(inout) :: bucket
    integer(ip),    intent(in)    :: k, item

    integer(ip) :: s, nextitem

    s = max(0, k + bucket%offset)
    s = min(s, bucket%maxbin)

    bucket%minbin = min(bucket%minbin, s)
    bucket%nobj   = bucket%nobj + 1
    bucket%key(item) = k

    ! insert at head of bin s
    nextitem = bucket%bin(s)
    bucket%bin(s) = item
    bucket%nextp(item) = nextitem
    bucket%lastp(item) = -1
    if (nextitem /= -1) bucket%lastp(nextitem) = item
  end subroutine insertBucket

  !===========================================================================
  subroutine removeBucket(bucket, item)
    type(bucket_t), intent(inout) :: bucket
    integer(ip),    intent(in)    :: item

    integer(ip) :: s, previtem, nextitem

    if (bucket%key(item) == MAX_INT) return  ! already not in bucket

    s = max(0, bucket%key(item) + bucket%offset)
    s = min(s, bucket%maxbin)

    previtem = bucket%lastp(item)
    nextitem = bucket%nextp(item)

    if (previtem == -1) then
      bucket%bin(s) = nextitem
    else
      bucket%nextp(previtem) = nextitem
    end if
    if (nextitem /= -1) bucket%lastp(nextitem) = previtem

    bucket%nextp(item) = -1
    bucket%lastp(item) = -1
    bucket%key(item)   = MAX_INT
    bucket%nobj        = bucket%nobj - 1
  end subroutine removeBucket

end module mod_monolis_pord_bucket
