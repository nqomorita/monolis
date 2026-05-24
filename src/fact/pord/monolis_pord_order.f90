!*****************************************************************************
! SPACE (SPArse Cholesky Elimination) Library: mod_monolis_pord_gelim.f90
!
! Fortran implementation of gelim.c
! Key structures: gelim_t (elimination graph)
!*****************************************************************************

module mod_monolis_pord_gelim
  use mod_monolis_pord_const
  use mod_monolis_pord_types
  use mod_monolis_pord_graph
  use mod_monolis_pord_tree
  implicit none
  private

  public :: setupElimGraph, freeElimGraph
  public :: crunchElimGraph
  public :: buildElement
  public :: updateAdjncy
  public :: findIndNodes
  public :: updateDegree
  public :: updateScore
  public :: extractElimTree

  ! Grow factor for adjncy array reallocation
  integer(ip), parameter :: GROW_FACTOR = 2

contains

  !===========================================================================
  subroutine freeElimGraph(Gelim)
    type(gelim_t), intent(inout) :: Gelim
    call freeGraph(Gelim%G)
    if (allocated(Gelim%len))    deallocate(Gelim%len)
    if (allocated(Gelim%elen))   deallocate(Gelim%elen)
    if (allocated(Gelim%parent)) deallocate(Gelim%parent)
    if (allocated(Gelim%degree)) deallocate(Gelim%degree)
    if (allocated(Gelim%score))  deallocate(Gelim%score)
    Gelim%maxedges = 0
  end subroutine freeElimGraph

  !===========================================================================
  subroutine setupElimGraph(Gelim, G)
    type(gelim_t), intent(out) :: Gelim
    type(graph_t), intent(in)  :: G

    integer(ip) :: nvtx, nedges, deg, u, i, istart, istop

    nvtx   = G%nvtx
    nedges = G%nedges

    ! allocate with extra space (nvtx more than needed)
    Gelim%maxedges = nedges + nvtx
    call newGraph(Gelim%G, nvtx, Gelim%maxedges)

    ! grow adjncy to maxedges size
    deallocate(Gelim%G%adjncy)
    allocate(Gelim%G%adjncy(0:Gelim%maxedges-1))
    Gelim%G%adjncy = 0

    allocate(Gelim%len(0:nvtx-1))
    allocate(Gelim%elen(0:nvtx-1))
    allocate(Gelim%parent(0:nvtx-1))
    allocate(Gelim%degree(0:nvtx-1))
    allocate(Gelim%score(0:nvtx-1))

    ! copy graph
    Gelim%G%gtype    = G%gtype
    Gelim%G%totvwght = G%totvwght
    Gelim%G%nedges   = nedges
    do u = 0, nvtx
      Gelim%G%xadj(u) = G%xadj(u)
    end do
    do u = 0, nvtx-1
      Gelim%G%vwght(u) = G%vwght(u)
    end do
    do i = 0, nedges-1
      Gelim%G%adjncy(i) = G%adjncy(i)
    end do

    ! initialize per-vertex data
    do u = 0, nvtx-1
      istart = G%xadj(u); istop = G%xadj(u+1)
      Gelim%len(u)    = istop - istart
      Gelim%elen(u)   = 0
      Gelim%parent(u) = -1
      Gelim%score(u)  = -1

      select case (Gelim%G%gtype)
        case (UNWEIGHTED)
          deg = Gelim%len(u)
        case (WEIGHTED)
          deg = 0
          do i = istart, istop-1
            deg = deg + G%vwght(G%adjncy(i))
          end do
        case default
          deg = Gelim%len(u)
      end select
      Gelim%degree(u) = deg

      if (Gelim%len(u) == 0) Gelim%G%xadj(u) = -1
    end do
  end subroutine setupElimGraph

  !===========================================================================
  ! crunchElimGraph: compacts the adjacency array in-place.
  ! Returns .true. if successful (i.e. freed space).
  !===========================================================================
  logical function crunchElimGraph(Gelim)
    type(gelim_t), intent(inout) :: Gelim

    integer(ip) :: nvtx, nedges_old, u, i, isrc, idest

    nvtx      = Gelim%G%nvtx
    nedges_old = Gelim%G%nedges

    ! mark heads of adjacency lists
    do u = 0, nvtx-1
      i = Gelim%G%xadj(u)
      if (i /= -1) then
        Gelim%G%xadj(u)   = Gelim%G%adjncy(i)
        Gelim%G%adjncy(i) = -(u+1)
      end if
    end do

    ! compact
    idest = 0; isrc = 0
    do while (isrc < Gelim%G%nedges)
      u = Gelim%G%adjncy(isrc)
      isrc = isrc + 1
      if (u < 0) then
        u = -u - 1
        Gelim%G%adjncy(idest) = Gelim%G%xadj(u)
        Gelim%G%xadj(u) = idest
        idest = idest + 1
        do i = 1, Gelim%len(u)-1
          Gelim%G%adjncy(idest) = Gelim%G%adjncy(isrc)
          idest = idest + 1
          isrc  = isrc + 1
        end do
      end if
    end do
    Gelim%G%nedges = idest
    crunchElimGraph = (idest < nedges_old)
  end function crunchElimGraph

  !===========================================================================
  ! buildElement: turns variable me into an element.
  !===========================================================================
  subroutine buildElement(Gelim, me)
    type(gelim_t), intent(inout) :: Gelim
    integer(ip),   intent(in)    :: me

    integer(ip) :: degme, elenme, vlenme, mesrcptr, medeststart, medeststart2
    integer(ip) :: medestptr, ln, p, i, j, v, e
    logical :: ok

    degme   = 0
    Gelim%G%totvwght = Gelim%G%totvwght - Gelim%G%vwght(me)
    Gelim%G%vwght(me) = -Gelim%G%vwght(me)
    Gelim%score(me) = -3

    elenme  = Gelim%elen(me)
    vlenme  = Gelim%len(me) - elenme
    mesrcptr = Gelim%G%xadj(me)

    if (elenme == 0) then
      ! leaf case: build in-place
      medeststart = Gelim%G%xadj(me)
      medestptr   = medeststart
      do i = 0, vlenme-1
        v = Gelim%G%adjncy(mesrcptr + i)
        if (Gelim%G%vwght(v) > 0) then
          degme = degme + Gelim%G%vwght(v)
          Gelim%G%vwght(v) = -Gelim%G%vwght(v)
          Gelim%G%adjncy(medestptr) = v
          medestptr = medestptr + 1
        end if
      end do
    else
      ! non-leaf: build in empty space
      medeststart = Gelim%G%nedges
      medestptr   = medeststart

      do i = 0, elenme    ! iterate over elements (0..elenme-1) then variables
        if (i < elenme) then
          Gelim%len(me) = Gelim%len(me) - 1
          e   = Gelim%G%adjncy(mesrcptr)
          mesrcptr = mesrcptr + 1
          p   = Gelim%G%xadj(e)
          ln  = Gelim%len(e)
        else
          e   = me
          p   = mesrcptr
          ln  = vlenme
        end if

        do j = 0, ln-1
          Gelim%len(e) = Gelim%len(e) - 1
          v = Gelim%G%adjncy(p)         ! read current position
          p = p + 1                     ! match C's "p++"
          if (Gelim%G%vwght(v) > 0) then
            degme = degme + Gelim%G%vwght(v)
            Gelim%G%vwght(v) = -Gelim%G%vwght(v)

            ! grow array if needed
            if (medestptr == Gelim%maxedges) then
              if (Gelim%len(me) == 0) then
                Gelim%G%xadj(me) = -1
              else
                Gelim%G%xadj(me) = mesrcptr
              end if
              if (Gelim%len(e) == 0) then
                Gelim%G%xadj(e) = -1
              else
                Gelim%G%xadj(e) = p
              end if

              ok = crunchElimGraph(Gelim)
              if (.not. ok) then
                call grow_adjncy(Gelim)
              end if

              ! re-position partially built element
              medeststart2 = Gelim%G%nedges
              do p = medeststart, medestptr-1
                Gelim%G%adjncy(Gelim%G%nedges) = Gelim%G%adjncy(p)
                Gelim%G%nedges = Gelim%G%nedges + 1
              end do
              medeststart = medeststart2
              medestptr   = Gelim%G%nedges

              mesrcptr = Gelim%G%xadj(me)
              p = Gelim%G%xadj(e)
            end if

            Gelim%G%adjncy(medestptr) = v
            medestptr = medestptr + 1
          end if
        end do

        if (e /= me) then
          Gelim%G%xadj(e) = -1
          Gelim%parent(e) = me
          Gelim%score(e) = -4
        end if
      end do
      Gelim%G%nedges = medestptr
    end if

    Gelim%degree(me) = degme
    Gelim%G%xadj(me) = medeststart
    Gelim%G%vwght(me) = -Gelim%G%vwght(me)
    Gelim%elen(me) = 0
    Gelim%len(me) = medestptr - medeststart
    if (Gelim%len(me) == 0) Gelim%G%xadj(me) = -1

    ! unmark
    if (Gelim%G%xadj(me) /= -1) then
      do i = 0, Gelim%len(me)-1
        v = Gelim%G%adjncy(Gelim%G%xadj(me) + i)
        Gelim%G%vwght(v) = -Gelim%G%vwght(v)
      end do
    end if
  end subroutine buildElement

  !===========================================================================
  ! Grow the adjacency array by GROW_FACTOR.
  ! IMPORTANT: must preserve the FULL contents of adjncy (not just the first
  ! nedges entries), because buildElement keeps a partially-constructed element
  ! at positions [medeststart..medestptr-1] which lie ABOVE nedges.
  !===========================================================================
  subroutine grow_adjncy(Gelim)
    type(gelim_t), intent(inout) :: Gelim
    integer(ip), allocatable :: tmp(:)
    integer(ip) :: new_size, old_size

    old_size = Gelim%maxedges
    new_size = old_size * GROW_FACTOR
    allocate(tmp(0:new_size-1))
    tmp(0:old_size-1) = Gelim%G%adjncy(0:old_size-1)
    call move_alloc(tmp, Gelim%G%adjncy)
    Gelim%maxedges = new_size
  end subroutine grow_adjncy

  !===========================================================================
  ! updateAdjncy: updates adjacency lists of all vertices in reachset.
  !===========================================================================
  subroutine updateAdjncy(Gelim, reachset, nreach, tmp, pflag)
    type(gelim_t), intent(inout) :: Gelim
    integer(ip),   intent(in)    :: nreach
    integer(ip),   intent(in)    :: reachset(0:nreach-1)
    integer(ip),   intent(inout) :: tmp(0:Gelim%G%nvtx-1)
    integer(ip),   intent(inout) :: pflag

    integer(ip) :: u, v, e, me
    integer(ip) :: i, j, jdest, jfirstolde, jfirstv, jstart, jstop
    integer(ip) :: covered, marku
    integer(ip) :: jj

    ! mark all variables in reachset
    do i = 0, nreach-1
      u = reachset(i)
      Gelim%G%vwght(u) = -Gelim%G%vwght(u)
    end do

    ! update element/variable lists
    do i = 0, nreach-1
      u = reachset(i)
      jstart     = Gelim%G%xadj(u)
      jstop      = jstart + Gelim%len(u)
      jdest      = jstart
      jfirstolde = jstart

      ! scan element list
      do j = jstart, jstart + Gelim%elen(u) - 1
        e = Gelim%G%adjncy(j)
        if (Gelim%score(e) == -4) then
          me = Gelim%parent(e)
          if (tmp(me) < pflag) then
            Gelim%G%adjncy(jdest) = Gelim%G%adjncy(jfirstolde)
            Gelim%G%adjncy(jfirstolde) = me
            jfirstolde = jfirstolde + 1
            jdest = jdest + 1
            tmp(me) = pflag
          end if
        else
          if (tmp(e) < pflag) then
            Gelim%G%adjncy(jdest) = e
            jdest = jdest + 1
            tmp(e) = pflag
          end if
        end if
      end do
      jfirstv = jdest

      ! scan variable list
      do j = jstart + Gelim%elen(u), jstop-1
        v = Gelim%G%adjncy(j)
        if (Gelim%score(v) == -3) then
          if (tmp(v) < pflag) then
            Gelim%G%adjncy(jdest) = Gelim%G%adjncy(jfirstv)
            Gelim%G%adjncy(jfirstv) = Gelim%G%adjncy(jfirstolde)
            Gelim%G%adjncy(jfirstolde) = v
            jfirstolde = jfirstolde + 1
            jfirstv    = jfirstv + 1
            jdest      = jdest + 1
            tmp(v) = pflag
          end if
        else
          Gelim%G%adjncy(jdest) = v
          jdest = jdest + 1
        end if
      end do

      Gelim%elen(u) = jfirstv - jstart
      Gelim%len(u)  = jdest - jstart
      pflag = pflag + 1
    end do

    ! remove covered edges
    do i = 0, nreach-1
      u = reachset(i)
      jstart = Gelim%G%xadj(u)
      jstop  = jstart + Gelim%len(u)
      marku  = 0

      jdest = jstart + Gelim%elen(u)
      do j = jstart + Gelim%elen(u), jstop-1
        v = Gelim%G%adjncy(j)
        if (Gelim%G%vwght(v) > 0) then
          Gelim%G%adjncy(jdest) = v
          jdest = jdest + 1
        else if (Gelim%G%vwght(v) < 0) then
          if (marku == 0) then
            do jj = jstart, jstart + Gelim%elen(u) - 1
              tmp(Gelim%G%adjncy(jj)) = pflag
            end do
            marku = 1
          end if
          covered = 0
          do jj = Gelim%G%xadj(v), Gelim%G%xadj(v) + Gelim%elen(v) - 1
            if (tmp(Gelim%G%adjncy(jj)) == pflag) then
              covered = 1
              exit
            end if
          end do
          if (covered == 0) then
            Gelim%G%adjncy(jdest) = v
            jdest = jdest + 1
          end if
        end if
      end do
      Gelim%len(u) = jdest - jstart
      pflag = pflag + 1
    end do

    ! unmark
    do i = 0, nreach-1
      u = reachset(i)
      Gelim%G%vwght(u) = -Gelim%G%vwght(u)
    end do
  end subroutine updateAdjncy

  !===========================================================================
  ! findIndNodes: detects indistinguishable vertices in reachset.
  !===========================================================================
  subroutine findIndNodes(Gelim, reachset, nreach, bin, nextn, tmp, pflag)
    type(gelim_t), intent(inout) :: Gelim
    integer(ip),   intent(in)    :: nreach
    integer(ip),   intent(in)    :: reachset(0:nreach-1)
    integer(ip),   intent(inout) :: bin(0:Gelim%G%nvtx-1)
    integer(ip),   intent(inout) :: nextn(0:Gelim%G%nvtx-1)
    integer(ip),   intent(inout) :: tmp(0:Gelim%G%nvtx-1)
    integer(ip),   intent(inout) :: pflag

    integer(ip) :: nvtx, chk, u, v, w, wlast, i, j, jstart, jstop
    integer(ip) :: jstep, jj, jjstop
    logical     :: keepon

    nvtx = Gelim%G%nvtx

    ! compute checksums
    do i = 0, nreach-1
      u = reachset(i)
      chk = 0
      jstart = Gelim%G%xadj(u)
      jstop  = jstart + Gelim%len(u)
      jstep  = max(1, 1000000000 / max(1, nvtx))
      j = jstart
      do while (j < jstop)
        jjstop = min(jstop, j + jstep)
        do jj = j, jjstop-1
          chk = chk + Gelim%G%adjncy(jj)
        end do
        chk = mod(chk, nvtx)
        j = j + jstep
      end do
      Gelim%parent(u) = chk
      nextn(u) = bin(chk)
      bin(chk) = u
    end do

    ! detect supervariables
    do i = 0, nreach-1
      u = reachset(i)
      if (Gelim%G%vwght(u) > 0) then
        chk = Gelim%parent(u)
        v = bin(chk)
        bin(chk) = -1   ! process each bin once
        do while (v /= -1)
          jstart = Gelim%G%xadj(v)
          jstop  = jstart + Gelim%len(v)
          do j = jstart, jstop-1
            tmp(Gelim%G%adjncy(j)) = pflag
          end do
          w     = nextn(v)
          wlast = v
          do while (w /= -1)
            keepon = .true.
            if ((Gelim%len(w) /= Gelim%len(v)) .or. &
                (Gelim%elen(w) /= Gelim%elen(v))) keepon = .false.
            if (keepon .and. (Gelim%score(w) < 0) .neqv. (Gelim%score(v) < 0)) &
              keepon = .false.
            if (keepon) then
              do jj = Gelim%G%xadj(w), Gelim%G%xadj(w) + Gelim%len(w) - 1
                if (tmp(Gelim%G%adjncy(jj)) < pflag) then
                  keepon = .false.; exit
                end if
              end do
            end if
            if (keepon) then
              Gelim%parent(w) = v
              Gelim%G%vwght(v) = Gelim%G%vwght(v) + Gelim%G%vwght(w)
              Gelim%G%vwght(w) = 0
              Gelim%G%xadj(w) = -1
              Gelim%score(w) = -2
              w = nextn(w)
              nextn(wlast) = w
            else
              wlast = w; w = nextn(w)
            end if
          end do
          v = nextn(v)
          pflag = pflag + 1
        end do
      end if
    end do

    ! reset parent for principal variables
    do i = 0, nreach-1
      u = reachset(i)
      if (Gelim%G%vwght(u) > 0) Gelim%parent(u) = -1
    end do
  end subroutine findIndNodes

  !===========================================================================
  ! updateDegree: approximate degree update.
  !===========================================================================
  subroutine updateDegree(Gelim, reachset, nreach, bin)
    type(gelim_t), intent(inout) :: Gelim
    integer(ip),   intent(in)    :: nreach
    integer(ip),   intent(in)    :: reachset(0:nreach-1)
    integer(ip),   intent(inout) :: bin(0:Gelim%G%nvtx-1)

    integer(ip) :: totvwght, deg, vwghtv, u, v, w, e, me, r
    integer(ip) :: i, istart, istop, j, jstart, jstop

    totvwght = Gelim%G%totvwght

    ! mark vertices adjacent to any element
    do r = 0, nreach-1
      u = reachset(r)
      if (Gelim%elen(u) > 0) bin(u) = 1
    end do

    do r = 0, nreach-1
      u = reachset(r)
      if (bin(u) == 1) then
        me = Gelim%G%adjncy(Gelim%G%xadj(u))  ! most recently formed element
        istart = Gelim%G%xadj(me)
        istop  = istart + Gelim%len(me)

        ! compute bin[e] = |Le \ Lme| for elements e != me
        do i = istart, istop-1
          v = Gelim%G%adjncy(i)
          vwghtv = Gelim%G%vwght(v)
          if (vwghtv > 0) then
            jstart = Gelim%G%xadj(v); jstop = jstart + Gelim%elen(v)
            do j = jstart, jstop-1
              e = Gelim%G%adjncy(j)
              if (e /= me) then
                if (bin(e) > 0) then
                  bin(e) = bin(e) - vwghtv
                else
                  bin(e) = Gelim%degree(e) - vwghtv
                end if
              end if
            end do
          end if
        end do

        ! update degree for principals in Lme
        do i = istart, istop-1
          v = Gelim%G%adjncy(i)
          vwghtv = Gelim%G%vwght(v)
          if (bin(v) == 1) then
            deg = 0
            jstart = Gelim%G%xadj(v); jstop = jstart + Gelim%len(v)
            do j = jstart, jstart + Gelim%elen(v) - 1
              e = Gelim%G%adjncy(j)
              if (e /= me) deg = deg + bin(e)
            end do
            do j = jstart + Gelim%elen(v), jstop-1
              w = Gelim%G%adjncy(j)
              deg = deg + Gelim%G%vwght(w)
            end do
            deg = min(Gelim%degree(v), deg)
            Gelim%degree(v) = max(1_ip, min(deg + Gelim%degree(me) - vwghtv, &
                                            totvwght - vwghtv))
            bin(v) = -1
          end if
        end do

        ! clear bin[e] for elements e != me
        do i = istart, istop-1
          v = Gelim%G%adjncy(i)
          if (Gelim%G%vwght(v) > 0) then
            jstart = Gelim%G%xadj(v); jstop = jstart + Gelim%elen(v)
            do j = jstart, jstop-1
              e = Gelim%G%adjncy(j)
              if (e /= me) bin(e) = -1
            end do
          end if
        end do
      end if
    end do
  end subroutine updateDegree

  !===========================================================================
  ! updateScore: updates score values.
  !===========================================================================
  subroutine updateScore(Gelim, reachset, nreach, scoretype, bin)
    type(gelim_t), intent(inout) :: Gelim
    integer(ip),   intent(in)    :: nreach, scoretype
    integer(ip),   intent(in)    :: reachset(0:nreach-1)
    integer(ip),   intent(inout) :: bin(0:Gelim%G%nvtx-1)

    integer(ip) :: vwghtv, deg, degme, u, v, me, r, i, istart, istop
    integer(ip) :: stype
    real(dp)    :: scr_dbl

    ! mark vertices adjacent to any element
    do r = 0, nreach-1
      u = reachset(r)
      if (Gelim%elen(u) > 0) bin(u) = 1
    end do

    stype = mod(scoretype, 10)

    do r = 0, nreach-1
      u = reachset(r)
      if (bin(u) == 1) then
        me = Gelim%G%adjncy(Gelim%G%xadj(u))
        istart = Gelim%G%xadj(me)
        istop  = istart + Gelim%len(me)

        do i = istart, istop-1
          v = Gelim%G%adjncy(i)
          if (bin(v) == 1) then
            vwghtv = Gelim%G%vwght(v)
            deg    = Gelim%degree(v)
            degme  = Gelim%degree(me) - vwghtv

            if ((deg > 40000) .or. (degme > 40000)) then
              select case (stype)
                case (AMD)
                  scr_dbl = real(deg, dp)
                case (AMF)
                  scr_dbl = real(deg,dp)*real(deg-1,dp)/2.0_dp &
                          - real(degme,dp)*real(degme-1,dp)/2.0_dp
                case (AMMF)
                  scr_dbl = (real(deg,dp)*real(deg-1,dp)/2.0_dp &
                           - real(degme,dp)*real(degme-1,dp)/2.0_dp) &
                           / real(vwghtv, dp)
                case (AMIND)
                  scr_dbl = max(0.0_dp, &
                    real(deg,dp)*real(deg-1,dp)/2.0_dp &
                  - real(degme,dp)*real(degme-1,dp)/2.0_dp &
                  - real(deg,dp)*real(vwghtv,dp))
                case default
                  scr_dbl = real(deg, dp)
              end select
              Gelim%score(v) = int(min(scr_dbl, real(MAX_INT - Gelim%G%nvtx, dp)))
            else
              select case (stype)
                case (AMD)
                  Gelim%score(v) = deg
                case (AMF)
                  Gelim%score(v) = deg*(deg-1)/2 - degme*(degme-1)/2
                case (AMMF)
                  Gelim%score(v) = (deg*(deg-1)/2 - degme*(degme-1)/2) / vwghtv
                case (AMIND)
                  Gelim%score(v) = max(0, (deg*(deg-1)/2 - degme*(degme-1)/2) &
                                        - deg*vwghtv)
                case default
                  Gelim%score(v) = deg
              end select
            end if
            if (Gelim%score(v) < 0) Gelim%score(v) = 0
            bin(v) = -1
          end if
        end do
      end if
    end do
  end subroutine updateScore

  !===========================================================================
  ! extractElimTree: builds the elimination tree from a completed Gelim.
  !===========================================================================
  subroutine extractElimTree(T, Gelim)
    type(elimtree_t), intent(out) :: T
    type(gelim_t),    intent(in)  :: Gelim

    integer(ip) :: nvtx, nfronts, root, u, v, front
    integer(ip), allocatable :: sib(:), fch(:)
    integer(ip), allocatable :: par(:)

    nvtx = Gelim%G%nvtx
    allocate(sib(0:nvtx-1), fch(0:nvtx-1), par(0:nvtx-1))
    sib = -1; fch = -1
    do u = 0, nvtx-1
      par(u) = Gelim%parent(u)
    end do

    ! build top-down tree
    nfronts = 0; root = -1
    do u = 0, nvtx-1
      select case (Gelim%score(u))
        case (-2)   ! nonprincipal: skip
        case (-3)   ! became element (root-level)
          sib(u) = root; root = u; nfronts = nfronts + 1
        case (-4)   ! absorbed element
          v = par(u)
          sib(u) = fch(v); fch(v) = u; nfronts = nfronts + 1
        case default
          write(*,'(A,I0,A,I0)') &
            "Error in extractElimTree: ordering not complete, score(", u, &
            ") = ", Gelim%score(u)
          stop
      end select
    end do

    call newElimTree(T, nvtx, nfronts)

    ! fill vtx2front in post-order
    nfronts = 0
    u = root
    do while (u /= -1)
      do while (fch(u) /= -1)
        u = fch(u)
      end do
      T%vtx2front(u) = nfronts; nfronts = nfronts + 1
      do while ((sib(u) == -1) .and. (par(u) /= -1))
        u = par(u)
        T%vtx2front(u) = nfronts; nfronts = nfronts + 1
      end do
      u = sib(u)
    end do

    ! fill nonprincipal vtx2front
    do u = 0, nvtx-1
      if (Gelim%score(u) == -2) then
        v = u
        do while ((par(v) /= -1) .and. (Gelim%score(v) == -2))
          v = par(v)
        end do
        T%vtx2front(u) = T%vtx2front(v)
      end if
    end do

    ! fill ncolfactor, ncolupdate, parent (assignment, matching C)
    do u = 0, nvtx-1
      front = T%vtx2front(u)
      if (Gelim%score(u) == -3) then
        T%parent(front)     = -1
        T%ncolfactor(front) = Gelim%G%vwght(u)
        T%ncolupdate(front) = Gelim%degree(u)
      else if (Gelim%score(u) == -4) then
        T%parent(front)     = T%vtx2front(par(u))
        T%ncolfactor(front) = Gelim%G%vwght(u)
        T%ncolupdate(front) = Gelim%degree(u)
      end if
    end do

    call initFchSilbRoot(T)
    deallocate(sib, fch, par)
  end subroutine extractElimTree

end module mod_monolis_pord_gelim
!*****************************************************************************
! SPACE (SPArse Cholesky Elimination) Library: mod_monolis_pord_minpriority.f90
!
! Fortran implementation of minpriority.c
!*****************************************************************************

module mod_monolis_pord_minpriority
  use mod_monolis_pord_const
  use mod_monolis_pord_types
  use mod_monolis_pord_bucket
  use mod_monolis_pord_gelim
  use mod_monolis_pord_tree
  implicit none
  private

  public :: setupMinPriority, freeMinPriority
  public :: orderMinPriority
  public :: stageinfo_sum

contains

  !===========================================================================
  subroutine freeMinPriority(minprior)
    type(minprior_t), intent(inout) :: minprior

    if (associated(minprior%Gelim)) then
      call freeElimGraph(minprior%Gelim)
      deallocate(minprior%Gelim)
    end if
    if (associated(minprior%bucket)) then
      call freeBucket(minprior%bucket)
      deallocate(minprior%bucket)
    end if
    if (allocated(minprior%stageinfo)) deallocate(minprior%stageinfo)
    if (allocated(minprior%reachset))  deallocate(minprior%reachset)
    if (allocated(minprior%auxaux))    deallocate(minprior%auxaux)
    if (allocated(minprior%auxbin))    deallocate(minprior%auxbin)
    if (allocated(minprior%auxtmp))    deallocate(minprior%auxtmp)
    minprior%ms => null()
  end subroutine freeMinPriority

  !===========================================================================
  subroutine setupMinPriority(minprior, ms)
    type(minprior_t),    intent(out) :: minprior
    type(multisector_t), target, intent(inout) :: ms

    integer(ip) :: nvtx, nstages, u

    nvtx    = ms%G%nvtx
    nstages = ms%nstages

    minprior%ms => ms
    allocate(minprior%Gelim)
    call setupElimGraph(minprior%Gelim, ms%G)
    allocate(minprior%bucket)
    call setupBucket(minprior%bucket, nvtx, nvtx, 0)

    allocate(minprior%stageinfo(0:nstages-1))
    allocate(minprior%reachset(0:nvtx-1))
    allocate(minprior%auxaux(0:nvtx-1))
    allocate(minprior%auxbin(0:nvtx-1))
    allocate(minprior%auxtmp(0:nvtx-1))

    do u = 0, nvtx-1
      minprior%auxbin(u) = -1
      minprior%auxtmp(u) = 0
    end do

    minprior%stageinfo%nstep = 0
    minprior%stageinfo%welim = 0
    minprior%stageinfo%nzf   = 0
    minprior%stageinfo%ops   = 0.0_dp
    minprior%nreach = 0
    minprior%flag   = 1
  end subroutine setupMinPriority

  !===========================================================================
  subroutine orderMinPriority(T, minprior, options, cpus)
    type(elimtree_t), intent(out)   :: T
    type(minprior_t), intent(inout) :: minprior
    integer(ip),      intent(in)    :: options(0:*)
    real(dp),         intent(inout) :: cpus(0:*)

    integer(ip) :: nstages, istage, scoretype, ordtype

    nstages   = minprior%ms%nstages
    ordtype   = options(OPTION_ORDTYPE)
    scoretype = options(OPTION_NODE_SELECTION1)  ! first stage

    if ((nstages < 1) .or. (nstages > minprior%Gelim%G%nvtx)) then
      write(*,'(A,I0)') "Error in orderMinPriority: invalid nstages = ", nstages
      stop
    end if

    ! stage 0: eliminate domain vertices
    call eliminateStage(minprior, 0, scoretype, cpus)

    ! remaining stages
    scoretype = options(OPTION_NODE_SELECTION2)
    select case (ordtype)
      case (MINIMUM_PRIORITY)
        ! nothing more to do
      case (INCOMPLETE_ND)
        do istage = 1, nstages-1
          call eliminateStage(minprior, istage, scoretype, cpus)
        end do
      case (MULTISECTION)
        call eliminateStage(minprior, nstages-1, scoretype, cpus)
      case (TRISTAGE_MULTISECTION)
        do istage = 1, nstages-1
          call eliminateStage(minprior, istage, scoretype, cpus)
        end do
      case default
        write(*,'(A,I0)') "Error in orderMinPriority: unknown ordtype ", ordtype
        stop
    end select

    if ((ordtype /= MINIMUM_PRIORITY) .and. (options(OPTION_MSGLVL) > 1)) then
      do istage = 0, nstages-1
        write(*,'(I4,A,I6,A,I6,A,I8,A,ES12.4)') &
          istage, '. stage: #steps ', minprior%stageinfo(istage)%nstep, &
          ', weight ', minprior%stageinfo(istage)%welim, &
          ', nzl ', minprior%stageinfo(istage)%nzf, &
          ', ops ', minprior%stageinfo(istage)%ops
      end do
    end if

    call extractElimTree(T, minprior%Gelim)
  end subroutine orderMinPriority

  !===========================================================================
  subroutine eliminateStage(minprior, istage, scoretype, cpus)
    type(minprior_t), intent(inout) :: minprior
    integer(ip),      intent(in)    :: istage, scoretype
    real(dp),         intent(inout) :: cpus(0:*)

    integer(ip) :: nvtx, nreach, r, u, i
    integer(ip), pointer :: score(:), degree(:)

    nvtx   = minprior%Gelim%G%nvtx
    score  => minprior%Gelim%score
    degree => minprior%Gelim%degree

    ! collect all principal vertices in stage <= istage
    nreach = 0
    do u = 0, nvtx-1
      if ((score(u) == -1) .and. (minprior%ms%stage(u) <= istage)) then
        minprior%reachset(nreach) = u
        nreach = nreach + 1
        score(u) = degree(u)
      end if
    end do

    ! initial update
    call pord_starttimer(cpus(TIME_UPDSCORE))
    call updateDegree(minprior%Gelim, minprior%reachset, nreach, minprior%auxbin)
    call updateScore(minprior%Gelim, minprior%reachset, nreach, scoretype, minprior%auxbin)
    call pord_stoptimer(cpus(TIME_UPDSCORE))

    do i = 0, nreach-1
      u = minprior%reachset(i)
      call insertBucket(minprior%bucket, score(u), u)
    end do

    ! main elimination loop
    do
      if (eliminateStep(minprior, istage, scoretype) == 0) exit
      nreach = minprior%nreach

      call pord_starttimer(cpus(TIME_UPDADJNCY))
      call updateAdjncy(minprior%Gelim, minprior%reachset, nreach, &
                        minprior%auxtmp, minprior%flag)
      call pord_stoptimer(cpus(TIME_UPDADJNCY))

      call pord_starttimer(cpus(TIME_FINDINODES))
      call findIndNodes(minprior%Gelim, minprior%reachset, nreach, &
                        minprior%auxbin, minprior%auxaux, minprior%auxtmp, minprior%flag)
      call pord_stoptimer(cpus(TIME_FINDINODES))

      ! remove nonprincipal from reachset
      r = 0
      do i = 0, nreach-1
        u = minprior%reachset(i)
        if (score(u) >= 0) then
          minprior%reachset(r) = u; r = r + 1
        end if
      end do
      nreach = r

      call pord_starttimer(cpus(TIME_UPDSCORE))
      call updateDegree(minprior%Gelim, minprior%reachset, nreach, minprior%auxbin)
      call updateScore(minprior%Gelim, minprior%reachset, nreach, scoretype, minprior%auxbin)
      call pord_stoptimer(cpus(TIME_UPDSCORE))

      do i = 0, nreach-1
        u = minprior%reachset(i)
        call insertBucket(minprior%bucket, score(u), u)
      end do

      minprior%stageinfo(istage)%nstep = minprior%stageinfo(istage)%nstep + 1
    end do
  end subroutine eliminateStage

  !===========================================================================
  integer(ip) function eliminateStep(minprior, istage, scoretype)
    type(minprior_t), intent(inout) :: minprior
    integer(ip),      intent(in)    :: istage, scoretype

    integer(ip) :: nelim, minscr, vwghtu, u, v, i, istart, istop
    real(dp)    :: tri, rec
    integer(ip), pointer :: xadj(:), adjncy(:), vwght(:), len(:), degree(:), score(:)
    integer(ip), pointer :: stage(:)

    xadj   => minprior%Gelim%G%xadj
    adjncy => minprior%Gelim%G%adjncy
    vwght  => minprior%Gelim%G%vwght
    len    => minprior%Gelim%len
    degree => minprior%Gelim%degree
    score  => minprior%Gelim%score
    stage  => minprior%ms%stage

    u = minBucket(minprior%bucket)
    if (u == -1) then; eliminateStep = 0; return; end if
    minscr = score(u)

    nelim = 0
    minprior%nreach = 0

    main_loop: do
      vwghtu = vwght(u)
      call removeBucket(minprior%bucket, u)
      minprior%stageinfo(istage)%welim = minprior%stageinfo(istage)%welim + vwghtu
      nelim = nelim + 1

      call buildElement(minprior%Gelim, u)
      istart = xadj(u); istop = istart + len(u)
      do i = istart, istop-1
        v = adjncy(i)
        if (minprior%auxtmp(v) < minprior%flag) then
          minprior%auxtmp(v) = minprior%flag
          if (stage(v) <= istage) call removeBucket(minprior%bucket, v)
          minprior%reachset(minprior%nreach) = v
          minprior%nreach = minprior%nreach + 1
        end if
      end do

      ! update statistics
      tri = real(vwghtu, dp)
      rec = real(degree(u), dp)
      minprior%stageinfo(istage)%nzf = minprior%stageinfo(istage)%nzf + &
        int(tri * (tri + 1.0_dp) / 2.0_dp) + int(tri * rec)
      minprior%stageinfo(istage)%ops = minprior%stageinfo(istage)%ops + &
        (tri**3) / 3.0_dp + (tri**2) / 2.0_dp - 5.0_dp*tri / 6.0_dp + &
        tri**2 * rec + rec * (rec + 1.0_dp) * tri

      ! check for multiple elimination
      ! C uses integer division: (scoretype / 10 == 0) i.e. scoretype < 10
      if (scoretype / 10 == 0) exit main_loop  ! no multiple elimination
      u = minBucket(minprior%bucket)
      if (u == -1) exit main_loop
      if (score(u) > minscr) exit main_loop
    end do main_loop

    minprior%flag = minprior%flag + 1
    eliminateStep = nelim
  end function eliminateStep

  !===========================================================================
  ! stageinfo_sum: sum statistics over all stages.
  !===========================================================================
  subroutine stageinfo_sum(minprior, nstages, totnstep, totnzf, totops)
    type(minprior_t), intent(in)  :: minprior
    integer(ip),      intent(in)  :: nstages
    integer(ip),      intent(out) :: totnstep, totnzf
    real(dp),         intent(out) :: totops

    integer(ip) :: i
    totnstep = 0; totnzf = 0; totops = 0.0_dp
    do i = 0, nstages-1
      totnstep = totnstep + minprior%stageinfo(i)%nstep
      totnzf   = totnzf   + minprior%stageinfo(i)%nzf
      totops   = totops   + minprior%stageinfo(i)%ops
    end do
  end subroutine stageinfo_sum

end module mod_monolis_pord_minpriority
