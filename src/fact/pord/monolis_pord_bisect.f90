!*****************************************************************************
! SPACE (SPArse Cholesky Elimination) Library: mod_monolis_pord_gbipart.f90
!
! Fortran implementation of gbipart.c
!*****************************************************************************

module mod_monolis_pord_gbipart
  use mod_monolis_pord_const
  use mod_monolis_pord_types
  use mod_monolis_pord_graph
  implicit none
  private

  public :: setupBipartiteGraph, freeBipartiteGraph
  public :: maximumMatching, maximumFlow
  public :: DMviaMatching, DMviaFlow

contains

  !===========================================================================
  subroutine freeBipartiteGraph(Gbipart)
    type(gbipart_t), intent(inout) :: Gbipart
    call freeGraph(Gbipart%G)
    Gbipart%nX = 0; Gbipart%nY = 0
  end subroutine freeBipartiteGraph

  !===========================================================================
  subroutine setupBipartiteGraph(Gbipart, G, bipartvertex, nX, nY, vtxmap)
    type(gbipart_t), intent(out)   :: Gbipart
    type(graph_t),   intent(in)    :: G
    integer(ip),     intent(in)    :: nX, nY
    integer(ip),     intent(inout) :: bipartvertex(0:nX+nY-1)
    integer(ip),     intent(inout) :: vtxmap(0:G%nvtx-1)

    integer(ip) :: nedgesGb, totvwght, u, x, y, i, j, jstart, jstop, ptr

    nedgesGb = 0
    do i = 0, nX+nY-1
      u = bipartvertex(i)
      jstart = G%xadj(u); jstop = G%xadj(u+1)
      do j = jstart, jstop-1
        vtxmap(G%adjncy(j)) = -1
      end do
      nedgesGb = nedgesGb + (jstop - jstart)
    end do
    do i = 0, nX+nY-1
      vtxmap(bipartvertex(i)) = i
    end do

    call newGraph(Gbipart%G, nX+nY, nedgesGb)
    Gbipart%nX = nX; Gbipart%nY = nY

    totvwght = 0; ptr = 0
    do i = 0, nX-1
      x = bipartvertex(i)
      Gbipart%G%xadj(i)  = ptr
      Gbipart%G%vwght(i) = G%vwght(x)
      totvwght = totvwght + G%vwght(x)
      jstart = G%xadj(x); jstop = G%xadj(x+1)
      do j = jstart, jstop-1
        y = G%adjncy(j)
        if (vtxmap(y) >= nX) then
          Gbipart%G%adjncy(ptr) = vtxmap(y); ptr = ptr + 1
        end if
      end do
    end do
    do i = nX, nX+nY-1
      y = bipartvertex(i)
      Gbipart%G%xadj(i)  = ptr
      Gbipart%G%vwght(i) = G%vwght(y)
      totvwght = totvwght + G%vwght(y)
      jstart = G%xadj(y); jstop = G%xadj(y+1)
      do j = jstart, jstop-1
        x = G%adjncy(j)
        if ((vtxmap(x) >= 0) .and. (vtxmap(x) < nX)) then
          Gbipart%G%adjncy(ptr) = vtxmap(x); ptr = ptr + 1
        end if
      end do
    end do
    Gbipart%G%xadj(nX+nY) = ptr
    Gbipart%G%gtype    = G%gtype
    Gbipart%G%totvwght = totvwght
  end subroutine setupBipartiteGraph

  !===========================================================================
  ! maximumMatching: Hopcroft-Karp BFS + DFS
  ! matching(0:nX+nY-1): matching(x)=y means x matched to y, -1 = free
  !===========================================================================
  subroutine maximumMatching(Gbipart, matching)
    type(gbipart_t), intent(in)  :: Gbipart
    integer(ip),     intent(out) :: matching(0:Gbipart%nX+Gbipart%nY-1)

    integer(ip) :: nX, nY, nXY, x, x2, y, y2, i, istart, istop
    integer(ip) :: qhead, qtail, max_level, top, top2
    integer(ip), allocatable :: level(:), marker(:), queue(:), stack(:)

    nX = Gbipart%nX; nY = Gbipart%nY; nXY = nX + nY
    allocate(level(0:nXY-1), marker(0:nXY-1))
    allocate(queue(0:nX-1), stack(0:nY-1))
    matching = -1   ! -1 = FREE_NODE

    ! greedy initial matching
    do x = 0, nX-1
      istart = Gbipart%G%xadj(x); istop = Gbipart%G%xadj(x+1)
      do i = istart, istop-1
        y = Gbipart%G%adjncy(i)
        if (matching(y) == -1) then
          matching(x) = y; matching(y) = x; exit
        end if
      end do
    end do

    ! augment with Hopcroft-Karp
    outer: do
      level  = -1; marker = -1
      qhead  = 0;  qtail  = 0
      do x = 0, nX-1
        if (matching(x) == -1) then
          queue(qtail) = x; qtail = qtail + 1; level(x) = 0
        end if
      end do

      top = 0; max_level = MAX_INT
      do while (qhead /= qtail)
        x = queue(qhead); qhead = qhead + 1
        if (level(x) >= max_level) cycle
        istart = Gbipart%G%xadj(x); istop = Gbipart%G%xadj(x+1)
        do i = istart, istop-1
          y = Gbipart%G%adjncy(i)
          if (level(y) /= -1) cycle
          level(y) = level(x) + 1
          if (matching(y) == -1) then
            max_level = level(y)
            stack(top) = y; top = top + 1
          else if (level(y) < max_level) then
            x2 = matching(y)
            level(x2) = level(y) + 1
            queue(qtail) = x2; qtail = qtail + 1
          end if
        end do
      end do

      if (top == 0) exit outer

      do while (top > 0)
        top = top - 1
        y = stack(top)
        marker(y) = Gbipart%G%xadj(y)
        top2 = top + 1

        do while (top2 > top)
          y = stack(top2-1)
          i = marker(y)
          marker(y) = marker(y) + 1
          if (i < Gbipart%G%xadj(y+1)) then
            x = Gbipart%G%adjncy(i)
            if ((marker(x) == -1) .and. (level(x) == level(y)-1)) then
              marker(x) = 0
              if (level(x) == 0) then
                ! augment: pop entire stack segment and update matching
                do while (top2 > top)
                  y2 = stack(top2-1); top2 = top2 - 1
                  x2 = matching(y2)
                  matching(x) = y2; matching(y2) = x
                  x = x2
                end do
              else
                y2 = matching(x)
                stack(top2) = y2; top2 = top2 + 1
                marker(y2) = Gbipart%G%xadj(y2)
              end if
            end if
          else
            top2 = top2 - 1
          end if
        end do
      end do
    end do outer

    deallocate(level, marker, queue, stack)
  end subroutine maximumMatching

  !===========================================================================
  ! maximumFlow: Edmonds-Karp on vertex-capacitated bipartite graph.
  ! flow(0:nedges-1): positive = forward x->y, negative = backward y->x
  ! rc(0:nXY-1): residual capacity
  !===========================================================================
  subroutine maximumFlow(Gbipart, flow, rc)
    type(gbipart_t), intent(in)  :: Gbipart
    integer(ip),     intent(out) :: flow(0:Gbipart%G%nedges-1)
    integer(ip),     intent(out) :: rc(0:Gbipart%nX+Gbipart%nY-1)

    integer(ip) :: nedges, nX, nY, nXY
    integer(ip) :: u, v, x, y, j, jj, i, istart, istop
    integer(ip) :: qhead, qtail, capacity
    integer(ip), allocatable :: parent(:), marker(:), queue(:)

    nedges = Gbipart%G%nedges
    nX = Gbipart%nX; nY = Gbipart%nY; nXY = nX + nY
    allocate(parent(0:nXY-1), marker(0:nXY-1), queue(0:nXY-1))

    ! initialise residual capacities and zero flow
    do u = 0, nXY-1
      rc(u) = Gbipart%G%vwght(u)
    end do
    flow = 0

    ! greedy initial flow
    do x = 0, nX-1
      if (rc(x) == 0) cycle
      istart = Gbipart%G%xadj(x); istop = Gbipart%G%xadj(x+1)
      do i = istart, istop-1
        y = Gbipart%G%adjncy(i)
        capacity = min(rc(x), rc(y))
        if (capacity > 0) then
          rc(x) = rc(x) - capacity; rc(y) = rc(y) - capacity
          flow(i) = capacity
          ! find reverse edge y->x
          jj = Gbipart%G%xadj(y)
          do while (Gbipart%G%adjncy(jj) /= x)
            jj = jj + 1
          end do
          flow(jj) = -capacity
        end if
        if (rc(x) == 0) exit
      end do
    end do

    ! Edmonds-Karp augmentation
    bfs_outer: do
      parent = -1; marker = -1
      qhead = 0; qtail = 0
      do x = 0, nX-1
        if (rc(x) > 0) then
          queue(qtail) = x; qtail = qtail + 1; parent(x) = x
        end if
      end do

      capacity = 0
      bfs_inner: do while (qhead /= qtail)
        u = queue(qhead); qhead = qhead + 1
        istart = Gbipart%G%xadj(u); istop = Gbipart%G%xadj(u+1)
        do i = istart, istop-1
          v = Gbipart%G%adjncy(i)
          ! traverse: X->Y forward (unlimited cap) or Y->X via backward edge
          if ((parent(v) /= -1) .or. &
              (.not. ((v >= nX) .or. (flow(i) < 0)))) cycle
          parent(v) = u; marker(v) = i
          queue(qtail) = v; qtail = qtail + 1

          if ((v >= nX) .and. (rc(v) > 0)) then
            ! augmenting path found: trace back to get min capacity
            u = v
            capacity = rc(u)
            do while (parent(u) /= u)
              j = marker(u); u = parent(u)
              if (u >= nX) capacity = min(capacity, -flow(j))
            end do
            capacity = min(capacity, rc(u))

            ! augment flow along the path
            rc(v) = rc(v) - capacity
            do while (parent(v) /= v)
              j = marker(v); u = parent(v)
              flow(j) = flow(j) + capacity
              ! find reverse edge u->v
              jj = Gbipart%G%xadj(v)
              do while (Gbipart%G%adjncy(jj) /= u)
                jj = jj + 1
              end do
              flow(jj) = -flow(j)
              v = u
            end do
            rc(v) = rc(v) - capacity
            qhead = qtail  ! break BFS
            exit bfs_inner
          end if
        end do
      end do bfs_inner

      if (capacity == 0) exit bfs_outer
    end do bfs_outer

    deallocate(parent, marker, queue)
  end subroutine maximumFlow

  !===========================================================================
  ! DMviaMatching
  !===========================================================================
  subroutine DMviaMatching(Gbipart, matching, dmflag, dmwght)
    type(gbipart_t), intent(in)  :: Gbipart
    integer(ip),     intent(in)  :: matching(0:Gbipart%nX+Gbipart%nY-1)
    integer(ip),     intent(out) :: dmflag(0:Gbipart%nX+Gbipart%nY-1)
    integer(ip),     intent(out) :: dmwght(0:5)

    integer(ip) :: nX, nY, nXY, u, x, y, i, istart, istop
    integer(ip) :: qhead, qtail
    integer(ip), allocatable :: queue(:)

    nX = Gbipart%nX; nY = Gbipart%nY; nXY = nX + nY
    allocate(queue(0:nXY-1))

    qhead = 0; qtail = 0
    do x = 0, nX-1
      if (matching(x) == -1) then
        queue(qtail) = x; qtail = qtail + 1; dmflag(x) = SI
      else
        dmflag(x) = SR
      end if
    end do
    do y = nX, nXY-1
      if (matching(y) == -1) then
        queue(qtail) = y; qtail = qtail + 1; dmflag(y) = BI
      else
        dmflag(y) = BR
      end if
    end do

    do while (qhead /= qtail)
      u = queue(qhead); qhead = qhead + 1
      istart = Gbipart%G%xadj(u); istop = Gbipart%G%xadj(u+1)
      select case (dmflag(u))
        case (SI)
          do i = istart, istop-1
            y = Gbipart%G%adjncy(i)
            if (dmflag(y) == BR) then
              dmflag(y) = BX; queue(qtail) = y; qtail = qtail + 1
            end if
          end do
        case (BX)
          x = matching(u)
          if (dmflag(x) /= SI) then
            dmflag(x) = SI; queue(qtail) = x; qtail = qtail + 1
          end if
        case (BI)
          do i = istart, istop-1
            x = Gbipart%G%adjncy(i)
            if (dmflag(x) == SR) then
              dmflag(x) = SX; queue(qtail) = x; qtail = qtail + 1
            end if
          end do
        case (SX)
          y = matching(u)
          if (dmflag(y) /= BI) then
            dmflag(y) = BI; queue(qtail) = y; qtail = qtail + 1
          end if
      end select
    end do

    dmwght = 0
    do u = 0, nXY-1
      dmwght(dmflag(u)) = dmwght(dmflag(u)) + Gbipart%G%vwght(u)
    end do
    deallocate(queue)
  end subroutine DMviaMatching

  !===========================================================================
  ! DMviaFlow
  ! Direct Fortran translation of gbipart.c's DMviaFlow.
  ! Internally uses C's SOURCE/SINK/FREE sentinels (-2/-3/-1) and only at the
  ! end converts to the Dulmage-Mendelsohn labels SI/SX/SR/BI/BX/BR.
  !===========================================================================
  subroutine DMviaFlow(Gbipart, flow, rc, dmflag, dmwght)
    type(gbipart_t), intent(in)  :: Gbipart
    integer(ip),     intent(in)  :: flow(0:Gbipart%G%nedges-1)
    integer(ip),     intent(in)  :: rc(0:Gbipart%nX+Gbipart%nY-1)
    integer(ip),     intent(out) :: dmflag(0:Gbipart%nX+Gbipart%nY-1)
    integer(ip),     intent(out) :: dmwght(0:5)

    integer(ip), parameter :: TAG_FREE   = -1
    integer(ip), parameter :: TAG_SOURCE = -2
    integer(ip), parameter :: TAG_SINK   = -3

    integer(ip) :: nX, nY, nXY, u, v, x, y, i, istart, istop
    integer(ip) :: qhead, qtail
    integer(ip), allocatable :: queue(:)

    nX = Gbipart%nX; nY = Gbipart%nY; nXY = nX + nY
    allocate(queue(0:nXY-1))

    ! mark all X/Y nodes reachable from source/sink
    qhead = 0; qtail = 0
    do x = 0, nX-1
      if (rc(x) > 0) then
        queue(qtail) = x; qtail = qtail + 1
        dmflag(x) = TAG_SOURCE
      else
        dmflag(x) = TAG_FREE
      end if
    end do
    do y = nX, nXY-1
      if (rc(y) > 0) then
        queue(qtail) = y; qtail = qtail + 1
        dmflag(y) = TAG_SINK
      else
        dmflag(y) = TAG_FREE
      end if
    end do

    do while (qhead /= qtail)
      u = queue(qhead); qhead = qhead + 1
      istart = Gbipart%G%xadj(u); istop = Gbipart%G%xadj(u+1)
      select case (dmflag(u))
        case (TAG_SOURCE)
          ! u is X reachable from source; follow forward edges u->v (v in Y)
          ! or backward edges u<-v (v in X) along negative flow
          do i = istart, istop-1
            v = Gbipart%G%adjncy(i)
            if ((dmflag(v) == TAG_FREE) .and. &
                ((v >= nX) .or. (flow(i) < 0))) then
              queue(qtail) = v; qtail = qtail + 1
              dmflag(v) = TAG_SOURCE
            end if
          end do
        case (TAG_SINK)
          ! u is Y reachable from sink; follow forward edges v->u (v in X)
          ! or backward edges v<-u (v in Y) along positive flow
          do i = istart, istop-1
            v = Gbipart%G%adjncy(i)
            if ((dmflag(v) == TAG_FREE) .and. &
                ((v < nX) .or. (flow(i) > 0))) then
              queue(qtail) = v; qtail = qtail + 1
              dmflag(v) = TAG_SINK
            end if
          end do
      end select
    end do

    ! convert tags to D-M labels and accumulate weights
    dmwght = 0
    do x = 0, nX-1
      select case (dmflag(x))
        case (TAG_SOURCE); dmflag(x) = SI; dmwght(SI) = dmwght(SI) + Gbipart%G%vwght(x)
        case (TAG_SINK);   dmflag(x) = SX; dmwght(SX) = dmwght(SX) + Gbipart%G%vwght(x)
        case default;      dmflag(x) = SR; dmwght(SR) = dmwght(SR) + Gbipart%G%vwght(x)
      end select
    end do
    do y = nX, nXY-1
      select case (dmflag(y))
        case (TAG_SOURCE); dmflag(y) = BX; dmwght(BX) = dmwght(BX) + Gbipart%G%vwght(y)
        case (TAG_SINK);   dmflag(y) = BI; dmwght(BI) = dmwght(BI) + Gbipart%G%vwght(y)
        case default;      dmflag(y) = BR; dmwght(BR) = dmwght(BR) + Gbipart%G%vwght(y)
      end select
    end do

    deallocate(queue)
  end subroutine DMviaFlow

end module mod_monolis_pord_gbipart
!*****************************************************************************
! SPACE (SPArse Cholesky Elimination) Library: mod_monolis_pord_ddbisect.f90
!
! Fortran implementation of ddbisect.c
! (initialDDSep, improveDDSep)
!*****************************************************************************

module mod_monolis_pord_ddbisect
  use mod_monolis_pord_const
  use mod_monolis_pord_types
  use mod_monolis_pord_bucket
  implicit none
  private

  public :: initialDDSep, improveDDSep

contains

  !===========================================================================
  ! findPseudoPeripheralDomain
  !===========================================================================
  integer(ip) function findPseudoPeripheralDomain(dd, domain_in)
    type(domdec_t), intent(in)  :: dd
    integer(ip),    intent(in)  :: domain_in

    integer(ip) :: nvtx, qhead, qtail, nlev, lastdomain, u, v, i, istart, istop
    integer(ip) :: domain
    integer(ip), allocatable :: level(:), queue(:)

    nvtx = dd%G%nvtx
    allocate(level(0:nvtx-1), queue(0:nvtx-1))

    domain = domain_in; nlev = 0; lastdomain = domain

    outer: do
      level = -1
      queue(0) = domain; level(domain) = 0
      qhead = 0; qtail = 1
      do while (qhead /= qtail)
        u = queue(qhead); qhead = qhead + 1
        if (dd%vtype(u) == 1) lastdomain = u
        istart = dd%G%xadj(u); istop = dd%G%xadj(u+1)
        do i = istart, istop-1
          v = dd%G%adjncy(i)
          if (level(v) == -1) then
            queue(qtail) = v; qtail = qtail + 1
            level(v) = level(u) + 1
          end if
        end do
      end do
      if (level(lastdomain) > nlev) then
        nlev = level(lastdomain); domain = lastdomain
      else
        exit outer
      end if
    end do outer

    findPseudoPeripheralDomain = domain
    deallocate(level, queue)
  end function findPseudoPeripheralDomain

  !===========================================================================
  ! constructLevelSep
  !===========================================================================
  subroutine constructLevelSep(dd, domain)
    type(domdec_t), intent(inout) :: dd
    integer(ip),    intent(in)    :: domain

    integer(ip) :: nvtx, bestvalue, weight, dS, dB, dW
    integer(ip) :: qhead, qtail, qopt, q, u, v, w
    integer(ip) :: i, istart, istop, j, jstart, jstop
    integer(ip), allocatable :: queue(:), deltaS(:), deltaB(:), deltaW(:)

    nvtx = dd%G%nvtx
    allocate(queue(0:nvtx-1), deltaS(0:nvtx-1), deltaB(0:nvtx-1), deltaW(0:nvtx-1))

    deltaS = 0; deltaB = 0; deltaW = 0
    do u = 0, nvtx-1
      if (dd%vtype(u) == 2) deltaW(u) = dd%G%xadj(u+1) - dd%G%xadj(u)
    end do

    qhead = 0; qtail = 1; queue(0) = domain
    dd%vtype(domain) = -1

    do while ((dd%cwght(BLACK) < dd%cwght(WHITE)) .and. (qhead /= qtail))
      qopt = qhead; bestvalue = MAX_INT

      do q = qhead, qtail-1
        u = queue(q)
        if (dd%vtype(u) == -1) then
          dB = dd%G%vwght(u); dW = -dB; dS = 0
          istart = dd%G%xadj(u); istop = dd%G%xadj(u+1)
          do i = istart, istop-1
            v = dd%G%adjncy(i)
            weight = dd%G%vwght(v)
            if (dd%color(v) == WHITE) then
              dW = dW - weight; dS = dS + weight
            else if (deltaW(v) == 1) then
              dB = dB + weight; dS = dS - weight
            end if
          end do
          deltaS(u) = dS; deltaB(u) = dB; deltaW(u) = dW
          dd%vtype(u) = -2
        end if
        if (dd%cwght(GRAY) + deltaS(u) < bestvalue) then
          bestvalue = dd%cwght(GRAY) + deltaS(u); qopt = q
        end if
      end do

      u = queue(qopt)
      ! swap to head
      queue(qopt) = queue(qhead); queue(qhead) = u
      qhead = qhead + 1

      dd%color(u) = BLACK
      dd%cwght(GRAY)  = dd%cwght(GRAY)  + deltaS(u)
      dd%cwght(BLACK) = dd%cwght(BLACK) + deltaB(u)
      dd%cwght(WHITE) = dd%cwght(WHITE) + deltaW(u)
      dd%vtype(u) = -3

      istart = dd%G%xadj(u); istop = dd%G%xadj(u+1)
      do i = istart, istop-1
        v = dd%G%adjncy(i)
        deltaB(v) = deltaB(v) + 1
        deltaW(v) = deltaW(v) - 1
        if (deltaW(v) == 0) then
          dd%color(v) = BLACK
        else if (deltaB(v) == 1) then
          dd%color(v) = GRAY
          jstart = dd%G%xadj(v); jstop = dd%G%xadj(v+1)
          do j = jstart, jstop-1
            w = dd%G%adjncy(j)
            if (dd%vtype(w) == 1) then
              queue(qtail) = w; qtail = qtail + 1; dd%vtype(w) = -1
            else if (dd%vtype(w) == -2) then
              dd%vtype(w) = -1
            end if
          end do
        else if (deltaW(v) == 1) then
          jstart = dd%G%xadj(v); jstop = dd%G%xadj(v+1)
          do j = jstart, jstop-1
            w = dd%G%adjncy(j)
            if (dd%vtype(w) == -2) dd%vtype(w) = -1
          end do
        end if
      end do
    end do

    do i = 0, qtail-1
      u = queue(i); dd%vtype(u) = 1
    end do
    deallocate(queue, deltaS, deltaB, deltaW)
  end subroutine constructLevelSep

  !===========================================================================
  ! initialDDSep
  !===========================================================================
  subroutine initialDDSep(dd)
    type(domdec_t), intent(inout) :: dd

    integer(ip) :: nvtx, totvwght, domain, u

    nvtx     = dd%G%nvtx
    totvwght = dd%G%totvwght

    dd%cwght(GRAY)  = 0
    dd%cwght(BLACK) = 0
    dd%cwght(WHITE) = totvwght
    do u = 0, nvtx-1
      dd%color(u) = WHITE
    end do

    do u = 0, nvtx-1
      if ((dd%vtype(u) == 1) .and. (dd%color(u) == WHITE)) then
        domain = findPseudoPeripheralDomain(dd, u)
        call constructLevelSep(dd, domain)
        if (dd%cwght(BLACK) >= dd%cwght(WHITE)) exit
      end if
    end do
  end subroutine initialDDSep

  !===========================================================================
  ! updateB2W / updateW2B helpers (inline into improveDDSep)
  !===========================================================================
  subroutine updateB2W(w_bucket, b_bucket, dd, domain, tmp_color, deltaW, deltaB, deltaS)
    type(bucket_t), intent(inout) :: w_bucket, b_bucket
    type(domdec_t), intent(inout) :: dd
    integer(ip),    intent(in)    :: domain
    integer(ip),    intent(inout) :: tmp_color(0:dd%G%nvtx-1)
    integer(ip),    intent(inout) :: deltaW(0:dd%G%nvtx-1)
    integer(ip),    intent(inout) :: deltaB(0:dd%G%nvtx-1)
    integer(ip),    intent(inout) :: deltaS(0:dd%G%nvtx-1)

    integer(ip) :: weight, u, v, i, istart, istop, j, jstart, jstop

    istart = dd%G%xadj(domain); istop = dd%G%xadj(domain+1)
    do i = istart, istop-1
      u = dd%G%adjncy(i); weight = dd%G%vwght(u)
      jstart = dd%G%xadj(u); jstop = dd%G%xadj(u+1)

      if (deltaW(u) < 0) then
        v = -(deltaW(u)+1); deltaW(u) = 1
        call removeBucket(w_bucket, v)
        deltaB(v) = deltaB(v) - weight; deltaS(v) = deltaS(v) + weight
        call insertBucket(w_bucket, deltaS(v), v)
      end if
      if (deltaW(u) == 0) then
        tmp_color(u) = GRAY
        do j = jstart, jstop-1
          v = dd%G%adjncy(j)
          if (dd%vtype(v) == 1) then
            call removeBucket(b_bucket, v)
            deltaB(v) = deltaB(v) + weight; deltaS(v) = deltaS(v) - weight
            call insertBucket(b_bucket, deltaS(v), v)
          end if
        end do
      end if
      if (deltaB(u) < 0) deltaB(u) = 1
      deltaB(u) = deltaB(u) - 1; deltaW(u) = deltaW(u) + 1
      if (deltaB(u) == 1) then
        do j = jstart, jstop-1
          v = dd%G%adjncy(j)
          if ((tmp_color(v) == BLACK) .and. (dd%vtype(v) == 1)) then
            call removeBucket(b_bucket, v)
            deltaW(v) = deltaW(v) + weight; deltaS(v) = deltaS(v) - weight
            deltaB(u) = -(v+1)
            call insertBucket(b_bucket, deltaS(v), v)
          end if
        end do
      end if
      if (deltaB(u) == 0) then
        tmp_color(u) = WHITE
        do j = jstart, jstop-1
          v = dd%G%adjncy(j)
          if (dd%vtype(v) == 1) then
            call removeBucket(w_bucket, v)
            deltaW(v) = deltaW(v) - weight; deltaS(v) = deltaS(v) + weight
            call insertBucket(w_bucket, deltaS(v), v)
          end if
        end do
      end if
    end do
  end subroutine updateB2W

  subroutine updateW2B(w_bucket, b_bucket, dd, domain, tmp_color, deltaW, deltaB, deltaS)
    type(bucket_t), intent(inout) :: w_bucket, b_bucket
    type(domdec_t), intent(inout) :: dd
    integer(ip),    intent(in)    :: domain
    integer(ip),    intent(inout) :: tmp_color(0:dd%G%nvtx-1)
    integer(ip),    intent(inout) :: deltaW(0:dd%G%nvtx-1)
    integer(ip),    intent(inout) :: deltaB(0:dd%G%nvtx-1)
    integer(ip),    intent(inout) :: deltaS(0:dd%G%nvtx-1)

    integer(ip) :: weight, u, v, i, istart, istop, j, jstart, jstop

    istart = dd%G%xadj(domain); istop = dd%G%xadj(domain+1)
    do i = istart, istop-1
      u = dd%G%adjncy(i); weight = dd%G%vwght(u)
      jstart = dd%G%xadj(u); jstop = dd%G%xadj(u+1)

      if (deltaB(u) < 0) then
        v = -(deltaB(u)+1); deltaB(u) = 1
        call removeBucket(b_bucket, v)
        deltaW(v) = deltaW(v) - weight; deltaS(v) = deltaS(v) + weight
        call insertBucket(b_bucket, deltaS(v), v)
      end if
      if (deltaB(u) == 0) then
        tmp_color(u) = GRAY
        do j = jstart, jstop-1
          v = dd%G%adjncy(j)
          if (dd%vtype(v) == 1) then
            call removeBucket(w_bucket, v)
            deltaW(v) = deltaW(v) + weight; deltaS(v) = deltaS(v) - weight
            call insertBucket(w_bucket, deltaS(v), v)
          end if
        end do
      end if
      if (deltaW(u) < 0) deltaW(u) = 1
      deltaB(u) = deltaB(u) + 1; deltaW(u) = deltaW(u) - 1
      if (deltaW(u) == 1) then
        do j = jstart, jstop-1
          v = dd%G%adjncy(j)
          if ((tmp_color(v) == WHITE) .and. (dd%vtype(v) == 1)) then
            call removeBucket(w_bucket, v)
            deltaB(v) = deltaB(v) + weight; deltaS(v) = deltaS(v) - weight
            deltaW(u) = -(v+1)
            call insertBucket(w_bucket, deltaS(v), v)
          end if
        end do
      end if
      if (deltaW(u) == 0) then
        tmp_color(u) = BLACK
        do j = jstart, jstop-1
          v = dd%G%adjncy(j)
          if (dd%vtype(v) == 1) then
            call removeBucket(b_bucket, v)
            deltaB(v) = deltaB(v) - weight; deltaS(v) = deltaS(v) + weight
            call insertBucket(b_bucket, deltaS(v), v)
          end if
        end do
      end if
    end do
  end subroutine updateW2B

  !===========================================================================
  ! improveDDSep: FM variant for domain decomposition separator improvement.
  !===========================================================================
  subroutine improveDDSep(dd)
    type(domdec_t), intent(inout) :: dd

    type(bucket_t) :: b_bucket, w_bucket
    integer(ip) :: nvtx, weight, tmp_S, tmp_B, tmp_W
    integer(ip) :: pos, bestglobalpos, badflips, b_domain, w_domain, domain, nxtdomain
    integer(ip) :: fhead, ftail, u, v, i, istart, istop
    real(dp)    :: bestglobalvalue, b_value, w_value, value
    integer(ip), allocatable :: tmp_color(:), deltaS(:), deltaB(:), deltaW(:)

    nvtx = dd%G%nvtx
    allocate(tmp_color(0:nvtx-1), deltaS(0:nvtx-1), deltaB(0:nvtx-1), deltaW(0:nvtx-1))

    outer_loop: do

      tmp_S = dd%cwght(GRAY); tmp_B = dd%cwght(BLACK); tmp_W = dd%cwght(WHITE)
      bestglobalpos = 0; badflips = 0
      bestglobalvalue = eval_sep(tmp_S, tmp_B, tmp_W)

      call setupBucket(b_bucket, nvtx, nvtx, nvtx/2)
      call setupBucket(w_bucket, nvtx, nvtx, nvtx/2)

      fhead = 0; ftail = -1; pos = 0

      ! initialise multisec colors and deltas
      do u = 0, nvtx-1
        if (dd%vtype(u) == 2) then
          deltaB(u) = 0; deltaW(u) = 0
          istart = dd%G%xadj(u); istop = dd%G%xadj(u+1)
          do i = istart, istop-1
            v = dd%G%adjncy(i)
            if (dd%color(v) == BLACK) then
              deltaB(u) = deltaB(u) + 1
            else
              deltaW(u) = deltaW(u) + 1
            end if
          end do
          if ((deltaB(u) > 0) .and. (deltaW(u) > 0)) then
            tmp_color(u) = GRAY
          else if (deltaB(u) > 0) then
            tmp_color(u) = BLACK
          else
            tmp_color(u) = WHITE
          end if
          dd%color(u) = tmp_color(u)
        end if
      end do

      ! initialise domain colors and fill buckets
      do u = 0, nvtx-1
        if (dd%vtype(u) == 1) then
          tmp_color(u) = dd%color(u)
          if (tmp_color(u) == BLACK) then
            deltaW(u) = dd%G%vwght(u); deltaB(u) = -deltaW(u); deltaS(u) = 0
            istart = dd%G%xadj(u); istop = dd%G%xadj(u+1)
            do i = istart, istop-1
              v = dd%G%adjncy(i); weight = dd%G%vwght(v)
              if (tmp_color(v) == BLACK) then
                deltaB(u) = deltaB(u) - weight; deltaS(u) = deltaS(u) + weight
              else if (deltaB(v) == 1) then
                deltaW(u) = deltaW(u) + weight; deltaS(u) = deltaS(u) - weight
                deltaB(v) = -(u+1)
              end if
            end do
            call insertBucket(b_bucket, deltaS(u), u)
          end if
          if (tmp_color(u) == WHITE) then
            deltaB(u) = dd%G%vwght(u); deltaW(u) = -deltaB(u); deltaS(u) = 0
            istart = dd%G%xadj(u); istop = dd%G%xadj(u+1)
            do i = istart, istop-1
              v = dd%G%adjncy(i); weight = dd%G%vwght(v)
              if (tmp_color(v) == WHITE) then
                deltaW(u) = deltaW(u) - weight; deltaS(u) = deltaS(u) + weight
              else if (deltaW(v) == 1) then
                deltaB(u) = deltaB(u) + weight; deltaS(u) = deltaS(u) - weight
                deltaW(v) = -(u+1)
              end if
            end do
            call insertBucket(w_bucket, deltaS(u), u)
          end if
        end if
      end do

      inner_loop: do

        b_value = MAX_FLOAT; w_value = MAX_FLOAT
        b_domain = minBucket(b_bucket)
        if (b_domain /= -1) &
          b_value = eval_sep(tmp_S+deltaS(b_domain), tmp_B+deltaB(b_domain), tmp_W+deltaW(b_domain))
        w_domain = minBucket(w_bucket)
        if (w_domain /= -1) &
          w_value = eval_sep(tmp_S+deltaS(w_domain), tmp_B+deltaB(w_domain), tmp_W+deltaW(w_domain))

        if ((b_domain == -1) .and. (w_domain == -1)) exit inner_loop

        if (b_value + EPS < w_value) then
          domain = b_domain; value = b_value
          call removeBucket(b_bucket, domain)
        else
          domain = w_domain; value = w_value
          call removeBucket(w_bucket, domain)
        end if

        if (ftail /= -1) then
          dd%vtype(ftail) = -(domain+1)
        else
          fhead = -(domain+1)
        end if
        dd%vtype(domain) = 0; ftail = domain

        if (tmp_color(domain) == BLACK) then
          tmp_color(domain) = WHITE
          call updateB2W(w_bucket, b_bucket, dd, domain, tmp_color, deltaW, deltaB, deltaS)
        else
          tmp_color(domain) = BLACK
          call updateW2B(w_bucket, b_bucket, dd, domain, tmp_color, deltaW, deltaB, deltaS)
        end if
        tmp_S = tmp_S + deltaS(domain)
        tmp_B = tmp_B + deltaB(domain)
        tmp_W = tmp_W + deltaW(domain)

        pos = pos + 1
        if (value + EPS < bestglobalvalue) then
          bestglobalvalue = value; bestglobalpos = pos; badflips = 0
        else
          badflips = badflips + 1
        end if
        if (badflips >= MAX_BAD_FLIPS) exit inner_loop

      end do inner_loop

      ! apply best flips
      pos = 0; nxtdomain = fhead
      do while (nxtdomain /= 0)
        domain = -nxtdomain - 1
        if (pos < bestglobalpos) then
          if (dd%color(domain) == BLACK) then
            dd%color(domain) = WHITE
          else
            dd%color(domain) = BLACK
          end if
          dd%cwght(GRAY)  = dd%cwght(GRAY)  + deltaS(domain)
          dd%cwght(BLACK) = dd%cwght(BLACK) + deltaB(domain)
          dd%cwght(WHITE) = dd%cwght(WHITE) + deltaW(domain)
          pos = pos + 1
        end if
        nxtdomain = dd%vtype(domain)
        dd%vtype(domain) = 1
      end do

      call freeBucket(b_bucket)
      call freeBucket(w_bucket)

      if (bestglobalpos > 0) cycle outer_loop
      exit outer_loop

    end do outer_loop

    deallocate(tmp_color, deltaS, deltaB, deltaW)
  end subroutine improveDDSep

end module mod_monolis_pord_ddbisect
!*****************************************************************************
! SPACE (SPArse Cholesky Elimination) Library: mod_monolis_pord_ddcreate.f90
!
! Fortran implementation of ddcreate.c
!*****************************************************************************

module mod_monolis_pord_ddcreate
  use mod_monolis_pord_const
  use mod_monolis_pord_types
  use mod_monolis_pord_graph
  use mod_monolis_pord_sort
  implicit none
  private

  public :: newDomainDecomposition, freeDomainDecomposition
  public :: constructDomainDecomposition
  public :: shrinkDomainDecomposition

contains

  !===========================================================================
  subroutine newDomainDecomposition(dd, nvtx, nedges)
    type(domdec_t), pointer, intent(inout) :: dd
    integer(ip),    intent(in)  :: nvtx, nedges

    call newGraph(dd%G, nvtx, nedges)
    dd%ndom = 0; dd%domwght = 0
    dd%cwght = 0
    allocate(dd%vtype(0:nvtx-1))
    allocate(dd%color(0:nvtx-1))
    allocate(dd%map_(0:nvtx-1))
    dd%vtype = 0; dd%color = -1; dd%map_ = -1
    dd%prev => null(); dd%next_ptr => null()
  end subroutine newDomainDecomposition

  !===========================================================================
  subroutine freeDomainDecomposition(dd)
    type(domdec_t), pointer, intent(inout) :: dd
    if (.not. associated(dd)) return
    call freeGraph(dd%G)
    if (allocated(dd%vtype)) deallocate(dd%vtype)
    if (allocated(dd%color)) deallocate(dd%color)
    if (allocated(dd%map_))  deallocate(dd%map_)
    deallocate(dd)
  end subroutine freeDomainDecomposition

  !===========================================================================
  subroutine buildInitialDomains(G, vtxlist, vtype, rep)
    type(graph_t), intent(in)    :: G
    integer(ip),   intent(in)    :: vtxlist(0:G%nvtx-1)
    integer(ip),   intent(inout) :: vtype(0:G%nvtx-1)
    integer(ip),   intent(inout) :: rep(0:G%nvtx-1)

    integer(ip) :: nvtx, u, v, w, i, j, jstart, jstop

    nvtx = G%nvtx

    do i = 0, nvtx-1
      u = vtxlist(i)
      if (vtype(u) == 0) then
        vtype(u) = 1
        jstart = G%xadj(u); jstop = G%xadj(u+1)
        do j = jstart, jstop-1
          vtype(G%adjncy(j)) = 2
        end do
      end if
    end do

    ! absorb multisecs adjacent to only one domain
    do i = 0, nvtx-1
      u = vtxlist(i)
      if (vtype(u) == 2) then
        v = -1
        jstart = G%xadj(u); jstop = G%xadj(u+1)
        do j = jstart, jstop-1
          w = G%adjncy(j)
          if (vtype(w) == 1) then
            if (v == -1) then
              v = rep(w)
            else if (v /= rep(w)) then
              v = -1; exit
            end if
          end if
        end do
        if (v /= -1) then
          vtype(u) = 1; rep(u) = v
        end if
      end if
    end do
  end subroutine buildInitialDomains

  !===========================================================================
  subroutine mergeMultisecs(G, vtype, rep)
    type(graph_t), intent(in)    :: G
    integer(ip),   intent(inout) :: vtype(0:G%nvtx-1)
    integer(ip),   intent(inout) :: rep(0:G%nvtx-1)

    integer(ip) :: nvtx, u, v, w, x, flag, keepon
    integer(ip) :: qhead, qtail, i, istart, istop, j, jstart, jstop
    integer(ip), allocatable :: tmp(:), queue(:)

    nvtx = G%nvtx
    allocate(tmp(0:nvtx-1), queue(0:nvtx-1))
    tmp = -1; flag = 1

    do u = 0, nvtx-1
      if (vtype(u) == 2) then
        qhead = 0; qtail = 1; queue(0) = u; vtype(u) = -2

        ! mark adjacent domains
        istart = G%xadj(u); istop = G%xadj(u+1)
        do i = istart, istop-1
          v = G%adjncy(i)
          if (vtype(v) == 1) tmp(rep(v)) = flag
        end do

        do while (qhead /= qtail)
          v = queue(qhead); qhead = qhead + 1
          istart = G%xadj(v); istop = G%xadj(v+1)
          do i = istart, istop-1
            w = G%adjncy(i)
            if (vtype(w) == 2) then
              keepon = 1
              jstart = G%xadj(w); jstop = G%xadj(w+1)
              do j = jstart, jstop-1
                x = G%adjncy(j)
                if ((vtype(x) == 1) .and. (tmp(rep(x)) == flag)) then
                  keepon = 0; exit
                end if
              end do
              if (keepon == 1) then
                do j = jstart, jstop-1
                  x = G%adjncy(j)
                  if (vtype(x) == 1) tmp(rep(x)) = flag
                end do
                queue(qtail) = w; qtail = qtail + 1
                rep(w) = u; vtype(w) = -2
              end if
            end if
          end do
        end do
        flag = flag + 1
      end if
    end do

    do u = 0, nvtx-1
      if (vtype(u) == -2) vtype(u) = 2
    end do
    deallocate(tmp, queue)
  end subroutine mergeMultisecs

  !===========================================================================
  subroutine initialDomainDecomposition(dd_out, G, map_, vtype, rep)
    type(domdec_t), pointer, intent(inout) :: dd_out
    type(graph_t),  intent(in)  :: G
    integer(ip),    intent(out) :: map_(0:G%nvtx-1)
    integer(ip),    intent(in)  :: vtype(0:G%nvtx-1)
    integer(ip),    intent(in)  :: rep(0:G%nvtx-1)

    integer(ip) :: nvtx, nedges, u, v, w
    integer(ip) :: nvtxdd, nedgesdd, ndom, domwght, flag
    integer(ip) :: i, j, jstart, jstop
    integer(ip), allocatable :: tmp(:), bin(:)

    nvtx   = G%nvtx
    nedges = G%nedges
    allocate(tmp(0:nvtx-1), bin(0:nvtx-1))
    tmp = -1; bin = -1

    call newDomainDecomposition(dd_out, nvtx, nedges)

    ! put all nodes u belonging to representative v in bin[v]
    do u = 0, nvtx-1
      v = rep(u)
      if (u /= v) then
        bin(u) = bin(v); bin(v) = u
      end if
    end do

    flag = 1; nedgesdd = 0; nvtxdd = 0; ndom = 0; domwght = 0
    do u = 0, nvtx-1
      if (rep(u) == u) then
        dd_out%G%xadj(nvtxdd)   = nedgesdd
        dd_out%vtype(nvtxdd)    = vtype(u)
        dd_out%G%vwght(nvtxdd)  = 0
        tmp(u) = flag

        v = u
        do
          map_(v) = nvtxdd
          dd_out%G%vwght(nvtxdd) = dd_out%G%vwght(nvtxdd) + G%vwght(v)
          jstart = G%xadj(v); jstop = G%xadj(v+1)
          do j = jstart, jstop-1
            w = G%adjncy(j)
            if ((vtype(w) /= vtype(u)) .and. (tmp(rep(w)) /= flag)) then
              tmp(rep(w)) = flag
              dd_out%G%adjncy(nedgesdd) = rep(w)
              nedgesdd = nedgesdd + 1
            end if
          end do
          v = bin(v)
          if (v == -1) exit
        end do

        if (dd_out%vtype(nvtxdd) == 1) then
          ndom    = ndom + 1
          domwght = domwght + dd_out%G%vwght(nvtxdd)
        end if
        nvtxdd = nvtxdd + 1
        flag = flag + 1
      end if
    end do

    dd_out%G%xadj(nvtxdd) = nedgesdd
    dd_out%G%nvtx   = nvtxdd
    dd_out%G%nedges = nedgesdd
    dd_out%G%gtype  = WEIGHTED
    dd_out%G%totvwght = G%totvwght
    do i = 0, nedgesdd-1
      dd_out%G%adjncy(i) = map_(dd_out%G%adjncy(i))
    end do
    dd_out%ndom    = ndom
    dd_out%domwght = domwght
    do i = 0, nvtxdd-1
      dd_out%color(i) = -1; dd_out%map_(i) = -1
    end do

    deallocate(tmp, bin)
  end subroutine initialDomainDecomposition

  !===========================================================================
  subroutine constructDomainDecomposition(dd_out, G, map_)
    type(domdec_t), pointer, intent(inout) :: dd_out
    type(graph_t),  intent(in)  :: G
    integer(ip),    intent(out) :: map_(0:G%nvtx-1)

    integer(ip) :: nvtx, deg, u, i, istart, istop
    integer(ip), allocatable :: vtxlist(:), key(:), vtype(:), rep(:)

    nvtx = G%nvtx
    allocate(vtxlist(0:nvtx-1), key(0:nvtx-1))
    allocate(vtype(0:nvtx-1), rep(0:nvtx-1))

    do u = 0, nvtx-1
      vtxlist(u) = u
      istart = G%xadj(u); istop = G%xadj(u+1)
      select case (G%gtype)
        case (UNWEIGHTED)
          deg = istop - istart
        case (WEIGHTED)
          deg = 0
          do i = istart, istop-1
            deg = deg + G%vwght(G%adjncy(i))
          end do
        case default
          deg = istop - istart
      end select
      key(u) = deg
    end do
    call distributionCounting(nvtx, vtxlist, key)
    deallocate(key)

    vtype = 0
    do u = 0, nvtx-1; rep(u) = u; end do

    call buildInitialDomains(G, vtxlist, vtype, rep)
    call mergeMultisecs(G, vtype, rep)
    deallocate(vtxlist)

    allocate(dd_out)
    call initialDomainDecomposition(dd_out, G, map_, vtype, rep)
    deallocate(vtype, rep)
  end subroutine constructDomainDecomposition

  !===========================================================================
  subroutine computePriorities(dd, msvtxlist, key, scoretype)
    type(domdec_t), intent(inout) :: dd
    integer(ip),    intent(in)    :: msvtxlist(0:*)
    integer(ip),    intent(out)   :: key(0:dd%G%nvtx-1)
    integer(ip),    intent(in)    :: scoretype

    integer(ip) :: nvtx, nlist, k, weight, deg, u, v, w
    integer(ip) :: i, istart, istop, j, jstart, jstop
    integer(ip), allocatable :: marker(:)

    nvtx  = dd%G%nvtx
    nlist = nvtx - dd%ndom
    allocate(marker(0:nvtx-1))

    select case (scoretype)
      case (QMRDV)
        do k = 0, nlist-1
          u = msvtxlist(k)
          weight = dd%G%vwght(u)
          istart = dd%G%xadj(u); istop = dd%G%xadj(u+1)
          do i = istart, istop-1
            weight = weight + dd%G%vwght(dd%G%adjncy(i))
          end do
          key(u) = weight / max(1, dd%G%vwght(u))
        end do

      case (QMD)
        marker = -1
        do k = 0, nlist-1; marker(msvtxlist(k)) = msvtxlist(k); end do
        do k = 0, nlist-1
          u = msvtxlist(k)
          marker(u) = u; deg = 0
          istart = dd%G%xadj(u); istop = dd%G%xadj(u+1)
          do i = istart, istop-1
            v = dd%G%adjncy(i)
            jstart = dd%G%xadj(v); jstop = dd%G%xadj(v+1)
            do j = jstart, jstop-1
              w = dd%G%adjncy(j)
              if (marker(w) /= u) then
                marker(w) = u; deg = deg + dd%G%vwght(w)
              end if
            end do
          end do
          key(u) = deg
        end do

      case (QRAND)
        do k = 0, nlist-1
          u = msvtxlist(k)
          call random_number_int(key(u), nvtx)
        end do

      case default
        key = 0
    end select
    deallocate(marker)
  end subroutine computePriorities

  subroutine random_number_int(k, range)
    integer(ip), intent(out) :: k
    integer(ip), intent(in)  :: range
    real :: r
    call random_number(r)
    k = int(r * real(range))
    if (k >= range) k = range - 1
    if (k < 0) k = 0
  end subroutine random_number_int

  !===========================================================================
  subroutine eliminateMultisecs(dd, msvtxlist, rep)
    type(domdec_t), intent(inout) :: dd
    integer(ip),    intent(in)    :: msvtxlist(0:*)
    integer(ip),    intent(inout) :: rep(0:dd%G%nvtx-1)

    integer(ip) :: nvtx, nlist, keepon, u, v, w, k, i, istart, istop

    nvtx  = dd%G%nvtx
    nlist = nvtx - dd%ndom

    do k = 0, nlist-1
      u = msvtxlist(k)
      istart = dd%G%xadj(u); istop = dd%G%xadj(u+1)
      keepon = 1
      do i = istart, istop-1
        v = dd%G%adjncy(i)
        if (rep(v) /= v) then; keepon = 0; exit; end if
      end do
      if (keepon == 1) then
        dd%vtype(u) = 3
        do i = istart, istop-1
          v = dd%G%adjncy(i); rep(v) = u
        end do
      end if
    end do

    ! absorb multisecs adjacent to only one (super)domain
    do k = 0, nlist-1
      u = msvtxlist(k)
      if (dd%vtype(u) == 2) then
        v = -1
        istart = dd%G%xadj(u); istop = dd%G%xadj(u+1)
        do i = istart, istop-1
          w = dd%G%adjncy(i)
          if (v == -1) then
            v = rep(w)
          else if (v /= rep(w)) then
            v = -1; exit
          end if
        end do
        if (v /= -1) then
          dd%vtype(u) = 4; rep(u) = v
        end if
      end if
    end do
  end subroutine eliminateMultisecs

  !===========================================================================
  subroutine findIndMultisecs(dd, msvtxlist, rep)
    type(domdec_t), intent(inout) :: dd
    integer(ip),    intent(in)    :: msvtxlist(0:*)
    integer(ip),    intent(inout) :: rep(0:dd%G%nvtx-1)

    integer(ip) :: nvtx, nlist, flag, keepon, deg, chk, u, v, k
    integer(ip) :: ulast, i, istart, istop
    integer(ip), allocatable :: tmp(:), bin(:), nextn(:), key(:)
    integer(ip), allocatable :: checksum(:)

    nvtx  = dd%G%nvtx
    nlist = nvtx - dd%ndom
    allocate(tmp(0:nvtx-1), bin(0:nvtx-1), nextn(0:nvtx-1), key(0:nvtx-1))
    allocate(checksum(0:nvtx-1))
    tmp = -1; bin = -1; flag = 1

    do k = 0, nlist-1
      u = msvtxlist(k)
      if (dd%vtype(u) == 2) then
        deg = 0; chk = 0
        istart = dd%G%xadj(u); istop = dd%G%xadj(u+1)
        do i = istart, istop-1
          v = dd%G%adjncy(i)
          if (tmp(rep(v)) /= flag) then
            tmp(rep(v)) = flag; chk = chk + rep(v); deg = deg + 1
          end if
        end do
        chk = mod(chk, nvtx)
        checksum(u) = chk; key(u) = deg
        nextn(u) = bin(chk); bin(chk) = u
        flag = flag + 1
      end if
    end do

    do k = 0, nlist-1
      u = msvtxlist(k)
      if (dd%vtype(u) == 2) then
        chk = checksum(u)
        v = bin(chk); bin(chk) = -1
        do while (v /= -1)
          istart = dd%G%xadj(v); istop = dd%G%xadj(v+1)
          do i = istart, istop-1
            tmp(rep(dd%G%adjncy(i))) = flag
          end do
          ulast = v; u = nextn(v)
          do while (u /= -1)
            keepon = 1
            if (key(u) /= key(v)) keepon = 0
            if (keepon == 1) then
              do i = dd%G%xadj(u), dd%G%xadj(u+1)-1
                if (tmp(rep(dd%G%adjncy(i))) /= flag) then
                  keepon = 0; exit
                end if
              end do
            end if
            if (keepon == 1) then
              rep(u) = v; dd%vtype(u) = 4
              u = nextn(u); nextn(ulast) = u
            else
              ulast = u; u = nextn(u)
            end if
          end do
          v = nextn(v); flag = flag + 1
        end do
      end if
    end do

    deallocate(tmp, bin, nextn, key, checksum)
  end subroutine findIndMultisecs

  !===========================================================================
  subroutine coarserDomainDecomposition(dd2_out, dd1, rep)
    type(domdec_t), pointer, intent(out)   :: dd2_out
    type(domdec_t), pointer, intent(inout) :: dd1
    integer(ip),    intent(in)             :: rep(0:dd1%G%nvtx-1)

    integer(ip) :: nvtxdd1, nedgesdd1, nvtxdd2, nedgesdd2
    integer(ip) :: ndom, domwght, flag, u, v, w, i, istart, istop
    integer(ip), allocatable :: tmp(:), bin(:)

    nvtxdd1   = dd1%G%nvtx
    nedgesdd1 = dd1%G%nedges

    allocate(tmp(0:nvtxdd1-1), bin(0:nvtxdd1-1))
    tmp = -1; bin = -1

    ! Allocate with upper bounds (as in C), shrink at end.
    allocate(dd2_out)
    call newDomainDecomposition(dd2_out, nvtxdd1, nedgesdd1)

    ! build bin chains: non-representatives link to their representative
    do u = 0, nvtxdd1-1
      v = rep(u)
      if (u /= v) then
        bin(u) = bin(v); bin(v) = u
      end if
    end do

    ! set up vertex mapping dd1 -> dd2
    nvtxdd2 = 0; nedgesdd2 = 0; ndom = 0; domwght = 0
    flag = 1
    do u = 0, nvtxdd1-1
      if (rep(u) == u) then
        if ((dd1%vtype(u) == 1) .or. (dd1%vtype(u) == 2) .or. (dd1%vtype(u) == 3)) then
          dd1%map_(u) = nvtxdd2
          dd2_out%G%xadj(nvtxdd2) = nedgesdd2
          dd2_out%G%vwght(nvtxdd2) = 0
          tmp(u) = flag

          ! collect all cluster members via bin chain
          v = u
          do
            dd1%map_(v) = nvtxdd2
            dd2_out%G%vwght(nvtxdd2) = dd2_out%G%vwght(nvtxdd2) + dd1%G%vwght(v)
            if ((dd1%vtype(v) == 1) .or. (dd1%vtype(v) == 2)) then
              istart = dd1%G%xadj(v); istop = dd1%G%xadj(v+1)
              do i = istart, istop-1
                w = rep(dd1%G%adjncy(i))
                if (tmp(w) /= flag) then
                  tmp(w) = flag
                  dd2_out%G%adjncy(nedgesdd2) = w
                  nedgesdd2 = nedgesdd2 + 1
                end if
              end do
            end if
            v = bin(v)
            if (v == -1) exit
          end do

          dd2_out%vtype(nvtxdd2) = dd1%vtype(u)
          if (dd2_out%vtype(nvtxdd2) == 3) dd2_out%vtype(nvtxdd2) = 1
          if (dd2_out%vtype(nvtxdd2) == 1) then
            ndom = ndom + 1
            domwght = domwght + dd2_out%G%vwght(nvtxdd2)
          end if
          nvtxdd2 = nvtxdd2 + 1; flag = flag + 1
        end if
      end if
    end do

    ! fix adjacency: rep-id -> local index
    do i = 0, nedgesdd2-1
      dd2_out%G%adjncy(i) = dd1%map_(dd2_out%G%adjncy(i))
    end do
    dd2_out%G%xadj(nvtxdd2) = nedgesdd2
    dd2_out%G%nvtx    = nvtxdd2
    dd2_out%G%nedges  = nedgesdd2
    dd2_out%G%gtype   = WEIGHTED
    dd2_out%G%totvwght = dd1%G%totvwght
    dd2_out%ndom      = ndom
    dd2_out%domwght   = domwght
    ! reset vtype==3/4 in dd1 back to 2 (as in C)
    do u = 0, nvtxdd1-1
      if ((dd1%vtype(u) == 3) .or. (dd1%vtype(u) == 4)) dd1%vtype(u) = 2
    end do
    do u = 0, nvtxdd2-1
      dd2_out%color(u) = -1; dd2_out%map_(u) = -1
    end do

    deallocate(tmp, bin)
  end subroutine coarserDomainDecomposition

  !===========================================================================
  subroutine shrinkDomainDecomposition(dd, scoretype)
    type(domdec_t), pointer, intent(inout) :: dd
    integer(ip),             intent(in)    :: scoretype

    integer(ip) :: nvtx, nlist, k, u
    integer(ip), allocatable :: msvtxlist(:), key(:), rep(:)
    type(domdec_t), pointer  :: dd2

    nvtx  = dd%G%nvtx
    ! count multisec nodes (vtype==2) explicitly
    nlist = 0
    do u = 0, nvtx-1
      if (dd%vtype(u) == 2) nlist = nlist + 1
    end do
    allocate(msvtxlist(0:nlist-1), key(0:nvtx-1), rep(0:nvtx-1))

    k = 0
    do u = 0, nvtx-1
      if (dd%vtype(u) == 2) then
        msvtxlist(k) = u; k = k + 1
      end if
    end do

    call computePriorities(dd, msvtxlist, key, scoretype)

    ! sort msvtxlist by priority key
    call sort_by_key(msvtxlist, key, nlist)

    do u = 0, nvtx-1; rep(u) = u; end do
    call eliminateMultisecs(dd, msvtxlist, rep)
    call findIndMultisecs(dd, msvtxlist, rep)
    call coarserDomainDecomposition(dd2, dd, rep)
    dd2%prev  => dd
    dd%next_ptr => dd2

    deallocate(msvtxlist, key, rep)
  end subroutine shrinkDomainDecomposition

  ! Simple insertion sort of msvtxlist(0:n-1) ascending by key(msvtxlist(i))
  subroutine sort_by_key(arr, key, n)
    integer(ip), intent(inout) :: arr(0:n-1)
    integer(ip), intent(in)    :: key(0:*)
    integer(ip), intent(in)    :: n
    integer(ip) :: i, j, e, ke

    do i = 1, n-1
      e = arr(i); ke = key(e); j = i
      do while (j > 0)
        if (key(arr(j-1)) <= ke) exit
        arr(j) = arr(j-1); j = j - 1
      end do
      arr(j) = e
    end do
  end subroutine sort_by_key

end module mod_monolis_pord_ddcreate
!*****************************************************************************
! SPACE (SPArse Cholesky Elimination) Library: mod_monolis_pord_gbisect.f90
!
! Fortran implementation of gbisect.c
! (constructSeparator, smoothSeparator, smoothBy2Layers)
!*****************************************************************************

module mod_monolis_pord_gbisect
  use mod_monolis_pord_const
  use mod_monolis_pord_types
  use mod_monolis_pord_graph
  use mod_monolis_pord_gbipart
  use mod_monolis_pord_ddcreate
  use mod_monolis_pord_ddbisect
  implicit none
  private

  public :: newGbisect, freeGbisect
  public :: constructSeparator, smoothSeparator

contains

  !===========================================================================
  subroutine newGbisect(Gbisect, G)
    type(gbisect_t), intent(out) :: Gbisect
    type(graph_t),   target, intent(inout) :: G

    Gbisect%G => G
    Gbisect%cwght = 0
    allocate(Gbisect%color(0:G%nvtx-1))
    Gbisect%color = WHITE
  end subroutine newGbisect

  !===========================================================================
  subroutine freeGbisect(Gbisect)
    type(gbisect_t), intent(inout) :: Gbisect
    if (allocated(Gbisect%color)) deallocate(Gbisect%color)
    Gbisect%G => null()
  end subroutine freeGbisect

  !===========================================================================
  ! constructSeparator: multilevel separator construction.
  !===========================================================================
  subroutine constructSeparator(Gbisect, options, cpus)
    type(gbisect_t), intent(inout) :: Gbisect
    integer(ip),     intent(in)    :: options(0:*)
    real(dp),        intent(inout) :: cpus(0:*)

    type(domdec_t), pointer :: dd_init
    type(domdec_t), pointer :: dd, dd2
    integer(ip), allocatable :: map_(:)
    integer(ip) :: nvtx, nstep, u

    nvtx = Gbisect%G%nvtx
    allocate(map_(0:nvtx-1))

    call pord_starttimer(cpus(TIME_INITDOMDEC))
    allocate(dd_init)
    call constructDomainDecomposition(dd_init, Gbisect%G, map_)
    dd => dd_init
    call pord_stoptimer(cpus(TIME_INITDOMDEC))

    call pord_starttimer(cpus(TIME_COARSEDOMDEC))
    nstep = 0
    do while ((dd%ndom > MIN_DOMAINS) .and. (nstep < MAX_COARSENING_STEPS) .and. &
              ((dd%G%nedges/2) > dd%G%nvtx))
      call shrinkDomainDecomposition(dd, options(OPTION_NODE_SELECTION3))
      dd => dd%next_ptr
      nstep = nstep + 1
    end do
    call pord_stoptimer(cpus(TIME_COARSEDOMDEC))

    call pord_starttimer(cpus(TIME_INITSEP))
    call initialDDSep(dd)
    if (dd%cwght(GRAY) > 0) call improveDDSep(dd)
    call pord_stoptimer(cpus(TIME_INITSEP))

    call pord_starttimer(cpus(TIME_REFINESEP))
    do while (associated(dd%prev))
      dd2 => dd%prev
      dd2%cwght(GRAY)  = dd%cwght(GRAY)
      dd2%cwght(BLACK) = dd%cwght(BLACK)
      dd2%cwght(WHITE) = dd%cwght(WHITE)
      do u = 0, dd2%G%nvtx-1
        dd2%color(u) = dd%color(dd2%map_(u))
      end do
      call freeDomainDecomposition(dd)
      if (dd2%cwght(GRAY) > 0) call improveDDSep(dd2)
      dd => dd2
    end do
    call pord_stoptimer(cpus(TIME_REFINESEP))

    ! copy coloring back to Gbisect
    Gbisect%cwght(GRAY)  = dd%cwght(GRAY)
    Gbisect%cwght(BLACK) = dd%cwght(BLACK)
    Gbisect%cwght(WHITE) = dd%cwght(WHITE)
    do u = 0, nvtx-1
      Gbisect%color(u) = dd%color(map_(u))
    end do

    call freeDomainDecomposition(dd)
    deallocate(map_)
  end subroutine constructSeparator

  !===========================================================================
  ! smoothBy2Layers: smooth separator by pairing with black or white side.
  ! Returns .true. if smoothing improved the separator.
  !===========================================================================
  logical function smoothBy2Layers(Gbisect, bipartvertex, pnX, black, white)
    type(gbisect_t), intent(inout) :: Gbisect
    integer(ip),     intent(inout) :: bipartvertex(0:*)
    integer(ip),     intent(inout) :: pnX
    integer(ip),     intent(in)    :: black, white

    type(gbipart_t) :: Gbipart
    integer(ip), allocatable :: map_(:), dmflag(:), matching(:), flow(:), rc(:)
    integer(ip) :: dmwght(0:5)
    integer(ip) :: nvtx, nX, nX2, nY, x, y, u, i, j, jstart, jstop

    nvtx = Gbisect%G%nvtx
    nX   = pnX
    allocate(map_(0:nvtx-1))

    ! build set Y
    nY = 0
    do i = 0, nX-1
      x = bipartvertex(i)
      jstart = Gbisect%G%xadj(x); jstop = Gbisect%G%xadj(x+1)
      do j = jstart, jstop-1
        y = Gbisect%G%adjncy(j)
        if (Gbisect%color(y) == black) then
          bipartvertex(nX+nY) = y; nY = nY + 1
          Gbisect%color(y) = GRAY
        end if
      end do
    end do
    do i = nX, nX+nY-1
      y = bipartvertex(i); Gbisect%color(y) = black
    end do

    call setupBipartiteGraph(Gbipart, Gbisect%G, bipartvertex, nX, nY, map_)

    allocate(dmflag(0:nX+nY-1))
    smoothBy2Layers = .false.

    select case (Gbipart%G%gtype)
      case (UNWEIGHTED)
        allocate(matching(0:nX+nY-1))
        call maximumMatching(Gbipart, matching)
        call DMviaMatching(Gbipart, matching, dmflag, dmwght)
        deallocate(matching)
      case (WEIGHTED)
        allocate(flow(0:Gbipart%G%nedges-1), rc(0:nX+nY-1))
        call maximumFlow(Gbipart, flow, rc)
        call DMviaFlow(Gbipart, flow, rc, dmflag, dmwght)
        deallocate(flow, rc)
    end select

    ! test 1: exchange SI with BX
    if (eval_sep(Gbisect%cwght(GRAY)-dmwght(SI)+dmwght(BX), &
                 Gbisect%cwght(black)-dmwght(BX), &
                 Gbisect%cwght(white)+dmwght(SI)) + EPS &
        < eval_sep(Gbisect%cwght(GRAY), Gbisect%cwght(black), Gbisect%cwght(white))) then
      smoothBy2Layers = .true.
      Gbisect%cwght(white) = Gbisect%cwght(white) + dmwght(SI)
      Gbisect%cwght(GRAY)  = Gbisect%cwght(GRAY)  - dmwght(SI)
      Gbisect%cwght(black) = Gbisect%cwght(black) - dmwght(BX)
      Gbisect%cwght(GRAY)  = Gbisect%cwght(GRAY)  + dmwght(BX)
      do i = 0, nX+nY-1
        u = bipartvertex(i)
        if (dmflag(map_(u)) == SI) Gbisect%color(u) = white
        if (dmflag(map_(u)) == BX) Gbisect%color(u) = GRAY
      end do
    end if

    ! test 2: exchange SR with BR (allowed if smoothed or SI==0)
    if ((eval_sep(Gbisect%cwght(GRAY)-dmwght(SR)+dmwght(BR), &
                  Gbisect%cwght(black)-dmwght(BR), &
                  Gbisect%cwght(white)+dmwght(SR)) + EPS &
        < eval_sep(Gbisect%cwght(GRAY), Gbisect%cwght(black), Gbisect%cwght(white))) &
        .and. (smoothBy2Layers .or. (dmwght(SI) == 0))) then
      smoothBy2Layers = .true.
      Gbisect%cwght(white) = Gbisect%cwght(white) + dmwght(SR)
      Gbisect%cwght(GRAY)  = Gbisect%cwght(GRAY)  - dmwght(SR)
      Gbisect%cwght(black) = Gbisect%cwght(black) - dmwght(BR)
      Gbisect%cwght(GRAY)  = Gbisect%cwght(GRAY)  + dmwght(BR)
      do i = 0, nX+nY-1
        u = bipartvertex(i)
        if (dmflag(map_(u)) == SR) Gbisect%color(u) = white
        if (dmflag(map_(u)) == BR) Gbisect%color(u) = GRAY
      end do
    end if

    nX2 = 0
    do i = 0, nX+nY-1
      u = bipartvertex(i)
      if (Gbisect%color(u) == GRAY) then
        bipartvertex(nX2) = u; nX2 = nX2 + 1
      end if
    end do
    pnX = nX2

    deallocate(map_, dmflag)
    call freeBipartiteGraph(Gbipart)
  end function smoothBy2Layers

  !===========================================================================
  ! smoothSeparator
  !===========================================================================
  subroutine smoothSeparator(Gbisect, options)
    type(gbisect_t), intent(inout) :: Gbisect
    integer(ip),     intent(in)    :: options(0:*)

    integer(ip) :: dummy_opt
    integer(ip) :: nvtx, nX, nX2, x, y, a, b, i, j, jstart, jstop
    integer(ip), allocatable :: bipartvertex(:)
    logical :: kept_on, tmp_a

    dummy_opt = options(0)  ! suppress unused-argument warning
    nvtx = Gbisect%G%nvtx
    allocate(bipartvertex(0:nvtx-1))

    nX = 0
    do x = 0, nvtx-1
      if (Gbisect%color(x) == GRAY) then
        bipartvertex(nX) = x; nX = nX + 1
      end if
    end do

    kept_on = .true.
    do while (kept_on)
      ! minimise separator
      Gbisect%cwght(GRAY) = 0; nX2 = 0
      do i = 0, nX-1
        x = bipartvertex(i); a = 0; b = 0
        jstart = Gbisect%G%xadj(x); jstop = Gbisect%G%xadj(x+1)
        do j = jstart, jstop-1
          y = Gbisect%G%adjncy(j)
          if (Gbisect%color(y) == WHITE) a = 1
          if (Gbisect%color(y) == BLACK) b = 1
        end do
        if ((a == 1) .and. (b == 0)) then
          Gbisect%color(x) = WHITE
          Gbisect%cwght(WHITE) = Gbisect%cwght(WHITE) + Gbisect%G%vwght(x)
        else if ((a == 0) .and. (b == 1)) then
          Gbisect%color(x) = BLACK
          Gbisect%cwght(BLACK) = Gbisect%cwght(BLACK) + Gbisect%G%vwght(x)
        else
          bipartvertex(nX2) = x; nX2 = nX2 + 1
          Gbisect%cwght(GRAY) = Gbisect%cwght(GRAY) + Gbisect%G%vwght(x)
        end if
      end do
      nX = nX2

      ! smooth via bipartite matching/flow
      if (Gbisect%cwght(BLACK) >= Gbisect%cwght(WHITE)) then
        tmp_a = smoothBy2Layers(Gbisect, bipartvertex, nX, BLACK, WHITE)
        if (.not. tmp_a) tmp_a = smoothBy2Layers(Gbisect, bipartvertex, nX, WHITE, BLACK)
      else
        tmp_a = smoothBy2Layers(Gbisect, bipartvertex, nX, WHITE, BLACK)
        if (.not. tmp_a) tmp_a = smoothBy2Layers(Gbisect, bipartvertex, nX, BLACK, WHITE)
      end if
      kept_on = tmp_a
    end do

    deallocate(bipartvertex)
  end subroutine smoothSeparator

end module mod_monolis_pord_gbisect
