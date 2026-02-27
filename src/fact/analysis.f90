module mod_monolis_fact_analysis
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc

  implicit none

contains

!==============================================================================
! Reverse Cuthill-McKee ordering
! O(N + nz*logN) — much better fill-reduction for FEM matrices
! than the simple greedy degree ordering.
!
! 1. Find pseudo-peripheral starting node (2-round BFS refinement)
! 2. BFS with degree-ordered neighbor insertion
! 3. Reverse the resulting ordering
!==============================================================================
subroutine reverse_cuthill_mckee_ordering(mat, lu)
  type(matrix_data), intent(inout) :: mat
  type(monolis_mat_lu), intent(inout) :: lu

  integer(kint) :: n, i, k, j, r, c
  integer(kint), allocatable :: degree(:), queue(:), nbr_buf(:)
  logical, allocatable :: visited(:)
  integer(kint) :: head, tail, node, nbr, start
  integer(kint) :: min_deg, nbr_count, tmp

  n = mat%n
  allocate(mat%perm(n), mat%invperm(n))
  allocate(degree(n), visited(n), queue(n), nbr_buf(n))

  ! Compute degrees (off-diagonal only)
  degree = 0
  do r = 1, n
    do k = mat%row_ptr(r)+1, mat%row_ptr(r+1)
      c = mat%col_ind(k)
      if (r /= c) degree(r) = degree(r) + 1
    end do
  end do

  ! Find pseudo-peripheral starting node
  ! Start with minimum degree node, then refine with 2 BFS rounds
  min_deg = n + 1
  start = 1
  do i = 1, n
    if (degree(i) < min_deg .and. degree(i) > 0) then
      min_deg = degree(i)
      start = i
    end if
  end do

  ! Refine: 2 rounds of BFS; last node of BFS is pseudo-peripheral
  do j = 1, 2
    visited = .false.
    head = 1; tail = 1
    queue(1) = start
    visited(start) = .true.

    do while (head <= tail)
      node = queue(head); head = head + 1
      do k = mat%row_ptr(node)+1, mat%row_ptr(node+1)
        nbr = mat%col_ind(k)
        if (nbr /= node .and. .not. visited(nbr)) then
          visited(nbr) = .true.
          tail = tail + 1
          queue(tail) = nbr
        end if
      end do
    end do

    start = queue(tail)   ! last node = pseudo-peripheral
  end do

  ! Main BFS with degree-sorted neighbor insertion (Cuthill-McKee)
  visited = .false.
  head = 1; tail = 1
  queue(1) = start
  visited(start) = .true.

  do while (head <= tail)
    node = queue(head); head = head + 1

    ! Collect unvisited neighbors
    nbr_count = 0
    do k = mat%row_ptr(node)+1, mat%row_ptr(node+1)
      nbr = mat%col_ind(k)
      if (nbr /= node .and. .not. visited(nbr)) then
        nbr_count = nbr_count + 1
        nbr_buf(nbr_count) = nbr
      end if
    end do

    ! Sort neighbors by non-decreasing degree (insertion sort — small lists)
    do i = 2, nbr_count
      tmp = nbr_buf(i)
      j = i - 1
      do while (j >= 1)
        if (degree(nbr_buf(j)) <= degree(tmp)) exit
        nbr_buf(j+1) = nbr_buf(j)
        j = j - 1
      end do
      nbr_buf(j+1) = tmp
    end do

    ! Add sorted neighbors to queue
    do i = 1, nbr_count
      nbr = nbr_buf(i)
      if (.not. visited(nbr)) then
        visited(nbr) = .true.
        tail = tail + 1
        queue(tail) = nbr
      end if
    end do
  end do

  ! Handle disconnected components
  do i = 1, n
    if (.not. visited(i)) then
      tail = tail + 1
      queue(tail) = i
    end if
  end do

  ! Reverse: RCM = reverse of CM ordering
  do i = 1, n
    mat%perm(n - i + 1) = queue(i)
    mat%invperm(queue(i)) = n - i + 1
  end do

  deallocate(degree, visited, queue, nbr_buf)
end subroutine

end module mod_monolis_fact_analysis
