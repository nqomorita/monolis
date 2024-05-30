module mod_monolis_fact_analysis
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc

  implicit none

contains

  subroutine monolis_matrix_get_super_node_information(monoTREE, n_super_node, super_node)
    implicit none
    type(monolis_mat) :: monoTREE
    integer(kint), allocatable :: super_node(:)
    integer(kint) :: n_super_node
    integer(kint) :: N, i, j, in, jn, iS, d1, d2

    N = monoTREE%N
    allocate(super_node(N), source = 0)

    jn = 1
    super_node(N) = 1
    do i = N-1, 1, -1
      j = monoTREE%SCSR%indexU(i) + 2
      in = monoTREE%SCSR%itemU(j)

      d1 = monoTREE%SCSR%indexU(i  + 1) - monoTREE%SCSR%indexU(i)
      d2 = monoTREE%SCSR%indexU(in + 1) - monoTREE%SCSR%indexU(in)

      !> check super node
      if(d1 == d2 + 1)then
        super_node(i) = jn
      else
        jn = jn + 1
        super_node(i) = jn
      endif
      !write(*,*)d1, d2
    enddo
    n_super_node = jn
    !write(*,*)"super_node", super_node
  end subroutine monolis_matrix_get_super_node_information

  subroutine monolis_matrix_get_factorize_order(monoTREE, fact_order)
    implicit none
    type(monolis_mat) :: monoTREE
    integer(kint), allocatable :: fact_level(:)
    integer(kint), allocatable :: fact_order(:)
    integer(kint) :: N, i, j, in, jn, iS

    N = monoTREE%N

    allocate(fact_level(N), source = 0)
    allocate(fact_order(N), source = 0)

    fact_level(N) = 1

    do i = N-1, 1, -1
      j = monoTREE%SCSR%indexU(i) + 2
      in = monoTREE%SCSR%itemU(j)
      fact_level(i) = fact_level(in) + 1
    enddo

    do i = 1, N
      fact_level(i) = N + 1 - fact_level(i)
    enddo

    do i = 1, N
      fact_order(i) = i
    enddo

!write(*,*)"fact_order", fact_order
!write(*,*)"fact_level", fact_level
    call monolis_qsort_I_2d(fact_level, fact_order, 1, N)
!write(*,*)"fact_order", fact_order
!write(*,*)"fact_level", fact_level
    !> 最適化必要

    iS = 1
    do i = 1, N
      in = fact_level(iS)
      jn = 0
      aa:do j = iS + 1, N
        if(in == fact_level(j))then
          jn = jn + 1
        else
          exit aa
        endif
      enddo aa
      call monolis_qsort_I_1d(fact_order, iS, iS + jn)
      iS = iS + jn + 1
      if(iS > N) exit
    enddo
!write(*,*)"fact_order", fact_order
!write(*,*)"fact_level", fact_level
  end subroutine monolis_matrix_get_factorize_order

  subroutine monolis_matrix_get_factorize_array(monoTREE, fact_order, n_fact_array, fact_array, fact_array_index)
    implicit none
    type(monolis_mat) :: monoTREE
    integer(kint) :: fact_order(:)
    integer(kint) :: n_fact_array
    real(kdouble), allocatable :: fact_array(:)
    integer(kint), allocatable :: fact_array_index(:)
    integer(kint) :: i, j, N, in, jn

    n_fact_array = 0

    N = monoTREE%N

    allocate(fact_array_index(N + 1), source = 0)

    do j = 1, N
      i = j
      !i = fact_order(j)
      in = monoTREE%SCSR%indexU(i + 1) - monoTREE%SCSR%indexU(i)
      jn = in*(in + 1)/2
      n_fact_array = n_fact_array + jn
      fact_array_index(j + 1) = jn
    enddo

    allocate(fact_array(n_fact_array), source = 0.0d0)

    do j = 1, N
      fact_array_index(j + 1) = fact_array_index(j + 1) + fact_array_index(j)
    enddo

!write(*,*)"n_fact_array", n_fact_array
!write(*,*)"fact_array_index", fact_array_index
  end subroutine monolis_matrix_get_factorize_array

  subroutine monolis_matrix_set_value_of_factorize_array(monoMAT, monoTREE, fact_order, n_fact_array, fact_array, &
    & fact_array_index)
    implicit none
    type(monolis_mat) :: monoMAT
    type(monolis_mat) :: monoTREE
    integer(kint) :: fact_order(:)
    integer(kint) :: n_fact_array
    real(kdouble) :: fact_array(:)
    integer(kint) :: fact_array_index(:)
    integer(kint) :: N
    integer(kint) :: i, j, k, l, iSorg, iS, iE, jS, jE, in, jn, ln

    N = monoTREE%N

    do k = 1, N
      iS = monoTREE%SCSR%indexU(k) + 1
      iE = monoTREE%SCSR%indexU(k + 1)
      iSorg = iS
      jS = monoMAT%CSR%index(k) + 1
      jE = monoMAT%CSR%index(k + 1)
      aa:do j = jS, jE
        jn = monoMAT%CSR%item(j)
        if(jn < k) cycle aa
        do i = iS, iE
          in = monoTREE%SCSR%itemU(i)
          if(jn == in)then
            ln = fact_array_index(k)
            !lS = NDOF2*(i-1)
            !lE = NDOF2*(j-1)
            !do l = 1, NDOF2
            fact_array(ln + i - iSorg + 1) = monoMAT%R%A(j)
            !enddo
            iS = i + 1
            cycle aa
          endif
        enddo
      enddo aa
    enddo
  end subroutine monolis_matrix_set_value_of_factorize_array

  subroutine monolis_matrix_get_add_location(monoTREE, fact_order, fact_array_index, add_location)
    implicit none
    type(monolis_mat) :: monoTREE
    integer(kint) :: fact_order(:)
    integer(kint) :: fact_array_index(:)
    integer(kint), allocatable :: add_location(:)
    integer(kint) :: N, NZ
    integer(kint) :: child_id, parent_id
    integer(kint) :: n_child, n_parent, frontal_size
    integer(kint) :: i, j, k, l, iSorg, iS, iE, jS, jE, in, jn, ln
    integer(kint), allocatable :: child_rows(:)
    integer(kint), allocatable :: parent_rows(:)
    logical, allocatable :: is_add(:)

    !> extended add 演算の fact_array における足し込み先 index

    N = monoTREE%N
    NZ = fact_array_index(N + 1)

    allocate(add_location(NZ), source = 0)

    !> Frontal 行列の最初の行は 0 をおく
    !> 短い方でループを回す
    !> 長い方を検索する
    do k = 1, N - 1
      child_id = k
      in = monoTREE%SCSR%indexU(k) + 1
      parent_id = monoTREE%SCSR%itemU(in + 1)

      n_child = monoTREE%SCSR%indexU(child_id + 1) - monoTREE%SCSR%indexU(child_id) - 1
      iS = monoTREE%SCSR%indexU(child_id) + 2
      iE = monoTREE%SCSR%indexU(child_id + 1)
      allocate(child_rows (n_child), source = 0)
      do i = 1, n_child
        child_rows(i) = monoTREE%SCSR%itemU(iS + i - 1)
      enddo

      n_parent = monoTREE%SCSR%indexU(parent_id + 1) - monoTREE%SCSR%indexU(parent_id)
      jS = monoTREE%SCSR%indexU(parent_id) + 1
      jE = monoTREE%SCSR%indexU(parent_id + 1)
      allocate(parent_rows(n_parent), source = 0)
      allocate(is_add(n_parent), source = .false.)
      do i = 1, n_parent
        parent_rows(i) = monoTREE%SCSR%itemU(jS + i - 1)
      enddo

      jS = 1
      do i = 1, n_child
        in = child_rows(i)
        do j = jS, n_parent
          jn = parent_rows(j)
          if(in == jn)then
            is_add(j) = .true.
            jS = j
          endif
        enddo
      enddo

      !write(*,*)"child_rows : ", child_rows
      !write(*,*)"parent_rows: ", parent_rows
      !write(*,*)"is_add     : ", is_add

      frontal_size = n_child + 1
      iS = fact_array_index(k) + frontal_size
      jS = fact_array_index(parent_id) 
      in = 0
      jn = 0
      do i = 1, n_parent
        aa:do j = i, n_parent
          jn = jn + 1
          if(is_add(i) .and. is_add(j))then
            in = in + 1
            add_location(iS + in) = jS + jn
          endif
        enddo aa
      enddo

      deallocate(child_rows)
      deallocate(parent_rows)
      deallocate(is_add)
    enddo

    do i = 1, N + 1
      if(fact_array_index(i) < 0) stop "minus fact_array_index"
    enddo

    do i = 1, NZ
      if(add_location(i) < 0) stop "minus add_location"
    enddo

!write(*,*)"add_location"
!write(*,"(20i4)")add_location
!call flush()
!call sleep(1)
  end subroutine monolis_matrix_get_add_location
end module mod_monolis_fact_analysis
