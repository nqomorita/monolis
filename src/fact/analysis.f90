module mod_monolis_fact_analysis
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc

  implicit none

contains

  subroutine monolis_matrix_get_factorize_order(monoTREE, fact_order)
    implicit none
    type(monolis_mat) :: monoTREE
    integer(kint), allocatable :: fact_level(:)
    integer(kint), allocatable :: fact_order(:)
    integer(kint) :: N, i, j, in

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

    call monolis_qsort_I_2d(fact_level, fact_order, 1, N)
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
      i = fact_order(j)
      in = monoTREE%SCSR%indexU(i + 1) - monoTREE%SCSR%indexU(i)
      jn = in*(in + 1)/2
      n_fact_array = n_fact_array + jn
      fact_array_index(j + 1) = jn
    enddo

    allocate(fact_array(n_fact_array), source = 0.0d0)

    do j = 1, N
      fact_array_index(j + 1) = fact_array_index(j + 1) + fact_array_index(j)
    enddo
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
      iS = fact_array_index(k) + 1 + frontal_size
      jS = fact_array_index(parent_id) + 1
      in = 0
      do i = 1, n_child
        if(.not. is_add(i)) cycle
        aa:do j = i, n_child
          if(.not. is_add(j)) cycle aa
          add_location(iS + in) = jS + in
          in = in + 1
        enddo aa
      enddo

      deallocate(child_rows)
      deallocate(parent_rows)
      deallocate(is_add)
    enddo
  end subroutine monolis_matrix_get_add_location
end module mod_monolis_fact_analysis
