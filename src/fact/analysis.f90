module mod_monolis_fact_analysis
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc

  implicit none

contains

  subroutine monolis_matrix_get_super_node_information(monoTREE, n_super_node, &
    & super_node_id, super_node_size, super_node_parent_id)
    implicit none
    type(monolis_mat) :: monoTREE
    integer(kint), allocatable :: super_node_level(:)
    integer(kint), allocatable :: super_node_id(:)
    integer(kint), allocatable :: super_node_size(:)
    integer(kint), allocatable :: super_node_parent_id(:)
    integer(kint), allocatable :: temp1(:), temp2(:)
    integer(kint) :: n_super_node
    integer(kint) :: N, i, j, k, in, jn, kn, d1, d2

    N = monoTREE%N
    allocate(super_node_level(N), source = 0)

    jn = 1
    super_node_level(N) = 1
    do i = N-1, 1, -1
      j = monoTREE%SCSR%indexU(i) + 2
      in = monoTREE%SCSR%itemU(j)

      d1 = monoTREE%SCSR%indexU(i  + 1) - monoTREE%SCSR%indexU(i)
      d2 = monoTREE%SCSR%indexU(in + 1) - monoTREE%SCSR%indexU(in)

      !> check super node
      if(d1 == d2 + 1)then
        super_node_level(i) = jn
      else
        jn = jn + 1
        super_node_level(i) = jn
      endif
      !write(*,*)d1, d2
    enddo
    n_super_node = jn

!write(*,*)"super_node_level", super_node_level

    allocate(super_node_id(n_super_node), source = 0)
    allocate(super_node_size(n_super_node), source = 0)
    allocate(temp1(n_super_node), source = 0)
    allocate(temp2(n_super_node), source = 0)
    temp1 = N + 1

    !> 逆順に並び替え
    do i = 1, N
      in = super_node_level(i)
      temp2(in) = temp2(in) + 1
      if(i < temp1(in)) temp1(in) = i
    enddo

    do i = 1, n_super_node
      super_node_id(i) = temp1(n_super_node - i + 1)
      super_node_size(i) = temp2(n_super_node - i + 1)
    enddo
!write(*,*)"super_node_id  ", super_node_id
!write(*,*)"super_node_size", super_node_size

    !> supernode の親ノードの取得
    allocate(super_node_parent_id(n_super_node), source = 0)
    do i = 1, n_super_node - 1
      kn = super_node_id(i)
      do j = 1, super_node_size(i)
        k = monoTREE%SCSR%indexU(kn) + 2 
        kn = monoTREE%SCSR%itemU(k)
      enddo
      super_node_parent_id(i) = kn
    enddo
!write(*,*)"super_node_parent_id", super_node_parent_id
  end subroutine monolis_matrix_get_super_node_information

  subroutine monolis_matrix_get_factorize_order(monoTREE, n_super_node, super_node_id, super_node_size)
    implicit none
    type(monolis_mat) :: monoTREE
    integer(kint) :: super_node_id(:)
    integer(kint) :: super_node_size(:)
    integer(kint), allocatable :: fact_level(:)
    integer(kint), allocatable :: fact_order(:)
    integer(kint) :: n_super_node
    integer(kint) :: N, i, j, in, jn, iS

    N = monoTREE%N

    allocate(fact_level(N), source = 0)

    fact_level(N) = 1

    do i = N-1, 1, -1
      j = monoTREE%SCSR%indexU(i) + 2
      in = monoTREE%SCSR%itemU(j)
      fact_level(i) = fact_level(in) + 1
    enddo

    do i = 1, N
      fact_level(i) = N + 1 - fact_level(i)
    enddo

    allocate(fact_order(n_super_node), source = 0)

    do i = 1, n_super_node
      in = super_node_id(i)
      fact_order(i) = fact_level(in)
    enddo

!write(*,*)"fact_level", fact_level
!write(*,*)"fact_order", fact_order

    call monolis_qsort_I_2d(fact_order, super_node_id, 1, n_super_node)

    do i = 1, n_super_node
      in = super_node_id(i)
      fact_order(i) = fact_level(in)
    enddo

    call monolis_qsort_I_2d(fact_order, super_node_size, 1, n_super_node)

!write(*,*)"super_node_id  ", super_node_id
!write(*,*)"super_node_size", super_node_size
  end subroutine monolis_matrix_get_factorize_order

  subroutine monolis_matrix_get_factorize_array(monoTREE, n_super_node, super_node_id, &
      & n_fact_array, fact_array, fact_array_index)
    implicit none
    type(monolis_mat) :: monoTREE
    integer(kint) :: n_super_node
    integer(kint) :: super_node_id(:)
    integer(kint) :: n_fact_array
    real(kdouble), allocatable :: fact_array(:)
    integer(kint), allocatable :: fact_array_index(:)
    integer(kint) :: i, j, in, jn

    n_fact_array = 0

    allocate(fact_array_index(n_super_node + 1), source = 0)

    do j = 1, n_super_node
      i = super_node_id(j)
      in = monoTREE%SCSR%indexU(i + 1) - monoTREE%SCSR%indexU(i)
      jn = in*(in + 1)/2
      n_fact_array = n_fact_array + jn
      fact_array_index(j + 1) = jn
    enddo

    allocate(fact_array(n_fact_array), source = 0.0d0)

    do j = 1, n_super_node
      fact_array_index(j + 1) = fact_array_index(j + 1) + fact_array_index(j)
    enddo

!write(*,*)"n_fact_array", n_fact_array
!write(*,*)"fact_array_index", fact_array_index
  end subroutine monolis_matrix_get_factorize_array

  subroutine monolis_matrix_set_value_of_factorize_array(monoMAT, monoTREE, &
    & n_super_node, super_node_id, super_node_size, fact_array, fact_array_index)
    implicit none
    type(monolis_mat) :: monoMAT
    type(monolis_mat) :: monoTREE
    integer(kint) :: n_super_node
    integer(kint) :: super_node_id(:)
    integer(kint) :: super_node_size(:)
    real(kdouble) :: fact_array(:)
    integer(kint) :: fact_array_index(:)
    integer(kint) :: i, j, k, l, m, iSorg, iS, iE, jS, jE, in, jn, ln, kn

    do k = 1, n_super_node
!write(*,*)"super_node_size(k)", super_node_size(k)
      ln = fact_array_index(k)
      do m = 1, super_node_size(k)
        if(m == 1)then
          i = super_node_id(k)
          kn = i
        else
          j = monoTREE%SCSR%indexU(kn) + 2 
          i = monoTREE%SCSR%itemU(j)
          kn = i
        endif
!write(*,*)"item", kn
        iS = monoTREE%SCSR%indexU(kn) + 1
        iE = monoTREE%SCSR%indexU(kn + 1)
        jS = monoMAT%CSR%index(kn) + 1
        jE = monoMAT%CSR%index(kn + 1)
        aa:do j = jS, jE
          jn = monoMAT%CSR%item(j)
          if(jn < kn) cycle aa
          do i = iS, iE
            in = monoTREE%SCSR%itemU(i)
            if(jn == in)then
              ln = ln + 1
!write(*,*)"ln", ln, monoMAT%R%A(j)
              fact_array(ln) = monoMAT%R%A(j)
              cycle aa
            endif
          enddo
        enddo aa
      enddo
    enddo

!write(*,*)"fact_array A"
!write(*,"(1p10e12.3)")fact_array
!call sleep(1)
  end subroutine monolis_matrix_set_value_of_factorize_array

  subroutine monolis_matrix_get_add_location(monoTREE, n_super_node, super_node_id, super_node_size, &
    & super_node_parent_id, fact_array_index, add_location)
    implicit none
    type(monolis_mat) :: monoTREE
    integer(kint) :: n_super_node
    integer(kint) :: super_node_id(:)
    integer(kint) :: super_node_size(:)
    integer(kint) :: super_node_parent_id(:)
    integer(kint) :: fact_array_index(:)
    integer(kint), allocatable :: add_location(:)
    integer(kint) :: NZ
    integer(kint) :: child_id, parent_id
    integer(kint) :: n_child, n_parent, frontal_size
    integer(kint) :: i, j, k, l, m, iSorg, iS, iE, jS, jE, in, jn, ln
    integer(kint), allocatable :: child_rows(:)
    integer(kint), allocatable :: parent_rows(:)
    logical, allocatable :: is_add(:)

    !> extended add 演算の fact_array における足し込み先 index

    NZ = fact_array_index(n_super_node + 1)
    allocate(add_location(NZ), source = 0)

    !> Frontal 行列の最初の行は 0 をおく
    !> 短い方でループを回す
    !> 長い方を検索する
    do m = 1, n_super_node - 1
      k = super_node_id(m)
      child_id = k + super_node_size(m) - 1
      n_child = monoTREE%SCSR%indexU(child_id + 1) - monoTREE%SCSR%indexU(child_id) - 1
      iS = monoTREE%SCSR%indexU(child_id) + 2
      iE = monoTREE%SCSR%indexU(child_id + 1)
      allocate(child_rows (n_child), source = 0)
      do i = 1, n_child
        child_rows(i) = monoTREE%SCSR%itemU(iS + i - 1)
      enddo

      parent_id = super_node_parent_id(m)
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
      
      iS = fact_array_index(k)
      do i = 1, super_node_size(m)
        iS = iS + n_child + super_node_size(m) + 1 - i 
      enddo

      !> 最適化必須
      in = 0
      do i = 1, n_super_node
        if(parent_id == super_node_id(i)) in = i
      enddo

      jS = fact_array_index(in) 
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

    do i = 1, n_super_node + 1
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
