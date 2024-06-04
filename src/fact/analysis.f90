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
    logical, allocatable :: is_used(:)

    N = monoTREE%N
    allocate(super_node_level(N), source = 0)
    allocate(is_used(N), source = .false.)

    jn = 1
    super_node_level(N) = 1
    do i = N-1, 1, -1
      j = monoTREE%SCSR%indexU(i) + 2
      in = monoTREE%SCSR%itemU(j)

      d1 = monoTREE%SCSR%indexU(i  + 1) - monoTREE%SCSR%indexU(i)
      d2 = monoTREE%SCSR%indexU(in + 1) - monoTREE%SCSR%indexU(in)

      !> check super node
      if(d1 == d2 + 1 .and. .not. is_used(in))then
        super_node_level(i) = jn
        is_used(in) = .true.
      else
        jn = jn + 1
        super_node_level(i) = jn
      endif
    enddo
    n_super_node = jn

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

    aa:do i = 1, n_super_node - 1
      in = 0
      do j = 1, n_super_node
        if(super_node_parent_id(i) == super_node_id(j))then
          cycle aa
        elseif(super_node_parent_id(i) > super_node_id(j))then
          in = super_node_id(j)
        endif
      enddo
      super_node_parent_id(i) = in
    enddo aa

write(*,*)"super_node_id"
write(*,"(20i4)")super_node_id
write(*,*)"super_node_size"
write(*,"(20i4)")super_node_size
write(*,*)"super_node_parent_id"
write(*,"(20i4)")super_node_parent_id
  end subroutine monolis_matrix_get_super_node_information

  subroutine monolis_matrix_get_factorize_array(monoTREE, n_super_node, super_node_id, &
      & n_fact_array, fact_array, fact_array_index, front_size)
    implicit none
    type(monolis_mat) :: monoTREE
    integer(kint) :: n_super_node
    integer(kint) :: super_node_id(:)
    integer(kint) :: n_fact_array
    real(kdouble), allocatable :: fact_array(:)
    integer(kint), allocatable :: fact_array_index(:)
    integer(kint), allocatable :: front_size(:)
    integer(kint) :: i, j, in, jn

    n_fact_array = 0

    allocate(fact_array_index(n_super_node + 1), source = 0)
    allocate(front_size(n_super_node), source = 0)

    do j = 1, n_super_node
      i = super_node_id(j)
      in = monoTREE%SCSR%indexU(i + 1) - monoTREE%SCSR%indexU(i)
      jn = in*(in + 1)/2
      n_fact_array = n_fact_array + jn
      fact_array_index(j + 1) = jn
      front_size(j) = in
    enddo

    allocate(fact_array(n_fact_array), source = 0.0d0)

    do j = 1, n_super_node
      fact_array_index(j + 1) = fact_array_index(j + 1) + fact_array_index(j)
    enddo
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
      ln = fact_array_index(k)
      do m = 1, super_node_size(k)
        if(m == 1)then
          kn = super_node_id(k)
        else
          j = monoTREE%SCSR%indexU(kn) + 2 
          kn = monoTREE%SCSR%itemU(j)
        endif
        iS = monoTREE%SCSR%indexU(kn) + 1
        iE = monoTREE%SCSR%indexU(kn + 1)
        jS = monoMAT%CSR%index(kn) + 1
        jE = monoMAT%CSR%index(kn + 1)
        aa:do i = iS, iE
          ln = ln + 1
          in = monoTREE%SCSR%itemU(i)
          if(in < kn) cycle aa
          do j = jS, jE
            jn = monoMAT%CSR%item(j)
            if(jn == in)then
              fact_array(ln) = monoMAT%R%A(j)
              jS = j + 1
              cycle aa
            endif
          enddo
        enddo aa
      enddo
    enddo
  end subroutine monolis_matrix_set_value_of_factorize_array

  subroutine monolis_matrix_get_add_location(monoTREE, n_super_node, super_node_id, super_node_size, &
    & super_node_parent_id, fact_array_index, front_size, add_location)
    implicit none
    type(monolis_mat) :: monoTREE
    integer(kint) :: n_super_node
    integer(kint) :: super_node_id(:)
    integer(kint) :: super_node_size(:)
    integer(kint) :: super_node_parent_id(:)
    integer(kint) :: fact_array_index(:)
    integer(kint) :: front_size(:)
    integer(kint), allocatable :: add_location(:)
    integer(kint), allocatable :: update_id(:)
    integer(kint) :: NZ
    integer(kint) :: child_id, parent_id, near_parent_id
    integer(kint) :: n_child, n_parent, frontal_size
    integer(kint) :: i, j, k, l, m, iSorg, iS, iE, jS, jE, in, jn, ln
    integer(kint), allocatable :: child_rows(:)
    integer(kint), allocatable :: parent_rows(:)
    logical, allocatable :: is_add(:)

    !> get update_id
    allocate(update_id(n_super_node), source = 0)

    do k = 1, n_super_node - 1
      do m = 1, super_node_size(k)
        if(m == 1)then
          in = super_node_id(k)
          update_id(k) = in
        else
          j = monoTREE%SCSR%indexU(in) + 2 
          in = monoTREE%SCSR%itemU(j)
          update_id(k) = in
        endif
      enddo
    enddo

    !> extended add 演算の fact_array における足し込み先 index

    NZ = fact_array_index(n_super_node + 1)
    allocate(add_location(NZ), source = 0)

    !> Frontal 行列の最初の行は 0 をおく
    !> 短い方でループを回す
    !> 長い方を検索する
    do m = 1, n_super_node - 1
      child_id = update_id(m)
      n_child = monoTREE%SCSR%indexU(child_id + 1) - monoTREE%SCSR%indexU(child_id) - 1
      iS = monoTREE%SCSR%indexU(child_id) + 2
      allocate(child_rows(n_child), source = 0)
      do i = 1, n_child
        child_rows(i) = monoTREE%SCSR%itemU(iS + i - 1)
      enddo

      parent_id = super_node_parent_id(m)
      n_parent = monoTREE%SCSR%indexU(parent_id + 1) - monoTREE%SCSR%indexU(parent_id)
      jS = monoTREE%SCSR%indexU(parent_id) + 1
      allocate(parent_rows(n_parent), source = 0)
      do i = 1, n_parent
        parent_rows(i) = monoTREE%SCSR%itemU(jS + i - 1)
      enddo

      allocate(is_add(n_parent), source = .false.)
      jS = 1
      do i = 1, n_child
        in = child_rows(i)
        do j = jS, n_parent
          jn = parent_rows(j)
          if(in == jn)then
            is_add(j) = .true.
            jS = j + 1
          endif
        enddo
      enddo

      iS = fact_array_index(m)
      do i = 1, super_node_size(m)
        iS = iS + front_size(m) + 1 - i 
      enddo

      near_parent_id = 0
      do i = 1, n_super_node
        if(parent_id == super_node_id(i))then
          near_parent_id = i
          exit
        endif
      enddo
      jS = fact_array_index(near_parent_id) 

      in = 0
      jn = 0
      do i = 1, n_parent
        if(.not. is_add(i))then
          jn = jn + n_parent - i + 1
          cycle
        endif
        do j = i, n_parent
          jn = jn + 1
          if(is_add(j))then
            in = in + 1
            add_location(iS + in) = jS + jn
          endif
        enddo 
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
  end subroutine monolis_matrix_get_add_location
end module mod_monolis_fact_analysis
