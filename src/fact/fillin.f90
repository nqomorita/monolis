module mod_monolis_fact_fillin
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc

  implicit none

  type monolis_fillin
    integer(kint) :: n_ancestor
    integer(kint), pointer :: ancestor(:)
  end type monolis_fillin

contains

  subroutine monolis_matrix_get_fillin(monoMAT, monoTREE, is_asym, is_fillin)
    implicit none
    type(monolis_mat) :: monoMAT
    type(monolis_mat) :: monoTREE
    type(monolis_fillin), allocatable:: tree(:)
    integer(kint), pointer :: array(:)
    integer(kint), allocatable :: fillin_mask(:)
    integer(kint), allocatable :: child_mask(:)
    integer(kint), allocatable :: parent_mask(:)
    integer(kint), allocatable :: parent_id(:)
    integer(kint), allocatable :: perm(:)
    integer(kint), allocatable :: temp(:)
    integer(kint), allocatable :: index(:)
    integer(kint), allocatable :: is_used(:)
    integer(kint) :: N, M, NPU
    integer(kint) :: i, j, k, p, jS, jE, in, c, d, s
    integer(kint) :: nbytes
    integer(kint) :: is, ie
    integer(kint) :: range, parent, child
    integer(kint) :: bit = kint*8
    integer(kint), allocatable :: count(:)
    logical :: is_asym, is_fillin, is_merge

    N = monoMAT%N
    allocate(tree(N))

    nbytes = N/bit+1
    allocate(child_mask (nbytes), source = 0)
    allocate(parent_mask(nbytes), source = 0)
    allocate(fillin_mask(nbytes), source = 0)

    !# fillin 決定
    do i = 1, N
      tree(i)%n_ancestor = 0

      in = 0
      jS = monoMAT%CSR%index(i) + 1
      jE = monoMAT%CSR%index(i + 1)
      do j = jS, jE
        if(i < monoMAT%CSR%item(j) .and. monoMAT%CSR%item(j) <= N)then
          in = in + 1
        endif
      enddo

      tree(i)%n_ancestor = in
      allocate(tree(i)%ancestor(in))

      in = 0
      do j = jS, jE
        if(i < monoMAT%CSR%item(j) .and. monoMAT%CSR%item(j) <= N)then
          in = in + 1
          tree(i)%ancestor(in) = monoMAT%CSR%item(j)
        endif
      enddo
    enddo

    if(is_fillin)then
      do i = 1, N
        if(tree(i)%n_ancestor < 2) cycle

        parent = tree(i)%ancestor(1)
        call merge_ancestor_union(tree(i), tree(parent), child_mask, parent_mask, fillin_mask, nbytes, bit)
      enddo
    endif

    !> 緩和スーパーノード決定
    if(.false.)then
      ! tree(i) と tree(parent) 両方の ancestor を和集合で更新
      is_merge = .true.
      do while(is_merge)
        do i = N - 1, 1, -1
          parent = tree(i)%ancestor(1)
          call merge_ancestor_union_both(tree, i, parent, child_mask, parent_mask, fillin_mask, &
            nbytes, bit, is_merge)
        enddo
      enddo
    endif

    !> SCSR allocation part
    monoTREE%N = monoMAT%N
    monoTREE%NDOF = monoMAT%NDOF

    allocate(monoTREE%SCSR%indexU(N+1))
    NPU = 0
    monoTREE%SCSR%indexU(1) = 0
    do i = 1, N
      monoTREE%SCSR%indexU(i+1) = monoTREE%SCSR%indexU(i) + tree(i)%n_ancestor + 1
      NPU = NPU + tree(i)%n_ancestor + 1
    enddo

    allocate(monoTREE%SCSR%itemU(NPU))
    in = 0
    do i = 1, N
      in = in + 1
      monoTREE%SCSR%itemU(in) = i
      do j = 1, tree(i)%n_ancestor
        in = in + 1
        monoTREE%SCSR%itemU(in) = tree(i)%ancestor(j)
      enddo
    enddo

    deallocate(child_mask )
    deallocate(parent_mask)
    deallocate(fillin_mask)
  end subroutine monolis_matrix_get_fillin

  subroutine merge_ancestor_union(tree_child, tree_parent, child_mask, parent_mask, fillin_mask, nbytes, bit)
    implicit none
    type(monolis_fillin), intent(in) :: tree_child
    type(monolis_fillin), intent(inout) :: tree_parent
    integer(kint), intent(inout) :: child_mask(:)
    integer(kint), intent(inout) :: parent_mask(:)
    integer(kint), intent(inout) :: fillin_mask(:)
    integer(kint), intent(in) :: nbytes, bit
    integer(kint), pointer :: array(:)
    integer(kint) :: i, j, k, c, in, is, ie, range, ancestor_id

    ! 初期化
    is = 1
    child_mask(is:nbytes) = 0
    parent_mask(is:nbytes) = 0
    
    ! 子ノードの祖先をビットマスクに設定（最初の祖先を除く）
    range = 0
    do j = 2, tree_child%n_ancestor
      ancestor_id = tree_child%ancestor(j)
      call set_bit(ancestor_id, child_mask, range, bit)
    enddo
    
    ! 親ノードの祖先をビットマスクに設定
    do j = 1, tree_parent%n_ancestor
      ancestor_id = tree_parent%ancestor(j)
      call set_bit(ancestor_id, parent_mask, range, bit)
    enddo
    
    ie = range/bit + 1
    
    ! 和集合を計算
    fillin_mask(is:ie) = ior(child_mask(is:ie), parent_mask(is:ie))
    
    ! ビット数をカウント
    c = 0
    do j = is, ie
      c = c + popcnt(fillin_mask(j))
    enddo
    
    ! 新しい祖先配列を作成
    if(0 < c)then
      allocate(array(c))
      in = 0
      do j = is, ie
        do k = 1, popcnt(fillin_mask(j))
          in = in + 1
          c = popcnt( iand(fillin_mask(j), - fillin_mask(j)) -1 )
          fillin_mask(j) = ibclr(fillin_mask(j),c)
          array(in) = bit*(j-1)+c
        enddo
      enddo
      
      ! 親ノードの祖先を更新
      deallocate(tree_parent%ancestor)
      tree_parent%ancestor => array
      tree_parent%n_ancestor = in
    endif
  end subroutine merge_ancestor_union

  subroutine merge_ancestor_union_both(tree, idx1, idx2, child_mask, parent_mask, fillin_mask, &
    nbytes, bit, is_merge)
    implicit none
    type(monolis_fillin), intent(inout) :: tree(:)
    integer(kint), intent(in) :: idx1, idx2
    integer(kint), intent(inout) :: child_mask(:)
    integer(kint), intent(inout) :: parent_mask(:)
    integer(kint), intent(inout) :: fillin_mask(:)
    integer(kint), intent(in) :: nbytes, bit
    logical, intent(inout) :: is_merge
    integer(kint), pointer :: array1(:), array2(:)
    integer(kint) :: i, j, k, c, num, i1, i2, is, ie, range, ancestor_id
    real(kdouble) :: z_ratio
    logical :: do_merge

    if(idx1 >= idx2) stop "merge_ancestor_union_both"

    !# 初期化
    is = 1
    child_mask(is:nbytes) = 0
    parent_mask(is:nbytes) = 0
    
    range = 0
    do j = 1, tree(idx1)%n_ancestor
      ancestor_id = tree(idx1)%ancestor(j)
      call set_bit(ancestor_id, child_mask, range, bit)
    enddo

    do j = 1, tree(idx2)%n_ancestor
      ancestor_id = tree(idx2)%ancestor(j)
      call set_bit(ancestor_id, parent_mask, range, bit)
    enddo
    
    ie = range/bit + 1
    
    ! 和集合を計算
    fillin_mask(is:ie) = ior(child_mask(is:ie), parent_mask(is:ie))
    
    ! ビット数をカウント
    num = 0
    do j = is, ie
      num = num + popcnt(fillin_mask(j))
    enddo

    z_ratio = dble(tree(idx1)%n_ancestor + tree(idx2)%n_ancestor) / dble(2*num)

    !# 新しい祖先配列を作成
    do_merge = .false.

    if(num <= 4) do_merge = .true.
    if(num <= 16 .and. z_ratio >= 0.5d0) do_merge = .true.
    if(num <= 48 .and. z_ratio >= 0.8d0) do_merge = .true.
    if(z_ratio >= 0.95d0) do_merge = .true.
    if(tree(idx1)%n_ancestor + tree(idx2)%n_ancestor == 2*num - 1) do_merge = .false.

    is_merge = .false.
    if(0 < num .and. do_merge)then
      allocate(array1(num))
      allocate(array2(num-1))

      i1 = 0
      i2 = 0
      do j = is, ie
        do k = 1, popcnt(fillin_mask(j))
          c = popcnt( iand(fillin_mask(j), - fillin_mask(j)) -1 )
          fillin_mask(j) = ibclr(fillin_mask(j),c)
          ancestor_id = bit*(j-1)+c
          i1 = i1 + 1
          array1(i1) = ancestor_id

          if(ancestor_id /= idx2)then
            i2 = i2 + 1
            array2(i2) = ancestor_id
          endif
        enddo
      enddo

      deallocate(tree(idx1)%ancestor)
      deallocate(tree(idx2)%ancestor)

      tree(idx1)%n_ancestor = num
      tree(idx1)%ancestor => array1

      tree(idx2)%n_ancestor = num - 1
      tree(idx2)%ancestor => array2
      
      ! マージが実行されたことを記録
      is_merge = .true.
    endif
  end subroutine merge_ancestor_union_both

  subroutine set_bit(in, parent_mask, range, bit)
    implicit none
    integer(kint) :: in, range, ie
    integer(kint) :: parent_mask(:)
    integer(kint) :: bit 
    ie = in/bit + 1
    parent_mask(ie) = ibset(parent_mask(ie), mod(in, bit))
    range = max(range, in)
  end subroutine set_bit

  subroutine get_true_bit_number(array, iS, iE, n)
    implicit none
    integer(kint) :: array(:)
    integer(kint) :: iS, iE, n, i
    n = 0
    do i = iS, iE
      n = n + popcnt(array(i))
    enddo
  end subroutine get_true_bit_number

  subroutine monolis_matrix_alloc_with_fillin(monoTREE, is_asym)
    implicit none
    type(monolis_mat) :: monoTREE
    logical :: is_asym
    integer(kint) :: N, NDOF, NZ

    NDOF = monoTREE%NDOF
    N = monoTREE%N
    NZ = monoTREE%SCSR%indexU(N + 1)

    call monolis_palloc_R_1d(monoTREE%R%A, NZ*NDOF*NDOF)
  end subroutine monolis_matrix_alloc_with_fillin
end module mod_monolis_fact_fillin
