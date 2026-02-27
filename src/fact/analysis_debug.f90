






  !> @ingroup fact_analysis
  !> フィルイン前の非ゼロパターンを可視化（デバッグ用）
  subroutine monolis_debug_print_original_pattern(monoMAT, filename)
    implicit none
    type(monolis_mat) :: monoMAT
    character(len=*) :: filename
    integer(kint) :: N, i, j, k, unit_num
    integer(kint) :: iS, iE, entry
    logical, allocatable :: pattern(:,:)
    
    N = monoMAT%N
    allocate(pattern(N, N))
    pattern = .false.
    
    ! 元の非ゼロパターンを記録
    do i = 1, N
      iS = monoMAT%CSR%index(i) + 1
      iE = monoMAT%CSR%index(i + 1)
      do j = iS, iE
        entry = monoMAT%CSR%item(j)
        if (entry >= 1 .and. entry <= N) then
          pattern(i, entry) = .true.
        endif
      enddo
    enddo
    
    ! ファイルに出力
    open(newunit=unit_num, file=trim(filename), status='replace')
    write(unit_num, '(a)') '# Original nonzero pattern (before fill-in)'
    write(unit_num, '(a,i0)') '# Matrix size: ', N
    write(unit_num, '(a)') '# Format: row col (1-based indexing)'
    
    do i = 1, N
      do j = 1, N
        if (pattern(i, j)) then
          write(unit_num, '(i0,1x,i0)') i, j
        endif
      enddo
    enddo
    
    close(unit_num)
    deallocate(pattern)
    
    write(*,'(a,a)') 'Original pattern saved to: ', trim(filename)
  end subroutine monolis_debug_print_original_pattern

  !> @ingroup fact_analysis
  !> フィルイン後の非ゼロパターンを可視化（デバッグ用）
  subroutine monolis_debug_print_fillin_pattern(monoTREE, filename)
    implicit none
    type(monolis_mat) :: monoTREE
    character(len=*) :: filename
    integer(kint) :: N, i, j, k, unit_num
    integer(kint) :: iS, iE, entry
    logical, allocatable :: pattern(:,:)
    
    N = monoTREE%N
    allocate(pattern(N, N))
    pattern = .false.
    
    ! フィルイン後の非ゼロパターンを記録（SCSR構造から）
    do i = 1, N
      iS = monoTREE%SCSR%indexU(i) + 1
      iE = monoTREE%SCSR%indexU(i + 1)
      do j = iS, iE
        entry = monoTREE%SCSR%itemU(j)
        if (entry >= 1 .and. entry <= N) then
          pattern(i, entry) = .true.
          pattern(entry, i) = .true.  ! 対称性を考慮
        endif
      enddo
    enddo
    
    ! ファイルに出力
    open(newunit=unit_num, file=trim(filename), status='replace')
    write(unit_num, '(a)') '# Fill-in nonzero pattern (after symbolic factorization)'
    write(unit_num, '(a,i0)') '# Matrix size: ', N
    write(unit_num, '(a)') '# Format: row col (1-based indexing)'
    
    do i = 1, N
      do j = 1, N
        if (pattern(i, j)) then
          write(unit_num, '(i0,1x,i0)') i, j
        endif
      enddo
    enddo
    
    close(unit_num)
    deallocate(pattern)
    
    write(*,'(a,a)') 'Fill-in pattern saved to: ', trim(filename)
  end subroutine monolis_debug_print_fillin_pattern

  !> @ingroup fact_analysis
  !> スーパーノード情報を可視化（デバッグ用）
  subroutine monolis_debug_print_supernode_info(n_super_node, super_node_id, &
      super_node_size, super_node_parent_id, filename)
    implicit none
    integer(kint) :: n_super_node
    integer(kint) :: super_node_id(:)
    integer(kint) :: super_node_size(:)
    integer(kint) :: super_node_parent_id(:)
    character(len=*) :: filename
    integer(kint) :: i, unit_num
    
    ! ファイルに出力
    open(newunit=unit_num, file=trim(filename), status='replace')
    write(unit_num, '(a)') '# Supernode information'
    write(unit_num, '(a,i0)') '# Number of supernodes: ', n_super_node
    write(unit_num, '(a)') '# Format: supernode_id representative_node size parent_supernode'
    
    do i = 1, n_super_node
      write(unit_num, '(i0,1x,i0,1x,i0,1x,i0)') i, super_node_id(i), &
        super_node_size(i), super_node_parent_id(i)
    enddo
    
    close(unit_num)
    
    write(*,'(a,a)') 'Supernode info saved to: ', trim(filename)
    write(*,'(a,i0,a)') 'Total supernodes: ', n_super_node, ' nodes'
  end subroutine monolis_debug_print_supernode_info

  !> @ingroup fact_analysis
  !> elimination treeを可視化（デバッグ用）
  subroutine monolis_debug_print_elimination_tree(monoTREE, filename)
    implicit none
    type(monolis_mat) :: monoTREE
    character(len=*) :: filename
    integer(kint) :: N, i, j, unit_num
    integer(kint) :: iS, iE, parent_node
    
    N = monoTREE%N
    
    ! ファイルに出力
    open(newunit=unit_num, file=trim(filename), status='replace')
    write(unit_num, '(a)') '# Elimination tree structure'
    write(unit_num, '(a,i0)') '# Matrix size: ', N
    write(unit_num, '(a)') '# Format: child_node parent_node'
    
    do i = 1, N
      iS = monoTREE%SCSR%indexU(i) + 1
      iE = monoTREE%SCSR%indexU(i + 1)
      
      parent_node = 0
      if (iE > iS + 1) then
        parent_node = monoTREE%SCSR%itemU(iS + 1)
        if (parent_node > i .and. parent_node <= N) then
          write(unit_num, '(i0,1x,i0)') i, parent_node
        endif
      endif
    enddo
    
    close(unit_num)
    
    write(*,'(a,a)') 'Elimination tree saved to: ', trim(filename)
  end subroutine monolis_debug_print_elimination_tree

  !> @ingroup fact_analysis
  !> 全ての可視化を一括実行（デバッグ用メイン関数）
  subroutine monolis_debug_visualize_all_patterns(monoMAT, monoTREE, n_super_node, &
      super_node_id, super_node_size, super_node_parent_id, prefix)
    implicit none
    type(monolis_mat) :: monoMAT
    type(monolis_mat) :: monoTREE
    integer(kint) :: n_super_node
    integer(kint) :: super_node_id(:)
    integer(kint) :: super_node_size(:)
    integer(kint) :: super_node_parent_id(:)
    character(len=*) :: prefix
    character(len=256) :: filename
    
    write(*,'(a)') 'Starting pattern visualization...'
    
    ! 元パターンの出力
    write(filename, '(a,a)') trim(prefix), '_original_pattern.dat'
    call monolis_debug_print_original_pattern(monoMAT, filename)
    
    ! フィルインパターンの出力
    write(filename, '(a,a)') trim(prefix), '_fillin_pattern.dat'
    call monolis_debug_print_fillin_pattern(monoTREE, filename)
    
    ! スーパーノード情報の出力
    write(filename, '(a,a)') trim(prefix), '_supernode_info.dat'
    call monolis_debug_print_supernode_info(n_super_node, super_node_id, &
      super_node_size, super_node_parent_id, filename)
    
    ! elimination treeの出力
    write(filename, '(a,a)') trim(prefix), '_elimination_tree.dat'
    call monolis_debug_print_elimination_tree(monoTREE, filename)
  end subroutine monolis_debug_visualize_all_patterns
