module mod_monolis_spmat_convert_sym
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc

  implicit none

contains

  !> @ingroup matrix
  !> 行列を対称行列に変換（1/2 (A + A^T)）
  subroutine monolis_matrix_convert_to_symmetric_R(monolis, com)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_structure), intent(inout) :: monolis
    !> [in,out] monolis 構造体
    type(monolis_COM), intent(inout) :: com

    call monolis_matrix_convert_to_symmetric_inner_R(monolis%MAT)
    call monolis_matrix_convert_to_symmetric_outer_R(monolis%MAT, com)
  end subroutine monolis_matrix_convert_to_symmetric_R

  !> @ingroup matrix
  !> 行列を対称行列に変換（内部自由度）
  subroutine monolis_matrix_convert_to_symmetric_inner_R(mat)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_mat), intent(inout) :: mat
    integer(kint) :: i, j, in, jn, n_dof, n_dof_2, jS, jE
    real(kdouble), allocatable :: A1(:,:), A2(:,:), A3(:,:)

    n_dof = mat%NDOF
    n_dof_2 = n_dof*n_dof

    call monolis_alloc_R_2d(A1, n_dof, n_dof)
    call monolis_alloc_R_2d(A2, n_dof, n_dof)
    call monolis_alloc_R_2d(A3, n_dof, n_dof)

    do i = 1, mat%N
      jS = mat%CSR%index(i) + 1
      jE = mat%CSR%index(i + 1)
      aa:do j = jS, jE
        in = mat%CSR%item(j)
        A1 = 0.0d0
        if(i < in)then
          cycle aa
        elseif(i == in)then
          call get_block_matrix(mat, n_dof, j, A1)
          call get_sym_matrix(n_dof, A1, transpose(A1), A2)
          call set_block_matrix(mat, n_dof, j, A2)
        elseif(i > in)then
          jn = mat%CSC%perm(j)
          call get_block_matrix(mat, n_dof, j, A1)
          call get_block_matrix(mat, n_dof, jn, A2)
          call get_sym_matrix(n_dof, A1, transpose(A2), A3)
          call set_block_matrix(mat, n_dof, j, A3)
          call set_block_matrix(mat, n_dof, jn, transpose(A3))
        endif
      enddo aa
    enddo
  end subroutine monolis_matrix_convert_to_symmetric_inner_R

  subroutine get_block_matrix(mat, n_dof, idx, A1)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_mat), intent(inout) :: mat
    !> [in] 行番号
    integer(kint) :: n_dof
    !> [in] 行番号
    integer(kint) :: idx
    !> [in] 行番号
    real(kdouble) :: A1(:,:)
    integer(kint) :: i, j, n_dof_2

    n_dof_2 = n_dof*n_dof

    do i = 1, n_dof
      do j = 1, n_dof
        A1(i,j) = mat%R%A(n_dof_2*(idx - 1) + n_dof*(i - 1) + j)
      enddo
    enddo
  end subroutine get_block_matrix

  subroutine get_sym_matrix(n_dof, A1, A2, A3)
    implicit none
    !> [in] 行番号
    integer(kint) :: n_dof
    !> [in] 行番号
    real(kdouble) :: A1(:,:)
    !> [in] 行番号
    real(kdouble) :: A2(:,:)
    !> [in] 行番号
    real(kdouble) :: A3(:,:)
    integer(kint) :: i, j

    do i = 1, n_dof
      do j = 1, n_dof
        A3(j,i) = 0.5d0*(A1(j,i) + A2(j,i))
      enddo
    enddo
  end subroutine get_sym_matrix

  subroutine set_block_matrix(mat, n_dof, idx, A1)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_mat), intent(inout) :: mat
    !> [in] 行番号
    integer(kint) :: n_dof
    !> [in] 行番号
    integer(kint) :: idx
    !> [in] 行番号
    real(kdouble) :: A1(:,:)
    integer(kint) :: i, j, n_dof_2

    n_dof_2 = n_dof*n_dof

    do i = 1, n_dof
      do j = 1, n_dof
        mat%R%A(n_dof_2*(idx - 1) + n_dof*(i - 1) + j) = A1(i,j)
      enddo
    enddo
  end subroutine set_block_matrix

  !> @ingroup matrix
  !> 行列を対称行列に変換（外部自由度、並列計算向け）
  subroutine monolis_matrix_convert_to_symmetric_outer_R(mat, com)
    implicit none
    !> [in,out] monolis 構造体
    type(monolis_mat), intent(inout) :: mat
    !> [in,out] monolis 構造体
    type(monolis_COM), intent(inout) :: com
    integer(kint) :: n_send, n_recv
    integer(kint) :: i, in, j, jn, jS, jE, k, kS, kE, kn, km, m, comm_size, n_dof
    integer(kint) :: mi, mj
    real(kdouble) :: tcomm
    integer(kint), allocatable :: vertex_id(:)
    integer(kint), allocatable :: n_send_list_all(:)
    integer(kint), allocatable :: n_recv_list_all(:)
    integer(kint), allocatable :: send_n_dof(:), send_n_dof_tmp(:)
    integer(kint), allocatable :: recv_n_dof(:), recv_n_dof_tmp(:)
    integer(kint), allocatable :: send_id(:), recv_id(:)
    real(kdouble), allocatable :: send_val(:), recv_val(:)

    !> 先にグローバル ID を振り直しておく
    n_dof = MAT%NDOF
    comm_size = monolis_mpi_get_local_comm_size(com%comm)
    call monolis_alloc_I_1d(vertex_id, mat%NP)
    call monolis_generate_global_vertex_id(mat%N, mat%NP, vertex_id, com)

    !> RECV の情報を見て、どの領域に送るかがわかる
    call monolis_alloc_I_1d(n_send_list_all, comm_size)
    call monolis_alloc_I_1d(n_recv_list_all, comm_size)

    n_send = 0
    do i = 1, com%recv_n_neib
      in = com%recv_neib_pe(i) + 1
      jS = com%recv_index(i) + 1
      jE = com%recv_index(i + 1)
      !> recv で受ける列を検索
      do j = jS, jE
        jn = com%recv_item(j)
        kS = mat%CSC%index(jn) + 1
        kE = mat%CSC%index(jn + 1)
        do k = kS, kE
          kn = mat%CSC%item(k)
          if(kn < mat%N)then
            n_send_list_all(in) = n_send_list_all(in) + 1
            n_send = n_send + 1
          endif
        enddo
      enddo
    enddo

    !> alltoall 交換
    n_recv_list_all = n_send_list_all
    call monolis_alltoall_I1(comm_size, n_recv_list_all, com%comm)

    !> recv 配列を作成
    n_recv = 0
    do i = 1, comm_size
      n_recv = n_recv + n_recv_list_all(i)
    enddo
    call monolis_alloc_I_1d(recv_id,  2*n_recv)
    call monolis_alloc_R_1d(recv_val, n_recv*n_dof*n_dof)

    !> send 配列を作成
    !> 新しいグローバル ID の (i,j) のデータを送る
    !> 行列の (i,j) のデータを送る(ndof * ndof のブロック）
    call monolis_alloc_I_1d(send_id,  2*n_send)
    call monolis_alloc_R_1d(send_val, n_send*n_dof*n_dof)

    n_send = 0
    do i = 1, com%recv_n_neib
      jS = com%recv_index(i) + 1
      jE = com%recv_index(i + 1)
      !> send で送る行を検索
      do j = jS, jE
        jn = com%recv_item(j)
        kS = mat%CSC%index(jn) + 1
        kE = mat%CSC%index(jn + 1)
        do k = kS, kE
          kn = mat%CSC%item(k)
          km = mat%CSC%perm(kn)
          if(kn < mat%N)then
            n_send = n_send + 1
            send_id(2*n_send - 1) = vertex_id(jn)
            send_id(2*n_send    ) = vertex_id(kn)
            do m = 1, n_dof*n_dof
              send_val(n_dof*n_dof*(n_send - 1) + m) = mat%R%A(n_dof*n_dof*(km - 1) + m)
            enddo
          endif
        enddo
      enddo
    enddo

    !> 行列の値を並列に更新する
    call monolis_alloc_I_1d(send_n_dof, com%send_n_neib)
    in = 0
    do i = 1, comm_size
      if(n_send_list_all(i) /= 0)then
        in = in + 1
        send_n_dof(in) = n_send_list_all(i)
      endif
    enddo

    call monolis_alloc_I_1d(recv_n_dof, com%recv_n_neib)
    in = 0
    do i = 1, comm_size
      if(n_recv_list_all(i) /= 0)then
        in = in + 1
        recv_n_dof(in) = n_recv_list_all(i)
      endif
    enddo

    !> データ通信パート ID
    call monolis_SendRecv_V_I1(com, 2*send_n_dof, 2*recv_n_dof, send_id, recv_id)

    !> データ通信パート VAL
    call monolis_SendRecv_V_R1(com, n_dof*n_dof*send_n_dof, n_dof*n_dof*recv_n_dof, send_val, recv_val)

    !> データ更新パート
    do i = 1, n_recv
      mi = recv_id(2*i-1)
      mj = recv_id(2*i)

      jS = mat%CSR%index(i) + 1
      jE = mat%CSR%index(i + 1)
      aa:do j = jS, jE
        in = mat%CSR%item(j)
        A1 = 0.0d0
        if(i < in)then
          cycle aa
        elseif(i == in)then
          call get_block_matrix(mat, n_dof, j, A1)
          call get_sym_matrix(n_dof, A1, transpose(A1), A2)
          call set_block_matrix(mat, n_dof, j, A2)
        elseif(i > in)then
          jn = mat%CSC%perm(j)
          call get_block_matrix(mat, n_dof, j, A1)
          call get_block_matrix(mat, n_dof, jn, A2)
          call get_sym_matrix(n_dof, A1, transpose(A2), A3)
          call set_block_matrix(mat, n_dof, j, A3)
          call set_block_matrix(mat, n_dof, jn, transpose(A3))
        endif
      enddo aa
    enddo

  end subroutine monolis_matrix_convert_to_symmetric_outer_R

end module mod_monolis_spmat_convert_sym


