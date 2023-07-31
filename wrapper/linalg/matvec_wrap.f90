!> 疎行列ベクトル積関数群
module mod_monolis_matvec_wrap
  use mod_monolis_utils
  use mod_monolis_def_struc
  use mod_monolis_matvec
  use iso_c_binding
  implicit none

contains

  !> @ingroup dev_linalg
  !> 疎行列ベクトル積（実数型）
  subroutine monolis_matvec_product_R_c(N, NP, NZ, NDOF, A, X, Y, index, item, &
    my_rank, comm, comm_size, &
    recv_n_neib, recv_nitem, recv_neib_pe, recv_index, recv_item, &
    send_n_neib, send_nitem, send_neib_pe, send_index, send_item) &
    & bind(c, name = "monolis_matvec_product_R_c_main")
    implicit none
    !> [in] 配列サイズ
    integer(c_int), intent(in), value :: N
    !> [in] 配列サイズ
    integer(c_int), intent(in), value :: NP
    !> [in] 配列サイズ
    integer(c_int), intent(in), value :: NZ
    !> [in] 自由度
    integer(c_int), intent(in), value :: NDOF
    !> [in] 自由度
    real(c_double), intent(in), target :: A(NDOF*NDOF*NZ)
    !> [in,out] 自由度
    real(c_double), intent(inout), target :: X(NDOF*NP)
    !> [out] 自由度
    real(c_double), intent(out), target :: Y(NDOF*NP)
    !> [in] 自由度
    integer(c_int), intent(in), target :: index(NP + 1)
    !> [in] 自由度
    integer(c_int), intent(in), target :: item(NZ)
    !> [in] MPI コミュニケータ
    integer(c_int), intent(in), value :: my_rank
    !> [in] MPI コミュニケータ
    integer(c_int), intent(in), value :: comm
    !> [in] MPI コミュニケータ
    integer(c_int), intent(in), value :: comm_size
    !> [in] recv する隣接領域数
    integer(c_int), intent(in), value :: recv_n_neib
    !> [in] recv の item 数
    integer(c_int), intent(in), value :: recv_nitem
    !> [in] recv する隣接領域 id
    integer(c_int), intent(in), target :: recv_neib_pe(recv_n_neib)
    !> [in] recv の index 配列
    integer(c_int), intent(in), target :: recv_index(recv_n_neib + 1)
    !> [in] recv の item 配列（受信する節点番号データ）
    integer(c_int), intent(in), target :: recv_item(recv_nitem)
    !> [in] send する隣接領域数
    integer(c_int), intent(in), value :: send_n_neib
    !> [in] send の item 数
    integer(c_int), intent(in), value :: send_nitem
    !> [in] send する隣接領域 id
    integer(c_int), intent(in), target :: send_neib_pe(send_n_neib)
    !> [in] send の index 配列
    integer(c_int), intent(in), target :: send_index(send_n_neib + 1)
    !> [in] send の item 配列（送信する節点番号データ）
    integer(c_int), intent(in), target :: send_item(send_nitem)
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    real(kdouble) :: tspmv, tcomm

    !> for monoMAT
    monoMAT%N = N
    monoMAT%NP = NP
    monoMAT%NDOF = NDOF
    monoMAT%R%A => A
    monoMAT%R%X => X
    monoMAT%R%B => Y
    monoMAT%CSR%index => index
    monoMAT%CSR%item => item

    !> for monoCOM
    monoCOM%my_rank = my_rank
    monoCOM%comm = comm
    monoCOM%comm_size = comm_size
    monoCOM%recv_n_neib = recv_n_neib
    monoCOM%recv_neib_pe => recv_neib_pe
    monoCOM%recv_index => recv_index
    monoCOM%recv_item => recv_item
    monoCOM%send_n_neib = send_n_neib
    monoCOM%send_neib_pe => send_neib_pe
    monoCOM%send_index => send_index
    monoCOM%send_item => send_item

    monoCOM%recv_item = monoCOM%recv_item + 1
    monoCOM%send_item = monoCOM%send_item + 1

    call monolis_matvec_product_main_R(monoCOM, monoMAT, X, Y, tspmv, tcomm)

    call monolis_mpi_update_R(monoCOM, monoMAT%NDOF, Y, tcomm)

    monoCOM%recv_item = monoCOM%recv_item - 1
    monoCOM%send_item = monoCOM%send_item - 1
  end subroutine monolis_matvec_product_R_c

  !> @ingroup dev_linalg
  !> 疎行列ベクトル積（複素数型）
  subroutine monolis_matvec_product_C_c(N, NP, NZ, NDOF, A, X, Y, index, item, &
    my_rank, comm, comm_size, &
    recv_n_neib, recv_nitem, recv_neib_pe, recv_index, recv_item, &
    send_n_neib, send_nitem, send_neib_pe, send_index, send_item) &
    & bind(c, name = "monolis_matvec_product_C_c_main")
    implicit none
    !> [in] 配列サイズ
    integer(c_int), intent(in), value :: N
    !> [in] 配列サイズ
    integer(c_int), intent(in), value :: NP
    !> [in] 配列サイズ
    integer(c_int), intent(in), value :: NZ
    !> [in] 自由度
    integer(c_int), intent(in), value :: NDOF
    !> [in] 自由度
    complex(c_double), intent(in), target :: A(NDOF*NDOF*NZ)
    !> [in,out] 自由度
    complex(c_double), intent(inout), target :: X(NDOF*NP)
    !> [out] 自由度
    complex(c_double), intent(out), target :: Y(NDOF*NP)
    !> [in] 自由度
    integer(c_int), intent(in), target :: index(NP + 1)
    !> [in] 自由度
    integer(c_int), intent(in), target :: item(NZ)
    !> [in] MPI コミュニケータ
    integer(c_int), intent(in), value :: my_rank
    !> [in] MPI コミュニケータ
    integer(c_int), intent(in), value :: comm
    !> [in] MPI コミュニケータ
    integer(c_int), intent(in), value :: comm_size
    !> [in] recv する隣接領域数
    integer(c_int), intent(in), value :: recv_n_neib
    !> [in] recv の item 数
    integer(c_int), intent(in), value :: recv_nitem
    !> [in] recv する隣接領域 id
    integer(c_int), intent(in), target :: recv_neib_pe(recv_n_neib)
    !> [in] recv の index 配列
    integer(c_int), intent(in), target :: recv_index(recv_n_neib + 1)
    !> [in] recv の item 配列（受信する節点番号データ）
    integer(c_int), intent(in), target :: recv_item(recv_nitem)
    !> [in] send する隣接領域数
    integer(c_int), intent(in), value :: send_n_neib
    !> [in] send の item 数
    integer(c_int), intent(in), value :: send_nitem
    !> [in] send する隣接領域 id
    integer(c_int), intent(in), target :: send_neib_pe(send_n_neib)
    !> [in] send の index 配列
    integer(c_int), intent(in), target :: send_index(send_n_neib + 1)
    !> [in] send の item 配列（送信する節点番号データ）
    integer(c_int), intent(in), target :: send_item(send_nitem)
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    real(kdouble) :: tspmv, tcomm

    !> for monoMAT
    monoMAT%N = N
    monoMAT%NP = NP
    monoMAT%NDOF = NDOF
    monoMAT%C%A => A
    monoMAT%C%X => X
    monoMAT%C%B => Y
    monoMAT%CSR%index => index
    monoMAT%CSR%item => item

    !> for monoCOM
    monoCOM%my_rank = my_rank
    monoCOM%comm = comm
    monoCOM%comm_size = comm_size
    monoCOM%recv_n_neib = recv_n_neib
    monoCOM%recv_neib_pe => recv_neib_pe
    monoCOM%recv_index => recv_index
    monoCOM%recv_item => recv_item
    monoCOM%send_n_neib = send_n_neib
    monoCOM%send_neib_pe => send_neib_pe
    monoCOM%send_index => send_index
    monoCOM%send_item => send_item

    monoCOM%recv_item = monoCOM%recv_item + 1
    monoCOM%send_item = monoCOM%send_item + 1

    call monolis_matvec_product_main_C(monoCOM, monoMAT, X, Y, tspmv, tcomm)

    call monolis_mpi_update_C(monoCOM, monoMAT%NDOF, Y, tcomm)

    monoCOM%recv_item = monoCOM%recv_item - 1
    monoCOM%send_item = monoCOM%send_item - 1
  end subroutine monolis_matvec_product_C_c
end module mod_monolis_matvec_wrap