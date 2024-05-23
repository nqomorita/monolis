module mod_monolis_spmat_reorder
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc
  implicit none

contains

  subroutine monolis_matrix_reordering_fw_R(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_com) :: monoCOM_temp
    type(monolis_mat) :: monoMAT
    type(monolis_mat) :: monoMAT_temp
    integer(kint), pointer ::  perm(:), iperm(:)
    real(kdouble) :: t1, t2

    !if(monoPRM%is_debug) call monolis_debug_header("monolis_matrix_reordering_fw_R")
    t1 = monolis_get_time()

    !if(monoPRM%is_reordering)then
      !call monolis_reorder_matrix_metis(monoMAT, monoMAT_reorder)
      !call monolis_restruct_matrix(monoMAT, monoMAT_reorder, perm, iperm)
      !call monolis_restruct_comm(monoCOM, monoCOM_reorder, iperm)
      !call monolis_reorder_vector_fw(monoMAT, monoMAT%NP, monoMAT%NDOF, monoMAT%B, monoMAT_reorder%B)
      !call monolis_reorder_vector_fw(monoMAT, monoMAT%NP, monoMAT%NDOF, monoMAT%X, monoMAT_reorder%X)
    !endif

    t2 = monolis_get_time()
    !monoPRM%tprep = monoPRM%tprep + t2 - t1
  end subroutine monolis_matrix_reordering_fw_R

  subroutine monolis_matrix_reordering_bk_R(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    type(monolis_mat) :: monoMAT_reorder
    real(kdouble) :: t1, t2

    !if(monoPRM%is_debug) call monolis_debug_header("monolis_reorder_matrix_bk")
    t1 = monolis_get_time()

    !if(monoPRM%is_reordering)then
      !call monolis_reorder_back_vector_bk(monoMAT, monoMAT%NP, monoMAT%NDOF, monoMAT_reorder%X, monoMAT%X)
    !endif

    t2 = monolis_get_time()
    !monoPRM%tprep = monoPRM%tprep + t2 - t1
  end subroutine monolis_matrix_reordering_bk_R

  subroutine monolis_reorder_vector_fw(monoMAT, N, NDOF, A, B)
    implicit none
    type(monolis_mat) :: monoMAT
    integer(kint) :: N, NDOF
    real(kdouble) :: A(:)
    real(kdouble) :: B(:)
    integer(kint) :: i, in, jn, jo, j
    !do i = 1, N
    !  in = monoMAT%iperm(i)
    !  jn = (in-1)*NDOF
    !  jo = (i -1)*NDOF
    !  do j = 1, NDOF
    !    B(jn + j) = A(jo + j)
    !  enddo
    !enddo
  end subroutine monolis_reorder_vector_fw

  subroutine monolis_reorder_back_vector_bk(monoMAT, N, NDOF, B, A)
    implicit none
    type(monolis_mat) :: monoMAT
    integer(kint) :: N, NDOF
    real(kdouble) :: B(:)
    real(kdouble) :: A(:)
    integer(kint) :: i, in, jn, jo, j
    !do i = 1, N
    !  in = monoMAT%perm(i)
    !  jn = (i -1)*NDOF
    !  jo = (in-1)*NDOF
    !  do j = 1, NDOF
    !    A(jo + j) = B(jn + j)
    !  enddo
    !enddo
  end subroutine monolis_reorder_back_vector_bk

end module mod_monolis_spmat_reorder