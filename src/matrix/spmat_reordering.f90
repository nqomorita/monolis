module mod_monolis_spmat_reorder
  use mod_monolis_utils
  use mod_monolis_spmat_nonzero_pattern_util
  use mod_monolis_def_mat
  use mod_monolis_def_struc
  use mod_monolis_spmat_copy
  use mod_gedatsu
  implicit none

contains

  subroutine monolis_matrix_reordering_fw_R(monoMAT, monoMAT_reorder)
    implicit none
    type(monolis_mat) :: monoMAT
    type(monolis_mat) :: monoMAT_reorder
    real(kdouble) :: t1, t2
    integer(kint) :: N
    integer(kint), allocatable :: index(:)
    integer(kint), allocatable :: item(:)

    !if(monoPRM%is_debug) call monolis_debug_header("monolis_matrix_reordering_fw_R")
    t1 = monolis_get_time()

    call monolis_get_nodal_graph_by_nonzero_pattern(monoMAT, N, index, item)

    call monolis_palloc_I_1d(monoMAT%REORDER%perm, N)
    call monolis_palloc_I_1d(monoMAT%REORDER%iperm, N)

    call gedatsu_part_graph_metis_reordering(N, index, item, monoMAT%REORDER%perm, monoMAT%REORDER%iperm)
    !call monolis_restruct_matrix(monoMAT, monoMAT_reorder, monoMAT%REORDER%perm, monoMAT%REORDER%iperm)

    call monolis_copy_mat_nonzero_pattern_main_R(monoMAT, monoMAT_reorder)
    call monolis_copy_mat_value_main_R(monoMAT, monoMAT_reorder)

    t2 = monolis_get_time()
    !monoPRM%tprep = monoPRM%tprep + t2 - t1
  end subroutine monolis_matrix_reordering_fw_R

  subroutine monolis_matrix_reordering_bk_R(monoMAT)
    implicit none
    type(monolis_mat) :: monoMAT
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
    do i = 1, N
      in = monoMAT%REORDER%iperm(i)
      jn = (in-1)*NDOF
      jo = (i -1)*NDOF
      do j = 1, NDOF
        B(jn + j) = A(jo + j)
      enddo
    enddo
  end subroutine monolis_reorder_vector_fw

  subroutine monolis_reorder_back_vector_bk(monoMAT, N, NDOF, B, A)
    implicit none
    type(monolis_mat) :: monoMAT
    integer(kint) :: N, NDOF
    real(kdouble) :: B(:)
    real(kdouble) :: A(:)
    integer(kint) :: i, in, jn, jo, j
    do i = 1, N
      in = monoMAT%REORDER%perm(i)
      jn = (i -1)*NDOF
      jo = (in-1)*NDOF
      do j = 1, NDOF
        A(jo + j) = B(jn + j)
      enddo
    enddo
  end subroutine monolis_reorder_back_vector_bk

end module mod_monolis_spmat_reorder