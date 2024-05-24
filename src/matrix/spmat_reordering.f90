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
    integer(kint) :: N, i
    integer(kint), allocatable :: index(:)
    integer(kint), allocatable :: item(:)

    !if(monoPRM%is_debug) call monolis_debug_header("monolis_matrix_reordering_fw_R")
    t1 = monolis_get_time()

    call monolis_get_nodal_graph_by_nonzero_pattern(monoMAT, N, index, item)

    call monolis_palloc_I_1d(monoMAT%REORDER%perm, N)
    call monolis_palloc_I_1d(monoMAT%REORDER%iperm, N)

    call gedatsu_part_graph_metis_reordering(N, index, item, monoMAT%REORDER%perm, monoMAT%REORDER%iperm)

    !call monolis_palloc_I_1d(monoMAT_reorder%REORDER%perm, N)
    !call monolis_palloc_I_1d(monoMAT_reorder%REORDER%iperm, N)
    !monoMAT_reorder%REORDER%perm  = monoMAT%REORDER%perm
    !monoMAT_reorder%REORDER%iperm = monoMAT%REORDER%iperm

    call monolis_restruct_matrix(monoMAT, monoMAT_reorder, monoMAT%REORDER%perm, monoMAT%REORDER%iperm)

!write(*,*)"monoMAT%N", monoMAT%N, monoMAT_reorder%N
!write(*,*)"monoMAT%NP", monoMAT%NP, monoMAT_reorder%NP
!write(*,*)"monoMAT%NDOF", monoMAT%NDOF, monoMAT_reorder%NDOF
!write(*,*)"monoMAT%CSR%index A", monoMAT%CSR%index
!write(*,*)"monoMAT%CSR%index B", monoMAT_reorder%CSR%index
!write(*,*)"monoMAT%CSR%item A", monoMAT%CSR%item
!write(*,*)"monoMAT%CSR%item B", monoMAT_reorder%CSR%item

    !call monolis_copy_mat_nonzero_pattern_main_R(monoMAT, monoMAT_reorder)
    !call monolis_copy_mat_value_main_R(monoMAT, monoMAT_reorder)

    t2 = monolis_get_time()
    !monoPRM%tprep = monoPRM%tprep + t2 - t1
  end subroutine monolis_matrix_reordering_fw_R

  subroutine monolis_matrix_reordering_bk_R(monoMAT)
    implicit none
    type(monolis_mat) :: monoMAT
    real(kdouble) :: t1, t2

    !if(monoPRM%is_debug) call monolis_debug_header("monolis_reorder_matrix_bk")
    t1 = monolis_get_time()

    !call monolis_reorder_back_vector_bk(monoMAT, monoMAT%NP, monoMAT%NDOF, monoMAT_reorder%X, monoMAT%X)

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

  subroutine monolis_restruct_matrix(monoMAT, monoMAT_reorder, perm, iperm)
    implicit none
    type(monolis_mat) :: monoMAT
    type(monolis_mat) :: monoMAT_reorder
    integer(kint) :: perm(:), iperm(:)
    integer(kint) :: N, NP, NZ, NDOF, NDOF2

    N = monoMAT%N
    NP = monoMAT%NP
    NZ = monoMAT%CSR%index(NP + 1)
    NDOF = monoMAT%NDOF
    NDOF2 = NDOF*NDOF

    monoMAT_reorder%N = N
    monoMAT_reorder%NP = NP
    monoMAT_reorder%NDOF = NDOF

    allocate(monoMAT_reorder%CSR%index(NP + 1), source = 0)
    allocate(monoMAT_reorder%CSR%item(NZ), source = 0)

    call monolis_restruct_matrix_profile(NP, perm, iperm, &
       & monoMAT%CSR%index, monoMAT%CSR%item, monoMAT_reorder%CSR%index, monoMAT_reorder%CSR%item)

    allocate(monoMAT_reorder%R%A(NDOF2*NZ))
    call monolis_restruct_matrix_values(NP, NDOF, perm, iperm, &
       & monoMAT%CSR%index, monoMAT%CSR%item, monoMAT%R%A, &
       & monoMAT_reorder%CSR%index, monoMAT_reorder%CSR%item, monoMAT_reorder%R%A)

    allocate(monoMAT_reorder%R%X(NDOF*NP))
    allocate(monoMAT_reorder%R%B(NDOF*NP))
  end subroutine monolis_restruct_matrix

  subroutine monolis_restruct_matrix_profile(N, perm, iperm, &
    & index, item, indexp, itemp)
    implicit none
    integer(kint) :: N
    integer(kint) :: perm(:), iperm(:)
    integer(kint) :: index(:), item(:)
    integer(kint) :: indexp(:), itemp(:)
    integer(kint) :: cnt, i, in, j, jo, jn

    cnt = 0
    indexp(1) = 0
    do i = 1, N
      in = perm(i)
      do j = index(in) + 1, index(in + 1)
        jo = item(j)
        jn = iperm(jo)
        cnt = cnt + 1
        itemp(cnt) = jn
      enddo
      indexp(i + 1) = cnt
      call monolis_qsort_I_1d(itemp, indexp(i) + 1, indexp(i + 1))
    enddo
  end subroutine monolis_restruct_matrix_profile

  subroutine monolis_restruct_matrix_values(N, NDOF, perm, iperm, index, item, A, &
      & indexp, itemp, Ap)
    implicit none
    integer(kint) :: N, NDOF
    integer(kint) :: perm(:), iperm(:)
    integer(kint) :: index(:), item(:)
    real(kdouble) :: A(:)
    integer(kint) :: indexp(:), itemp(:)
    real(kdouble) :: Ap(:)
    integer(kint) :: NDOF2, in, i
    integer(kint) :: jSn, jEn
    integer(kint) :: jo, ko, kn, jn, lo, ln, l

    NDOF2 = NDOF*NDOF
    do i = 1, N
      in = iperm(i)
      jSn = indexp(in) + 1
      jEn = indexp(in + 1)
      do jo = index(i) + 1, index(i + 1)
        ko = item(jo)
        kn = iperm(ko)
        call monolis_bsearch_I(itemp, jSn, jEn, kn, jn)
        lo = (jo - 1)*NDOF2
        ln = (jn - 1)*NDOF2
        do l = 1, NDOF2
          Ap(ln + l) = A(lo + l)
        enddo
      enddo
    enddo
  end subroutine monolis_restruct_matrix_values
end module mod_monolis_spmat_reorder