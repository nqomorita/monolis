module mod_monolis_convert_coo
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat

  implicit none

contains

  subroutine monolis_convert_coo_get_size(N, NZ, indexI, indexJ, NPU, NPL)
    implicit none
    integer(kind=kint), pointer :: indexI(:)
    integer(kind=kint), pointer :: indexJ(:)
    integer(kind=kint) :: N, NZ
    integer(kind=kint) :: NPU, NPL
    integer(kind=kint) :: i, ni, nj

    if(N < 1 .or. NZ < 1 &
      & .or. (.not. associated(indexI)) .or. (.not. associated(indexJ)))then
      stop "  ** monolis error: monolis_convert_coo_get_size"
    endif

    NPL = 0
    NPU = 0
    do i = 1, NZ
      ni = indexI(i)
      nj = indexJ(i)
      if(nj < ni)then
        NPL = NPL + 1
      endif
      if(ni < nj)then
        NPU = NPU + 1
      endif
    enddo
  end subroutine monolis_convert_coo_get_size

  subroutine monolis_convert_coo_get_size_c(N_c, NZ_c, indexI_c, indexJ_c, &
    & NPU_c, NPL_c) bind(c, name="monolis_convert_coo_get_size")
    use iso_c_binding
    implicit none
    integer(c_int), target :: indexI_c(NZ_c)
    integer(c_int), target :: indexJ_c(NZ_c)
    integer(c_int), value  :: N_c, NZ_c
    integer(c_int) :: NPU_c, NPL_c
    !> for fortran
    integer(kind=kint), pointer :: indexI(:) => null()
    integer(kind=kint), pointer :: indexJ(:) => null()
    integer(kind=kint) :: N, NZ, NDOF, NPU, NPL, NDOF2

    indexI => indexI_c
    indexJ => indexJ_c
    N = N_c
    NZ = NZ_c

    call monolis_convert_coo_get_size(N, NZ, indexI, indexJ, NPU, NPL)

    NPU_c = NPU
    NPL_c = NPL
  end subroutine monolis_convert_coo_get_size_c

  subroutine monolis_convert_coo_get_index(N, NZ, indexI, indexJ, NPU, NPL, indexU, itemU, indexL, itemL)
    implicit none
    integer(kind=kint), pointer :: indexI(:)
    integer(kind=kint), pointer :: indexJ(:)
    integer(kind=kint), pointer :: indexU(:)
    integer(kind=kint), pointer :: indexL(:)
    integer(kind=kint), pointer :: itemU(:)
    integer(kind=kint), pointer :: itemL(:)
    integer(kind=kint) :: N, NZ, NDOF, NPU, NPL
    integer(kind=kint) :: i, j, k, jS, jE, in, ni, nj, id, iu, il

    if(N < 1 .or. NZ < 1 &
      & .or. (.not. associated(indexI)) .or. (.not. associated(indexJ)))then
      stop "  ** monolis error: monolis_convert_coo_get_index"
    endif

    indexU = 0
    indexL = 0

    do i = 1, NZ
      ni = indexI(i)
      nj = indexJ(i)
      if(nj < ni)then
        indexL(ni) = indexL(ni) + 1
      endif
      if(ni < nj)then
        indexU(ni) = indexU(ni) + 1
      endif
    enddo

    do i = 1, N
      indexL(i) = indexL(i-1) + indexL(i)
      indexU(i) = indexU(i-1) + indexU(i)
    enddo

    itemL = 0
    itemU = 0

    il = 0
    id = 0
    iu = 0
    do i = 1, NZ
      ni = indexI(i)
      nj = indexJ(i)
      if(nj < ni)then
        il = il + 1
        itemL(il) = nj
      endif
      if(ni < nj)then
        iu = iu + 1
        itemU(iu) = nj
      endif
    enddo
  end subroutine monolis_convert_coo_get_index

  subroutine monolis_convert_coo_get_index_c(N_c, NZf_c, indexI_c, indexJ_c, NPU_c, NPL_c, &
    & indexU_c, itemU_c, indexL_c, itemL_c) bind(c, name="monolis_convert_coo_get_index")
    use iso_c_binding
    implicit none
    integer(c_int), target :: indexI_c(NZf_c)
    integer(c_int), target :: indexJ_c(NZf_c)
    integer(c_int), target :: indexU_c(0:N_c)
    integer(c_int), target :: indexL_c(0:N_c)
    integer(c_int), target :: itemU_c(NPU_c)
    integer(c_int), target :: itemL_c(NPL_c)
    integer(c_int), value  :: N_c, NZf_c, NPU_c, NPL_c
    !> for fortran
    integer(kind=kint), pointer :: indexI(:) => null()
    integer(kind=kint), pointer :: indexJ(:) => null()
    integer(kind=kint), pointer :: indexU(:) => null()
    integer(kind=kint), pointer :: indexL(:) => null()
    integer(kind=kint), pointer :: itemU(:) => null()
    integer(kind=kint), pointer :: itemL(:) => null()
    integer(kind=kint) :: N, NZ, NPU, NPL

    N = N_c
    NZ = NZf_c
    NPU = NPU_c
    NPL = NPL_c

    indexI => indexI_c
    indexJ => indexJ_c
    indexU => indexU_c
    indexL => indexL_c
    itemU => itemU_c
    itemL => itemL_c

    call monolis_convert_coo_get_index(N, NZ, indexI, indexJ, NPU, NPL, &
         & indexU, itemU, indexL, itemL)
  end subroutine monolis_convert_coo_get_index_c

  subroutine monolis_convert_coo_update_matrix_entry(N, NZ, NDOF, A, indexI, indexJ, NPU, NPL, &
    & D, AU, AL, indexU, itemU, indexL, itemL)
    implicit none
    real(kind=kdouble), pointer :: A(:)
    real(kind=kdouble), pointer :: D(:), AU(:), AL(:)
    integer(kind=kint), pointer :: indexI(:)
    integer(kind=kint), pointer :: indexJ(:)
    integer(kind=kint), pointer :: indexU(:)
    integer(kind=kint), pointer :: indexL(:)
    integer(kind=kint), pointer :: itemU(:)
    integer(kind=kint), pointer :: itemL(:)
    integer(kind=kint) :: N, NZ, NDOF, NPU, NPL, NDOF2
    integer(kind=kint) :: i, j, k, jS, jE, in, ni, nj, id, iu, il

    if(N < 1 .or. NZ < 1 .or. NDOF < 1 .or. (.not. associated(A)) &
      & .or. (.not. associated(indexI)) .or. (.not. associated(indexJ)))then
      stop "  ** monolis error: monolis_convert_coo_update_matrix_entry"
    endif

    NDOF2 = NDOF*NDOF

    D = 0.0d0
    AL = 0.0d0
    AU = 0.0d0

    il = 0
    id = 0
    iu = 0
    do i = 1, NZ
      ni = indexI(i)
      nj = indexJ(i)
      if(nj < ni)then
        il = il + 1
        do k = 1, NDOF2
          AL(NDOF2*(il-1) + k) = dabs(A(NDOF2*(i-1) + k))
        enddo
      endif
      if(ni == nj)then
        id = id + 1
        do k = 1, NDOF2
          D (NDOF2*(id-1) + k) = dabs(A(NDOF2*(i-1) + k))
        enddo
      endif
      if(ni < nj)then
        iu = iu + 1
        do k = 1, NDOF2
          AU(NDOF2*(iu-1) + k) =  dabs(A(NDOF2*(i-1) + k))
        enddo
      endif
    enddo
  end subroutine monolis_convert_coo_update_matrix_entry

  subroutine monolis_convert_coo_update_matrix_entry_c(N_c, NZf_c, NDOF_c, A_c, indexI_c, indexJ_c, NPU_c, NPL_c, &
    & D_c, AU_c, AL_c, indexU_c, itemU_c, indexL_c, itemL_c) bind(c, name="monolis_convert_coo_update_matrix_entry")
    use iso_c_binding
    implicit none
    real(c_double), target :: A_c(NZf_c)
    integer(c_int), target :: indexI_c(NZf_c)
    integer(c_int), target :: indexJ_c(NZf_c)
    real(c_double), target :: D_c (N_c  *NDOF_c*NDOF_c)
    real(c_double), target :: AU_c(NPU_c*NDOF_c*NDOF_c)
    real(c_double), target :: AL_c(NPL_c*NDOF_c*NDOF_c)
    integer(c_int), target :: indexU_c(0:N_c)
    integer(c_int), target :: indexL_c(0:N_c)
    integer(c_int), target :: itemU_c(NPU_c)
    integer(c_int), target :: itemL_c(NPL_c)
    integer(c_int), value  :: N_c, NZf_c, NDOF_c, NPU_c, NPL_c
    !> for fortran
    real(kind=kdouble), pointer :: A(:) => null()
    real(kind=kdouble), pointer :: D(:) => null()
    real(kind=kdouble), pointer :: AU(:) => null()
    real(kind=kdouble), pointer :: AL(:) => null()
    integer(kind=kint), pointer :: indexI(:) => null()
    integer(kind=kint), pointer :: indexJ(:) => null()
    integer(kind=kint), pointer :: indexU(:) => null()
    integer(kind=kint), pointer :: indexL(:) => null()
    integer(kind=kint), pointer :: itemU(:) => null()
    integer(kind=kint), pointer :: itemL(:) => null()
    integer(kind=kint) :: N, NZ, NDOF, NPU, NPL, NDOF2

    N = N_c
    NZ = NZf_c
    NDOF = NDOF_c
    NPU = NPU_c
    NPL = NPL_c

    A => A_c
    indexI => indexI_c
    indexJ => indexJ_c
    D  => D_c
    AU => AU_c
    AL => AL_c
    indexU => indexU_c
    indexL => indexL_c
    itemU => itemU_c
    itemL => itemL_c

    call monolis_convert_coo_update_matrix_entry(N, NZ, NDOF, A, indexI, indexJ, NPU, NPL, &
    & D, AU, AL, indexU, itemU, indexL, itemL)
  end subroutine monolis_convert_coo_update_matrix_entry_c
end module mod_monolis_convert_coo