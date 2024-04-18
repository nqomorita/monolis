
subroutine monolis_ML_get_nlocal(id, nlocal, nlocal_allcolumns, ierr)
  use mod_monolis_utils
  implicit none
  integer(kint), intent(in) :: id
  integer(kint), intent(out) :: nlocal
  integer(kint), intent(out) :: nlocal_allcolumns
  integer(kint), intent(out) :: ierr
  !type(monolisST_matrix), pointer :: hecMAT
  
  !call monolis_mat_id_get(id, hecMAT, hecMESH)
  !nlocal = hecMAT%N * hecMAT%NDOF
  !nlocal_allcolumns = hecMAT%NP * hecMAT%NDOF
  !ierr = 0
end subroutine monolis_ML_get_nlocal

subroutine monolis_ML_get_coord(id, x, y, z, ierr)
  use mod_monolis_utils
  implicit none
  integer(kint), intent(in) :: id
  real(kdouble), intent(out) :: x(*), y(*), z(*)
  integer(kint), intent(out) :: ierr
  !type(monolisST_matrix), pointer :: hecMAT
  integer(kint) :: offset, i

  !call monolis_mat_id_get(id, hecMAT, hecMESH)
  !offset = 0
  !do i = 1, hecMESH%nn_internal
  !  x(i) = hecMESH%node(offset+1)
  !  y(i) = hecMESH%node(offset+2)
  !  z(i) = hecMESH%node(offset+3)
  !  offset = offset + 3
  !enddo
  !ierr = 0
end subroutine monolis_ML_get_coord

subroutine monolis_ML_get_loglevel(id, level)
  use mod_monolis_utils
  implicit none
  integer(kint), intent(in) :: id
  integer(kint), intent(out) :: level
  !type(monolisST_matrix), pointer :: hecMAT

  !call monolis_mat_id_get(id, hecMAT, hecMESH)
  !level = monolis_mat_get_timelog(hecMAT)
end subroutine monolis_ML_get_loglevel

subroutine monolis_ML_get_opt(id, opt, ierr)
  use mod_monolis_utils
  implicit none
  integer(kint), intent(in) :: id
  integer(kint), intent(out) :: opt(*)
  integer(kint), intent(out) :: ierr
  !type(monolisST_matrix), pointer :: hecMAT
  integer(kint) :: iopt(10)

  !call monolis_mat_id_get(id, hecMAT, hecMESH)
  !call monolis_mat_get_solver_opt(hecMAT, iopt)
  !opt(1:10) = iopt(1:10)
  !ierr = 0
end subroutine monolis_ML_get_opt

subroutine monolis_ML_set_opt(id, opt, ierr)
  use mod_monolis_utils
  implicit none
  integer(kint), intent(in) :: id
  integer(kint), intent(in) :: opt(*)
  integer(kint), intent(out) :: ierr
  !type(monolisST_matrix), pointer :: hecMAT
  integer(kint) :: iopt(10)
  
  !call monolis_mat_id_get(id, hecMAT, hecMESH)
  !iopt(1:10) = opt(1:10)
  !call monolis_mat_set_solver_opt(hecMAT, iopt)
  !ierr = 0
end subroutine monolis_ML_set_opt

subroutine monolis_ML_getrow_nn(id, n_requested_rows, requested_rows, &
    allocated_space, cols, values, row_lengths, ierr)
  use mod_monolis_utils
  implicit none
  integer(kint), intent(in) :: id
  integer(kint), intent(in) :: n_requested_rows
  integer(kint), intent(in) :: requested_rows(n_requested_rows)
  integer(kint), intent(in) :: allocated_space
  integer(kint), intent(out) :: cols(allocated_space)
  real(kdouble), intent(out) :: values(allocated_space)
  integer(kint), intent(out) :: row_lengths(n_requested_rows)
  integer(kint), intent(out) :: ierr
  !type(monolisST_matrix), pointer :: hecMAT
  integer(kint) :: m, i, row, inod, idof, nl, nd, nu, js, je, j, jj, jdof, start, ndof

!  !call monolis_mat_id_get(id, hecMAT, hecMESH)
!  ndof = hecMAT%NDOF
!  m = 1
!  do i = 1, n_requested_rows
!    row = requested_rows(i) + 1 ! '+1' for Fortran-numbering
!    inod = (row-1)/ndof + 1
!    idof = row - (inod-1)*ndof
!    nl = (hecMAT%indexL(inod) - hecMAT%indexL(inod-1)) * ndof
!    nd = ndof
!    nu = (hecMAT%indexU(inod) - hecMAT%indexU(inod-1)) * ndof
!    if (allocated_space < m + nl + nd + nu) return
!    start = m
!    js = hecMAT%indexL(inod-1)+1
!    je = hecMAT%indexL(inod)
!    do j = js, je
!      jj = hecMAT%itemL(j)
!      do jdof = 1, ndof
!        cols(m) = (jj-1)*ndof + jdof - 1 ! '-1' for C-numbering
!        values(m) = hecMAT%AL((j-1)*ndof*ndof + (idof-1)*ndof + jdof)
!        m = m+1
!      enddo
!    enddo
!    do jdof = 1, ndof
!      cols(m) = (inod-1)*ndof + jdof - 1 ! '-1' for C-numbering
!      values(m) = hecMAT%D((inod-1)*ndof*ndof + (idof-1)*ndof + jdof)
!      m = m+1
!    enddo
!    js = hecMAT%indexU(inod-1)+1
!    je = hecMAT%indexU(inod)
!    do j = js, je
!      jj = hecMAT%itemU(j)
!      do jdof = 1, ndof
!        cols(m) = (jj-1)*ndof + jdof - 1 ! '-1' for C-numbering
!        values(m) = hecMAT%AU((j-1)*ndof*ndof + (idof-1)*ndof + jdof)
!        m = m+1
!      enddo
!    enddo
!    row_lengths(i) = m - start
!  enddo
!  ierr = 1
end subroutine monolis_ML_getrow_nn

subroutine monolis_ML_matvec_nn(id, in_length, p, out_length, ap, ierr)
  use mod_monolis_utils
  implicit none
  integer(kint), intent(in) :: id
  integer(kint), intent(in) :: in_length
  real(kdouble), intent(in) :: p(in_length)
  integer(kint), intent(in) :: out_length
  real(kdouble), intent(out) :: ap(out_length)
  integer(kint), intent(out) :: ierr
  !type(monolisST_matrix), pointer :: hecMAT
  real(kdouble), allocatable :: w(:)
  integer(kint) :: i
  !call monolis_mat_id_get(id, hecMAT, hecMESH)
  !allocate(w(hecMAT%NP*hecMAT%NDOF))
  !do i = 1, hecMAT%N*hecMAT%NDOF
  !  w(i) = p(i)
  !enddo
  !call monolis_matvec(hecMESH, hecMAT, w, ap)
  !deallocate(w)
  !ierr = 0
end subroutine monolis_ML_matvec_nn

subroutine monolis_ML_comm_nn(id, x, ierr)
  use mod_monolis_utils
  implicit none
  integer(kint), intent(in) :: id
  real(kdouble), intent(inout) :: x(*)
  integer(kint), intent(out) :: ierr
  !type(monolisST_matrix), pointer :: hecMAT
  !call monolis_mat_id_get(id, hecMAT, hecMESH)
  !call monolis_update_R (hecMESH, x, hecMAT%NP, hecMAT%NDOF)
  ierr = 0
end subroutine monolis_ML_comm_nn

