interface
  subroutine monolis(N, NP, NDOF, NPU, NPL, D, AU, AL, X, B, &
    & indexU, itemU, indexL, itemL, myrank, comm, commsize, n_neib, &
    & neib_pe, recv_index, send_index, recv_item, send_item, &
    & method, precond, maxiter, tol, is_scaling)
    implicit none
    !> for monoMAT
    integer, parameter :: kint = 4
    integer, parameter :: kdouble = 8
    integer(kind=kint), intent(in) :: N, NP, NDOF, NPU, NPL
    integer(kind=kint), intent(in), pointer :: indexU(:)
    integer(kind=kint), intent(in), pointer :: indexL(:)
    integer(kind=kint), intent(in), pointer :: itemU(:)
    integer(kind=kint), intent(in), pointer :: itemL(:)
    real(kind=kdouble), intent(in), pointer :: D(:)
    real(kind=kdouble), intent(in), pointer :: AU(:)
    real(kind=kdouble), intent(in), pointer :: AL(:)
    real(kind=kdouble), intent(in), pointer :: B(:)
    real(kind=kdouble), intent(out),pointer :: X(:)
    !> for monoCOM
    integer(kind=kint), intent(in)          :: myrank
    integer(kind=kint), intent(in)          :: comm
    integer(kind=kint), intent(in)          :: commsize
    integer(kind=kint), intent(in)          :: n_neib
    integer(kind=kint), intent(in), pointer :: neib_pe(:)
    integer(kind=kint), intent(in), pointer :: recv_index(:)
    integer(kind=kint), intent(in), pointer :: recv_item(:)
    integer(kind=kint), intent(in), pointer :: send_index(:)
    integer(kind=kint), intent(in), pointer :: send_item(:)
    !> for monoPRM
    integer(kind=kint), intent(in) :: method
    integer(kind=kint), intent(in) :: precond
    integer(kind=kint), intent(in) :: maxiter
    real(kind=kdouble), intent(in) :: tol
    logical, intent(in) :: is_scaling
  end subroutine monolis
end interface

interface
  subroutine monolis_serial(N, NDOF, NPU, NPL, D, AU, AL, X, B, &
    & indexU, itemU, indexL, itemL, &
    & method, precond, maxiter, tol, is_scaling)
    use mod_monolis_prm
    use mod_monolis_com
    use mod_monolis_mat
    use mod_monolis_solve
    implicit none
    type(monolis_prm), save :: monoPRM
    type(monolis_com), save :: monoCOM
    type(monolis_mat), save :: monoMAT
    !> for monoMAT
    integer(kind=kint), intent(in) :: N, NDOF, NPU, NPL
    integer(kind=kint), intent(in), pointer :: indexU(:)
    integer(kind=kint), intent(in), pointer :: indexL(:)
    integer(kind=kint), intent(in), pointer :: itemU(:)
    integer(kind=kint), intent(in), pointer :: itemL(:)
    real(kind=kdouble), intent(in), pointer :: D(:)
    real(kind=kdouble), intent(in), pointer :: AU(:)
    real(kind=kdouble), intent(in), pointer :: AL(:)
    real(kind=kdouble), intent(in), pointer :: B(:)
    real(kind=kdouble), intent(out),pointer :: X(:)
    !> for monoPRM
    integer(kind=kint), intent(in) :: method
    integer(kind=kint), intent(in) :: precond
    integer(kind=kint), intent(in) :: maxiter
    real(kind=kdouble), intent(in) :: tol
    logical, intent(in) :: is_scaling
  end subroutine monolis_wrapper
end interface

interface
  subroutine monolis_convert_full_matrix(Nf, NDOFf, Af, thresh, &
    & N, NDOF, NPU, NPL, D, AU, AL, indexU, indexL, itemU, itemL)
    use mod_monolis_prm
    use mod_monolis_convert
    implicit none
    real(kind=kdouble), pointer :: Af(:)
    real(kind=kdouble), pointer :: D(:), AU(:), AL(:)
    integer(kind=kint), pointer :: indexU(:)
    integer(kind=kint), pointer :: indexL(:)
    integer(kind=kint), pointer :: itemU(:)
    integer(kind=kint), pointer :: itemL(:)
    integer(kind=kint) :: Nf, NDOFf
    integer(kind=kint) :: N, NDOF, NPU, NPL
    integer(kind=kint) :: i, j, k, jS, jE, in
    real(kind=kdouble) :: thresh

    call monolis_convert_full_matrix_main(Nf, NDOFf, Af, thresh, &
    & N, NDOF, NPU, NPL, D, AU, AL, indexU, indexL, itemU, itemL)
  end subroutine monolis_convert_full_matrix
end interface

interface
  subroutine monolis_convert_coo_matrix(Nf, NZf, NDOFf, Af, indexI, indexJ, &
    & N, NDOF, NPU, NPL, D, AU, AL, indexU, indexL, itemU, itemL)
    use mod_monolis_prm
    use mod_monolis_convert
    implicit none
    real(kind=kdouble), pointer :: Af(:)
    real(kind=kdouble), pointer :: D(:), AU(:), AL(:)
    integer(kind=kint), pointer :: indexI(:)
    integer(kind=kint), pointer :: indexJ(:)
    integer(kind=kint), pointer :: indexU(:)
    integer(kind=kint), pointer :: indexL(:)
    integer(kind=kint), pointer :: itemU(:)
    integer(kind=kint), pointer :: itemL(:)
    integer(kind=kint) :: Nf, NZf, NDOFf
    integer(kind=kint) :: N, NDOF, NPU, NPL
    integer(kind=kint) :: i, j, k, jS, jE, in

    call  monolis_convert_coo_matrix_main(Nf, NZf, NDOFf, Af, indexI, indexJ, &
    & N, NDOF, NPU, NPL, D, AU, AL, indexU, indexL, itemU, itemL)
  end subroutine monolis_convert_coo_matrix
end interface

interface
  subroutine monolis_convert_csr_matrix(Nf, NDOFf, Af, index, item, &
    & N, NDOF, NPU, NPL, D, AU, AL, indexU, indexL, itemU, itemL)
    use mod_monolis_prm
    use mod_monolis_convert
    implicit none
    real(kind=kdouble), pointer :: Af(:)
    real(kind=kdouble), pointer :: D(:), AU(:), AL(:)
    integer(kind=kint), pointer :: indexU(:)
    integer(kind=kint), pointer :: indexL(:)
    integer(kind=kint), pointer :: index(:)
    integer(kind=kint), pointer :: itemU(:)
    integer(kind=kint), pointer :: itemL(:)
    integer(kind=kint), pointer :: item(:)
    integer(kind=kint) :: Nf, NDOFf
    integer(kind=kint) :: N, NDOF, NPU, NPL
    integer(kind=kint) :: i, j, k, jS, jE, in

    call monolis_convert_csr_matrix_main(Nf, NDOFf, Af, index, item, &
    & N, NDOF, NPU, NPL, D, AU, AL, indexU, indexL, itemU, itemL)
  end subroutine monolis_convert_csr_matrix
end interface
