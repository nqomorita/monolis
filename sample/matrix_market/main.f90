program main
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_util
  use mod_monolis_solve
  use mod_monolis_convert
  implicit none
  type(monolis_prm) :: monoPRM
  type(monolis_com) :: monoCOM
  type(monolis_mat) :: monoMAT
  integer(kind=kint) :: i, j, N, NDOF, NPU, NPL, Nf, NZf, NDOFf
  integer(kind=kint), pointer :: indexI(:)
  integer(kind=kint), pointer :: indexJ(:)
  integer(kind=kint), pointer :: indexU(:)
  integer(kind=kint), pointer :: indexL(:)
  integer(kind=kint), pointer :: itemU(:)
  integer(kind=kint), pointer :: itemL(:)
  real(kind=kdouble), pointer :: Af(:)
  real(kind=kdouble), pointer :: D(:), AU(:), AL(:)
  character :: ctemp*10

  !> reference: matrix market
  !> https://math.nist.gov/MatrixMarket/data/Harwell-Boeing/bcsstruc2/bcsstk14.html
  open(20, file="./bcsstk14.mtx", status="old")
    read(20,*)ctemp
    read(20,*)Nf, Nf, NZf
    allocate(indexI(NZf))
    allocate(indexJ(NZf))
    allocate(Af(NZf))
    do i = 1, NZf
      read(20,*)indexI(i), indexJ(i), Af(i)
    enddo
  close(20)

  NDOFf = 1
  call monolis_convert_coo_matrix_main(Nf, NZf, NDOFf, Af, indexI, indexJ, &
    & N, NDOF, NPU, NPL, D, AU, AL, indexU, indexL, itemU, itemL)
end program main
