module mod_monolis_precond_ROM
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_util

  implicit none

  private

  integer(kint), save :: n_get_eigen
  real(kdouble), save, allocatable :: val(:), vec(:,:)

  public :: monolis_precond_ROM_setup
  public :: monolis_precond_ROM_apply
  public :: monolis_precond_ROM_clear

contains

  subroutine monolis_precond_ROM_setup(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kint) :: maxiter, NPNDOF, n_bc
    real(kdouble) :: ths
    logical, allocatable :: is_bc(:)

    allocate(is_bc(monoMAT%NP*monoMAT%NDOF), source = .false.)
    call get_is_bc(monoPRM, monoCOM, monoMAT, is_bc, n_bc)
    !write(*,*)"n_bc", n_bc

    maxiter = monoMAT%NP*monoMAT%NDOF !min(100, monoMAT%NP*monoMAT%NDOF)
    n_get_eigen = monoMAT%NP*monoMAT%NDOF - n_bc !min(10, monoMAT%NP*monoMAT%NDOF)
    ths = 1.0d-4
    allocate(val(n_get_eigen), source = 0.0d0)
    allocate(vec(monoMAT%NP*monoMAT%NDOF, n_get_eigen), source = 0.0d0)

    monoPRM%precond = monolis_prec_SOR
    monoPRM%show_summary = .false.
    monoPRM%show_time = .false.
    monoPRM%show_iterlog = .false.

    NPNDOF = monoMAT%NP*monoMAT%NDOF
    call interface_monolis_eigen_inverted_standard_lanczos_ &
    & (monoPRM, monoCOM, monoMAT, NPNDOF, n_get_eigen, ths, maxiter, val, vec, is_bc)

    monoPRM%show_summary = .true.
    monoPRM%show_time = .true.
    monoPRM%show_iterlog = .true.

    monoPRM%precond = monolis_prec_ROM
  end subroutine monolis_precond_ROM_setup

  subroutine get_is_bc(monoPRM, monoCOM, monoMAT, is_bc, n_bc)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kint) :: i, j, jS, jE, k, l, in, jn, NDOF, NDOF2, n_bc
    integer(kint), allocatable :: flag(:)
    logical :: is_bc(:)

    NDOF  = monoMAT%NDOF
    NDOF2 = NDOF*NDOF
    allocate(flag(NDOF))

    n_bc = 0
    do i = 1, monoMAT%N
      jS = monoMAT%index(i-1) + 1
      jE = monoMAT%index(i  )
      flag = 0
      do j = jS, jE
        in = monoMAT%item(j)
        do k = 1, NDOF
          do l = 1, NDOF
            jn = NDOF2*(j-1)+NDOF*(k-1)+l
            if(monoMAT%A(jn) /= 0.0d0) flag(k) = flag(k) + 1
          enddo
        enddo
      enddo
      do j = 1, NDOF
        if(flag(j) == 1)then
          is_bc(NDOF*(i-1)+j) = .true.
          n_bc = n_bc + 1
        endif
      enddo
    enddo
  end subroutine get_is_bc

  subroutine monolis_precond_ROM_apply(monoPRM, monoCOM, monoMAT, X, Y)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kint) :: i, j
    real(kdouble) :: X(:), Y(:), phi
    real(kdouble), allocatable :: coef(:)

    allocate(coef(n_get_eigen), source = 0.0d0)

    do i = 1, n_get_eigen
      phi = 0.0d0
      do j = 1, monoMAT%NP*monoMAT%NDOF
        phi = phi + vec(j,i)*X(j)
      enddo
      coef(i) = phi/val(i)
    enddo

    Y = 0.0d0
    do i = 1, n_get_eigen
      do j = 1, monoMAT%NP*monoMAT%NDOF
        Y(j) = Y(j) + vec(j,i)*coef(i)
      enddo
    enddo
  end subroutine monolis_precond_ROM_apply

  subroutine monolis_precond_ROM_clear(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT

  end subroutine monolis_precond_ROM_clear
end module mod_monolis_precond_ROM