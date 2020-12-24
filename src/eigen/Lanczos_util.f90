module mod_monolis_eigen_lanczos_util
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat

  implicit none

contains

  subroutine lanczos_initialze(n, q)
    implicit none
    integer(kint) :: i, n
    real(kdouble) :: q(:), norm

    norm = 0.0d0
    do i = 1, n
      q(i) = 1.0d0
      norm = norm + q(i)*q(i)
    enddo

    norm = 1.0d0/dsqrt(norm)
    do i = 1, n
      q(i) = q(i)*norm
    enddo
  end subroutine lanczos_initialze

  subroutine monolis_get_eigen_pair_from_tridiag(iter, alpha_t, beta_t, q, e_value, e_mode)
    implicit none
    integer(kint) :: iter, i, n, m, iu, il, ldz, info, liwork, lwork
    real(kdouble) :: alpha_t(:), beta_t(:), q(:,:), e_value(:), e_mode(:,:)
    real(kdouble) :: vl, vu, abstol
    integer(kint), allocatable :: isuppz(:), idum(:)
    real(kdouble), allocatable :: alpha(:), beta(:), rdum(:)

    !> DSTEVR
    if(iter <= 1) return

    allocate(alpha(iter), source = 0.0d0)
    allocate(beta (iter-1), source = 0.0d0)
    allocate(isuppz(2*iter), source = 0)
    allocate(idum(10*iter), source = 0)
    allocate(rdum(20*iter), source = 0.0d0)

    alpha = alpha_t(1:iter)
    beta = beta_t(1:iter-1)

    vl = 0.0d0
    vu = 0.0d0
    il = 0
    iu = 0
    abstol = 1.0d-8
    n = iter
    m = iter
    ldz = iter
    lwork = 20*iter
    liwork = 10*iter

    call dstevr("V", "A", n, alpha, beta, vl, vu, il, iu, abstol, m, e_value, &
        e_mode, ldz, isuppz, rdum, lwork, idum, liwork, info)

    write(*,*)"e_value"
    write(*,"(1pe12.5)")e_value(1:iter)

    deallocate(alpha)
    deallocate(beta)
    deallocate(isuppz)
    deallocate(idum)
    deallocate(rdum)
  end subroutine monolis_get_eigen_pair_from_tridiag
end module mod_monolis_eigen_lanczos_util
