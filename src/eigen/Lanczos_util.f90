module mod_monolis_eigen_lanczos_util
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_linalg

  implicit none

contains

  subroutine lanczos_initialze(monoCOM, N, NDOF, q, is_bc, beta)
    use mod_monolis_linalg_util
    implicit none
    type(monolis_com) :: monoCOM
    integer(kint) :: i, N, NDOF
    real(kdouble) :: q(:), norm, t1, t2, beta
    logical :: is_bc(:)

    call get_rundom_number(N*NDOF, q, monolis_global_myrank())

    do i = 1, N*NDOF
      if(is_bc(i)) q(i) = 0.0d0
    enddo

    call monolis_update_R(monoCOM, NDOF, q, t1)

    call monolis_inner_product_R(monoCOM, N, NDOF, q, q, norm, t1, t2)

    beta = dsqrt(norm)
    norm = 1.0d0/dsqrt(norm)
    do i = 1, N*NDOF
      q(i) = q(i)*norm
    enddo
  end subroutine lanczos_initialze

  subroutine get_rundom_number(N, X, SHIFT)
    implicit none
    real(kdouble) :: X(N), INVM
    integer(kint), parameter :: MM = 1664501
    integer(kint), parameter :: LAMBDA = 1229
    integer(kint), parameter :: MU = 351750
    integer(kint) :: i, N, IR, SHIFT

    IR = 0
    INVM = 1.0D0 / MM
    do I = 1, SHIFT
      IR = mod( LAMBDA * IR + MU, MM)
    enddo
    do I = SHIFT+1, SHIFT+N
      IR = mod( LAMBDA * IR + MU, MM)
      X(I-SHIFT) = INVM * IR
    enddo
  end subroutine get_rundom_number

  subroutine monolis_gram_schmidt(monoPRM, monoCOM, monoMAT, iter, q, p)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kint) :: i, j, iter, N, NDOF
    real(kdouble) :: q(:,0:), p(:), norm
    real(kdouble), allocatable :: t(:)

    N    = monoMAT%N
    NDOF = monoMAT%NDOF

    do i = 1, iter
      call monolis_inner_product_R(monoCOM, N, NDOF, p, q(:,i), norm, monoPRM%tdotp, monoPRM%tcomm_dotp)
      do j = 1, N*NDOF
        p(j) = p(j) - norm*q(j,i)
      enddo
    enddo
  end subroutine monolis_gram_schmidt

  subroutine monolis_lobpcg_initialze(q)
    implicit none
    real(kdouble) :: q(:,:)
    call random_number(q)
  end subroutine monolis_lobpcg_initialze

  subroutine monolis_gram_schmidt_lobpcg(monoPRM, monoCOM, monoMAT, q, prm, ng)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kint) :: i, j, k, ng, N, NDOF, prm
    real(kdouble) :: q(:,:), norm

    N    = monoMAT%N
    NDOF = monoMAT%NDOF

    do i = prm*ng + 1, prm*ng
      do j = 1, i - 1
        call monolis_inner_product_R(monoCOM, N, NDOF, q(:,i), q(:,j), norm, monoPRM%tdotp, monoPRM%tcomm_dotp)
        do k = 1, N*NDOF
          q(k,i) = q(k,i) - norm*q(k,j)
        enddo
      enddo

      norm = 0.0d0
      do j = 1, N*NDOF
        norm = norm + q(j,i)*q(j,i)
      enddo
      if(norm == 0.0d0) cycle

      norm = 1.0d0/dsqrt(norm)
      q(:,i) = q(:,i)*norm
    enddo
  end subroutine monolis_gram_schmidt_lobpcg

  subroutine monolis_get_smallest_eigen_pair_m(dof, Sa, vec, val)
    implicit none
    integer(kint) :: lda, info, iw, N, dof
    real(kdouble) :: Sa(:,:), vec(:,:), val(:)
    real(kdouble) :: ev(dof)
    real(kdouble), allocatable :: rw(:)

    N = dof
    lda = dof
    iw = 3*dof-1
    allocate(rw(iw), source = 0.0d0)
    call dsyev("V", "U", N, Sa, lda, ev, rw, iw, info)
    deallocate(rw)

    vec = Sa
    val = ev
  end subroutine monolis_get_smallest_eigen_pair_m

  subroutine monolis_get_smallest_eigen_pair_3m(iter, NG, Sa, Sb, lambda, coef)
    implicit none
    integer(kint) :: iter, it, lda, ldb, info, iw, N, dof, i, j, NG
    real(kdouble) :: Sa(:,:), Sb(:,:), lambda(:), coef(:,:)
    real(kdouble) :: e_value(3*NG)
    real(kdouble), allocatable :: rw(:)

    dof = 3*NG
    if(iter == 1) dof = 2*NG

    it = 1 !> A x = B lambda x
    N = dof
    lda = dof
    ldb = dof
    iw = 3*dof-1
    e_value = 0.0d0
    allocate(rw(iw), source = 0.0d0)

    call dsygv(it, "V", "L", N, Sa(1:dof,1:dof), lda, Sb(1:dof,1:dof), ldb, e_value(1:dof), rw, iw, info)

    lambda = 0.0d0
    do i = 1, NG
      lambda(i) = e_value(i)
    enddo

    coef = 0.0d0
    do j = 1, NG
      do i = 1, dof
        coef(i,j) = Sa(i,j)
      enddo
    enddo

    deallocate(rw)
  end subroutine monolis_get_smallest_eigen_pair_3m

  subroutine monolis_get_smallest_eigen_pair_from_3x3(iter, Sa, Sb, lambda, coef)
    implicit none
    integer(kint) :: iter, it, lda, ldb, info, iw, N, dof, i
    real(kdouble) :: Sa(3,3), Sb(3,3), lambda, coef(3), e_value(3), rw(8)

    dof = 3
    if(iter == 1) dof = 2

    it = 1
    N = dof
    lda = dof
    ldb = dof
    iw = 8
    call dsygv(it, "V", "L", N, Sa(1:dof,1:dof), lda, Sb(1:dof,1:dof), ldb, e_value(1:dof), rw, iw, info)

    lambda = e_value(1)

    coef = 0.0d0
    do i = 1, dof
      coef(i) = Sa(i,1)
    enddo
  end subroutine monolis_get_smallest_eigen_pair_from_3x3

  subroutine monolis_get_eigen_pair_from_tridiag(iter, n_get_eigen, &
    & alpha_t, beta_t, q, e_value, e_mode, norm)
    implicit none
    integer(kint), intent(in) :: iter
    integer(kint) :: i, n, ldz, info, n_get_eigen
    real(kdouble) :: alpha_t(:), beta_t(:), q(:,0:), e_value(:), e_mode(:,:), norm
    integer(kint), allocatable :: isuppz(:), idum(:)
    real(kdouble), allocatable :: alpha(:), beta(:), rdum(:), e_mode_t(:,:)

    !> DSTEVR
    allocate(alpha(iter), source = 0.0d0)
    allocate(beta (max(1,iter-1)), source = 0.0d0)
    allocate(rdum(20*iter), source = 0.0d0)
    allocate(e_mode_t(iter,iter), source = 0.0d0)

    alpha = alpha_t(1:iter)
    beta = beta_t(2:max(1,iter-1)+1)
    if(iter == 1) beta = 0.0d0

    n = iter
    ldz = iter

    call dstev("V", n, alpha, beta, e_mode_t, ldz, rdum, info)

    norm = 0.0d0
    do i = 1, min(iter, n_get_eigen)
      e_value(i) = 1.0d0/alpha(iter - i +1)
      e_mode(:,i) = matmul(q(:,1:iter), e_mode_t(1:iter,iter - i + 1))
      if(norm < sqrt(e_mode_t(iter,iter - i + 1)**2)*beta_t(iter+1))then
        norm = sqrt(e_mode_t(iter,iter - i + 1)**2)*beta_t(iter+1)
      endif
    enddo

    if(info /= 0) stop "monolis_get_eigen_pair_from_tridiag"

    deallocate(alpha)
    deallocate(beta)
    deallocate(e_mode_t)
    deallocate(rdum)
  end subroutine monolis_get_eigen_pair_from_tridiag

  subroutine monolis_get_normarize_vectors(vec, N, M)
    implicit none
    integer(kint) :: N, M, i, j
    real(kdouble) :: vec(:,:), norm

    do i = 1, M
      norm = 0.0d0
      do j = 1, N
        norm = norm + vec(j,i)*vec(j,i)
      enddo
      if(norm == 0.0d0) cycle
      norm = 1.0d0/dsqrt(norm)
      vec(:,i) = vec(:,i)*norm
    enddo
  end subroutine monolis_get_normarize_vectors
end module mod_monolis_eigen_lanczos_util
