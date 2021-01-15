module mod_monolis_eigen_LOBPCG
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_solve
  use mod_monolis_converge
  use mod_monolis_eigen_lanczos_util

  implicit none

contains

  subroutine monolis_eigen_inverted_lobpcg(monolis, n_get_eigen, ths, vec)
    implicit none
    type(monolis_structure) :: monolis
    integer(kint) :: n_get_eigen
    real(kdouble) :: ths
    real(kdouble) :: vec(:,:)

    call monolis_eigen_inverted_lobpcg_mat (monolis%PRM, monolis%COM, monolis%MAT, n_get_eigen, ths, vec)
  end subroutine monolis_eigen_inverted_lobpcg

  subroutine monolis_eigen_inverted_lobpcg_mat(monoPRM, monoCOM, monoMAT, n_get_eigen, ths, vec)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kint), parameter :: XX = 0 !> A
    integer(kint), parameter :: RR = 1 !> A
    integer(kint), parameter :: PP = 2 !> A
    integer(kint) :: N, NP, NDOF, NG, total_dof
    integer(kint) :: i, j, iter, maxiter,  n_get_eigen
    real(kdouble) :: ths, mu
    real(kdouble) :: vec(:,:)
    real(kdouble), allocatable :: X(:,:), R(:,:), P(:,:), XAX(:,:), XBX(:,:), evec(:,:), eval(:), T(:)
    real(kdouble), allocatable :: A(:,:), B(:,:), R0(:), R2(:), resid(:)

    if(monoPRM%is_debug) call monolis_debug_header("monolis_eigen_inverted_lobpcg_mat")

    N    = monoMAT%N
    NP   = monoMAT%NP
    NDOF = monoMAT%NDOF
    NG   = n_get_eigen

    total_dof = N*NDOF
    call monolis_allreduce_I1(total_dof, monolis_sum, monoCOM%comm)
    maxiter = 10*total_dof

    if(3*NG > total_dof)then
      stop "monolis_eigen_inverted_lobpcg_mat: 3*NG > DOF"
    endif

    allocate(A(NP*NDOF,3*NG), source = 0.0d0)
    allocate(B(NP*NDOF,3*NG), source = 0.0d0)
    allocate(X(NP*NDOF,NG), source = 0.0d0)
    allocate(R(NP*NDOF,NG), source = 0.0d0)
    allocate(P(NP*NDOF,NG), source = 0.0d0)
    allocate(T(NP*NDOF), source = 0.0d0)
    allocate(XAX(NG,NG), source = 0.0d0)
    allocate(evec(NG,NG), source = 0.0d0)
    allocate(eval(NG), source = 0.0d0)
    allocate(R0(NG), source = 0.0d0)
    allocate(R2(NG), source = 0.0d0)
    allocate(resid(NG), source = 0.0d0)

    !> 0 step
    call monolis_precond_setup(monoPRM, monoCOM, monoMAT)

    call random_number(X)
    call monolis_gram_schmidt_lobpcg(monoPRM, monoCOM, monoMAT, X, 0, NG)

    do i = 1, NG
      call monolis_matvec(monoCOM, monoMAT, X(:,i), R(:,i), monoPRM%tspmv, monoPRM%tcomm_spmv)
    enddo
    XAX = matmul(transpose(X), R)
    call monolis_get_smallest_eigen_pair_m(NG, XAX, evec, eval)

    X(:,1:NG) = matmul(X, evec)
    call monolis_get_normarize_vectors(X, NP*NDOF, NG)
    do i = 1, NG
      call monolis_matvec(monoCOM, monoMAT, X(:,i), T, monoPRM%tspmv, monoPRM%tcomm_spmv)
      R(:,i) = T - XAX(i,i)*X(:,i)
    enddo

    !> 1 step
    deallocate(XAX)
    deallocate(evec)
    deallocate(eval)
    allocate(XAX(2*NG,2*NG), source = 0.0d0)
    allocate(evec(2*NG,2*NG), source = 0.0d0)
    allocate(eval(2*NG), source = 0.0d0)

    do i = 1, NG
      A(:,NG*XX+i) = X(:,i)
      A(:,NG*RR+i) = R(:,i)
    enddo
    call monolis_gram_schmidt_lobpcg(monoPRM, monoCOM, monoMAT, A, 0, 2*NG)

    do i = 1, NG
      call monolis_matvec(monoCOM, monoMAT, A(:,NG*XX+i), B(:,NG*XX+i), monoPRM%tspmv, monoPRM%tcomm_spmv)
      call monolis_matvec(monoCOM, monoMAT, A(:,NG*RR+i), B(:,NG*RR+i), monoPRM%tspmv, monoPRM%tcomm_spmv)
    enddo

    XAX = matmul(transpose(A(:,1:2*NG)), B(:,1:2*NG))
    XBX = matmul(transpose(A(:,1:2*NG)), A(:,1:2*NG))
    call monolis_get_smallest_eigen_pair_3m(1, NG, XAX, XBX, eval, evec)

    X(:,1:NG) = matmul(A(:,1:2*NG), evec(1:2*NG,1:NG))
    call monolis_get_normarize_vectors(X, NP*NDOF, NG)
    do i = 1, NG
      call monolis_matvec(monoCOM, monoMAT, X(:,i), T, monoPRM%tspmv, monoPRM%tcomm_spmv)
      R(:,i) = T - eval(i)*X(:,i)
      call monolis_precond_apply(monoPRM, monoCOM, monoMAT, R(:,i), R(:,i))
      call monolis_inner_product_R(monoCOM, monoMAT%N, monoMAT%NDOF, R(:,i), R(:,i), R0(i), monoPRM%tdotp, monoPRM%tcomm_dotp)
    enddo

    P(:,1:NG) = matmul(A(:,NG+1:2*NG), evec(NG+1:2*NG,1:NG))

    !> 2 step
    deallocate(XAX)
    deallocate(evec)
    deallocate(eval)
    allocate(XAX(3*NG,3*NG), source = 0.0d0)
    allocate(evec(3*NG,3*NG), source = 0.0d0)
    allocate(eval(3*NG), source = 0.0d0)

    do iter = 1, maxiter
      do i = 1, NG
        A(:,NG*XX+i) = X(:,i)
        A(:,NG*RR+i) = R(:,i)
        A(:,NG*PP+i) = P(:,i)
      enddo
      call monolis_gram_schmidt_lobpcg(monoPRM, monoCOM, monoMAT, A, 0, 3*NG)

      do i = 1, NG
        call monolis_matvec(monoCOM, monoMAT, A(:,NG*XX+i), B(:,NG*XX+i), monoPRM%tspmv, monoPRM%tcomm_spmv)
        call monolis_matvec(monoCOM, monoMAT, A(:,NG*RR+i), B(:,NG*RR+i), monoPRM%tspmv, monoPRM%tcomm_spmv)
        call monolis_matvec(monoCOM, monoMAT, A(:,NG*PP+i), B(:,NG*PP+i), monoPRM%tspmv, monoPRM%tcomm_spmv)
      enddo

      XAX = matmul(transpose(A), B)
      XBX = matmul(transpose(A), A)
      call monolis_get_smallest_eigen_pair_3m(iter+1, NG, XAX, XBX, eval, evec)

      X(:,1:NG) = matmul(A, evec(:,1:NG))
      P(:,1:NG) = matmul(A(:,NG+1:3*NG), evec(NG+1:3*NG,1:NG))

      do i = 1, NG
        call monolis_matvec(monoCOM, monoMAT, X(:,i), T, monoPRM%tspmv, monoPRM%tcomm_spmv)
        R(:,i) = T - eval(i)*X(:,i)
        call monolis_precond_apply(monoPRM, monoCOM, monoMAT, R(:,i), R(:,i))
        call monolis_inner_product_R(monoCOM, monoMAT%N, monoMAT%NDOF, R(:,i), R(:,i), R2(i), monoPRM%tdotp, monoPRM%tcomm_dotp)
        resid(i) = dsqrt(R2(i)/R0(i))
      enddo
      write (*,"(i7, 1pe16.6)") iter, maxval(resid)

      if(maxval(resid) < ths)then
!write(*,*)"eigen_value"
!write(*,"(1pe12.5)")eval(1:NG)
!write(*,*)"e_mode"
!write(*,"(1p6e12.5)")X(:,1:NG)
        vec = X
        exit
      endif
    enddo

    !call monolis_precond_clear(monoPRM, monoCOM, monoMAT)
  end subroutine monolis_eigen_inverted_lobpcg_mat

  subroutine monolis_eigen_inverted_lobpcg_(monoPRM, monoCOM, monoMAT, n_get_eigen, ths)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kint), parameter :: W = 1 !> A
    integer(kint), parameter :: X = 2 !> A
    integer(kint), parameter :: P = 3 !> A
    integer(kint), parameter :: V = 1 !> B
    integer(kint), parameter :: Y = 2 !> B
    integer(kint), parameter :: Q = 3 !> B
    integer(kint) :: N, NP, NDOF, total_dof
    integer(kint) :: i, iter, n_get_eigen
    real(kdouble) :: ths, mu, Sa(3,3), Sb(3,3), lambda, coef(3), norm_x, norm_p, R0, R2, resid
    real(kdouble), allocatable :: A(:,:), B(:,:)

    if(monoPRM%is_debug) call monolis_debug_header("monolis_eigen_inverted_lobpcg_")

    N     = monoMAT%N
    NP    = monoMAT%NP
    NDOF  = monoMAT%NDOF

    total_dof = N*NDOF
    call monolis_allreduce_I1(total_dof, monolis_sum, monoCOM%comm)

    allocate(A(NP*NDOF,3), source = 0.0d0)
    allocate(B(NP*NDOF,3), source = 0.0d0)

    call monolis_precond_setup(monoPRM, monoCOM, monoMAT)

    call monolis_lobpcg_initialze(A(:,X:X))
    call monolis_gram_schmidt_lobpcg(monoPRM, monoCOM, monoMAT, A, X, 1)

    call monolis_matvec(monoCOM, monoMAT, A(:,X), B(:,Y), monoPRM%tspmv, monoPRM%tcomm_spmv)
    mu = dot_product(A(:,X), B(:,Y)) !> if block matrix, calculate eigen value

    do iter = 1, 2*total_dof
      A(:,W) = B(:,Y) - mu*A(:,X)
      call monolis_precond_apply(monoPRM, monoCOM, monoMAT, A(:,W), A(:,W))
      call monolis_matvec(monoCOM, monoMAT, A(:,W), B(:,V), monoPRM%tspmv, monoPRM%tcomm_spmv)

      call monolis_inner_product_R(monoCOM, monoMAT%N, monoMAT%NDOF, A(:,W), A(:,W), R2, monoPRM%tdotp, monoPRM%tcomm_dotp)
      if(iter == 1) R0 = R2
      resid = dsqrt(R2/R0)
      write (*,"(i7, 1pe16.6)") iter, resid

      Sa = matmul(transpose(A), B)
      Sb = matmul(transpose(A), A)

      call monolis_get_smallest_eigen_pair_from_3x3(iter, Sa, Sb, lambda, coef)

      mu = 0.5d0*(lambda + dot_product(A(:,X), B(:,Y)))
      A(:,X) = coef(1)*A(:,W) + coef(2)*A(:,X) + coef(3)*A(:,P)
      A(:,P) = coef(1)*A(:,W)                  + coef(3)*A(:,P)
      B(:,Y) = coef(1)*B(:,V) + coef(2)*B(:,Y) + coef(3)*B(:,Q)
      B(:,Q) = coef(1)*B(:,V)                  + coef(3)*B(:,Q)

      norm_x = monolis_get_l2_norm(N*NDOF, A(:,X))
      norm_p = monolis_get_l2_norm(N*NDOF, A(:,P))
      if(norm_x == 0.0d0) stop "monolis_eigen_inverted_lobpcg_ norm_x"
      if(norm_p == 0.0d0) stop "monolis_eigen_inverted_lobpcg_ norm_p"

      A(:,X) = A(:,X)/norm_x
      A(:,P) = A(:,P)/norm_p
      B(:,Y) = B(:,Y)/norm_x
      B(:,Q) = B(:,Q)/norm_p

      if(resid < ths .or. iter == 2*total_dof)then
write(*,*)"eigen_value"
write(*,"(1p2e12.5)")lambda

write(*,*)"e_mode"
write(*,"(1p10e12.5)")A(:,X)
        exit
      endif
    enddo

    !call monolis_precond_clear(monoPRM, monoCOM, monoMAT)
  end subroutine monolis_eigen_inverted_lobpcg_

end module mod_monolis_eigen_LOBPCG
