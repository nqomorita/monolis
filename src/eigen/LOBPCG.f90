module mod_monolis_eigen_LOBPCG
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_solve
  use mod_monolis_converge
  use mod_monolis_eigen_lanczos_util

  implicit none

contains

  subroutine monolis_eigen_inverted_lobpcg(monolis, n_get_eigen, ths)
    implicit none
    type(monolis_structure) :: monolis
    integer(kint) :: n_get_eigen
    real(kdouble) :: ths
    call monolis_eigen_inverted_lobpcg_mat(monolis%PRM, monolis%COM, monolis%MAT, n_get_eigen, ths)
  end subroutine monolis_eigen_inverted_lobpcg

  subroutine monolis_eigen_inverted_lobpcg_mat(monoPRM, monoCOM, monoMAT, n_get_eigen, ths)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kint), parameter :: W = 0 !> A
    integer(kint), parameter :: X = 1 !> A
    integer(kint), parameter :: P = 2 !> A
    integer(kint), parameter :: V = 0 !> B
    integer(kint), parameter :: Y = 1 !> B
    integer(kint), parameter :: Q = 2 !> B
    integer(kint) :: N, NP, NDOF, NG, total_dof
    integer(kint) :: i, j, iter, n_get_eigen
    real(kdouble) :: ths
    real(kdouble), allocatable :: mu(:), Sa(:,:), Sb(:,:)
    real(kdouble), allocatable :: A(:,:), B(:,:)
    real(kdouble), allocatable :: norm_x(:), norm_p(:), R0(:), R2(:), resid(:), lambda(:), coef(:,:)

    if(monoPRM%is_debug) call monolis_debug_header("monolis_eigen_inverted_lobpcg_mat")

    N     = monoMAT%N
    NP    = monoMAT%NP
    NDOF  = monoMAT%NDOF
    NG    = n_get_eigen

write(*,"(a,i8)")"n_get_eigen: ", n_get_eigen
write(*,"(a,1pe12.5)")"ths:     ", ths

    total_dof = N*NDOF
    call monolis_allreduce_I1(total_dof, monolis_sum, monoCOM%comm)

    allocate(A(NP*NDOF,3*NG), source = 0.0d0)
    allocate(B(NP*NDOF,3*NG), source = 0.0d0)
    allocate(mu(NG), source = 0.0d0)
    allocate(Sa(3*NG,3*NG), source = 0.0d0)
    allocate(Sb(3*NG,3*NG), source = 0.0d0)
    allocate(coef(3*NG,NG), source = 0.0d0)
    allocate(R0(NG), source = 0.0d0)
    allocate(R2(NG), source = 0.0d0)
    allocate(resid(NG), source = 0.0d0)
    allocate(norm_x(NG), source = 0.0d0)
    allocate(norm_p(NG), source = 0.0d0)
    allocate(lambda(NG), source = 0.0d0)

    call monolis_precond_setup(monoPRM, monoCOM, monoMAT)

    do i = 1, NG
      call lanczos_initialze(N*NDOF, A(:,NG*X+i))
      call monolis_matvec(monoCOM, monoMAT, A(:,NG*X+i), B(:,NG*Y+i), monoPRM%tspmv, monoPRM%tcomm_spmv)

      mu(i) = dot_product(A(:,NG*X+i), B(:,NG*Y+i))
      A(:,NG*W+i) = B(:,NG*Y+i) - mu(i)*A(:,NG*X+i)
      call monolis_precond_apply(monoPRM, monoCOM, monoMAT, A(:,NG*W+i), A(:,NG*W+i))

      call monolis_inner_product_R(monoCOM, monoMAT%N, monoMAT%NDOF, A(:,NG*W+i), A(:,NG*W+i), R0(i), &
        & monoPRM%tdotp, monoPRM%tcomm_dotp)
    enddo

    do iter = 1, 2*total_dof
      do i = 1, NG
        call monolis_matvec(monoCOM, monoMAT, A(:,NG*W+i), B(:,NG*V+i), monoPRM%tspmv, monoPRM%tcomm_spmv)
      enddo
      Sa = matmul(transpose(A), B)
      Sb = matmul(transpose(A), A)
      call monolis_get_smallest_eigen_pair_from_ng(iter, NG, Sa, Sb, lambda, coef)

      do i = 1, NG
        mu(i) = 0.5d0*(lambda(i) + dot_product(A(:,NG*X+i), B(:,NG*Y+i)))

        do j = 1, NG
          A(:,NG*X+i) = coef(i,j)*A(:,NG*W+i) + coef(NG+i,j)*A(:,NG*X+i) + coef(2*NG+i,j)*A(:,NG*P+i)
          A(:,NG*P+i) = coef(i,j)*A(:,NG*W+i)                            + coef(2*NG+i,j)*A(:,NG*P+i)
          B(:,NG*Y+i) = coef(i,j)*B(:,NG*V+i) + coef(NG+i,j)*B(:,NG*Y+i) + coef(2*NG+i,j)*B(:,NG*Q+i)
          B(:,NG*Q+i) = coef(i,j)*B(:,NG*V+i)                            + coef(2*NG+i,j)*B(:,NG*Q+i)
        enddo

        norm_x(i) = monolis_get_l2_norm(N*NDOF, A(:,NG*X+i))
        if(norm_x(i) == 0.0d0) stop "monolis_eigen_inverted_lobpcg_mat norm_x"

        norm_p(i) = monolis_get_l2_norm(N*NDOF, A(:,NG*P+i))
        if(norm_p(i) == 0.0d0) stop "monolis_eigen_inverted_lobpcg_mat norm_p"

        A(:,NG*X+i) = A(:,NG*X+i)/norm_x(i)
        A(:,NG*P+i) = A(:,NG*P+i)/norm_p(i)
        B(:,NG*Y+i) = B(:,NG*Y+i)/norm_x(i)
        B(:,NG*Q+i) = B(:,NG*Q+i)/norm_p(i)

        A(:,NG*W+i) = B(:,NG*Y+i) - mu(i)*A(:,NG*X+i)

        call monolis_inner_product_R(monoCOM, monoMAT%N, monoMAT%NDOF, A(:,NG*W+i), A(:,NG*W+i), R2(i), &
        &  monoPRM%tdotp, monoPRM%tcomm_dotp)

        resid(i) = dsqrt(R2(i)/R0(i))
      enddo
      write (*,"(i7, 1pe16.6)") iter, maxval(resid)

      if(maxval(resid) < ths)then
write(*,*)"eigen_value"
write(*,"(1p2e12.5)")lambda

write(*,*)"e_mode"
write(*,"(1p10e12.5)")A(:,X)
        exit
      endif

      do i = 1, NG
        call monolis_precond_apply(monoPRM, monoCOM, monoMAT, A(:,NG*W+i), A(:,NG*W+i))

        norm_x(i) = monolis_get_l2_norm(N*NDOF, A(:,NG*W+i))
        if(norm_x(i) == 0.0d0) stop "monolis_eigen_inverted_lobpcg_mat norm_w"

        A(:,NG*W+i) = A(:,NG*W+i)/norm_x(i)
      enddo
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

    call lanczos_initialze(N*NDOF, A(:,X))
    call monolis_matvec(monoCOM, monoMAT, A(:,X), B(:,Y), monoPRM%tspmv, monoPRM%tcomm_spmv)
    mu = dot_product(A(:,X), B(:,Y))

    A(:,W) = B(:,Y) - mu*A(:,X)
    call monolis_precond_apply(monoPRM, monoCOM, monoMAT, A(:,W), A(:,W))

    call monolis_inner_product_R(monoCOM, monoMAT%N, monoMAT%NDOF, A(:,W), A(:,W), R0, monoPRM%tdotp, monoPRM%tcomm_dotp)

    do iter = 1, 10*total_dof
      call monolis_matvec(monoCOM, monoMAT, A(:,W), B(:,V), monoPRM%tspmv, monoPRM%tcomm_spmv)
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

      A(:,W) = B(:,Y) - mu*A(:,X)
      call monolis_inner_product_R(monoCOM, monoMAT%N, monoMAT%NDOF, A(:,W), A(:,W), R2, monoPRM%tdotp, monoPRM%tcomm_dotp)

      resid = dsqrt(R2/R0)
      write (*,"(i7, 1pe16.6)") iter, resid

      if(resid < ths)then
write(*,*)"eigen_value"
write(*,"(1p2e12.5)")lambda

write(*,*)"e_mode"
write(*,"(1p10e12.5)")A(:,X)
        exit
      endif

      call monolis_precond_apply(monoPRM, monoCOM, monoMAT, A(:,W), A(:,W))
      norm_x = monolis_get_l2_norm(N*NDOF, A(:,W))
      if(norm_x == 0.0d0) stop "monolis_eigen_inverted_lobpcg_ norm_w"

      A(:,W) = A(:,W)/norm_x
    enddo

    !call monolis_precond_clear(monoPRM, monoCOM, monoMAT)
  end subroutine monolis_eigen_inverted_lobpcg_

end module mod_monolis_eigen_LOBPCG
