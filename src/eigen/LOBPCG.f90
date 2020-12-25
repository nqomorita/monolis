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
    call monolis_eigen_inverted_lobpcg_(monolis%PRM, monolis%COM, monolis%MAT, n_get_eigen, ths)
  end subroutine monolis_eigen_inverted_lobpcg

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
    real(kdouble) :: ths, mu, Sa(3,3), Sb(3,3), lambda, coef(3), norm_x, norm_p, R0
    real(kdouble), allocatable :: A(:,:), B(:,:)
    logical :: is_converge

    if(monoPRM%is_debug) call monolis_debug_header("monolis_eigen_inverted_lobpcg_")

    N     = monoMAT%N
    NP    = monoMAT%NP
    NDOF  = monoMAT%NDOF

    total_dof = N*NDOF
    call monolis_allreduce_I1(total_dof, monolis_sum, monoCOM%comm)

    allocate(A(NP*NDOF,3), source = 0.0d0)
    allocate(B(NP*NDOF,3), source = 0.0d0)

    call lanczos_initialze(N*NDOF, A(:,X))
    call monolis_matvec(monoCOM, monoMAT, A(:,X), B(:,Y), monoPRM%tspmv, monoPRM%tcomm_spmv)
    mu = dot_product(A(:,X), B(:,Y))

    A(:,W) = B(:,Y) - mu*A(:,X)
    call monolis_set_converge(monoPRM, monoCOM, monoMAT, A(:,W), R0, is_converge, monoPRM%tdotp, monoPRM%tcomm_dotp)

    do iter = 1, total_dof
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
      call monolis_check_converge(monoPRM, monoCOM, monoMAT, A(:,W), R0, iter, is_converge, monoPRM%tdotp, monoPRM%tcomm_dotp)
      if(is_converge) exit

      !call monolis_precond_apply(monoPRM, monoCOM, monoMAT, A(:,W), A(:,W))
      norm_x = monolis_get_l2_norm(N*NDOF, A(:,W))
      if(norm_x == 0.0d0) stop "monolis_eigen_inverted_lobpcg_ norm_w"

      A(:,W) = A(:,W)/norm_x
    enddo
  end subroutine monolis_eigen_inverted_lobpcg_

end module mod_monolis_eigen_LOBPCG
