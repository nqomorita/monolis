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
    call monolis_eigen_inverted_lobpcg_    (monolis%PRM, monolis%COM, monolis%MAT, n_get_eigen, ths)
    call monolis_eigen_inverted_lobpcg_mat(monolis%PRM, monolis%COM, monolis%MAT, n_get_eigen, ths)
  end subroutine monolis_eigen_inverted_lobpcg

  subroutine monolis_eigen_inverted_lobpcg_mat(monoPRM, monoCOM, monoMAT, n_get_eigen, ths)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kint) :: N, NP, NDOF, NG, total_dof
    integer(kint) :: i, j, iter, n_get_eigen
    real(kdouble) :: ths
    real(kdouble), allocatable :: X(:,:), W(:,:), P(:,:), AX(:,:), AW(:,:), AP(:,:)
    real(kdouble), allocatable :: XAW(:,:), XAP(:,:), WAW(:,:), WAP(:,:), PAP(:,:), XW(:,:), XP(:,:), WP(:,:)
    real(kdouble), allocatable :: lt(:,:), lv(:,:), Sa(:,:), Sb(:,:), St(:,:), eval(:), evec(:,:)
    real(kdouble), allocatable :: R0(:), R2(:), resid(:), norm_x(:), norm_p(:)
    real(kdouble), allocatable :: A(:,:), B(:,:)

    if(monoPRM%is_debug) call monolis_debug_header("monolis_eigen_inverted_lobpcg_mat")

    N    = monoMAT%N
    NP   = monoMAT%NP
    NDOF = monoMAT%NDOF
    NG   = n_get_eigen

write(*,"(a,i8)")"n_get_eigen: ", n_get_eigen
write(*,"(a,1pe12.5)")"ths:     ", ths

    total_dof = N*NDOF
    call monolis_allreduce_I1(total_dof, monolis_sum, monoCOM%comm)

    allocate(A(NP*NDOF,3*NG), source = 0.0d0)
    allocate(B(NP*NDOF,3*NG), source = 0.0d0)

    allocate(X(NP*NDOF,NG), source = 0.0d0)
    allocate(W(NP*NDOF,NG), source = 0.0d0)
    allocate(P(NP*NDOF,NG), source = 0.0d0)
    allocate(AX(NP*NDOF,NG), source = 0.0d0)
    allocate(AW(NP*NDOF,NG), source = 0.0d0)
    allocate(AP(NP*NDOF,NG), source = 0.0d0)
    allocate(XAW(NG,NG), source = 0.0d0)
    allocate(XAP(NG,NG), source = 0.0d0)
    allocate(WAW(NG,NG), source = 0.0d0)
    allocate(WAP(NG,NG), source = 0.0d0)
    allocate(PAP(NG,NG), source = 0.0d0)
    allocate(XW(NG,NG), source = 0.0d0)
    allocate(XP(NG,NG), source = 0.0d0)
    allocate(WP(NG,NG), source = 0.0d0)
    allocate(St(NG,NG), source = 0.0d0)
    allocate(Sa(3*NG,3*NG), source = 0.0d0)
    allocate(Sb(3*NG,3*NG), source = 0.0d0)
    allocate(evec(3*NG,NG), source = 0.0d0)
    allocate(R0(NG), source = 0.0d0)
    allocate(R2(NG), source = 0.0d0)
    allocate(resid(NG), source = 0.0d0)
    allocate(norm_x(NG), source = 0.0d0)
    allocate(norm_p(NG), source = 0.0d0)
    allocate(eval(NG), source = 0.0d0)
    allocate(lt(NG,NG), source = 0.0d0)
    allocate(lv(NG,NG), source = 0.0d0)

    call monolis_precond_setup(monoPRM, monoCOM, monoMAT)
    call monolis_lobpcg_initialze(X)
    call monolis_gram_schmidt_lobpcg(monoPRM, monoCOM, monoMAT, X, 0, NG)
    do i = 1, NG
      call monolis_matvec(monoCOM, monoMAT, X(:,i), AX(:,i), monoPRM%tspmv, monoPRM%tcomm_spmv)
    enddo

    St = matmul(transpose(X), AX)
    call monolis_get_smallest_eigen_pair(NG, St, lt, lv)
     X = matmul( X, lv)
    AX = matmul(AX, lv)

    do i = 1, NG
      eval(i) = lt(i,i)
    enddo

    !write(*,"(1pe12.5)")eval

    do iter = 1, 2*total_dof
      do i = 1, NG
        W(:,i) = AX(:,i) - eval(i)*X(:,i)
        call monolis_precond_apply(monoPRM, monoCOM, monoMAT, W(:,i), W(:,i))

        call monolis_matvec(monoCOM, monoMAT, W(:,i), AW(:,i), monoPRM%tspmv, monoPRM%tcomm_spmv)
        !call monolis_matvec(monoCOM, monoMAT, P(:,i), AP(:,i), monoPRM%tspmv, monoPRM%tcomm_spmv)

        call monolis_inner_product_R(monoCOM, monoMAT%N, monoMAT%NDOF, W(:,i), W(:,i), R2(i), &
           & monoPRM%tdotp, monoPRM%tcomm_dotp)
        if(iter == 1) R0(i) = R2(i)
        resid(i) = dsqrt(R2(i)/R0(i))
      enddo

      write (*,"(i7, 1pe16.6)") iter, maxval(resid)

      XAW = matmul(transpose(X), AW)
      XAP = matmul(transpose(X), AP)
      WAW = matmul(transpose(W), AW)
      WAP = matmul(transpose(W), AP)
      PAP = matmul(transpose(P), AP)
      XW  = 0.0d0 !matmul(transpose(X), W)
      XP  = matmul(transpose(X), P)
      WP  = matmul(transpose(W), P)

      Sa = 0.0d0
      Sb = 0.0d0
      do i = 1, NG
        Sa(i,i) = eval(i)
      enddo
      do i = 1, 3*NG
        Sb(i,i) = 1.0d0
      enddo

      Sa(     1:  NG,  NG+1:2*NG) = XAW
      Sa(     1:  NG,2*NG+1:3*NG) = XAP
      Sa(  NG+1:2*NG,     1:  NG) = transpose(XAW)
      Sa(  NG+1:2*NG,  NG+1:2*NG) = WAW
      Sa(  NG+1:2*NG,2*NG+1:3*NG) = WAP
      Sa(2*NG+1:3*NG,     1:  NG) = transpose(XAP)
      Sa(2*NG+1:3*NG,  NG+1:2*NG) = transpose(WAP)
      Sa(2*NG+1:3*NG,2*NG+1:3*NG) = PAP

      Sb(     1:1*NG,  NG+1:2*NG) = XW
      Sb(     1:1*NG,2*NG+1:3*NG) = XP
      Sb(  NG+1:2*NG,     1:  NG) = transpose(XW)
      Sb(  NG+1:2*NG,2*NG+1:3*NG) = WP
      Sb(2*NG+1:3*NG,     1:  NG) = transpose(XP)
      Sb(2*NG+1:3*NG,  NG+1:2*NG) = transpose(WP)

      call monolis_get_smallest_eigen_pair_from_ng(iter, NG, Sa, Sb, eval, evec)

!write(*,*)"eval"
!write(*,"(1p5e12.5)")eval
!write(*,*)"evec"
!write(*,"(1p5e12.5)")evec

       P = matmul( W, evec(NG+1:2*NG,:)) + matmul( P, evec(2*NG+1:3*NG,:))
      AP = matmul(AW, evec(NG+1:2*NG,:)) + matmul(AP, evec(2*NG+1:3*NG,:))

       X = matmul( X, evec(   1:  NG,:)) +  P
      AX = matmul(AX, evec(   1:  NG,:)) + AP

      ! P = matmul( W, evec(NG+1:2*NG,:)) + matmul( P, evec(2*NG+1:3*NG,:))
      !AP = matmul(AW, evec(NG+1:2*NG,:)) + matmul(AP, evec(2*NG+1:3*NG,:))

      ! X = matmul( X, evec(   1:  NG,:)) +  P
      !AX = matmul(AX, evec(   1:  NG,:)) + AP

      !call monolis_gram_schmidt_lobpcg(monoPRM, monoCOM, monoMAT, X, 0, NG)
      !call monolis_gram_schmidt_lobpcg(monoPRM, monoCOM, monoMAT, P, 0, NG)

      do i = 1, NG
        norm_x(i) = monolis_get_l2_norm(N*NDOF, X(:,i))
        norm_p(i) = monolis_get_l2_norm(N*NDOF, P(:,i))
        if(norm_x(i) == 0.0d0) stop "monolis_eigen_inverted_lobpcg_mat norm_x"
        if(norm_p(i) == 0.0d0) stop "monolis_eigen_inverted_lobpcg_mat norm_p"

         X(:,i) =  X(:,i)/norm_x(i)
         P(:,i) =  P(:,i)/norm_p(i)
        AX(:,i) = AX(:,i)/norm_x(i)
        AP(:,i) = AP(:,i)/norm_p(i)
      enddo

      if(maxval(resid) < ths .or. iter == 2*total_dof)then
write(*,*)"eigen_value"
write(*,"(1p2e12.5)")eval(1)

write(*,*)"e_mode"
write(*,"(1p10e12.5)")X(:,1)
        exit
      endif
    enddo
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
