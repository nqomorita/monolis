module mod_monolis_eigen_lanczos
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_solve
  use mod_monolis_eigen_lanczos_util
  use mod_monolis_util_debug

  implicit none

contains

  subroutine monolis_eigen_inverted_standard_lanczos &
    & (monolis, n_get_eigen, ths, maxiter, val, vec, is_bc)
    implicit none
    type(monolis_structure) :: monolis
    integer(kint) :: n_get_eigen, maxiter
    real(kdouble) :: ths, val(:), vec(:,:)
    logical :: is_bc(:)

    call monolis_eigen_inverted_standard_lanczos_(monolis%PRM, monolis%COM, monolis%MAT, &
      & n_get_eigen, ths, maxiter, val, vec, is_bc)
  end subroutine monolis_eigen_inverted_standard_lanczos

  subroutine monolis_eigen_standard_lanczos &
    & (monolis, n_get_eigen, ths, maxiter, val, vec, is_bc)
    implicit none
    type(monolis_structure) :: monolis
    integer(kint) :: n_get_eigen, maxiter
    real(kdouble) :: ths, val(:), vec(:,:)
    logical :: is_bc(:)

    call monolis_eigen_standard_lanczos_(monolis%PRM, monolis%COM, monolis%MAT, &
      & n_get_eigen, ths, maxiter, val, vec, is_bc)
  end subroutine monolis_eigen_standard_lanczos

  subroutine monolis_eigen_inverted_standard_lanczos_ &
    & (monoPRM, monoCOM, monoMAT, n_get_eigen, ths, maxiter, val, vec, is_bc)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kint) :: N, NP, NDOF, total_dof, j, k
    integer(kint) :: i, iter, maxiter, n_get_eigen, n_bc
    real(kdouble) :: beta_t, ths, norm, tmp
    real(kdouble) :: vec(:,:), val(:)
    real(kdouble), allocatable :: p(:), q(:,:), alpha(:), beta(:), eigen_value(:), eigen_mode(:,:), prev(:)
    logical :: is_converge
    logical :: is_bc(:)

    if(monoPRM%is_debug) call monolis_debug_header("monolis_eigen_inverted_standard_lanczos_")

    call monolis_set_initial_comm(monoCOM, monoMAT)

    N     = monoMAT%N
    NP    = monoMAT%NP
    NDOF  = monoMAT%NDOF
    norm = 0.0d0
    is_converge = .false.

    total_dof = N*NDOF

    n_bc = 0
    do i = 1, N*NDOF
      if(is_bc(i)) n_bc = n_bc + 1
    enddo

    call monolis_allreduce_I1(total_dof, monolis_sum, monoCOM%comm)
    call monolis_allreduce_I1(n_bc, monolis_sum, monoCOM%comm)
    total_dof = total_dof - n_bc

    if(n_get_eigen > total_dof) n_get_eigen = total_dof

    allocate(alpha(maxiter), source = 0.0d0)
    allocate(beta(maxiter+1), source = 0.0d0)
    allocate(eigen_value(maxiter), source = 0.0d0)
    allocate(prev(maxiter), source = 0.0d0)
    allocate(q(NP*NDOF,0:maxiter+1), source = 0.0d0)
    allocate(p(NP*NDOF), source = 0.0d0)
    allocate(eigen_mode(NP*NDOF,n_get_eigen), source = 0.0d0)

    call lanczos_initialze(monoCOM, N, NDOF, q(:,1), is_bc, beta(1))

    do iter = 1, maxiter
      call monolis_set_RHS(monoMAT, q(:,iter))

      call monolis_solve_(monoPRM, monoCOM, monoMAT)

      do i = 1, N*NDOF
        if(is_bc(i)) monoMAT%X(i) = 0.0d0
      enddo

      call monolis_vec_AXPY(N, NDOF, -beta(iter), q(:,iter-1), monoMAT%X, p)

      call monolis_inner_product_R(monoCOM, N, NDOF, p, q(:,iter), alpha(iter))

      call monolis_vec_AXPY(N, NDOF, -alpha(iter), q(:,iter), p, p)

      call monolis_gram_schmidt(monoPRM, monoCOM, monoMAT, iter, q, p)

      call monolis_inner_product_R(monoCOM, N, NDOF, p, p, beta_t)

      beta(iter+1) = dsqrt(beta_t)
      beta_t = 1.0d0/beta(iter+1)
      do i = 1, NP*NDOF
        q(i,iter+1) = p(i)*beta_t
      enddo

      call monolis_get_eigen_pair_from_tridiag(iter, n_get_eigen, &
        & alpha, beta, q, eigen_value, eigen_mode, norm)

      if(norm < ths) is_converge = .true.

if(monolis_global_myrank() == 0)then
  write(*,"(a,i6,a,1p2e12.4)")"iter: ", iter, ", ths: ", norm
endif

      if(is_converge .or. iter >= total_dof .or. iter == maxiter)then
        do i = 1, n_get_eigen
          val(i) = eigen_value(i)
          do j = 1, NP*NDOF
            vec(j,i) = eigen_mode(j,i)
          enddo
        enddo
        exit
      endif
    enddo
  end subroutine monolis_eigen_inverted_standard_lanczos_

  subroutine monolis_eigen_standard_lanczos_ &
    & (monoPRM, monoCOM, monoMAT, n_get_eigen, ths, maxiter, val, vec, is_bc)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kint) :: N, NP, NDOF, total_dof, n_bc, j, k
    integer(kint) :: i, iter, maxiter, n_get_eigen
    real(kdouble) :: beta_t, ths, norm, tmp
    real(kdouble) :: vec(:,:), val(:)
    real(kdouble), allocatable :: p(:), q(:,:), alpha(:), beta(:), eigen_value(:), eigen_mode(:,:), prev(:)
    logical :: is_converge
    logical :: is_bc(:)

    call monolis_set_initial_comm(monoCOM, monoMAT)

    if(monoPRM%is_debug) call monolis_debug_header("monolis_eigen_standard_lanczos_")

    N     = monoMAT%N
    NP    = monoMAT%NP
    NDOF  = monoMAT%NDOF
    norm = 0.0d0
    is_converge = .false.

    total_dof = N*NDOF
    write(*,*)N, NP, NDOF

    n_bc = 0
    do i = 1, N*NDOF
      if(is_bc(i)) n_bc = n_bc + 1
    enddo

    call monolis_allreduce_I1(total_dof, monolis_sum, monoCOM%comm)
    call monolis_allreduce_I1(n_bc, monolis_sum, monoCOM%comm)
    total_dof = total_dof - n_bc

    if(n_get_eigen > total_dof) n_get_eigen = total_dof

    allocate(alpha(maxiter), source = 0.0d0)
    allocate(beta(maxiter+1), source = 0.0d0)
    allocate(eigen_value(maxiter), source = 0.0d0)
    allocate(prev(maxiter), source = 0.0d0)
    allocate(q(NP*NDOF,0:maxiter+1), source = 0.0d0)
    allocate(p(NP*NDOF), source = 0.0d0)
    allocate(eigen_mode(NP*NDOF,n_get_eigen), source = 0.0d0)

    call lanczos_initialze(monoCOM, N, NDOF, q(:,1), is_bc, beta(1))

    do iter = 1, maxiter
      call monolis_matvec(monoCOM, monoMAT, q(:,iter), monoMAT%X, monoPRM%tspmv, monoPRM%tcomm_spmv)

      do i = 1, N*NDOF
        if(is_bc(i)) monoMAT%X(i) = 0.0d0
      enddo

      call monolis_vec_AXPY(N, NDOF, -beta(iter), q(:,iter-1), monoMAT%X, p)

      call monolis_inner_product_R(monoCOM, N, NDOF, p, q(:,iter), alpha(iter))

      call monolis_vec_AXPY(N, NDOF, -alpha(iter), q(:,iter), p, p)

      call monolis_gram_schmidt(monoPRM, monoCOM, monoMAT, iter, q, p)

      call monolis_inner_product_R(monoCOM, N, NDOF, p, p, beta_t)

      beta(iter+1) = dsqrt(beta_t)
      beta_t = 1.0d0/beta(iter+1)
      do i = 1, NP*NDOF
        q(i,iter+1) = p(i)*beta_t
      enddo

      call monolis_get_eigen_pair_from_tridiag(iter, n_get_eigen, &
        & alpha, beta, q, eigen_value, eigen_mode, norm)

      do i = 1, min(iter, n_get_eigen)
        eigen_value(i) = 1.0d0/eigen_value(i)
      enddo

      if(norm < ths) is_converge = .true.

if(monolis_global_myrank() == 0)then
  write(*,"(a,i6,a,1p2e12.4)")"iter: ", iter, ", ths: ", norm
endif

      if(is_converge .or. iter >= total_dof .or. iter == maxiter)then
        do i = 1, n_get_eigen
          val(i) = eigen_value(i)
          do j = 1, NP*NDOF
            vec(j,i) = eigen_mode(j,i)
          enddo
        enddo
        exit
      endif
    enddo
  end subroutine monolis_eigen_standard_lanczos_

  !> c interface
  subroutine monolis_eigen_inverted_standard_lanczos_c(N, NP, NZ, NDOF, A, index, item, &
    myrank, comm, commsize, &
    recv_n_neib, recv_nitem, recv_neib_pe, recv_index, recv_item, &
    send_n_neib, send_nitem, send_neib_pe, send_index, send_item, &
    method, precond, maxiter, tol, &
    iterlog, timelog_statistics, timelog, summary, &
    is_check_diag, is_measurement, is_init_x, curiter, curresid, time, &
    n_get_eigen, ths, eigen_maxiter, eigen_value, eigen_mode_tmp, is_bc_int) &
    & bind(c, name = "monolis_eigen_inverted_standard_lanczos_c_main")
    implicit none
    type(monolis_structure) :: monolis
    integer(c_int), intent(in), value :: N, NP, NZ, NDOF
    integer(c_int), intent(in), value :: myrank, comm, commsize
    integer(c_int), intent(in), value :: recv_n_neib, send_n_neib, recv_nitem, send_nitem
    integer(c_int), intent(in), value :: method, precond, maxiter
    integer(c_int), intent(in), value :: iterlog, timelog, timelog_statistics, summary
    integer(c_int), intent(in), value :: is_check_diag, is_measurement, is_init_x
    integer(c_int), intent(in), value :: n_get_eigen
    integer(c_int), intent(in), value :: eigen_maxiter
    integer(c_int), intent(in), target :: index(0:NP)
    integer(c_int), intent(in), target :: item(NZ)
    integer(c_int), intent(in), target :: recv_neib_pe(recv_n_neib)
    integer(c_int), intent(in), target :: recv_index(0:recv_n_neib), recv_item(recv_nitem)
    integer(c_int), intent(in), target :: send_neib_pe(send_n_neib)
    integer(c_int), intent(in), target :: send_index(0:send_n_neib), send_item(send_nitem)
    integer(c_int), intent(in), target :: is_bc_int(NDOF*NP)
    integer(c_int), intent(out), target :: curiter
    real(c_double), intent(in), value :: tol, ths
    real(c_double), intent(in), target :: A(NDOF*NDOF*NZ)
    real(c_double), intent(inout), target :: eigen_value(NDOF*NP)
    real(c_double), intent(inout), target :: eigen_mode_tmp(NDOF*NP*n_get_eigen)
    real(c_double), intent(out), target :: curresid
    real(c_double), intent(out), target :: time(7)
    integer(kint) :: i, j
    real(kdouble), allocatable :: vec(:,:), val(:)
    logical, allocatable :: is_bc(:)

    allocate(val(NDOF*NP), source = 0.0d0)
    allocate(vec(NDOF*NP,n_get_eigen), source = 0.0d0)
    allocate(is_bc(NDOF*NP), source = .false.)

    do i = 1, NDOF*NP
      if(is_bc_int(i) == 1) is_bc(i) = .true.
    enddo

    !> for monoMAT
    monolis%MAT%N = N
    monolis%MAT%NP = NP
    monolis%MAT%NZ = NZ
    monolis%MAT%NDOF = NDOF
    monolis%MAT%A => A
    !monolis%MAT%X => X
    !monolis%MAT%B => B
    allocate(monolis%MAT%X(NDOF*NP), source = 0.0d0)
    allocate(monolis%MAT%B(NDOF*NP), source = 0.0d0)
    monolis%MAT%index => index
    monolis%MAT%item => item

    !> for monoCOM
    monoliS%COM%internal_nnode = N
    monoliS%COM%myrank = myrank
    monoliS%COM%comm = comm
    monoliS%COM%commsize = commsize
    monoliS%COM%recv_n_neib = recv_n_neib
    monoliS%COM%recv_neib_pe => recv_neib_pe
    monoliS%COM%recv_index => recv_index
    monoliS%COM%recv_item => recv_item
    monoliS%COM%send_n_neib = send_n_neib
    monoliS%COM%send_neib_pe => send_neib_pe
    monoliS%COM%send_index => send_index
    monoliS%COM%send_item => send_item

    !> for monoPRM
    monolis%PRM%method = method
    monolis%PRM%precond = precond
    monolis%PRM%maxiter = maxiter
    monolis%PRM%tol = tol

    monolis%PRM%is_scaling    = .false.
    monolis%PRM%is_reordering = .false.
    monolis%PRM%is_init_x     = .true.
    monolis%PRM%is_debug      = .false.
    monolis%PRM%is_check_diag = .true.
    monolis%PRM%is_measurement= .false.
    monolis%PRM%show_iterlog  = .false.
    monolis%PRM%show_time     = .false.
    monolis%PRM%show_summary  = .false.

    if(iterlog == 1) monolis%PRM%show_iterlog  = .true.
    if(timelog == 1) monolis%PRM%show_time     = .true.
    if(timelog_statistics == 1) monolis%PRM%show_time_statistics = .true.
    if(summary == 1) monolis%PRM%show_summary  = .true.
    if(is_check_diag == 0) monolis%PRM%is_check_diag= .false.
    if(is_init_x == 0) monolis%PRM%is_init_x= .false.
    if(is_measurement == 1) monolis%PRM%is_measurement= .true.

    call monolis_eigen_inverted_standard_lanczos(monolis, n_get_eigen, ths, eigen_maxiter, val, vec, is_bc)

    do i = 1, n_get_eigen
      eigen_value(i) = val(i)
    enddo

    do j = 1, n_get_eigen
      do i = 1, NDOF*NP
        eigen_mode_tmp(NDOF*NP*(j-1) + i) = vec(i,j)
      enddo
    enddo

    curiter = monolis%PRM%curiter
    curresid = monolis%PRM%curresid
    time(1) = monolis%PRM%tsol
    time(2) = monolis%PRM%tprep
    time(3) = monolis%PRM%tspmv
    time(4) = monolis%PRM%tdotp
    time(5) = monolis%PRM%tprec
    time(6) = monolis%PRM%tcomm_dotp
    time(7) = monolis%PRM%tcomm_spmv

    deallocate(monolis%MAT%X)
    deallocate(monolis%MAT%B)
  end subroutine monolis_eigen_inverted_standard_lanczos_c

  !> c interface
  subroutine monolis_eigen_standard_lanczos_c(N, NP, NZ, NDOF, A, index, item, &
    myrank, comm, commsize, &
    recv_n_neib, recv_nitem, recv_neib_pe, recv_index, recv_item, &
    send_n_neib, send_nitem, send_neib_pe, send_index, send_item, &
    method, precond, maxiter, tol, &
    iterlog, timelog_statistics, timelog, summary, &
    is_check_diag, is_measurement, is_init_x, curiter, curresid, time, &
    n_get_eigen, ths, eigen_maxiter, eigen_value, eigen_mode_tmp, is_bc_int) &
    & bind(c, name = "monolis_eigen_standard_lanczos_c_main")
    implicit none
    type(monolis_structure) :: monolis
    integer(c_int), intent(in), value :: N, NP, NZ, NDOF
    integer(c_int), intent(in), value :: myrank, comm, commsize
    integer(c_int), intent(in), value :: recv_n_neib, send_n_neib, recv_nitem, send_nitem
    integer(c_int), intent(in), value :: method, precond, maxiter
    integer(c_int), intent(in), value :: iterlog, timelog, timelog_statistics, summary
    integer(c_int), intent(in), value :: is_check_diag, is_measurement, is_init_x
    integer(c_int), intent(in), value :: n_get_eigen
    integer(c_int), intent(in), value :: eigen_maxiter
    integer(c_int), intent(in), target :: index(0:NP)
    integer(c_int), intent(in), target :: item(NZ)
    integer(c_int), intent(in), target :: recv_neib_pe(recv_n_neib)
    integer(c_int), intent(in), target :: recv_index(0:recv_n_neib), recv_item(recv_nitem)
    integer(c_int), intent(in), target :: send_neib_pe(send_n_neib)
    integer(c_int), intent(in), target :: send_index(0:send_n_neib), send_item(send_nitem)
    integer(c_int), intent(in), target :: is_bc_int(NDOF*NP)
    integer(c_int), intent(out), target :: curiter
    real(c_double), intent(in), value :: tol, ths
    real(c_double), intent(in), target :: A(NDOF*NDOF*NZ)
    real(c_double), intent(inout), target :: eigen_value(NDOF*NP)
    real(c_double), intent(inout), target :: eigen_mode_tmp(NDOF*NP*n_get_eigen)
    real(c_double), intent(out), target :: curresid
    real(c_double), intent(out), target :: time(7)
    integer(kint) :: i, j
    real(kdouble), allocatable :: vec(:,:), val(:)
    logical, allocatable :: is_bc(:)

    allocate(val(NDOF*NP), source = 0.0d0)
    allocate(vec(NDOF*NP,n_get_eigen), source = 0.0d0)
    allocate(is_bc(NDOF*NP), source = .false.)

    do i = 1, NDOF*NP
      if(is_bc_int(i) == 1) is_bc(i) = .true.
    enddo

    !> for monoMAT
    monolis%MAT%N = N
    monolis%MAT%NP = NP
    monolis%MAT%NZ = NZ
    monolis%MAT%NDOF = NDOF
    monolis%MAT%A => A
    !monolis%MAT%X => X
    !monolis%MAT%B => B
    allocate(monolis%MAT%X(NDOF*NP), source = 0.0d0)
    allocate(monolis%MAT%B(NDOF*NP), source = 0.0d0)
    monolis%MAT%index => index
    monolis%MAT%item => item

    !> for monoCOM
    monoliS%COM%internal_nnode = N
    monoliS%COM%myrank = myrank
    monoliS%COM%comm = comm
    monoliS%COM%commsize = commsize
    monoliS%COM%recv_n_neib = recv_n_neib
    monoliS%COM%recv_neib_pe => recv_neib_pe
    monoliS%COM%recv_index => recv_index
    monoliS%COM%recv_item => recv_item
    monoliS%COM%send_n_neib = send_n_neib
    monoliS%COM%send_neib_pe => send_neib_pe
    monoliS%COM%send_index => send_index
    monoliS%COM%send_item => send_item

    !> for monoPRM
    monolis%PRM%method = method
    monolis%PRM%precond = precond
    monolis%PRM%maxiter = maxiter
    monolis%PRM%tol = tol

    monolis%PRM%is_scaling    = .false.
    monolis%PRM%is_reordering = .false.
    monolis%PRM%is_init_x     = .true.
    monolis%PRM%is_debug      = .false.
    monolis%PRM%is_check_diag = .true.
    monolis%PRM%is_measurement= .false.
    monolis%PRM%show_iterlog  = .false.
    monolis%PRM%show_time     = .false.
    monolis%PRM%show_summary  = .false.

    if(iterlog == 1) monolis%PRM%show_iterlog  = .true.
    if(timelog == 1) monolis%PRM%show_time     = .true.
    if(timelog_statistics == 1) monolis%PRM%show_time_statistics = .true.
    if(summary == 1) monolis%PRM%show_summary  = .true.
    if(is_check_diag == 0) monolis%PRM%is_check_diag= .false.
    if(is_init_x == 0) monolis%PRM%is_init_x= .false.
    if(is_measurement == 1) monolis%PRM%is_measurement= .true.

    call monolis_eigen_standard_lanczos(monolis, n_get_eigen, ths, eigen_maxiter, val, vec, is_bc)

    do i = 1, n_get_eigen
      eigen_value(i) = val(i)
    enddo

    do j = 1, n_get_eigen
      do i = 1, NDOF*NP
        eigen_mode_tmp(NDOF*NP*(j-1) + i) = vec(i,j)
      enddo
    enddo

    curiter = monolis%PRM%curiter
    curresid = monolis%PRM%curresid
    time(1) = monolis%PRM%tsol
    time(2) = monolis%PRM%tprep
    time(3) = monolis%PRM%tspmv
    time(4) = monolis%PRM%tdotp
    time(5) = monolis%PRM%tprec
    time(6) = monolis%PRM%tcomm_dotp
    time(7) = monolis%PRM%tcomm_spmv

    deallocate(monolis%MAT%X)
    deallocate(monolis%MAT%B)
  end subroutine monolis_eigen_standard_lanczos_c

end module mod_monolis_eigen_lanczos

  subroutine interface_monolis_eigen_inverted_standard_lanczos_ &
    & (monoPRM, monoCOM, monoMAT, NPNDOF, &
    & n_get_eigen, ths, maxiter, val, vec, is_bc)
    use mod_monolis_prm
    use mod_monolis_com
    use mod_monolis_mat
    use mod_monolis_eigen_lanczos
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kint) :: maxiter, n_get_eigen, NPNDOF
    real(kdouble) :: ths
    real(kdouble) :: vec(NPNDOF,n_get_eigen), val(n_get_eigen)
    logical :: is_bc(NPNDOF)

    call monolis_eigen_inverted_standard_lanczos_(monoPRM, monoCOM, monoMAT, &
      & n_get_eigen, ths, maxiter, val, vec, is_bc)
  end subroutine interface_monolis_eigen_inverted_standard_lanczos_
