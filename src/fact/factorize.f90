module mod_monolis_fact_factorize
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc

  implicit none

contains

  subroutine monolis_matrix_factorize_mf(monoTREE, n_super_node, super_node_id, super_node_size, &
    & fact_array, fact_array_index, front_size, add_location)
    implicit none
    type(monolis_mat) :: monoTREE
    integer(kint) :: n_super_node
    integer(kint) :: super_node_id(:)
    integer(kint) :: super_node_size(:)
    real(kdouble) :: fact_array(:)
    integer(kint) :: fact_array_index(:)
    integer(kint) :: front_size(:)
    integer(kint) :: add_location(:)
    integer(kint) :: N
    integer(kint) :: i, iS, iE, in, jn, j, k, n_fact
    integer(kint) :: uS
    integer(kint) :: update_size, n_update_size, snode_size

!write(*,*)"monoTREE%CSR%indexU", monoTREE%SCSR%indexU
!write(*,*)"monoTREE%CSR%itemU ", monoTREE%SCSR%itemU
!call sleep(1)
!write(*,*)"fact_array_index   ", fact_array_index
!write(*,*)"n_super_node   ", n_super_node
!write(*,*)"super_node_id  ", super_node_id
!write(*,*)"fact_array_index  ", fact_array_index

    do k = 1, n_super_node
!write(*,*)i
      !> factorization
      iS = fact_array_index(k) + 1
      iE = fact_array_index(k + 1)
      n_fact = front_size(k)
      snode_size = super_node_size(k)
      !n_fact = iE - iS + 1
!write(*,*)iS, iE, n_fact, snode_size
      !call monolis_matrix_update_LDLt(n_fact, fact_array(iS:iE))
!write(*,*)"k", k, snode_size
!write(*,"(1p10e12.3)")fact_array(iS:iE)
!write(*,*)"start", snode_size, n_fact, iS, iE
!write(*,"(1p10e12.3)")fact_array(iS:iE)
      call monolis_matrix_n_step_update_LDLt(snode_size, n_fact, fact_array(iS:iE))
!write(*,*)"end"

      !> update
      uS = fact_array_index(k)
      do i = 1, snode_size
        uS = uS + n_fact + 1 - i
      enddo

      n_update_size = 0
      do i = 1, n_fact
        if(i <= snode_size) cycle
        n_update_size = n_update_size + n_fact + 1 - i
      enddo

!write(*,*)"fact_array_index(k)", fact_array_index(k)
!write(*,*)"snode_size", snode_size
!write(*,*)"n_fact", n_fact
!write(*,*)"uS", uS
!write(*,*)"n_update_size", n_update_size
!write(*,*)"add_location"
!write(*,"(20i4)")add_location

!write(*,*)"n_update_size", n_update_size
      do j = 1, n_update_size
        in = uS + j
        jn = add_location(in)
!write(*,*)jn, in
        fact_array(jn) = fact_array(jn) + fact_array(in)
      enddo
    enddo

!write(*,*)"fact_array"
!write(*,"(1p10e12.3)")fact_array
!call sleep(1)
  end subroutine monolis_matrix_factorize_mf

  subroutine monolis_matrix_update_LDLt(N, A)
    implicit none
    integer(kint) :: N, j, k, in, next
    real(kdouble) :: A(:)
    real(kdouble) :: inv, al, au

    inv = 1.0d0/A(1)
    A(1) = inv
    next = N + 1
    in = 0
    do j = 1, N - 1
      al = A(1 + j)
      do k = j, N - 1
        au = A(1 + k)
        A(next + in) = A(next + in) - al*au*inv
        in = in + 1
      enddo
    enddo
  end subroutine monolis_matrix_update_LDLt

  subroutine monolis_matrix_n_step_update_LDLt(nstep, N, A)
    implicit none
    integer(kint) :: nstep
    integer(kint) :: N, i, j, k, in, next, index
    real(kdouble) :: A(:)
    real(kdouble) :: inv, al, au

!write(*,*)"nstep", nstep, N
    index = 1
    do i = 1, nstep
!write(*,*)"index", index, A(index)
      inv = 1.0d0/A(index)
      A(index) = inv
      next = index + N - i + 1
      in = 0
      do j = 1, N - i
        al = A(index + j)
        do k = j, N - i
          au = A(index + k)
!write(*,*)next + in, index + k, index + j
          A(next + in) = A(next + in) - al*au*inv
          in = in + 1
        enddo
      enddo
      index = next
    enddo
  end subroutine monolis_matrix_n_step_update_LDLt

  subroutine monolis_matrix_copy_lu_factor(monoTREE, n_super_node, super_node_id, super_node_size, &
    & fact_array, fact_array_index)
    implicit none
    type(monolis_mat) :: monoTREE
    integer(kint) :: n_super_node
    integer(kint) :: super_node_id(:)
    integer(kint) :: super_node_size(:)
    real(kdouble) :: fact_array(:)
    integer(kint) :: fact_array_index(:)
    integer(kint) :: N
    integer(kint) :: i, j, k, m, iS, in, jn, kn, n_fact

    in = 1
    do i = 1, n_super_node
      iS = fact_array_index(i)
      m = 1
      do j = 1, super_node_size(i)
        if(j == 1)then
          jn = super_node_id(i)
        else
          k = monoTREE%SCSR%indexU(jn) + 2 
          jn = monoTREE%SCSR%itemU(k)
        endif

        n_fact = monoTREE%SCSR%indexU(jn + 1) - monoTREE%SCSR%indexU(jn)
        do k = 1, n_fact
          monoTREE%R%A(in) = fact_array(iS + m)
          m = m + 1
          in = in + 1
        enddo
      enddo
    enddo
!write(*,*)"monoTREE%R%A"
!write(*,"(1p10e12.3)")monoTREE%R%A
!call sleep(1)
  end subroutine monolis_matrix_copy_lu_factor
end module mod_monolis_fact_factorize
