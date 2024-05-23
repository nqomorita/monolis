module mod_monolis_fact_factorize
  use mod_monolis_utils
  use mod_monolis_def_mat
  use mod_monolis_def_struc

  implicit none

contains

  subroutine monolis_matrix_factorize_mf(monoTREE, fact_array, fact_array_index, add_location)
    implicit none
    type(monolis_mat) :: monoTREE
    real(kdouble) :: fact_array(:)
    integer(kint) :: fact_array_index(:)
    integer(kint) :: add_location(:)
    integer(kint) :: N
    integer(kint) :: i, iS, iE, in, jn, j, n_fact

    N = monoTREE%N

    do i = 1, N
      !> factorization
      iS = fact_array_index(i) + 1
      iE = fact_array_index(i + 1)
      n_fact = monoTREE%SCSR%indexU(i + 1) - monoTREE%SCSR%indexU(i)
      call monolis_matrix_update_LDLt(n_fact, fact_array(iS:iE))

      !> update
      do j = 1, (n_fact - 1)*(n_fact - 2)/2
        in = iS + n_fact + j - 1
        jn = add_location(in)
        fact_array(jn) = fact_array(jn) + fact_array(in)
      enddo
    enddo
  end subroutine monolis_matrix_factorize_mf

  subroutine monolis_matrix_update_LDLt(N, A)
    implicit none
    integer(kint) :: N, j, k, in, next
    real(kdouble) :: A(:)
    real(kdouble) :: inv, al, au

    inv = 1.0d0/dsqrt(A(1))
    !inv = 1.0d0/A(1)
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

  subroutine monolis_solve_dense_LDLt(N, A)
    implicit none
    integer(kint) :: N, i, j, k, in, index, next
    real(kdouble) :: A(:)
    real(kdouble) :: inv, al, au

    index = 1
    do i = 1, N
      inv = 1.0d0/dsqrt(A(index))
      !inv = 1.0d0/A(index)
      A(index) = inv
      next = index + N - i + 1
      in = 0
      do j = 1, N - i
        al = A(index + j)
        do k = j, N - i
          au = A(index + k)
          A(next + in) = A(next + in) - al*au*inv
          in = in + 1
        enddo
      enddo
      index = next
    enddo
  end subroutine monolis_solve_dense_LDLt

  subroutine monolis_matrix_copy_lu_factor(monoTREE, fact_array, fact_array_index)
    implicit none
    type(monolis_mat) :: monoTREE
    real(kdouble) :: fact_array(:)
    integer(kint) :: fact_array_index(:)
    integer(kint) :: N
    integer(kint) :: i, iS, in, j, n_fact

    N = monoTREE%N

    in = 1
    do i = 1, N
      n_fact = monoTREE%SCSR%indexU(i + 1) - monoTREE%SCSR%indexU(i)
      iS = fact_array_index(i)
      do j = 1, n_fact
        monoTREE%R%A(in) = fact_array(iS + j)
        in = in + 1
      enddo
    enddo
  end subroutine monolis_matrix_copy_lu_factor
end module mod_monolis_fact_factorize
