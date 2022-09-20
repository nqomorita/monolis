module mod_monolis_matmat
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_linalg_util
  use mod_monolis_linalg_com
  use mod_monolis_util_debug
  implicit none

contains

  subroutine monolis_matmat(monoCOM, monoMAT, X, Y, tspmv, tcomm)
    implicit none
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    real(kdouble) :: X(:), Y(:)
    real(kdouble) :: t1, t2
    real(kdouble) :: tspmv, tcomm

#ifdef DEBUG
    call monolis_debug_header("monolis_matmat")
#endif
    t1 = monolis_get_time()

    !call monolis_update_pre_R(monoCOM, monoMAT%NDOF, X, tcomm)

    call monolis_matmat_nn(monoCOM, monoMAT, X, Y, monoMAT%NDOF)

    !call monolis_update_post_R(monoCOM, monoMAT%NDOF, X, tcomm)

    t2 = monolis_get_time()
    tspmv = tspmv + t2 - t1
  end subroutine monolis_matmat

  subroutine monolis_matmat_nn(monoCOM, monoMAT, X, Y, NDOF)
    implicit none
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kint) :: i, j, k, l, in, N, NDOF, NDOF2, jS, jE
    integer(kint), pointer :: index(:), item(:)
    real(kdouble) :: X(:), Y(:), XT(NDOF), YT(NDOF)
    real(kdouble), pointer :: A(:)

    N = monoMAT%N
    NDOF2 = NDOF*NDOF
    A => monoMAT%A
    index => monoMAT%index
    item  => monoMAT%item

    do i = 1, N
      YT = 0.0d0
      jS = index(i-1) + 1
      jE = index(i  )
      do j = jS, jE
        in = item(j)
        do k = 1, NDOF
          XT(k) = X(NDOF*(in-1)+k)
        enddo
        do k = 1, NDOF
          do l = 1, NDOF
            YT(k) = YT(k) + A(NDOF2*(j-1)+NDOF*(k-1)+l) * XT(l)
          enddo
        enddo
      enddo
      do k = 1, NDOF
        Y(NDOF*(i-1)+k) = YT(k)
      enddo
    enddo
  end subroutine monolis_matmat_nn

end module mod_monolis_matmat