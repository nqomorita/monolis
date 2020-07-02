module mod_monolis_matvec
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_linalg_util
  use mod_monolis_linalg_com
  use mod_monolis_util
  implicit none

contains

  subroutine monolis_residual(monoCOM, monoMAT, X, B, R, tcomm)
    implicit none
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kind=kint) :: i
    real(kind=kdouble) :: X(:), B(:), R(:)
    real(kind=kdouble), optional :: tcomm

    call monolis_matvec(monoCOM, monoMAT, X, R, tcomm)

    do i = 1, monoMAT%N*monoMAT%NDOF
      R(i) = B(i) - R(i)
    enddo
  end subroutine monolis_residual

  subroutine monolis_matvec_product(monolis, X, Y)
    implicit none
    type(monolis_structure) :: monolis
    real(kdouble) :: X(:), Y(:)

    call monolis_matvec(monolis%COM, monolis%MAT, X, Y)
  end subroutine monolis_matvec_product

  subroutine monolis_matvec(monoCOM, monoMAT, X, Y, tspmv)
    implicit none
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    real(kind=kdouble) :: X(:), Y(:)
    real(kind=kdouble) :: t1, t2
    real(kind=kdouble), optional :: tspmv

#ifdef DEBUG
    call monolis_debug_header("monolis_matvec")
#endif
    t1 = monolis_get_time()

    call monolis_update_pre_R(monoCOM, monoMAT%NDOF, X)

    if(monoMAT%NDOF == 3)then
      call monolis_matvec_33(monoCOM, monoMAT, X, Y)
    elseif(monoMAT%NDOF == 1)then
      call monolis_matvec_11(monoCOM, monoMAT, X, Y)
    else
      call monolis_matvec_nn(monoCOM, monoMAT, X, Y, monoMAT%NDOF)
    endif

    call monolis_update_post_R(monoCOM, monoMAT%NDOF, X)

    t2 = monolis_get_time()
    if(present(tspmv)) tspmv = tspmv + t2 - t1
  end subroutine monolis_matvec

  subroutine monolis_matvec_nn(monoCOM, monoMAT, X, Y, NDOF, tcomm)
    implicit none
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kind=kint) :: i, j, k, l, in, N, NDOF, NDOF2, jS, jE
    integer(kind=kint), pointer :: index(:), item(:)
    real(kind=kdouble) :: X(:), Y(:), XT(NDOF), YT(NDOF)
    real(kind=kdouble), pointer :: A(:)
    real(kind=kdouble) :: t1, t2
    real(kind=kdouble), optional :: tcomm

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
  end subroutine monolis_matvec_nn

  subroutine monolis_matvec_11(monoCOM, monoMAT, X, Y, tcomm)
    implicit none
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kind=kint) :: i, j, in, N, jS, jE
    integer(kind=kint), pointer :: index(:), item(:)
    real(kind=kdouble) :: Y1
    real(kind=kdouble) :: X(:), Y(:)
    real(kind=kdouble), pointer :: A(:)
    real(kind=kdouble) :: t1, t2
    real(kind=kdouble), optional :: tcomm

    N = monoMAT%N
    A => monoMAT%A
    index => monoMAT%index
    item  => monoMAT%item

    do i = 1, N
      Y1 = 0.0d0
      jS = index(i-1) + 1
      jE = index(i  )
      do j = jS, jE
        in = item(j)
        Y1 = Y1 + A(j)*X(in)
      enddo
      Y(i) = Y1
    enddo
  end subroutine monolis_matvec_11

  subroutine monolis_matvec_33(monoCOM, monoMAT, X, Y, tcomm)
    implicit none
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kind=kint) :: i, j, in, N, jS, jE
    integer(kind=kint), pointer :: index(:), item(:)
    real(kind=kdouble) :: X1, X2, X3, Y1, Y2, Y3
    real(kind=kdouble) :: X(:), Y(:)
    real(kind=kdouble), pointer :: A(:)
    real(kind=kdouble) :: t1, t2
    real(kind=kdouble), optional :: tcomm

    N = monoMAT%N
    A => monoMAT%A
    index => monoMAT%index
    item  => monoMAT%item

    do i = 1, N
      Y1 = 0.0d0
      Y2 = 0.0d0
      Y3 = 0.0d0
      jS = index(i-1) + 1
      jE = index(i  )
      do j = jS, jE
        in = item(j)
        X1 = X(3*in-2)
        X2 = X(3*in-1)
        X3 = X(3*in  )
        Y1 = Y1 + A(9*j-8)*X1 + A(9*j-7)*X2 + A(9*j-6)*X3
        Y2 = Y2 + A(9*j-5)*X1 + A(9*j-4)*X2 + A(9*j-3)*X3
        Y3 = Y3 + A(9*j-2)*X1 + A(9*j-1)*X2 + A(9*j  )*X3
      enddo
      Y(3*i-2) = Y1
      Y(3*i-1) = Y2
      Y(3*i  ) = Y3
    enddo
  end subroutine monolis_matvec_33

end module mod_monolis_matvec