module mod_monolis_matvec
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat
  use mod_monolis_linalg_util
  use mod_monolis_linalg_com
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

    do i=1,monoMAT%N*monoMAT%NDOF
      R(i) = B(i) - R(i)
    enddo
  end subroutine monolis_residual

  subroutine monolis_matvec(monoCOM, monoMAT, X, Y, tcomm)
    implicit none
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    real(kind=kdouble) :: X(:), Y(:)
    real(kind=kdouble), optional :: tcomm

    if(monoMAT%NDOF == 3)then
      call monolis_matvec_33(monoCOM, monoMAT, X, Y, tcomm)
    else
      call monolis_matvec_nn(monoCOM, monoMAT, X, Y, monoMAT%NDOF, tcomm)
    endif
  end subroutine monolis_matvec

  subroutine monolis_matvec_nn(monoCOM, monoMAT, X, Y, NDOF, tcomm)
    implicit none
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kind=kint) :: i, j, k, l, in, N, NDOF, NDOF2, jS, jE
    integer(kind=kint), pointer :: indexL(:), itemL(:)
    integer(kind=kint), pointer :: indexU(:), itemU(:)
    real(kind=kdouble) :: X(:), Y(:), XT(NDOF), YT(NDOF)
    real(kind=kdouble), pointer :: D(:), AU(:), AL(:)
    real(kind=kdouble) :: t1, t2
    real(kind=kdouble), optional :: tcomm

    N = monoMAT%N
    NDOF2 = NDOF*NDOF
    D  => monoMAT%D
    AU => monoMAT%AU
    AL => monoMAT%AL
    indexU => monoMAT%indexU
    indexL => monoMAT%indexL
    itemU  => monoMAT%itemU
    itemL  => monoMAT%itemL

    call monolis_update_R(monoCOM, NDOF, X, tcomm)

    do i = 1, N
      YT = 0.0d0
      do j = 1, NDOF
        XT(j) = X(NDOF*(i-1)+j)
      enddo
      do j = 1, NDOF
        do k = 1, NDOF
          YT(j) = YT(j) + D(NDOF2*(i-1)+NDOF*(j-1)+k) * XT(k)
        enddo
      enddo

      jS = indexL(i-1) + 1
      jE = indexL(i  )
      do j = jS, jE
        in = itemL(j)
        do k = 1, NDOF
          XT(k) = X(NDOF*(in-1)+k)
        enddo
        do k = 1, NDOF
          do l = 1, NDOF
            YT(k) = YT(k) + AL(NDOF2*(j-1)+NDOF*(k-1)+l) * XT(l)
          enddo
        enddo
      enddo

      jS = indexU(i-1) + 1
      jE = indexU(i  )
      do j = jS, jE
        in = itemU(j)
        do k = 1, NDOF
          XT(k) = X(NDOF*(in-1)+k)
        enddo
        do k = 1, NDOF
          do l = 1, NDOF
            YT(k) = YT(k) + AU(NDOF2*(j-1)+(k-1)*NDOF+l)*XT(l)
          enddo
        enddo
      enddo

      do k = 1,NDOF
        Y(NDOF*(i-1)+k) = YT(k)
      enddo
    enddo
  end subroutine monolis_matvec_nn

  subroutine monolis_matvec_33(monoCOM, monoMAT, X, Y, tcomm)
    implicit none
    type(monolis_com) :: monoCOM
    type(monolis_mat) :: monoMAT
    integer(kind=kint) :: i, j, in, N, jS, jE
    integer(kind=kint), pointer :: indexL(:), itemL(:)
    integer(kind=kint), pointer :: indexU(:), itemU(:)
    real(kind=kdouble) :: X1, X2, X3, Y1, Y2, Y3
    real(kind=kdouble) :: X(:), Y(:)
    real(kind=kdouble), pointer :: D(:), AU(:), AL(:)
    real(kind=kdouble) :: t1, t2
    real(kind=kdouble), optional :: tcomm

    N = monoMAT%N
    D  => monoMAT%D
    AU => monoMAT%AU
    AL => monoMAT%AL
    indexU => monoMAT%indexU
    indexL => monoMAT%indexL
    itemU  => monoMAT%itemU
    itemL  => monoMAT%itemL

    call monolis_update_R(monoCOM, monoMAT%NDOF, X, tcomm)

    do i = 1, N
      X1 = X(3*i-2)
      X2 = X(3*i-1)
      X3 = X(3*i  )
      Y1 = D(9*i-8)*X1 + D(9*i-7)*X2 + D(9*i-6)*X3
      Y2 = D(9*i-5)*X1 + D(9*i-4)*X2 + D(9*i-3)*X3
      Y3 = D(9*i-2)*X1 + D(9*i-1)*X2 + D(9*i  )*X3

      jS = indexL(i-1) + 1
      jE = indexL(i  )
      do j = jS, jE
        in = itemL(j)
        X1 = X(3*in-2)
        X2 = X(3*in-1)
        X3 = X(3*in  )
        Y1 = Y1 + AL(9*j-8)*X1 + AL(9*j-7)*X2 + AL(9*j-6)*X3
        Y2 = Y2 + AL(9*j-5)*X1 + AL(9*j-4)*X2 + AL(9*j-3)*X3
        Y3 = Y3 + AL(9*j-2)*X1 + AL(9*j-1)*X2 + AL(9*j  )*X3
      enddo

      jS = indexU(i-1) + 1
      jE = indexU(i  )
      do j = jS, jE
        in = itemU(j)
        X1 = X(3*in-2)
        X2 = X(3*in-1)
        X3 = X(3*in  )
        Y1 = Y1 + AU(9*j-8)*X1 + AU(9*j-7)*X2 + AU(9*j-6)*X3
        Y2 = Y2 + AU(9*j-5)*X1 + AU(9*j-4)*X2 + AU(9*j-3)*X3
        Y3 = Y3 + AU(9*j-2)*X1 + AU(9*j-1)*X2 + AU(9*j  )*X3
      enddo

      Y(3*i-2) = Y1
      Y(3*i-1) = Y2
      Y(3*i  ) = Y3
    enddo
  end subroutine monolis_matvec_33

end module mod_monolis_matvec