module mod_monolis_fact_LU_nn
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_mat

  implicit none
  private
  public :: monolis_init_LU_inner_nn
  public :: monolis_fact_LU_inner_nn
  public :: monolis_solv_LU_inner_nn
  public :: monolis_clear_LU_inner_nn

  !integer(kint) :: N
  real(kdouble), allocatable :: Dlu0(:)
  real(kdouble), allocatable :: ALlu0(:)
  real(kdouble), allocatable :: AUlu0(:)
  integer(kint), allocatable :: inumFI1L(:)
  integer(kint), allocatable :: inumFI1U(:)
  integer(kint), allocatable :: FI1L(:)
  integer(kint), allocatable :: FI1U(:)

  logical, save :: INITIALIZED = .false.

contains

  subroutine monolis_init_LU_inner_nn(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat_LDU) :: monoMAT

  end subroutine monolis_init_LU_inner_nn

  subroutine monolis_clear_LU_inner_nn(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat_LDU) :: monoMAT

  end subroutine monolis_clear_LU_inner_nn

  subroutine monolis_fact_LU_inner_nn(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat_LDU) :: monoMAT

    call FORM_ILU0_nn &
      & (monoMAT%N, monoMAT%NDOF, monoMAT%N, monoMAT%NPL, monoMAT%NPU, &
      &  monoMAT%D, monoMAT%L, monoMAT%indexL, monoMAT%itemL, &
      &  monoMAT%U, monoMAT%indexU, monoMAT%itemU, 1.0d0, 1.0d0)
  end subroutine monolis_fact_LU_inner_nn

  subroutine monolis_solv_LU_inner_nn(monoPRM, monoCOM, monoMAT)
    implicit none
    type(monolis_prm) :: monoPRM
    type(monolis_com) :: monoCOM
    type(monolis_mat_LDU) :: monoMAT

    monoMAT%X = monoMAT%B
    call BILU_nn_apply(monoMAT%X, monoMAT%N, monoMAT%ndof)
  end subroutine monolis_solv_LU_inner_nn

  subroutine BILU_nn_apply(WW, N, NDOF)
    implicit none
    real(kdouble), intent(inout) :: WW(:)
    integer(kint) :: i, j, ii, ij, isL, ieL, isU, ieU, k, NDOF, N
    real(kdouble) :: SW(NDOF), X(NDOF)

    do i= 1, N
      do ii = 1, NDOF
        SW(ii)= WW(NDOF*(i-1)+ii)
      end do
      isL= inumFI1L(i-1)+1
      ieL= inumFI1L(i)
      do j= isL, ieL
        k= FI1L(j)
        do ii = 1, NDOF
          X(ii)= WW(NDOF*(k-1)+ii)
        end do
        do ii = 1, NDOF
          do ij = 1, NDOF
            SW(ii)= SW(ii) - ALlu0(NDOF*NDOF*(j-1)+NDOF*(ii-1)+ij)*X(ij)
          end do
        end do
      enddo

      X= SW
      do ii=2,NDOF
        do ij = 1,ii-1
          X(ii)=X(ii)-Dlu0(NDOF*NDOF*(i-1)+NDOF*(ii-1)+ij )*X(ij)
        end do
      end do
      do ii=NDOF,1,-1
        do ij = NDOF,ii+1,-1
          X(ii)=X(ii)-Dlu0(NDOF*NDOF*(i-1)+NDOF*(ii-1)+ij )*X(ij)
        end do
        X(ii)=Dlu0(NDOF*NDOF*(i-1)+(NDOF+1)*(ii-1)+1 )*X(ii)
      end do
      do ii = 1, NDOF
        WW(NDOF*(i-1)+ii)=X(ii)
      end do
    enddo

    do i= N, 1, -1
      isU= inumFI1U(i-1) + 1
      ieU= inumFI1U(i)
      SW= 0.d0

      do j= ieU, isU, -1
        k= FI1U(j)
        do ii = 1, NDOF
          X(ii)= WW(NDOF*(k-1)+ii)
        end do
        do ii = 1, NDOF
          do ij = 1, NDOF
            SW(ii)= SW(ii) + AUlu0(NDOF*NDOF*(j-1)+NDOF*(ii-1)+ij)*X(ij)
          end do
        end do
      enddo
      X= SW
      do ii=2,NDOF
        do ij = 1,ii-1
          X(ii)=X(ii)-Dlu0(NDOF*NDOF*(i-1)+NDOF*(ii-1)+ij )*X(ij)
        end do
      end do
      do ii=NDOF,1,-1
        do ij = NDOF,ii+1,-1
          X(ii)=X(ii)-Dlu0(NDOF*NDOF*(i-1)+NDOF*(ii-1)+ij )*X(ij)
        end do
        X(ii)=Dlu0(NDOF*NDOF*(i-1)+(NDOF+1)*(ii-1)+1 )*X(ii)
      end do
      do ii = 1, NDOF
        WW(NDOF*(i-1)+ii)= WW(NDOF*(i-1)+ii)-X(ii)
      end do
    enddo
  end subroutine BILU_nn_apply

  subroutine FORM_ILU0_nn                                   &
      &   (N, NDOF, NP, NPL, NPU, D, AL, INL, IAL, AU, INU, IAU, &
      &    SIGMA, SIGMA_DIAG)
    implicit none
    integer(kint ), intent(in):: N, NDOF, NP, NPU, NPL
    real   (kdouble), intent(in):: SIGMA, SIGMA_DIAG

    real(kdouble), dimension(NDOF*NDOF*NPL), intent(in):: AL
    real(kdouble), dimension(NDOF*NDOF*NPU), intent(in):: AU
    real(kdouble), dimension(NDOF*NDOF*NP ), intent(in):: D

    integer(kint), dimension(0:NP) ,intent(in) :: INU, INL
    integer(kint), dimension(  NPL),intent(in) :: IAL
    integer(kint), dimension(  NPU),intent(in) :: IAU

    integer(kint), dimension(:), allocatable :: IW1, IW2
    real (kdouble),  dimension(NDOF,NDOF) :: RHS_Aij, DkINV, Aik, Akj
    integer(kint) :: i,jj,jj1,ij0,kk,kk1,NDOF2
    integer(kint) :: j,k,ii,ij

    NDOF2=NDOF*NDOF
    allocate(IW1(NP), IW2(NP))
    allocate(Dlu0(NDOF2*NP), ALlu0(NDOF2*NPL), AUlu0(NDOF2*NPU))
    allocate(inumFI1L(0:NP), inumFI1U(0:NP), FI1L(NPL), FI1U(NPU))

    do i=1,NDOF2*NP
      Dlu0(i) = D(i)
    end do
    do i=1,NDOF2*NPL
      ALlu0(i) = AL(i)
    end do
    do i=1,NDOF2*NPU
      AUlu0(i) = AU(i)
    end do
    do i=0,NP
      inumFI1L(i) = INL(i)
      inumFI1U(i) = INU(i)
    end do
    do i=1,NPL
      FI1L(i) = IAL(i)
    end do
    do i=1,NPU
      FI1U(i) = IAU(i)
    end do

    do i=1,NP
      do ii=1,NDOF
        Dlu0(NDOF2*(i-1)+(ii-1)*NDOF+ii)=Dlu0(NDOF2*(i-1)+(ii-1)*NDOF+ii)*SIGMA_DIAG
      end do
    enddo

    i = 1
    call ILU1aNN (DkINV, Dlu0(NDOF2*(i-1)+1:NDOF2*NDOF2),NDOF)
    do ii=1,NDOF
      do ij=1,NDOF
        Dlu0(NDOF2*(i-1)+NDOF*(ii-1)+ij)= DkINV(ii,ij)
      end do
    end do

    do i= 2, NP
      IW1= 0
      IW2= 0

      do k= inumFI1L(i-1)+1, inumFI1L(i)
        IW1(FI1L(k))= k
      enddo

      do k= inumFI1U(i-1)+1, inumFI1U(i)
        IW2(FI1U(k))= k
      enddo

      do kk= INL(i-1)+1, INL(i)
        k= IAL(kk)
        do ii=1,NDOF
          do ij=1,NDOF
            DkINV(ii,ij) = Dlu0(NDOF2*(k-1)+NDOF*(ii-1)+ij)
          end do
        end do
        do ii=1,NDOF
          do ij=1,NDOF
            Aik(ii,ij) = ALlu0(NDOF2*(kk-1)+NDOF*(ii-1)+ij)
          end do
        end do

        do jj= INU(k-1)+1, INU(k)
          j= IAU(jj)
          if (IW1(j).eq.0.and.IW2(j).eq.0) cycle
          do ii=1,NDOF
            do ij=1,NDOF
              Akj(ii,ij) = AUlu0(NDOF2*(jj-1)+NDOF*(ii-1)+ij)
            end do
          end do

          call ILU1bNN (RHS_Aij, DkINV, Aik, Akj,NDOF)

          if (j.eq.i) then
            do ii=1,NDOF
              do ij=1,NDOF
                Dlu0(NDOF2*(i-1)+NDOF*(ii-1)+ij) = Dlu0(NDOF2*(i-1)+NDOF*(ii-1)+ij) - RHS_Aij(ii,ij)
              end do
            end do
          endif

          if (j.lt.i) then
            ij0= IW1(j)
            do ii=1,NDOF
              do ij=1,NDOF
                ALlu0(NDOF2*(ij0-1)+NDOF*(ii-1)+ij) = ALlu0(NDOF2*(ij0-1)+NDOF*(ii-1)+ij) - RHS_Aij(ii,ij)
              end do
            end do
          endif

          if (j.gt.i) then
            ij0= IW2(j)
            do ii=1,NDOF
              do ij=1,NDOF
                AUlu0(NDOF2*(ij0-1)+NDOF*(ii-1)+ij) = AUlu0(NDOF2*(ij0-1)+NDOF*(ii-1)+ij) - RHS_Aij(ii,ij)
              end do
            end do
          endif

        enddo
      enddo
      call ILU1aNN (DkINV, Dlu0(NDOF2*(i-1)+1:NDOF2*NDOF2),NDOF)

      do ii=1,NDOF
        do ij=1,NDOF
          Dlu0(NDOF2*(i-1)+NDOF*(ii-1)+ij) = DkINV(ii,ij)
        end do
      end do
    enddo

    deallocate (IW1, IW2)
  end subroutine FORM_ILU0_nn

  subroutine ILU1aNN (ALU, D, NDOF)
    implicit none
    real(kdouble) :: ALU(NDOF,NDOF), D(NDOF*NDOF), PW(NDOF)
    real(kdouble) :: D11,D12,D13,D21,D22,D23,D31,D32,D33
    integer(kint) :: NDOF, i,j,k

    do i = 1, NDOF
      do j = 1, NDOF
        ALU(i,j) = D(NDOF*(i-1)+j)
      end do
    end do

    do k= 1, NDOF
      if (ALU(k,k) == 0.d0) then
        stop 'ERROR: Divide by zero in ILU setup'
      endif
      ALU(k,k)= 1.d0/ALU(k,k)
      do i= k+1, NDOF
        ALU(i,k)= ALU(i,k) * ALU(k,k)
        do j= k+1, NDOF
          PW(j)= ALU(i,j) - ALU(i,k)*ALU(k,j)
        enddo
        do j= k+1, NDOF
          ALU(i,j)= PW(j)
        enddo
      enddo
    enddo

    return
  end subroutine ILU1aNN

  subroutine ILU1bNN (RHS_Aij, DkINV, Aik, Akj, NDOF)
    implicit none
    real(kdouble) :: RHS_Aij(NDOF,NDOF), DkINV(NDOF,NDOF), Aik(NDOF,NDOF), Akj(NDOF,NDOF)
    real(kdouble) :: X(NDOF)
    integer(kint) :: NDOF,i,j,k

    do k=1,NDOF
      X(1:NDOF)= Akj(1:NDOF,k)
      do i=2,NDOF
        do j = 1,i-1
          X(i)=X(i)-DkINV(i,j)*X(j)
        end do
      end do
      do i=NDOF,1,-1
        do j = NDOF,i+1,-1
          X(i)=X(i)-DkINV(i,j)*X(j)
        end do
        X(i)=DkINV(i,i)*X(i)
      end do
      RHS_Aij(:,k)=0.0d0
      do i=1,NDOF
        do j=1,NDOF
          RHS_Aij(i,k) = RHS_Aij(i,k)+Aik(i,j)*X(j)
        end do
      end do
    end do
    return
  end subroutine ILU1bNN

end module mod_monolis_fact_LU_nn
