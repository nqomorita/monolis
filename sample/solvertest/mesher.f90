program main
  use mod_monolis
  implicit none
  integer(4) :: nnode, nelem, i
  integer(4), allocatable :: elem(:,:)
  real(kdouble), allocatable :: coef(:)
  character :: filename*128

  call monolis_global_initialize()

  call monolis_set_debug(.true.)

  call monolis_get_mtx_arg(filename)
  write(*,*) "filename: ", trim(filename)

  call monolis_input_mtx(filename, nnode, nelem, elem, coef)
  !call monolis_input_mtx_complex(filename, nnode, nelem, elem, coef)

  open(20, file = "node.dat", status = "replace")
    write(20,"(i8,i8)")nnode
    do i = 1, nnode
      write(20,*)"0., 0., 0."
    enddo
  close(20)

  open(20, file = "elem.dat", status = "replace")
    write(20,"(i8,i8)")nelem, 2
    do i = 1, nelem
      write(20,"(2i10)")elem(1,i), elem(2,i)
    enddo
  close(20)

  open(20, file = "coef.dat", status = "replace")
    write(20,"(i8,i8)")nelem
    do i = 1, nelem
      write(20,"(1p2e23.15)")coef(i)!, coef(2,i)
    enddo
  close(20)

  !call check_symmetric_C(nelem, elem, coef, 1)

  call monolis_global_finalize()

  contains

  subroutine check_symmetric_C(nelem, elem, coef, fix)  !非対称行列だった場合、fix = 1 入力で対角行列に修正
    implicit none
    integer(4) :: i, j, nelem, elem(:,:), fix, count_real, count_imag
    real(kdouble) :: coef(:,:)
    character :: onoff*3

    if(fix == 1)then
      onoff = "ON"
    else
      onoff = "OFF"
    endif
    count_real = 0
    count_imag = 0
    write(*,*)"check_symmetric  fix = ", trim(onoff)
    do i = 1, nelem
      do j = i+1, nelem 
        if(elem(1,i) == elem(2,i))then  !対角成分
          exit
        elseif(elem(1,i) == elem(2,j) .AND. elem(2,i) == elem(1,j))then
          if(coef(1,i) /= coef(1,j))then
            count_real = count_real + 1
            if(fix == 1)then
              ! coef(1,i)=(coef(1,i)+coef(1,j))/2
              coef(1,j) = coef(1,i)
            else
              write(*,*)"real part", elem(1,i), elem(2,i), "data", coef(1,i), coef(1,j)
            endif
          endif
          if(coef(2,i) /= coef(2,j))then
            count_imag = count_imag + 1
            if(fix == 1)then
              ! coef(2,i)=(coef(2,i)+coef(2,j))/2
              coef(2,j) = coef(2,i)
            else
              write(*,*)"imag part", elem(1,i), elem(2,i), "data", coef(2,i), coef(2,j)
            endif
          endif
        endif
      enddo
    enddo
    if(count_real > 0 .OR. count_imag > 0)then
      write(*,*)"mismatch num (real, imag)", count_real, count_imag
    else
      write(*,*)"input is symmetric."
    endif
  end subroutine check_symmetric_C

end program main
