program main
  use mod_monolis
  implicit none
  integer(kint) :: max_iter
  real(kdouble) :: tol, residual
  real(kdouble) :: A(8,5)
  real(kdouble) :: b(8)
  real(kdouble) :: x(5)

  A = 0.0d0

  A(1,1) = 1.0d0; A(1,2) = 0.0d0; A(1,3) = 1.0d0; A(1,4) = 0.0d0; A(1,5) = 0.0d0; 
  A(2,1) = 1.0d0; A(2,2) = 2.0d0; A(2,3) = 0.0d0; A(2,4) = 0.0d0; A(2,5) = 0.0d0; 
  A(3,1) = 2.0d0; A(3,2) = 0.0d0; A(3,3) = 0.0d0; A(3,4) = 1.0d0; A(3,5) = 0.0d0; 
  A(4,1) = 1.0d0; A(4,2) = 0.0d0; A(4,3) = 4.0d0; A(4,4) = 0.0d0; A(4,5) = 0.0d0; 
  A(5,1) = 3.0d0; A(5,2) = 0.0d0; A(5,3) = 0.0d0; A(5,4) = 0.0d0; A(5,5) = 1.0d0; 
  A(6,1) = 1.0d0; A(6,2) = 3.0d0; A(6,3) = 0.0d0; A(6,4) = 0.0d0; A(6,5) = 0.0d0; 
  A(7,1) = 4.0d0; A(7,2) = 0.0d0; A(7,3) = 1.0d0; A(7,4) = 2.0d0; A(7,5) = 0.0d0; 
  A(8,1) = 0.0d0; A(8,2) = 1.0d0; A(8,3) = 0.0d0; A(8,4) = 1.0d0; A(8,5) = 1.0d0; 

  b(1) =-1.0d0
  b(2) =-1.0d0
  b(3) = 2.0d0
  b(4) = 2.0d0
  b(5) =-3.0d0
  b(6) = 3.0d0
  b(7) = 4.0d0
  b(8) = 4.0d0

  max_iter = 10
  tol = 1.0d-6

  call monolis_optimize_nnls(A, b, x, 8, 5, max_iter, tol, residual)

  write(*,*)"x", x
  write(*,*)"residual", residual
end program main
