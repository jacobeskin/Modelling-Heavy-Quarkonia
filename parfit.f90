program parfit
  use qqModule
  
  implicit none

  double precision, allocatable :: u(:), E_m(:), E_c(:), x(:)
  double precision :: alpha, b, c, sigma, m, F
  integer :: k = 2500, i

  alpha = 0.1789d+0  !0.1951d+0 
  b = 0.2152d+0 !0.2569 
  c =  1.0d+0 !0.1d+0
  m = 4.18d+0  !1.27 
  sigma = 4.2d+0 !1.2 

  allocate(u(k), E_m(6), E_c(6), x(6))

  ! Measured energylevels
  E_m(1) = 9.3990d+0 !2.9834d+0
  E_m(2) = 9.4603d+0 !3.0969d+0
  E_m(3) = 9.8993d+0 !3.5252d+0
  E_m(4) = 9.8594d+0 !3.4148d+0
  E_m(5) = 9.8928d+0 !3.5107d+0
  E_m(6) = 9.9122d+0 !3.5562d+0

  call ln(alpha, b, c, m, sigma, k, E_m, E_c, F, x)

  print *,
  print *, "DONE!"
  print *,  F
  do i=1,6
     print *, E_c(i)
  end do
  print *,

end program parfit
