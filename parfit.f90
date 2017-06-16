program parfit
  use qqModule
  
  implicit none

  ! Fits the parameters to data

  double precision, allocatable :: u(:), E_m(:), E_c(:), x(:)
  double precision :: alpha, b, c, sigma, m, F
  integer :: k = 2500, i

  ! Starting values for Bottomonium
  !alpha = 0.1789d+0   
  !b = 0.2152d+0 
  !c =  1.0d+0 
  !m = 4.18d+0   
  !sigma = 4.2d+0  

  ! Starting values for Charmonium
  alpha = 0.1951d+0 
  b = 0.2569d+0 
  c =  0.1d+0
  m = 1.27d+0 
  sigma = 1.2d+0 

  allocate(u(k), E_m(6), E_c(6), x(5))

  ! Measured energylevels

  ! For Bottomonium
  !E_m(1) = 9.3990d+0 
  !E_m(2) = 9.4603d+0 
  !E_m(3) = 9.8993d+0 
  !E_m(4) = 9.8594d+0 
  !E_m(5) = 9.8928d+0 
  !E_m(6) = 9.9122d+0 

  ! For Charmonium
  E_m(1) = 2.9834d+0
  E_m(2) = 3.0969d+0
  E_m(3) = 3.5252d+0
  E_m(4) = 3.4148d+0
  E_m(5) = 3.5107d+0
  E_m(6) = 3.5562d+0

  call ln(alpha, b, c, m, sigma, k, E_m, E_c, F, x)

  print *,
  print *, "DONE!"
  print *,  F
  do i=1,6
     print *, E_c(i)
  end do
  print *, 
  print *, "Paramters:"
  do i=1,5
     print *, x(i)
  end do
  print *,

end program parfit
