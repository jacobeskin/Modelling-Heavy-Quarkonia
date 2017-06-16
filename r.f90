program r
  use qqModule

  implicit none

  ! Program for calculating the sqrt(<r**2>) and <v**2> for quarkonium state.
  ! Input parameters are:
  !
  ! file        filename of the wavefunction
  ! k           length of the wavefunction
  ! alpha       strong coupling constant
  ! b           linear constant
  ! m           quark mass
  ! sigma       self explanatory
  ! SS          s*s expectation value
  ! 
  ! Prints out results to the screen.

  integer :: i, k, ios, iarg
  double precision, allocatable :: u1(:), u2(:), u3(:), dV(:)
  double precision :: exp_r2, exp_v2, j
  double precision :: hc, alpha, b, m, sigma, x1, x2, x3, pi, SS
  character(len=80) :: arg, file

  ! Test if enough command line arguments
  iarg = command_argument_count()
  if (iarg/=7) then
     print *, "Need 7 arguments!"
     print *, "filename k alpha b m sigma SS"
     stop
  end if

  ! Read from command line

  ! Read filename
  call get_command_argument(1,arg)
  read(arg,*) file

  ! Read k
  call get_command_argument(2,arg)
  read(arg,*) k

  ! Read alpha
  call get_command_argument(3,arg)
  read(arg,*) alpha

  ! Read b
  call get_command_argument(4,arg)
  read(arg,*) b

  ! Read m
  call get_command_argument(5,arg)
  read(arg,*) m

  ! Read sigma
  call get_command_argument(6,arg)
  read(arg,*) sigma

  ! Read SS
  call get_command_argument(7,arg)
  read(arg,*) SS

  allocate(u1(k), u2(k), u3(k), dV(k))

  pi = 4d+0*atan(1d+0)
  hc = 197.327d+0
  x1 = 4d+0*alpha*hc/3
  x2 = b/hc
  x3 = 64*alpha*(sigma**5)/(9*sqrt(pi)*(m**2)*(hc**2))*SS

  ! Read wavefunction from file and create the u's
  open(unit=1, file=file, iostat=ios, status='old', action='read')
  if (ios/=0) then
     print *, "Problem opening file."
     stop
  end if
  do i=1,k
     read(1,*,iostat=ios) u1(i)
     j = i*1d+0
     u2(i) = j*j*u1(i)
  end do
  close(1)
  
  do i=1,k

     j = i*1d+0
     dV(i) = x1/j+x2*j-x3*j*j*exp(-sigma*sigma*j*j/(hc*hc))

     u3(i) = dV(i)*u1(i)
     
  end do

  ! Calculate <r>
  call simpson(k, exp_r2, u1, u2)
  call simpson(k, exp_v2, u1, u3)

  print *, 
  print *, "sqrt(<r**2>) is:"
  print *, sqrt(exp_r2)
  print *,
  print *, "<v**2> is;"
  print *, (1d+0/(2*m))*exp_v2

end program r

  
     

  
