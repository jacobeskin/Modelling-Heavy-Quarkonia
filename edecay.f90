program edecay
  use qqModule
  
  implicit none

  ! Program for calculating the E1 transition of quarkonia. Commad line
  ! arguments needed are
  !
  ! 'file1'        filename of the initial state wavefunction
  ! mass1          mass of the initial state
  ! 'file2'        filename of the final state wavefunction
  ! mass2          mass of the final state
  ! k              number of elements in the shorter wavefunction file
  ! c_fi           angluar matrix element
  !
  ! Prints out the 
  
  integer :: i, ios1, ios2, k, iarg
  double precision :: mass1, mass2, c_fi, X, Y, L
  double precision, allocatable :: u1(:), u2(:)
  character(len=80) :: file1, file2 ,arg

  ! Test if there is enough arguments on command line
  iarg = command_argument_count()
  if (iarg/=6) then
     print *, "Need 6 arguments on command line!"
     print *, "file1 mass1 file2 mass2 k c_fi"
     stop
  end if

  ! Read variables from command line

  ! Read 'file1'
  call get_command_argument(1,arg)
  read(arg,*) file1

  ! Read mass1
  call get_command_argument(2,arg)
  read(arg,*) mass1

  ! Read 'file2'
  call get_command_argument(3,arg)
  read(arg,*) file2

  ! Read mass2
  call get_command_argument(4,arg)
  read(arg,*) mass2

  ! Read k
  call get_command_argument(5,arg)
  read(arg,*) k

  ! Read c_fi
  call get_command_argument(6,arg)
  read(arg,*) c_fi

  ! Allocate wavefunctions
  allocate(u1(k), u2(k))

  ! Open files for reading and read in the wavefunctions

  open(unit=1, file='file1', iostat=ios1, status='old')
  if (ios1/=0) then
     print *, "Problem opening file1"
     stop
  end if

  open(unit=2, file='file2', iostat=ios2, status='old')
  if (ios2/=0) then
     print *, "Problem opening file2"
     stop
  end if

  do i=1,k
     read(1,*,iostat=ios1) u1(i)
     read(2,*,iostat=ios2) u2(i)
     u2(i) = i*u2(i)
  end do

  ! Calculate the transition L
  X = ((16d+0*c_fi)/(27d+0*137.036d+0))*((mass1-mass2)**3) ! Number part 
  call simpson(k, Y, u1, u2)                           ! <u1|r|u2>
  Y = Y*Y
  L = X*Y

  ! Print results
  print *,
  print *, "E1 decay width is:"
  print *, L
  print *,

end program edecay

  
