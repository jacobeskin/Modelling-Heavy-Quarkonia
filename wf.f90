program wf
  use qqModule

  implicit none

  ! Program for calcualting quarkonium masses. Parameters are read from command
  ! line. Prints out the calculated wafefunction to file and calculated mass
  ! to screen.
  !
  ! Example how to run:
  ! $ ./wf.exe 0.48682d+0 0.17518d+0 0.8743d+0 0.86352d+0 0.78861d+0 0d+0 1
  ! -0.75d+0 0d+0 0d+0 'P-Singlet-1' 3000

  !
  ! Parameters:
  !
  ! alpha             strong coupling constant
  ! b,c               parameters for linear term
  ! m                 quark mass
  ! sigma             sigma term of spin-spin interaction
  ! E_0               starting energy
  ! l                 the orbital angular momentum quantum number
  ! SS                <S_1*S_2>
  ! LS                <L*S>
  ! S_12              <S_12>
  ! filename          name of the file that has the wavefunction
  ! k                 distance of integration, in attometers

  ! From command line
  integer :: l, k
  double precision :: alpha, b, c, m, sigma, E_0, SS, LS, S_12
  character(len=80) :: filename

  integer :: i, iarg
  double precision :: E, MesonMass, str, nd, j, hc, hc3
  double precision :: E_LS, E_S12, r_1, r_3 ! For mass shifts
  double precision, allocatable :: u(:), u_1(:), u_3(:)
  character(len=80) :: arg

  call cpu_time(str)

  ! Test if there is enough arguments on command line
  iarg = command_argument_count()
  if (iarg/=12) then
     print *, "Need 12 arguments from command line"
     print *, "alpha b c m sigma E_0 l SS LS S_12 'filename' k" 
     stop
  end if

  ! --------------------------------------------------------------------
  ! Get variables from command line, and define some other constants.
  ! --------------------------------------------------------------------

  ! Read alpha
  call get_command_argument(1,arg)
  read(arg,*) alpha

  !Read b
  call get_command_argument(2,arg)
  read(arg,*) b

  ! Read c
  call get_command_argument(3,arg)
  read(arg,*) c

  ! Read m
  call get_command_argument(4,arg)
  read(arg,*) m

  ! Read sigma
  call get_command_argument(5,arg)
  read(arg,*) sigma

  ! Read E_0
  call get_command_argument(6,arg)
  read(arg,*) E_0

  ! Read l
  call get_command_argument(7,arg)
  read(arg,*) l

  ! Read <S*S>
  call get_command_argument(8,arg)
  read(arg,*) SS

  ! Read <L*S>
  call get_command_argument(9,arg)
  read(arg,*) LS

  ! Read <S_12>
  call get_command_argument(10,arg)
  read(arg,*) S_12

  ! Read filename
  call get_command_argument(11,arg)
  read(arg,*) filename

  ! Read k
  call get_command_argument(12,arg)
  read(arg,*) k
  

  ! Some constants
  hc = 197.327d+0 ! Reduced Planks constant times speed of light, units GeV am
  hc3 = hc**3

  ! Allocate memory for the wave function
  allocate(u(k))

  ! Call subroutine MPM from qqModule in order to get the wavefunction and E
  call shoot(alpha, b, c, m, sigma, SS, l, E_0, u, k, E)
  open(unit=99, file=filename, status='unknown')
  do i=1,k
     write(99,*) u(i)
  end do
  close(99)
  
  ! ----------------------------------
  ! Calculate the possible mass shifts
  ! ----------------------------------
  
  if ((l>0).and.(LS/=0d+0).and.(S_12/=0d+0)) then ! Higher triplet states
     
     ! Like above, form functions for calculating expectation values
     allocate(u_1(k), u_3(k))
     do i=1,k
        j = i*1d+0
        
        u_1(i) = u(i)/j
        
        u_3(i) = u(i)/(j**3)
        
     end do

     ! Calculate mass shifts due to spin-orbit and tensor interactions
     call SIMPSON(k, r_1, u, u_1)  ! <r**-1>
     call SIMPSON(k, r_3, u, u_3)  ! <r**-3>
     E_LS = LS*((2d+0*hc3*alpha*r_3)/(m*m) - (b*hc*r_1)/(2d+0*m*m))
     E_S12 = S_12*((alpha*hc3*r_3)/(3d+0*m*m))
     
     ! Meson mass
     MesonMass = E + 2*m + E_LS + E_S12

     

  else

     MesonMass = E + 2*m
     
  end if

  
  call cpu_time(nd)

  print *, 
  print *,
  print *, "Done!!"
  print *,
  print *, "Calculated energy level:", E
  print *,
  print *, "Calculated meson mass:", MesonMass
  print *,
  print *, "Code runtime:", nd-str
  print *, 

end program wf

  

     
