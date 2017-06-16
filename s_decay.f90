program s_decay
  use qqModule

  implicit none

  ! Program for calculating all the S-state charmonium decays needed for the
  ! project.
  ! Command line arguments are:
  !
  ! 'file1'    1S-singlet wavefunction
  ! mass1      1S-singlet mass
  ! 'file2'    2S-singlet wavefunction
  ! mass2      2S-singlet mass
  ! 'file3'    1S-triplet wavefunction
  ! mass3      1S-triplet mass
  ! 'file4'    2S-triplet wavefuncton
  ! mass4      2S-triplet mass
  ! 'file5'    3S-triplet wavefunction
  ! mass5      3S-triplet mass
  ! 'file6'    4S-triplet wavefunction
  ! mass6      4S-triplet mass
  !
  ! Prints out results on screen. 

  integer :: i, ios1, ios2, ios3, ios4, ios5, ios66, k, iarg
  double precision :: mass1, mass2, mass3, mass4, mass5, mass6
  double precision :: alpha, b ,sigma, m, j, pi, hc, alpha_cs, alpha_em, q
  double precision :: x1, x2, x3
  double precision :: R_1S, R_2S, R_1T, R_2T, R_3T, R_4T, SS1, SS3
  double precision :: lep_1, lep_2, lep_3, lep_4
  double precision :: rad_1, rad_2, rad_3, ggg_1, ggg_2, pgg_1, pgg_2
  double precision, allocatable :: u1(:), u2(:), u3(:), u4(:), u5(:), u6(:)
  double precision, allocatable :: dV1(:), dV3(:)
  character(len=80) :: file1, file2, file3, file4, file5, file6, arg

  ! Test if there is enough arguments on command line
  iarg = command_argument_count()
  if (iarg/=12) then
     print *, "Need 12 arguments on command line!"
     print *, "file1 mass1 file2 mass2 file3 mass3 file4 mass4 file5 mass5 file6 mass6"
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

  ! Read 'file3'
  call get_command_argument(5,arg)
  read(arg,*) file3

  ! Read mass3
  call get_command_argument(6,arg)
  read(arg,*) mass3

  ! Read 'file4'
  call get_command_argument(7,arg)
  read(arg,*) file4

  ! Read mass4
  call get_command_argument(8,arg)
  read(arg,*) mass4

  ! Read 'file5'
  call get_command_argument(9,arg)
  read(arg,*) file5

  ! Read mass5
  call get_command_argument(10,arg)
  read(arg,*) mass5

  ! Read 'file6'
  call get_command_argument(11,arg)
  read(arg,*) file6

  ! Read mass6
  call get_command_argument(12,arg)
  read(arg,*) mass6

  ! Create functions
  k = 3000
  allocate(u1(k), u2(k), u3(k), u4(k), u5(k), u6(k), dV1(k), dV3(k))
  
  ! Read in the wavefunctions and create everything needed for the integrals
  open(unit=1, file=file1, iostat=ios1, status='old', action='read')
  if (ios1/=0) then
     print *, "Problem opening file1!"
     stop
  end if
  do i=1,k  
     read(1,*,iostat=ios1) u1(i)
     u1(i) = u1(i)*u1(i)
  end do
  close(1)

  open(unit=2, file=file2, iostat=ios2, status='old', action='read')
  if (ios2/=0) then
     print *, "Problem opening file2!"
     stop
  end if
  do i =1,k
     read(2,*,iostat=ios2) u2(i)
     u2(i) = u2(i)*u2(i)
  end do
  close(2)

  open(unit=3, file=file3, iostat=ios3, status='old', action='read')
  if (ios3/=0) then
     print *, "Problem opening file3!"
     stop
  end if
  do i=1,k
     read(3,*,iostat=ios3) u3(i)
     u3(i) = u3(i)*u3(i)
  end do
  close(3)

  open(unit=4, file=file4, iostat=ios4, status='old', action='read')
  if (ios4/=0) then
     print *, "Problem opening file4!"
     stop
  end if
  do i=1,k
     read(4,*,iostat=ios4) u4(i)
     u4(i) = u4(i)*u4(i)
  end do
  close(4)

  open(unit=5, file=file5, iostat=ios5, status='old', action='read')
  if (ios5/=0) then
     print *, "Problem opening file5!"
     stop
  end if
  do i=1,k
     read(5,*,iostat=ios5) u5(i)
     u5(i) = u5(i)*u5(i)
  end do
  close(5)
  
  open(unit=66, file=file6, iostat=ios66, status='old', action='read')
  if (ios66/=0) then
     print *, "Problem opening file6"
     stop
  end if
  do i=1,k
     read(66,*,iostat=ios66) u6(i)
     u6(i) = u6(i)*u6(i)
  end do
  close(66)

  ! Some constants
  alpha = 0.54174d+0
  b = 0.15026d+0
  m = 1.27836d+0
  sigma = 1.0243
  alpha_cs = 0.1951d+0
  q = 2d+0/3
  alpha_em = 1d+0/137.036
  pi = 4d+0*atan(1d+0)
  hc = 197.327d+0
  SS1 = -0.75d+0
  SS3 = 0.25d+0
  x1 = 4d+0*alpha*hc/3
  x2 = b/hc
  x3 = 64d+0*alpha*(sigma**5)/(9d+0*sqrt(pi)*(hc**2)*(m**2))
  
  do i=1,k
     
     j = i*1d+0
     
     dV1(i) = x1/(j*j)+x2-x3*j*exp((-sigma*sigma*j*j)/(hc*hc))*SS1
     
     dV3(i) = x1/(j*j)+x2-x3*j*exp((-sigma*sigma*j*j)/(hc*hc))*SS3
    
  end do
  
  ! Calculate |R(0)|**2
  call simpson(k, R_1S, u1, dV1)
  R_1S = R_1S*m
  call simpson(k, R_2S, u2, dV1)
  R_2S = R_2S*m
  call simpson(k, R_1T, u3, dV3)
  R_1T = R_1T*m
  call simpson(k, R_2T, u4, dV3)
  R_2T = R_2T*m
  call simpson(k, R_3T, u5, dV3)
  R_3T = R_3T*m
  call simpson(k, R_4T, u6, dV3)
  R_4T = R_4T*m
  
  ! Leptonic widths
  lep_1 = ((4d+0*((alpha_em*q)**2)*R_1T)/(mass3**2))*hc
  lep_2 = ((4d+0*((alpha_em*q)**2)*R_2T)/(mass4**2))*hc
  lep_3 = ((4d+0*((alpha_em*q)**2)*R_3T)/(mass5**2))*hc
  lep_4 = ((4d+0*((alpha_em*q)**2)*R_4T)/(mass6**2))*hc

  ! 2 photon and 2 gluon widths
  rad_1 = ((12d+0*(alpha_em**2)*(q**4)*R_1S)/(mass1**2))*hc
  rad_2 = ((12d+0*(alpha_em**2)*(q**4)*R_2S)/(mass2**2))*hc
  
  ! 3 photon and 3 gluon widths
  rad_3 = (16d+0*(pi*pi-9d+0)*(q**6)*(alpha_em**2)*R_1T)/(3d+0*pi*mass3*mass3)&
       *hc
  ggg_1 = ((40d+0*(pi*pi-9d+0)*(alpha_cs**3)*R_1T)/(81d+0*pi*(mass3**2)))*hc
  ggg_2 = ((40d+0*(pi*pi-9d+0)*(alpha_cs**3)*R_2T)/(81d+0*pi*(mass4**2)))*hc

  ! photon-gluon-gluon widths
  pgg_1 = ((32d+0*(pi*pi-9d+0)*alpha_em*(alpha_cs**2)*(q**2)*R_1T)/&
       (9d+0*pi*mass3*mass3))*hc
  pgg_2 = ((32d+0*(pi*pi-9d+0)*alpha_em*(alpha_cs**2)*(q**2)*R_2T)/&
       (9d+0*pi*mass4*mass4))*hc

  ! Print out the results
  print *,
  print *, "Results for S-state Charmonium:"
  print *,
  print *, "Leptonic widths of triplet states:"
  print *, lep_1
  print *, lep_2
  print *, lep_3
  print *, lep_4
  print *,
  print *, "Two photon S-singlet state decay widths:"
  print *, rad_1
  print *, rad_2
  print *,
  print *, "3 photon and the two 3 gluon decay widths:"
  print *, rad_3
  print *, ggg_1
  print *, ggg_2
  print *,
  print *, "Photon-gluon-gluon widths:"
  print *, pgg_1
  print *, pgg_2
  print *,
  
end program s_decay
