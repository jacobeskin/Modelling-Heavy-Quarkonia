module qqModule

  ! This module contains subroutines for calculating parameters and bound state
  ! wavefunctions for different quarkonium states. Subroutines are the
  ! Midpoint method, Simpson integral and Newtons method. Even though Newtons
  ! method is primarily used for finding foots of functions, here it is used to
  ! find a minimum of a function.

contains

  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  !--------------------------- Shooting method --------------------------------
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------

  ! Calculating the wave function and energy level with midpoint method and
  ! bisection method.
  ! Parameters are:
  !
  ! alpha        strong coupling constant
  ! b, c         terms for the linear part of Cornell potential
  ! m            mass of the particle
  ! sigma        sigma term for spin-spin interaction
  ! SS           expectation value of spin-spin operator
  ! l            angular momentum quantum number
  ! E_0          starting value for searching the energy
  ! u            the wave function (as a vector array)
  ! k            number of elements in u
  ! E            calculated energy level
  
  subroutine shoot(alpha, b, c, m, sigma, SS, l, E_0, u, k, E)
    implicit none
    
    integer, intent(in) :: l, k
    double precision, intent(in) :: alpha, b, c, m, sigma, SS, E_0 
    double precision, intent(inout) :: u(k)
    double precision, intent(out) :: E
    double precision, allocatable :: u_l(:)
    
    integer :: i, j, counter
    double precision :: r, d_E, pi, hc, hc2
    double precision :: f, f_0, f_1, N
    double precision :: c0, c1, c2, c3, c4, X

    ! Calculate pi
    pi = 4.0*atan(1.0)
    pi = sqrt(pi)

    r = 1.0d+0 
    
    ! Define some of the parameter values (for speed)
    hc = 197.327d+0
    hc2 = hc**2
    X = m/hc2                                
    c0 = l*(l+1)                                     ! From l-quantum number
    c1 = X*(4.0d+0*alpha*hc)/3                       ! Coulombic term
    c2 = X*(b/hc)                                    ! Linear confining term
    c3 = X*(32*alpha*(sigma**3)/(9*pi*(m**2)))*SS    ! Spin-spin term
    c4 = X*c
    

    ! Build the two functions that are built from shooting form left and right
    allocate(u_l(k))

    ! Set rest of the parameters and boundary conditions before starting
    ! the iteration
    E = E_0                        ! Set the energylevel to starting energy
    d_E = 1.0d-3                    ! Starting energy increment
    u_l(1) = 0.0d+0                   ! Boundary condition in origin
    u_l(2) = r**(1.0d+0*l+1.0d+0)  ! From r->0 asymptotic behaviour
    
    
    ! Bisection method with Midpoint
    counter = 0
    do while (counter<10000)

       ! Construct the wave function 
       do i=3,k
          
          j = i-1

          f = X*E-c0/(j*j*r)+c1/(j*r)-c2*r*j-c4-c3*exp(-sigma*sigma*j*j/(hc2))
          
          ! The wave function is now by Midpoint/Verlet integration
          u_l(i) = (2-f)*u_l(i-1)-u_l(i-2) 
          
       end do
       
       f_0 = u_l(k)
       
       ! Bisection
       if (counter==0) then
          f_1 = f_0
       else if (counter>0) then
          if ((abs(f_0) < abs(f_1)).and.(abs(f_0)>1.0d-4)) then
             E = E+d_E
             f_1 = f_0
          else if ((abs(f_0) < abs(f_1)).and.(abs(f_0)<=1.0d-4)) then
             counter = counter+1
             exit
          else if (abs(f_0) >= abs(f_1)) then
             if (abs(f_0)>=1.0d-4) then
                d_E = -d_E/2
                E = E+d_E
                f_1 = f_0
             else if (abs(f_0)<1.0d-4) then
                counter = counter+1
                exit
             end if
          end if
       end if
       counter = counter + 1
       
    end do

    ! Form the normalized wavefunction, call Simpson()
    call simpson(k, N, u_l, u_l)  ! Calculate the normalization constant
    u = u_l/sqrt(N)               ! Make the normalized wavefunction

    E = E   ! Energylevel

    deallocate(u_l)

  end subroutine shoot


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  !------------------------ Simpson integration --------------------------------
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------

  ! Subroutine for preforming Simpson integration for the purpose of normalizing
  ! the wavefunctons and calculating expectation values.
  ! Parameters are:
  !
  ! k          length of the wavefunction
  ! N          Result of the integration
  ! u_1, u_2   the wavefunctions
  
  subroutine simpson(k, N, u_1, u_2)
    implicit none 
    
    integer, intent(in) :: k
    double precision, intent(out) :: N
    double precision, dimension(k), intent(in) :: u_1, u_2
    double precision :: r, pp1, pp2, pp3
    integer :: i

    r = 1.0d+0 ! Lengthscale
    
    N = 0.0d+0
    
    do i=1,k-2
       
       ! Calculate the integrand in advance
       pp1 = u_1(i)*u_2(i)
       pp2 = u_1(i+1)*u_2(i+1)
       pp3 = u_1(i+2)*u_2(i+2)
       
       !Calculating the integral
       N = N+r*((pp1/3)+(4.0*pp2/3)+(pp3/3))
       
    end do
    
  end subroutine simpson

  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------
  !-----------------------  Line Search --------------------------------------
  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------

  subroutine ln(alpha, b, c, m, sigma, k, E_m, E_c, F0, x1)
    implicit none

    double precision, intent(inout) :: alpha, b, c, m, sigma, x1(5)
    double precision, intent(in) :: E_m(6)
    double precision, intent(out) :: E_c(6), F0
    integer, intent(in) :: k

    ! Rest
    double precision, allocatable :: x0(:)
    double precision :: F1
    double precision :: d_alpha, d_b, d_c, d_m, d_sigma
    integer :: n, i

    allocate(x0(5))
    
    !---------- Initial values ----------

    x0(1) = alpha
    X0(2) = b
    x0(3) = c
    x0(4) = m
    x0(5) = sigma
    
    !---------------------------------------------------------------
    ! Begin the routine for finding the minimum by bisection method,
    ! going through each parameter one at a time.
    !---------------------------------------------------------------
    print *,
    n = 0
    call mss(x0(1), x0(2), x0(3), x0(4), x0(5), k, E_m, E_c, F0)
    do 
       
       ! ---------- alpha ----------

       ! Set initial values
       i = 0
       d_alpha = 1d-2
       x0(1) = x0(1) + d_alpha
       
       ! Start the bisection method
       do
          
          call mss(x0(1), x0(2), x0(3), x0(4), x0(5), k, E_m, E_c, F1)
          
          if (F1<=F0) then
             if (abs(d_alpha)<1d-4) then
                exit
             else
                x0(1) = x0(1) + d_alpha
                F0 = F1
             end if
          else
             d_alpha = (-0.5d+0)*d_alpha
             x0(1) = x0(1) + d_alpha
             F0 = F1
          end if
          i = i+1
          print *,
          print *, "alpha", F0, d_alpha, i
          print *, x0(1)
          print *,
       end do

       F0 = F1
       if (F0<1d-3) then
          exit
       end if
       
       
       ! ---------- b ----------

       ! Set initial values
       i = 0
       d_b = 1d-2
       x0(2) = x0(2) + d_b

       do
          
          call mss(x0(1), x0(2), x0(3), x0(4), x0(5), k, E_m, E_c, F1)
          
          if (F1<=F0) then
             if (abs(d_b)<1d-4) then
                exit
             else
                x0(2) = x0(2) + d_b
                F0 = F1
             end if
          else
             d_b = (-0.5d+0)*d_b
             x0(2) = x0(2) + d_b
             F0 = F1
          end if
          i = i+1
          print *,
          print *, "b", F0, d_b, i
          print *, x0(2)
          print *,
       end do

       F0 = F1
       if (F0<1d-3) then
          exit
       end if

       ! ---------- c ----------

       ! Set initial values
       i = 0
       d_c = 1d-2
       x0(3) = x0(3) + d_c

       do
          
          call mss(x0(1), x0(2), x0(3), x0(4), x0(5), k, E_m, E_c, F1)
          
          if (F1<=F0) then
             if (abs(d_c)<1d-4) then
                exit
             else
                x0(3) = x0(3) + d_c
                F0 = F1
             end if
          else
             d_c = (-0.5d+0)*d_c
             x0(3) = x0(3) + d_c
             F0 = F1
          end if
          i = i+1
          print *,
          print *, "c", F0, d_c, i
          print *, x0(3)
          print *,
       end do

       F0 = F1
       if (F0<1d-3) then
          exit
       end if
       
       ! ---------- m ----------

       ! Set initial values
       i = 0
       d_m = 1d-2
       x0(4) = x0(4) + d_m

       do
          
          call mss(x0(1), x0(2), x0(3), x0(4), x0(5), k, E_m, E_c, F1)
          
          if (F1<=F0) then
             if (abs(d_m)<1d-4) then
                exit
             else
                x0(4) = x0(4) + d_m
                F0 = F1
             end if
          else
             d_m = (-0.5d+0)*d_m
             x0(4) = x0(4) + d_m
             F0 = F1
          end if
          i = i+1
          print *,
          print *, "m", F0, d_m, i
          print *, x0(4)
          print *,
       end do

       F0 = F1
       if (F0<1d-3) then
          exit
       end if

       ! ---------- sigma ----------

       ! Set initial values
       i = 0
       d_sigma = 1d-2
       x0(5) = x0(5) + d_sigma

       do
          
          call mss(x0(1), x0(2), x0(3), x0(4), x0(5), k, E_m, E_c, F1)
          
          if (F1<=F0) then
             if (abs(d_sigma)<=1d-4) then
                exit
             else 
                x0(5) = x0(5) + d_sigma
                F0 = F1
             end if
          else 
             d_sigma = (-0.5d+0)*d_sigma
             x0(5) = x0(5) + d_sigma
             F0 = F1
          end if
          i = i+1
          print *,
          print *, "sigma", F0, d_sigma, i
          print *, x0(5)
          print *,
       end do

       F0 = F1
       if (F0<1d-3) then
          exit
       end if
       
       n = n + 1

       print *,
       print *, "Steps and current F:"
       print *, n, F0
       print *,
       do i=1,5
          print *, x0(i)
       end do
       print *,
             
    end do

    x1 = x0

    print *, " Done!" 
    print *, "----------------------------------------------"
    do i=1,5
       print *, x1(i)
    end do
    print *,
    print *, n, F0
    print *, "----------------------------------------------"

    deallocate(x0)
    
  end subroutine ln
  
  !----------------------------------------------------------------------
  !----------------------------------------------------------------------
  !------------------------ Masses --------------------------------------
  !----------------------------------------------------------------------
  !----------------------------------------------------------------------

  !
  
  subroutine mss(alpha, b, c, m, sigma, k, E_m, E_c, F)
    use omp_lib  ! Library for OpenMP
    implicit none
    
    double precision, intent(inout) :: E_c(6)
    double precision, intent(in) :: E_m(6), alpha, b, c, m, sigma  
    integer, intent(in) :: k
    
    double precision, allocatable :: u(:,:)             ! WF's  
    double precision, allocatable :: u3P(:,:), exp2(:)  ! For 3P exp vals
    
    ! Rest
    double precision, allocatable :: x0(:)
    double precision :: j, jj, pi, hc, hc2, hc3, F
    double precision :: E_1S, E_3S, E_1P, E_3P, V_LS, V_S12
    integer :: i, core
    
    allocate(u(k,4), u3P(k,2), x0(5), exp2(2))
    
    ! Calculate pi and sqrt(pi)
    pi = 4.0d+0*atan(1.0d+0)
    pi = sqrt(pi)
    
    !---------- Initial values ----------
    
    x0(1) = alpha
    X0(2) = b
    x0(3) = c
    x0(4) = m
    x0(5) = sigma
    
    hc = 197.327d+0 ! Reduced Planks constant times speed of light, for scaling
    hc2 = hc**2
    hc3 = hc**3
    
    
    call omp_set_num_threads(4)
    
    ! ===========
    ! Go parallel
    ! ===========
    
    !$omp parallel private(core, i, j , jj) 
    core = omp_get_thread_num() ! Core numbers
    
    ! 1S-singlet
    if (core==0) then
       
       ! Calculate the wavefunction and energylevel
       
       ! Wavefunction
       call shoot(x0(1), x0(2), x0(3), x0(4), x0(5), -0.75d+0, 0, 0d+0, &
            u(:,1), k, E_1S)
       
       ! Meson Mass
       E_c(1) = E_1S + 2*x0(4)
       
    end if

    ! 1S-triplet, similarly as above
    if (core==1) then
       
       call shoot(x0(1), x0(2), x0(3), x0(4), x0(5), 0.25d+0, 0, 0d+0, &
            u(:,2), k, E_3S)
       
       E_c(2) = E_3S + 2*x0(4)
       
    end if
    
    ! 1P-singlet, again as above
    if (core==2) then
       
       call shoot(x0(1), x0(2), x0(3), x0(4), x0(5), -0.75d+0, 1, 0d+0, &
            u(:,3), k, E_1P)
       
       ! Meson Mass
       E_c(3) = E_1P + 2*x0(4)
       
    end if

    ! 1P-triplet, like above except with additional terms 
    if (core==3) then ! 1P-triplet
       
       call shoot(x0(1), x0(2), x0(3), x0(4), x0(5), 0.25d+0, 1, 0d+0, &
            u(:,4), k, E_3P)
       
       ! Calculate again functions needed for dF and mass shifts
       do i=1,k
          
          j = i*1d+0
          jj = j**3
          
          u3P(i,1) = u(i,4)/j  ! For <r**-1>
          u3P(i,2) = u(i,4)/jj ! For <r**-3>
          
       end do
       
       call simpson(k, exp2(1), u(:,4), u3P(:,1)) ! <r**-1> 
       call simpson(k, exp2(2), u(:,4), u3P(:,2)) ! <r**-3>
       
       ! Calculate the mass shifts and energies

      ! Mass for 3P0
       V_LS = -2d+0*(2*x0(1)*exp2(2)*hc3/(x0(4)*x0(4)) - &
            (x0(2)*hc*exp2(1)/(2*x0(4)*x0(4))))
       V_S12 = -4d+0*(x0(1)*exp2(2)*hc3/(3*x0(4)*x0(4)))
       E_c(4) = E_3P + 2*x0(4) + V_LS + V_S12
       
       ! Mass for 3P1
       V_LS = -1d+0*(2*x0(1)*exp2(2)*hc3/(x0(4)*x0(4)) - &
            (x0(2)*hc*exp2(1)/(2*x0(4)*x0(4))))
       V_S12 = 2d+0*(x0(1)*exp2(2)*hc3/(3*x0(4)*x0(4)))
       E_c(5) = E_3P + 2*x0(4) + V_LS + V_S12
       
       ! Mass for 3P2
       V_LS = 1d+0*(2*x0(1)*exp2(2)*hc3/(x0(4)*x0(4)) - &
            (x0(2)*hc*exp2(1)/(2*x0(4)*x0(4))))
       V_S12 = (-2d+0/5)*(x0(1)*exp2(2)*hc3/(3*x0(4)*x0(4)))
       E_c(6) = E_3P + 2*x0(4) + V_LS + V_S12
       
       
    end if
    
    !$omp end parallel
       
    !===============
    ! Parallel ended
    !===============
       
    ! Make F
    F = norm2(E_m-E_c) 

    deallocate(u, u3P, x0, exp2)
    
  end subroutine mss
     
     
end module qqModule
         


  
       
