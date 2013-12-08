module mgh_eval

  implicit none

  ! PARAMETERS
  integer, parameter, private :: dp = kind(0.0d0)

contains

  subroutine initpt(n,x,nprob,factor)

    ! SCALAR ARGUMENTS
    integer,       intent(in) :: n,nprob
    real(kind=dp), intent(in) :: factor

    ! ARRAY ARGUMENT
    real(kind=dp), intent(out) :: x(n)

    ! LOCAL SCALARS
    integer       :: ivar,j
    real(kind=dp) :: dfloat
    real(kind=dp) :: c1,c2,c3,c4,five,h,half,one,ten,three, &
                     twenty,twntf,two,zero

    data zero,half,one,two,three,five,ten,twenty,twntf &
         /0.0d0,0.5d0,1.0d0,2.0d0,3.0d0,5.0d0,1.0d1,2.0d1,2.5d1/
    data c1,c2,c3,c4 /4.0d-1,2.5d0,1.5d-1,1.2d0/

    dfloat(ivar) = ivar

    ! SELECTION OF INITIAL POINT.

    select case ( nprob )

    case (1) 
       ! HELICAL VALLEY FUNCTION.
       x(1) = -one
       x(2) = zero
       x(3) = zero

    case (2)
       ! BIGGS EXP6 FUNCTION.
       x(1) = one
       x(2) = two
       x(3) = one
       x(4) = one
       x(5) = one
       x(6) = one

    case (3)
       ! GAUSSIAN FUNCTION.
       x(1) = c1
       x(2) = one
       x(3) = zero

    case (4)
       ! POWELL BADLY SCALED FUNCTION.
       x(1) = zero
       x(2) = one

    case (5)
       ! BOX 3-DIMENSIONAL FUNCTION.
       x(1) = zero
       x(2) = ten
       x(3) = twenty

    case (6)
       ! VARIABLY DIMENSIONED FUNCTION.
       h = one/dfloat(n)
       do j = 1, n
          x(j) = one - dfloat(j)*h
       end do

    case (7)
       ! WATSON FUNCTION.
       do j = 1, n
          x(j) = zero
       end do

    case (8)
       ! PENALTY FUNCTION I.
       do j = 1, n
          x(j) = dfloat(j)
       end do

    case (9)
       ! PENALTY FUNCTION II.
       do j = 1, n
          x(j) = half
       end do

    case (10)
       ! BROWN BADLY SCALED FUNCTION.
       x(1) = one
       x(2) = one

    case (11)
       ! BROWN AND DENNIS FUNCTION.
       x(1) = twntf
       x(2) = five
       x(3) = -five
       x(4) = -one

    case (12)
       ! GULF RESEARCH AND DEVELOPMENT FUNCTION.
       x(1) = five
       x(2) = c2
       x(3) = c3

    case (13)
       ! TRIGONOMETRIC FUNCTION.
       h = one/dfloat(n)
       do j = 1, n
          x(j) = h
       end do

    case (14)
       ! EXTENDED ROSENBROCK FUNCTION.
       do j = 1, n, 2
          x(j) = -c4
          x(j+1) = one
       end do

    case (15)
       ! EXTENDED POWELL SINGULAR FUNCTION.
       do j = 1, n, 4
          x(j) = three
          x(j+1) = -one
          x(j+2) = zero
          x(j+3) = one
       end do

    case (16)
       ! BEALE FUNCTION.
       x(1) = one
       x(2) = one

    case (17)
       ! WOOD FUNCTION.
       x(1) = -three
       x(2) = -one
       x(3) = -three
       x(4) = -one

    case (18)
       ! CHEBYQUAD FUNCTION.
       h = one/dfloat(n+1)
       do j = 1, n
          x(j) = dfloat(j)*h
       end do

    end select

    ! COMPUTE MULTIPLE OF INITIAL POINT.
    if ( factor .ne. one ) then
       if ( nprob .eq. 7 ) then
          do j = 1, n
             x(j) = factor
          end do
       else
          do j = 1, n
             x(j) = factor*x(j)
          end do
       end if
    end if

  end subroutine initpt

  ! ******************************************************************
  ! ******************************************************************

  subroutine objfcn(n,x,f,nprob)

    ! SCALAR ARGUMENTS
    integer,       intent(in)  :: n,nprob
    real(kind=dp), intent(out) :: f

    ! ARRAY ARGUMENTS
    real(kind=dp), intent(in) :: x(n)
    
    ! LOCAL SCALARS
    integer       :: i,iev,ivar,j
    real(kind=dp) :: dfloat
    real(kind=dp) :: ap,arg,bp,c2pdm6,cp0001,cp1,cp2,cp25,cp5,c1p5, &
                     c2p25,c2p625,c3p5,c25,c29,c90,c100,c10000,     &
                     c1pd6,d1,d2,eight,fifty,five,four,one,r,s1,s2, &
                     s3,t,t1,t2,t3,ten,th,three,tpi,two,zero
    ! LOCAL ARRAYS
    real(kind=dp) :: fvec(50),y(15)

    data zero,one,two,three,four,five,eight,ten,fifty &
         /0.0d0,1.0d0,2.0d0,3.0d0,4.0d0,5.0d0,8.0d0,1.0d1,5.0d1/
    data c2pdm6,cp0001,cp1,cp2,cp25,cp5,c1p5,c2p25,c2p625,c3p5,c25, &
         c29,c90,c100,c10000,c1pd6                                  &
         /2.0d-6,1.0d-4,1.0d-1,2.0d-1,2.5d-1,5.0d-1,1.5d0,2.25d0,   &
          2.625d0,3.5d0,2.5d1,2.9d1,9.0d1,1.0d2,1.0d4,1.0d6/
    data ap,bp /1.0d-5,1.0d0/
    data y(1),y(2),y(3),y(4),y(5),y(6),y(7),y(8),y(9),y(10),y(11),  &
         y(12),y(13),y(14),y(15)                                    &
         /9.0d-4,4.4d-3,1.75d-2,5.4d-2,1.295d-1,2.42d-1,3.521d-1,   &
          3.989d-1,3.521d-1,2.42d-1,1.295d-1,5.4d-2,1.75d-2,4.4d-3, &
          9.0d-4/

    dfloat(ivar) = ivar

    ! FUNCTION ROUTINE SELECTOR.

    select case ( nprob )

    case (1)
       ! HELICAL VALLEY FUNCTION.
       tpi = eight*datan(one)
       th = dsign(cp25,x(2))
       if (x(1) .gt. zero) th = datan(x(2)/x(1))/tpi
       if (x(1) .lt. zero) th = datan(x(2)/x(1))/tpi + cp5
       arg = x(1)**2 + x(2)**2
       r = dsqrt(arg)
       t = x(3) - ten*th
       f = c100*(t**2 + (r - one)**2) + x(3)**2

    case (2)
       ! BIGGS EXP6 FUNCTION.
       f = zero
       do i = 1, 13
          d1 = dfloat(i)/ten
          d2 = dexp(-d1) - five*dexp(-ten*d1) + three*dexp(-four*d1)
          s1 = dexp(-d1*x(1))
          s2 = dexp(-d1*x(2))
          s3 = dexp(-d1*x(5))
          t = x(3)*s1 - x(4)*s2 + x(6)*s3 - d2
          f = f + t**2
       end do

    case (3)
       ! GAUSSIAN FUNCTION.
       f = zero
       do i = 1, 15
          d1 = cp5*dfloat(i-1)
          d2 = c3p5 - d1 - x(3)
          arg = -cp5*x(2)*d2**2
          r = dexp(arg)
          t = x(1)*r - y(i)
          f = f + t**2
       end do

    case (4)
       ! POWELL BADLY SCALED FUNCTION.
       T1 = C10000*X(1)*X(2) - ONE
       S1 = DEXP(-X(1))
       S2 = DEXP(-X(2))
       T2 = S1 + S2 - ONE - CP0001
       F = T1**2 + T2**2

    case (5)
       ! BOX 3-DIMENSIONAL FUNCTION.
       f = zero
       do i = 1, 10
          d1 = dfloat(i)
          d2 = d1/ten
          s1 = dexp(-d2*x(1))
          s2 = dexp(-d2*x(2))
          s3 = dexp(-d2) - dexp(-d1)
          t = s1 - s2 - s3*x(3)
          f = f + t**2
       end do

    case (6)
       ! VARIABLY DIMENSIONED FUNCTION.
       t1 = zero
       t2 = zero
       do j = 1, n
          t1 = t1 + dfloat(j)*(x(j) - one)
          t2 = t2 + (x(j) - one)**2
       end do
       f = t2 + t1**2*(one + t1**2)

    case (7)
       ! WATSON FUNCTION.
       f = zero
       do i = 1, 29
          d1 = dfloat(i)/c29
          s1 = zero
          d2 = one
          do j = 2, n
             s1 = s1 + dfloat(j-1)*d2*x(j)
             d2 = d1*d2
          end do
          s2 = zero
          d2 = one
          do j = 1, n
             s2 = s2 + d2*x(j)
             d2 = d1*d2
          end do
          t = s1 - s2**2 - one
          f = f + t**2
       end do
       t1 = x(2) - x(1)**2 - one
       f = f + x(1)**2 + t1**2
       
    case (8)
       ! PENALTY FUNCTION I.
       t1 = -cp25
       t2 = zero
       do j = 1, n
          t1 = t1 + x(j)**2
          t2 = t2 + (x(j) - one)**2
       end do
       f = ap*t2 + bp*t1**2
       
    case (9)
       ! PENALTY FUNCTION II.
       t1 = -one
       t2 = zero
       t3 = zero
       d1 = dexp(cp1)
       d2 = one
       do j = 1, n
          t1 = t1 + dfloat(n-j+1)*x(j)**2
          s1 = dexp(x(j)/ten)
          if (j .ne. 1) then
             s3 = s1 + s2 - d2*(d1 + one)
             t2 = t2 + s3**2
             t3 = t3 + (s1 - one/d1)**2
          end if
          s2 = s1
          d2 = d1*d2
       end do
       f = ap*(t2 + t3) + bp*(t1**2 + (x(1) - cp2)**2)

    case (10)
       ! BROWN BADLY SCALED FUNCTION.
       t1 = x(1) - c1pd6
       t2 = x(2) - c2pdm6
       t3 = x(1)*x(2) - two
       f = t1**2 + t2**2 + t3**2

    case (11)
       ! BROWN AND DENNIS FUNCTION.
       f = zero
       do i = 1, 20
          d1 = dfloat(i)/five
          d2 = dsin(d1)
          t1 = x(1) + d1*x(2) - dexp(d1)
          t2 = x(3) + d2*x(4) - dcos(d1)
          t = t1**2 + t2**2
          f = f + t**2
       end do
       
    case (12)
       ! GULF RESEARCH AND DEVELOPMENT FUNCTION.
       f = zero
       d1 = two/three
       do i = 1, 99
          arg = dfloat(i)/c100
          r = dabs((-fifty*dlog(arg))**d1+c25-x(2))
          t1 = r**x(3)/x(1)
          t2 = dexp(-t1)
          t = t2 - arg
          f = f + t**2
       end do
       
    case (13)
       ! TRIGONOMETRIC FUNCTION.
       s1 = zero
       do j = 1, n
          s1 = s1 + dcos(x(j))
       end do
       f = zero
       do j = 1, n
          t = dfloat(n+j) - dsin(x(j)) - s1 - dfloat(j)*dcos(x(j))
          f = f + t**2
       end do
       
    case (14)
       ! EXTENDED ROSENBROCK FUNCTION.
       f = zero
       do j = 1, n, 2
          t1 = one - x(j)
          t2 = ten*(x(j+1) - x(j)**2)
          f = f + t1**2 + t2**2
       end do
       
    case (15)
       ! EXTENDED POWELL FUNCTION.
       f = zero
       do j = 1, n, 4
          t = x(j) + ten*x(j+1)
          t1 = x(j+2) - x(j+3)
          s1 = five*t1
          t2 = x(j+1) - two*x(j+2)
          s2 = t2**3
          t3 = x(j) - x(j+3)
          s3 = ten*t3**3
          f = f + t**2 + s1*t1 + s2*t2 + s3*t3
       end do
       
    case (16)
       ! BEALE FUNCTION.
       s1 = one - x(2)
       t1 = c1p5 - x(1)*s1
       s2 = one - x(2)**2
       t2 = c2p25 - x(1)*s2
       s3 = one - x(2)**3
       t3 = c2p625 - x(1)*s3
       f = t1**2 + t2**2 + t3**2

    case (17)
       ! WOOD FUNCTION.
       s1 = x(2) - x(1)**2
       s2 = one - x(1)
       s3 = x(2) - one
       t1 = x(4) - x(3)**2
       t2 = one - x(3)
       t3 = x(4) - one
       f = c100*s1**2 + s2**2 + c90*t1**2 + t2**2 + &
            ten*(s3 + t3)**2 + (s3 - t3)**2/ten
       
    case (18)
       ! CHEBYQUAD FUNCTION.
       do i = 1, n
          fvec(i) = zero
       end do
       do j = 1, n
          t1 = one
          t2 = two*x(j) - one
          t = two*t2
          do i = 1, n
             fvec(i) = fvec(i) + t2
             th = t*t2 - t1
             t1 = t2
             t2 = th
          end do
       end do
       f = zero
       d1 = one/dfloat(n)
       iev = -1
       do i = 1, n
          t = d1*fvec(i)
          if (iev .gt. 0) t = t + one/(dfloat(i)**2 - one)
          f = f + t**2
          iev = -iev
       end do

    end select

  end subroutine objfcn

  ! ******************************************************************
  ! ******************************************************************

  subroutine grdfcn(n,x,g,nprob)

    ! SCALAR ARGUMENTS
    integer, intent(in) :: n,nprob

    ! ARRAY ARGUMENTS
    real(kind=dp), intent(out) :: g(n)
    real(kind=dp), intent(in)  :: x(n)

    ! LOCAL SCALARS
    integer :: i,iev,ivar,j
    real(kind=dp) :: dfloat
    real(kind=dp) :: ap,arg,bp,c2pdm6,cp0001,cp1,cp2,cp25,cp5,c1p5, &
                     c2p25,c2p625,c3p5,c19p8,c20p2,c25,c29,c100,    &
                     c180,c200,c10000,c1pd6,d1,d2,eight,fifty,five, &
                     four,one,r,s1,s2,s3,t,t1,t2,t3,ten,th,three,   &
                     tpi,twenty,two,zero

    ! LOCAL ARRAYS
    real(kind=dp) :: fvec(50),y(15)
    data zero,one,two,three,four,five,eight,ten,twenty,fifty     &
         /0.0d0,1.0d0,2.0d0,3.0d0,4.0d0,5.0d0,8.0d0,1.0d1,2.0d1, &
          5.0d1/
    data c2pdm6,cp0001,cp1,cp2,cp25,cp5,c1p5,c2p25,c2p625,c3p5,   &
         c19p8,c20p2,c25,c29,c100,c180,c200,c10000,c1pd6          &
         /2.0d-6,1.0d-4,1.0d-1,2.0d-1,2.5d-1,5.0d-1,1.5d0,2.25d0, &
          2.625d0,3.5d0,1.98d1,2.02d1,2.5d1,2.9d1,1.0d2,1.8d2,    &
          2.0d2,1.0d4,1.0d6/
    data ap,bp /1.0d-5,1.0d0/
    data y(1),y(2),y(3),y(4),y(5),y(6),y(7),y(8),y(9),y(10),y(11),  &
         y(12),y(13),y(14),y(15)                                    &
         /9.0d-4,4.4d-3,1.75d-2,5.4d-2,1.295d-1,2.42d-1,3.521d-1,   &
          3.989d-1,3.521d-1,2.42d-1,1.295d-1,5.4d-2,1.75d-2,4.4d-3, &
          9.0d-4/

    dfloat(ivar) = ivar

    ! GRADIENT ROUTINE SELECTOR.

    select case ( nprob )
       
    case (1)
       ! HELICAL VALLEY FUNCTION.
       tpi = eight*datan(one)
       th = dsign(cp25,x(2))
       if (x(1) .gt. zero) th = datan(x(2)/x(1))/tpi
       if (x(1) .lt. zero) th = datan(x(2)/x(1))/tpi + cp5
       arg = x(1)**2 + x(2)**2
       r = dsqrt(arg)
       t = x(3) - ten*th
       s1 = ten*t/(tpi*arg)
       g(1) = c200*(x(1) - x(1)/r + x(2)*s1)
       g(2) = c200*(x(2) - x(2)/r - x(1)*s1)
       g(3) = two*(c100*t + x(3))
       
    case (2)
       ! BIGGS EXP6 FUNCTION.
       g(1:n) = zero
       do i = 1, 13
          d1 = dfloat(i)/ten
          d2 = dexp(-d1) - five*dexp(-ten*d1) + three*dexp(-four*d1)
          s1 = dexp(-d1*x(1))
          s2 = dexp(-d1*x(2))
          s3 = dexp(-d1*x(5))
          t = x(3)*s1 - x(4)*s2 + x(6)*s3 - d2
          th = d1*t
          g(1) = g(1) - s1*th
          g(2) = g(2) + s2*th
          g(3) = g(3) + s1*t
          g(4) = g(4) - s2*t
          g(5) = g(5) - s3*th
          g(6) = g(6) + s3*t
       end do
       g(1) = two*x(3)*g(1)
       g(2) = two*x(4)*g(2)
       g(3) = two*g(3)
       g(4) = two*g(4)
       g(5) = two*x(6)*g(5)
       g(6) = two*g(6)
       
    case (3)
       ! GAUSSIAN FUNCTION.
       g(1) = zero
       g(2) = zero
       g(3) = zero
       do i = 1, 15
          d1 = cp5*dfloat(i-1)
          d2 = c3p5 - d1 - x(3)
          arg = -cp5*x(2)*d2**2
          r = dexp(arg)
          t = x(1)*r - y(i)
          s1 = r*t
          s2 = d2*s1
          g(1) = g(1) + s1
          g(2) = g(2) - d2*s2
          g(3) = g(3) + s2
       end do
       g(1) = two*g(1)
       g(2) = x(1)*g(2)
       g(3) = two*x(1)*x(2)*g(3)

    case (4)
       ! POWELL BADLY SCALED FUNCTION.
       t1 = c10000*x(1)*x(2) - one
       s1 = dexp(-x(1))
       s2 = dexp(-x(2))
       t2 = s1 + s2 - one - cp0001
       g(1) = two*(c10000*x(2)*t1 - s1*t2)
       g(2) = two*(c10000*x(1)*t1 - s2*t2)
       
    case (5)
       ! BOX 3-DIMENSIONAL FUNCTION.
       g(1) = zero
       g(2) = zero
       g(3) = zero
       do i = 1, 10
          d1 = dfloat(i)
          d2 = d1/ten
          s1 = dexp(-d2*x(1))
          s2 = dexp(-d2*x(2))
          s3 = dexp(-d2) - dexp(-d1)
          t = s1 - s2 - s3*x(3)
          th = d2*t
          g(1) = g(1) - s1*th
          g(2) = g(2) + s2*th
          g(3) = g(3) - s3*t
       end do
       g(1) = two*g(1)
       g(2) = two*g(2)
       g(3) = two*g(3)
       
    case (6)
       ! VARIABLY DIMENSIONED FUNCTION.
       t1 = zero
       do j = 1, n
          t1 = t1 + dfloat(j)*(x(j) - one)
       end do
       t = t1*(one + two*t1**2)
       do j = 1, n
          g(j) = two*(x(j) - one + dfloat(j)*t)
       end do
       
    case (7)
       ! WATSON FUNCTION.
       g(1:n) = zero
       do i = 1, 29
          d1 = dfloat(i)/c29
          s1 = zero
          d2 = one
          do j = 2, n
             s1 = s1 + dfloat(j-1)*d2*x(j)
             d2 = d1*d2
          end do
          s2 = zero
          d2 = one
          do j = 1, n
             s2 = s2 + d2*x(j)
             d2 = d1*d2
          end do
          t = s1 - s2**2 - one
          s3 = two*d1*s2
          d2 = two/d1
          do j = 1, n
             g(j) = g(j) + d2*(dfloat(j-1) - s3)*t
             d2 = d1*d2
          end do
       end do
       t1 = x(2) - x(1)**2 - one
       g(1) = g(1) + x(1)*(two - four*t1)
       g(2) = g(2) + two*t1
       
    case (8)
       ! PENALTY FUNCTION I.
       t1 = -cp25
       do j = 1, n
          t1 = t1 + x(j)**2
       end do
       d1 = two*ap
       th = four*bp*t1
       do j = 1, n
          g(j) = d1*(x(j) - one) + x(j)*th
       end do
       
    case (9)
       ! PENALTY FUNCTION II.
       t1 = -one
       do j = 1, n
          t1 = t1 + dfloat(n-j+1)*x(j)**2
       end do
       d1 = dexp(cp1)
       d2 = one
       th = four*bp*t1
       do j = 1, n
          g(j) = dfloat(n-j+1)*x(j)*th
          s1 = dexp(x(j)/ten)
          if (j .ne. 1) then
             s3 = s1 + s2 - d2*(d1 + one)
             g(j) = g(j) + ap*s1*(s3 + s1 - one/d1)/five
             g(j-1) = g(j-1) + ap*s2*s3/five
          end if
          s2 = s1
          d2 = d1*d2
       end do
       g(1) = g(1) + two*bp*(x(1) - cp2)
       
    case (10)
       ! BROWN BADLY SCALED FUNCTION.
       t1 = x(1) - c1pd6
       t2 = x(2) - c2pdm6
       t3 = x(1)*x(2) - two
       g(1) = two*(t1 + x(2)*t3)
       g(2) = two*(t2 + x(1)*t3)
       
    case (11)
       ! BROWN AND DENNIS FUNCTION.
       g(1) = zero
       g(2) = zero
       g(3) = zero
       g(4) = zero
       do i = 1, 20
          d1 = dfloat(i)/five
          d2 = dsin(d1)
          t1 = x(1) + d1*x(2) - dexp(d1)
          t2 = x(3) + d2*x(4) - dcos(d1)
          t = t1**2 + t2**2
          s1 = t1*t
          s2 = t2*t
          g(1) = g(1) + s1
          g(2) = g(2) + d1*s1
          g(3) = g(3) + s2
          g(4) = g(4) + d2*s2
       end do
       g(1) = four*g(1)
       g(2) = four*g(2)
       g(3) = four*g(3)
       g(4) = four*g(4)
       
    case (12)
       ! GULF RESEARCH AND DEVELOPMENT FUNCTION.
       g(1) = zero
       g(2) = zero
       g(3) = zero
       d1 = two/three
       do i = 1, 99
          arg = dfloat(i)/c100
          r = dabs((-fifty*dlog(arg))**d1+c25-x(2))
          t1 = r**x(3)/x(1)
          t2 = dexp(-t1)
          t = t2 - arg
          s1 = t1*t2*t
          g(1) = g(1) + s1
          g(2) = g(2) + s1/r
          g(3) = g(3) - s1*dlog(r)
       end do
       g(1) = two*g(1)/x(1)
       g(2) = two*x(3)*g(2)
       g(3) = two*g(3)
       
    case (13)
       ! TRIGONOMETRIC FUNCTION.
       s1 = zero
       do j = 1, n
          g(j) = dcos(x(j))
          s1 = s1 + g(j)
       end do
       s2 = zero
       do j = 1, n
          th = dsin(x(j))
          t = dfloat(n+j) - th - s1 - dfloat(j)*g(j)
          s2 = s2 + t
          g(j) = (dfloat(j)*th - g(j))*t
       end do
       do j = 1, n
          g(j) = two*(g(j) + dsin(x(j))*s2)
       end do

    case (14)
       ! EXTENDED ROSENBROCK FUNCTION.
       do j = 1, n, 2
          t1 = one - x(j)
          g(j+1) = c200*(x(j+1) - x(j)**2)
          g(j) = -two*(x(j)*g(j+1) + t1)
       end do
       
    case (15)
       ! EXTENDED POWELL FUNCTION.
       do j = 1, n, 4
          t = x(j) + ten*x(j+1)
          t1 = x(j+2) - x(j+3)
          s1 = five*t1
          t2 = x(j+1) - two*x(j+2)
          s2 = four*t2**3
          t3 = x(j) - x(j+3)
          s3 = twenty*t3**3
          g(j) = two*(t + s3)
          g(j+1) = twenty*t + s2
          g(j+2) = two*(s1 - s2)
          g(j+3) = -two*(s1 + s3)
       end do

    case (16)
       ! BEALE FUNCTION.
       s1 = one - x(2)
       t1 = c1p5 - x(1)*s1
       s2 = one - x(2)**2
       t2 = c2p25 - x(1)*s2
       s3 = one - x(2)**3
       t3 = c2p625 - x(1)*s3
       g(1) = -two*(s1*t1 + s2*t2 + s3*t3)
       g(2) = two*x(1)*(t1 + x(2)*(two*t2 + three*x(2)*t3))
       
    case (17)
       ! WOOD FUNCTION.
       s1 = x(2) - x(1)**2
       s2 = one - x(1)
       s3 = x(2) - one
       t1 = x(4) - x(3)**2
       t2 = one - x(3)
       t3 = x(4) - one
       g(1) = -two*(c200*x(1)*s1 + s2)
       g(2) = c200*s1 + c20p2*s3 + c19p8*t3
       g(3) = -two*(c180*x(3)*t1 + t2)
       g(4) = c180*t1 + c20p2*t3 + c19p8*s3
       
    case (18)
       ! CHEBYQUAD FUNCTION.
       fvec(1:n) = zero
       do j = 1, n
          t1 = one
          t2 = two*x(j) - one
          t = two*t2
          do i = 1, n
             fvec(i) = fvec(i) + t2
             th = t*t2 - t1
             t1 = t2
             t2 = th
          end do
       end do
       d1 = one/dfloat(n)
       iev = -1
       do i = 1, n
          fvec(i) = d1*fvec(i)
          if (iev .gt. 0) fvec(i) = fvec(i) + one/(dfloat(i)**2 - one)
          iev = -iev
       end do
       do j = 1, n
          g(j) = zero
          t1 = one
          t2 = two*x(j) - one
          t = two*t2
          s1 = zero
          s2 = two
          do i = 1, n
             g(j) = g(j) + fvec(i)*s2
             th = four*t2 + t*s2 - s1
             s1 = s2
             s2 = th
             th = t*t2 - t1
             t1 = t2
             t2 = th
          end do
       end do
       d2 = two*d1
       do j = 1, n
          g(j) = d2*g(j)
       end do

    end select

  end subroutine grdfcn

  ! ******************************************************************
  ! ******************************************************************

  subroutine hesfcn(n,x,hesd,hesl,nprob)

    ! SCALAR ARGUMENTS
    integer, intent(in) :: n,nprob

    ! ARRAY ARGUMENTS
    real(kind=dp), intent(in)  :: x(n)
    real(kind=dp), intent(out) :: hesd(n),hesl(:)

    ! LOCAL SCALARS
    integer       :: i,j,k,m,ii,jj,ix,ivar
    logical       :: iev
    real(kind=dp) :: dfloat
    real(kind=dp) :: arg,d1,d2,d3,logr,p1,p2,piarg,piarg2,    &
                     r,r3inv,s1,s2,s3,s1s2,s1s3,s2s3,ss1,ss2, &
                     t,t1,t2,t3,th,tt,tt1,tt2,tth

    real(kind=dp), parameter :: zero=0.0d0,one=1.0d0,two=2.0d0,        &
         three=3.0d0,four=4.0d0,five=5.0d0,six=6.0d0,eight=8.0d0,      &
         nine=9.0d0,ten=1.0d1,fifty=5.0d1,cp0001=1.0d-4,cp1=1.0d-1,    &
         cp2=2.0d-1,cp25=2.5d-1,cp5=5.0d-1,c1p5=1.5d0,c2p25=2.25d0,    &
         c2p625=2.625d0,c3p5=3.5d0,c12=1.2d1,c19p8=1.98d1,c25=2.5d1,   &
         c29=2.9d1,c50=5.0d1,c90=9.0d1,c100=1.0d2,c120=1.2d2,          &
         c180=1.8d2,c200=2.0d2,c200p2=2.002d2,c202=2.02d2,             &
         c220p2=2.202d2,c360=3.6d2,c400=4.0d2,c1000=1.0d3,             &
         c1080=1.08d3,c1200=1.2d3,c2000=2.0d3,c20000=2.0d4,c2e8=2.0d8, &
         c4e8=4.0d8,ap=1.0d-5,bp=one,pi=3.141592653589793d0

    ! local arrays
    real(kind=dp) :: fvec(50),gvec(50),y(15)

    data y /9.0d-4,4.4d-3,1.75d-2,5.4d-2,1.295d-1,2.42d-1,      &
            3.521d-1,3.989d-1,3.521d-1,2.42d-1,1.295d-1,5.4d-2, &
            1.75d-2,4.4d-3,9.0d-4/

    ix(ii,jj)=(ii-1)*(ii-2)/2+jj
    dfloat(ivar) = ivar

    ! HESSIAN ROUTINE SELECTOR.

    select case ( nprob )
       
    case (1)
       ! HELICAL VALLEY FUNCTION.
       if (x(1) .eq. zero) then
          th = sign(cp25,x(2))
       else
          th = atan(x(2)/x(1)) / (two*pi)
          if (x(1) .lt. zero) th = th + cp5
       end if
       arg = x(1)**2 + x(2)**2
       piarg = pi * arg
       piarg2 = piarg * arg
       r3inv = one / sqrt(arg)**3
       t = x(3) - ten*th
       s1 = five*t / piarg
       p1 = c2000*x(1)*x(2)*t / piarg2
       p2 = (five/piarg)**2
       hesd(1) = c200 - c200*(r3inv-p2)*x(2)**2 - p1
       hesd(2) = c200 - c200*(r3inv-p2)*x(1)**2 + p1
       hesd(3) = c202
       hesl(1) = c200*x(1)*x(2)*r3inv + &
            c1000/piarg2 * ( t*(x(1)**2-x(2)**2) - five*x(1)*x(2)/pi )
       hesl(2) =  c1000*x(2) / piarg
       hesl(3) = -c1000*x(1) / piarg
       
    case (2)
       ! BIGGS EXP6 FUNCTION.
       hesd(1:6)  = zero
       hesl(1:15) = zero
       do i = 1, 13
          d1 = dfloat(i)/ten
          d2 = exp(-d1) - five*exp(-ten*d1) + three*exp(-four*d1)
          s1 = exp(-d1*x(1))
          s2 = exp(-d1*x(2))
          s3 = exp(-d1*x(5))
          t = x(3)*s1 - x(4)*s2 + x(6)*s3 - d2
          d2 = d1**2
          s1s2 = s1 * s2
          s1s3 = s1 * s3
          s2s3 = s2 * s3
          hesd(1) = hesd(1) + d2*s1*(t+x(3)*s1)
          hesd(2) = hesd(2) - d2*s2*(t-x(4)*s2)
          hesd(3) = hesd(3) + s1**2
          hesd(4) = hesd(4) + s2**2
          hesd(5) = hesd(5) + d2*s3*(t+x(6)*s3)
          hesd(6) = hesd(6) + s3**2
          hesl(1) = hesl(1) - d2*s1s2
          hesl(2) = hesl(2) - d1*s1*(t+x(3)*s1)
          hesl(3) = hesl(3) + d1*s1s2
          hesl(4) = hesl(4) + d1*s1s2
          hesl(5) = hesl(5) + d1*s2*(t-x(4)*s2)
          hesl(6) = hesl(6) - s1s2
          hesl(7) = hesl(7) + d2*s1s3
          hesl(8) = hesl(8) - d2*s2s3
          hesl(9) = hesl(9) - d1*s1s3
          hesl(10) = hesl(10) + d1*s2s3
          hesl(11) = hesl(11) - d1*s1s3
          hesl(12) = hesl(12) + d1*s2s3
          hesl(13) = hesl(13) + s1s3
          hesl(14) = hesl(14) - s2s3
          hesl(15) = hesl(15) - d1*s3*(t+x(6)*s3)
       end do
       hesd(1) = x(3)*hesd(1)
       hesd(2) = x(4)*hesd(2)
       hesd(5) = x(6)*hesd(5)
       hesl(1) = x(3)*x(4)*hesl(1)
       hesl(3) = x(4)*hesl(3)
       hesl(4) = x(3)*hesl(4)
       hesl(7) = x(3)*x(6)*hesl(7)
       hesl(8) = x(4)*x(6)*hesl(8)
       hesl(9) = x(6)*hesl(9)
       hesl(10) = x(6)*hesl(10)
       hesl(11) = x(3)*hesl(11)
       hesl(12) = x(4)*hesl(12)
       hesd(1:6)  = two*hesd(1:6)
       hesl(1:15) = two*hesl(1:15)
       
    case (3)
       ! GAUSSIAN FUNCTION.
       hesd(1) = zero
       hesd(2) = zero
       hesd(3) = zero
       hesl(1) = zero
       hesl(2) = zero
       hesl(3) = zero
       do i = 1, 15
          d1 = cp5*dfloat(i-1)
          d2 = c3p5 - d1 - x(3)
          arg = cp5*x(2)*d2**2
          r = exp(-arg)
          t = x(1)*r - y(i)
          t1 = two*x(1)*r - y(i)
          hesd(1) = hesd(1) + r**2
          hesd(2) = hesd(2) + r*t1*d2**4
          hesd(3) = hesd(3) + r*(x(2)*t1*d2**2-t)
          hesl(1) = hesl(1) - r*t1*d2**2
          hesl(2) = hesl(2) + d2*r*t1
          hesl(3) = hesl(3) + d2*r*(t-arg*t1)
       end do
       hesd(1) = two*hesd(1)
       hesd(2) = cp5*x(1)*hesd(2)
       hesd(3) = two*x(1)*x(2)*hesd(3)
       hesl(2) = two*x(2)*hesl(2)
       hesl(3) = two*x(1)*hesl(3)

    case (4)
       ! POWELL BADLY SCALED FUNCTION.
       s1 = exp(-x(1))
       s2 = exp(-x(2))
       t2 = s1 + s2 - one - cp0001
       hesd(1) = c2e8*x(2)**2 + two*s1*(s1+t2)
       hesd(2) = c2e8*x(1)**2 + two*s2*(s2+t2)
       hesl(1) = c4e8*x(1)*x(2) + two*s1*s2 - c20000
       
    case (5)
       ! BOX 3-DIMENSIONAL FUNCTION.
       hesd(1) = zero
       hesd(2) = zero
       hesd(3) = zero
       hesl(1) = zero
       hesl(2) = zero
       hesl(3) = zero
       do i = 1, 10
          d1 = dfloat(i)
          d2 = d1/ten
          s1 = exp(-d2*x(1))
          s2 = exp(-d2*x(2))
          s3 = exp(-d2) - exp(-d1)
          t = s1 - s2 - s3*x(3)
          th = t*d2**2
          hesd(1) = hesd(1) + th*s1 + (d2*s1)**2
          hesd(2) = hesd(2) - th*s2 + (d2*s2)**2
          hesd(3) = hesd(3) + s3**2
          hesl(1) = hesl(1) - s1*s2*d2**2
          hesl(2) = hesl(2) + d2*s1*s3
          hesl(3) = hesl(3) - d2*s2*s3
       end do
       hesd(1) = two*hesd(1)
       hesd(2) = two*hesd(2)
       hesd(3) = two*hesd(3)
       hesl(1) = two*hesl(1)
       hesl(2) = two*hesl(2)
       hesl(3) = two*hesl(3)
       
    case (6)
       ! VARIABLY DIMENSIONED FUNCTION.
       t1 = zero
       do j = 1, n
          t1 = t1 + dfloat(j)*(x(j)-one)
       end do
       t = one + six*t1**2
       m = 0
       do j = 1, n
          hesd(j) = two + two*t*dfloat(j)**2
          do k = 1, j-1
             m = m + 1
             hesl(m) = two*t*dfloat(j*k)
          end do
       end do
       
    case (7)
       ! WATSON FUNCTION.
       hesd(1:n)           = zero
       hesl(1:(n*(n-1))/2) = zero
       do i = 1, 29
          d1 = dfloat(i)/c29
          d2 = one
          s1 = zero
          s2 = x(1)
          do j = 2, n
             s1 = s1 + dfloat(j-1)*d2*x(j)
             d2 = d1*d2
             s2 = s2 + d2*x(j)
          end do
          t = two * (s1-s2**2-one) * d1**2
          s3 = two*d1*s2
          d2 = one/d1
          m = 0
          do j = 1, n
             t1 = dfloat(j-1) - s3
             hesd(j) = hesd(j) + (t1**2-t)*d2**2
             d3 = one/d1
             do k = 1, j-1
                m = m + 1
                hesl(m) = hesl(m) + (t1*(dfloat(k-1)-s3) - t) * d2*d3
                d3 = d1*d3
             end do
             d2 = d1*d2
          end do
       end do
       t3 = x(2) - x(1)**2 - one
       hesd(1) = hesd(1) + one - two*(t3-two*x(1)**2)
       hesd(2) = hesd(2) + one
       hesl(1) = hesl(1) - two*x(1)
       hesd(1:n) = two * hesd(1:n)
       hesl(1:(n*(n-1))/2) = two * hesl(1:(n*(n-1))/2)
       
    case (8)
       ! PENALTY FUNCTION I.
       t1 = -cp25
       do j = 1, n
          t1 = t1 + x(j)**2
       end do
       d1 = two*ap
       th = four*bp*t1
       m = 0
       do j = 1, n
          hesd(j) = d1 + th + eight*x(j)**2
          do k = 1, j-1
             m = m + 1
             hesl(m) = eight*x(j)*x(k)
          end do
       end do
       
    case (9)
       ! PENALTY FUNCTION II.
       t1 = -one
       do j = 1, n
          t1 = t1 + dfloat(n-j+1)*x(j)**2
       end do
       d1 = exp(cp1)
       d2 = one
       th = four*bp*t1
       m = 0
       do j = 1, n
          hesd(j) = eight*bp*(dfloat(n-j+1)*x(j))**2 + dfloat(n-j+1)*th
          s1 = exp(x(j)/ten)
          if (j .gt. 1) then
             s3 = s1 + s2 - d2*(d1 + one)
             hesd(j) = hesd(j) + ap*s1*(s3 + s1 - one/d1 + two*s1)/c50
             hesd(j-1) = hesd(j-1) + ap*s2*(s2+s3)/c50
             do k = 1, j-1
                m = m + 1
                t1 = exp(dfloat(k)/ten)
                hesl(m) = eight*dfloat(n-j+1)*dfloat(n-k+1)*x(j)*x(k)
             end do
             hesl(m) = hesl(m) + ap*s1*s2/c50
          end if
          s2 = s1
          d2 = d1*d2
       end do
       hesd(1) = hesd(1) + two*bp
       
    case (10)
       ! BROWN BADLY SCALED FUNCTION.
       hesd(1) = two + two*x(2)**2
       hesd(2) = two + two*x(1)**2
       hesl(1) = four*x(1)*x(2) - four
       
    case (11)
       ! BROWN AND DENNIS FUNCTION.
       hesd(1:4) = zero
       hesl(1:6) = zero
       do i = 1, 20
          d1 = dfloat(i)/five
          d2 = sin(d1)
          t1 = x(1) + d1*x(2) - exp(d1)
          t2 = x(3) + d2*x(4) - cos(d1)
          t = eight * t1 * t2
          s1 = c12*t1**2 + four*t2**2
          s2 = c12*t2**2 + four*t1**2
          hesd(1) = hesd(1) + s1
          hesd(2) = hesd(2) + s1*d1**2
          hesd(3) = hesd(3) + s2
          hesd(4) = hesd(4) + s2*d2**2
          hesl(1) = hesl(1) + s1*d1
          hesl(2) = hesl(2) + t
          hesl(4) = hesl(4) + t*d2
          hesl(3) = hesl(3) + t*d1
          hesl(5) = hesl(5) + t*d1*d2
          hesl(6) = hesl(6) + s2*d2
       end do
       
    case (12)
       ! GULF RESEARCH AND DEVELOPMENT FUNCTION.
       do i = 1, 3
          hesd(i) = zero
          hesl(i) = zero
       end do
       d1 = two/three
       do i = 1, 99
          arg = dfloat(i)/c100
          r = (-fifty*log(arg))**d1+c25-x(2)
          t1 = abs(r)**x(3)/x(1)
          t2 = exp(-t1)
          t3 = t1 * t2 * (t1*t2+(t1-one)*(t2-arg))
          t = t1 * t2 * (t2-arg)
          logr = log(abs(r))
          hesd(1) = hesd(1) + t3 - t
          hesd(2) = hesd(2) + (t+x(3)*t3)/r**2
          hesd(3) = hesd(3) + t3*logr**2
          hesl(1) = hesl(1) + t3/r
          hesl(2) = hesl(2) - t3*logr
          hesl(3) = hesl(3) + (t-x(3)*t3*logr)/r
       end do
       hesd(1) = hesd(1) / x(1)**2
       hesd(2) = hesd(2) * x(3)
       hesl(1) = hesl(1) * x(3)/x(1)
       hesl(2) = hesl(2) / x(1)
       do i = 1, 3
          hesd(i) = two * hesd(i)
          hesl(i) = two * hesl(i)
       end do
       
    case (13)
       ! TRIGONOMETRIC FUNCTION.
       s1 = zero
       do j = 1, n
          hesd(j) = sin(x(j))
          s1 = s1 + cos(x(j))
       end do
       s2 = zero
       m = 0
       do j = 1, n
          th = cos(x(j))
          t = dfloat(n+j) - hesd(j) - s1 - dfloat(j)*th
          s2 = s2 + t
          do k = 1, j-1
             m = m + 1
             hesl(m) = sin(x(k))*(dfloat(n+j+k)*hesd(j)-th) - &
                  hesd(j)*cos(x(k))
             hesl(m) = two*hesl(m)
          end do
          hesd(j) = dfloat(j*(j+2)+n)*hesd(j)**2 + &
               th*(th-dfloat(2*j+2)*hesd(j)) + t*(dfloat(j)*th+hesd(j))
       end do
       do j = 1, n
          hesd(j) = two*(hesd(j) + cos(x(j))*s2)
       end do
       
    case (14)
       ! EXTENDED ROSENBROCK FUNCTION.
       hesl(1:(n*(n-1))/2) = zero
       do j = 1, n, 2
          hesd(j+1) = c200
          hesd(j) = c1200*x(j)**2 - c400*x(j+1) + two
          hesl(ix(j+1,j)) = -c400*x(j)
       end do
       
    case (15)
       ! EXTENDED POWELL FUNCTION.
       hesl(1:(n*(n-1))/2) = zero
       do j = 1, n, 4
          t2 = x(j+1) - two*x(j+2)
          t3 = x(j) - x(j+3)
          s1 = c12 * t2**2
          s2 = c120 * t3**2
          hesd(j) = two + s2
          hesd(j+1) = c200 + s1
          hesd(j+2) = ten + four*s1
          hesd(j+3) = ten + s2
          hesl(ix(j+1,j)) = two*ten
          hesl(ix(j+2,j)) = zero
          hesl(ix(j+2,j+1)) = -two*s1
          hesl(ix(j+3,j)) = -s2
          hesl(ix(j+3,j+1)) = zero
          hesl(ix(j+3,j+2)) = -ten
       end do
       
    case (16)
       ! BEALE FUNCTION.
       s1 = one - x(2)
       t1 = c1p5 - x(1)*s1
       s2 = one - x(2)**2
       t2 = c2p25 - x(1)*s2
       s3 = one - x(2)**3
       t3 = c2p625 - x(1)*s3
       hesd(1) = two * (s1**2 + s2**2 + s3**2)
       hesd(2) = two*x(1) * (x(1) + two*t2 + four*x(1)*x(2)**2 + &
            six*x(2)*t3 + nine*x(1)*x(2)**4)
       hesl(1) = two*(t1-x(1)*s1) + four*x(2)*(t2-x(1)*s2) + &
            six*(t3-x(1)*s3)*x(2)**2
       
    case (17)
       ! WOOD FUNCTION.
       hesd(1) = c1200*x(1)**2 - c400*x(2) + two
       hesd(2) = c220p2
       hesd(3) = c1080*x(3)**2 - c360*x(4) + two
       hesd(4) = c200p2
       hesl(1) = -c400*x(1)
       hesl(2) = zero
       hesl(3) = zero
       hesl(4) = zero
       hesl(5) = c19p8
       hesl(6) = -c360*x(3)
       
    case (18)
       ! CHEBYQUAD FUNCTION.
       fvec(1:n) = zero
       do j = 1, n
          t1 = one
          t2 = two*x(j) - one
          t = two*t2
          do i = 1, n
             fvec(i) = fvec(i) + t2
             th = t*t2 - t1
             t1 = t2
             t2 = th
          end do
       end do
       d1 = one/float(n)
       iev = .false.
       do i = 1, n
          fvec(i) = d1*fvec(i)
          if (iev) fvec(i) = fvec(i) + one/(dfloat(i)**2 - one)
          iev = .not. iev
       end do
       d2 = two*d1
       m = 0
       do j = 1, n
          hesd(j) = four*d1
          t1 = one
          t2 = two*x(j) - one
          t = two*t2
          s1 = zero
          s2 = two
          p1 = zero
          p2 = zero
          gvec(1) = s2
          do i = 2, n
             th = four*t2 + t*s2 - s1
             s1 = s2
             s2 = th
             th = t*t2 - t1
             t1 = t2
             t2 = th
             th = eight*s1 + t*p2 - p1
             p1 = p2
             p2 = th
             gvec(i) = s2
             hesd(j) = hesd(j) + fvec(i)*th + d1*s2**2
          end do
          hesd(j) = d2*hesd(j)
          do k = 1, j-1
             m = m + 1
             hesl(m) = zero
             tt1 = one
             tt2 = two*x(k) - one
             tt = two*tt2
             ss1 = zero
             ss2 = two
             do i = 1, n
                hesl(m) = hesl(m) + ss2*gvec(i)
                tth = four*tt2 + tt*ss2 - ss1
                ss1 = ss2
                ss2 = tth
                tth = tt*tt2 - tt1
                tt1 = tt2
                tt2 = tth
             end do
             hesl(m) = d2*d1*hesl(m)
          end do
       end do

    end select

  end subroutine hesfcn

end module mgh_eval
