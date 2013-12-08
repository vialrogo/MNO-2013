!*****************************************************************************
!*****************************************************************************
program main

    implicit none

    ! Parameters
    integer, parameter :: nmax=1000

    ! Variables
    integer :: n,m,iter,i
    real*8  :: ftf
    real*8  :: x(nmax), f(nmax), g(nmax), fj(nmax,nmax)

    ! Call the function
    call getfun (x,n,f,m,ftf,fj,nmax,g,-1)

    ! Set initial point
    write(*,"('Insert the initial point (' I2.2 ')')"), n 

    do i=1,n
        read (*,*) x(i)
    end do

    ! Call the function
    call linearSearch(n,m,x,iter)
    write(*,"('The solutin is:' F10.6  F10.6 ' in ' I6, ' iterations')"), x(1), x(2), iter 

end program main
    
!*****************************************************************************
!*****************************************************************************
subroutine linearSearch(n,m,x,iter)
 
    ! Parameters
    real*8,  parameter :: gamma = 1.0d-04
    real*8,  parameter :: epsg  = 1.0d-010
    
    ! Arguments
    integer, intent(in) :: n,m
    integer, intent(out) :: iter
    real*8, intent(inout) :: x(n)

    ! Variables
    real*8  :: alpha,gnorm,gtd, ftf, ftfnew
    real*8  :: f(n), g(n), fj(m,n), d(n), xnew(n)

    ! Call the function
    call getfun(x,n,f,m,ftf,fj,m,g,-1)
    call getfun(x,n,f,m,ftf,fj,m,g,1111)

    gnorm = maxval(abs(g))
    iter = 0

    do while ( gnorm >= epsg )

        ! Calculate a new direction
        d = -g

        ! Calculate a new point
        alpha = 1.0d0
        xnew = x + alpha * d
        gtd = dot_product(g,d)
        call getfun (xnew,n,f,m,ftfnew,fj,m,g,1111)

        ! Confirm if it's a usefull point
        do while ( ftfnew > (ftf + alpha * gamma * gtd ) )
            
            alpha = 0.5d0 * alpha
            xnew = x + alpha * d
            call getfun (xnew,n,f,m,ftfnew,fj,m,g,1111)
        
        end do

        ! Update the point
        ftf = ftfnew
        x = xnew
        gnorm = maxval(abs(g))
        iter = iter + 1

    end do
   
end subroutine linearSearch

