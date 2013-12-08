!*****************************************************************************
!*****************************************************************************
program main

    implicit none

    ! Variables
    integer, parameter :: n=2
    real*8  :: xsol(n), xini(n)

    ! Set initial point
    xini(1) = 1.2d0
    xini(2) = 1.2d0

    ! Call the function
    call linearSearch(n,xini,xsol)

    ! Print solution
    write(*,*) 'The solutin is: ', xsol(1), xsol(2) 

end program main
    
!*****************************************************************************
!*****************************************************************************
subroutine linearSearch(n,xini,xnew)
 
    ! Input arguments
    integer, intent(in) :: n
    real*8, intent(in) :: xini(n)

    ! Output arguments
    real*8,  intent(out) :: xnew(n)
    
    ! Parameters
    integer, parameter :: nmax  = 1000
    real*8,  parameter :: gamma = 1.0d-04
    real*8,  parameter :: epsg  = 1.0d-010

    ! Variables
    integer :: iter
    real*8  :: alpha,f,fnew,gnorm,gtd
    real*8  :: g(n),d(n),x(n)

    ! Eval function
    x = xini
    call evalf(n,x,f)
    call evalg(n,x,g)

    ! Calculate gradient norm
    gnorm = maxval(abs(g))
    iter = 0

    do while ( gnorm >= epsg )

        ! Calculate a new direction
        d = -g 

        ! Calculate a new point
        alpha = 1.0d0
        xnew = x + alpha * d
        call evalf(n,xnew,fnew)

        ! Confirm if it's a usefull point
        gtd = dot_product(g,d)
        do while ( fnew > (f + alpha * gamma * gtd ) )
            
            alpha = 0.5d0 * alpha
            xnew = x + alpha * d
            call evalf(n,xnew,fnew)
        
        end do

        ! Update the point
        f = fnew
        x = xnew
        call evalg(n,x,g)
        gnorm = maxval(abs(g))
        
        ! Output
        iter = iter + 1
        write(*,"('iter: ' I8 '  f: ',E14.6,'  gnorm: ',E14.6)") iter, f ,gnorm

    end do
   
end subroutine linearSearch

!*****************************************************************************
!*****************************************************************************
subroutine evalf(n,x,f)

    ! Input arguments
    integer, intent(in) :: n
    real*8, intent(in) :: x(n)
    
    ! Output arguments
    real*8, intent(out) :: f

    f = 100.0d0 * ( x(2) - x(1) ** 2 ) ** 2 + ( 1.0d0 - x(1) ) ** 2

end subroutine evalf

!*****************************************************************************
!*****************************************************************************
subroutine evalg(n,x,g)

    ! Input arguments
    integer, intent(in) :: n
    real*8,  intent(in) :: x(n)

    ! Output arguments
    real*8,  intent(out) :: g(n)

    g(1) = - 400.0d0 * ( x(2) - x(1) ** 2 ) * x(1) - 2.0d0 * ( 1.0d0 - x(1) )
    g(2) =   200.0d0 * ( x(2) - x(1) ** 2 )

end subroutine evalg
