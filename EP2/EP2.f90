!*****************************************************************************
!*****************************************************************************
program main

    ! Module used
    use mgh_eval

    implicit none

    ! Problems parameters
    integer,parameter::  n = 2
    integer,parameter::  p = 14
    integer :: prob

    ! Variables
    real*8  :: x(n)
    
    do prob=p, p, 1
        
        call initpt(n, x, prob, 1.0d0) ! Set initial point
        call linearSearch(n,x,prob) ! Call the subrutine
    
    end do

end program main
    
!*****************************************************************************
!*****************************************************************************
subroutine linearSearch(n,x,prob)
 
    ! Arguments
    integer, intent(in) :: n
    integer, intent(in) :: prob
    real*8, intent(inout) :: x(n)

    ! Parameters
    integer, parameter :: nmax  = 1000
    real*8,  parameter :: gamma = 1.0d-04
    real*8,  parameter :: epsg  = 1.0d-010

    ! Variables
    integer :: iter,i
    real*8  :: alpha,f,fnew,gnorm,gtd
    real*8  :: g(n),d(n),xnew(n)

    ! Eval function
    call evalf(n,x,f,prob)
    call evalg(n,x,g,prob)

    ! Calculate gradient norm
    gnorm = maxval(abs(g))
    iter = 0

    do while ( gnorm >= epsg )

        ! Calculate a new direction
        call calculateDMaxDesc(n,x,d,prob) 

        ! Calculate a new point
        alpha = 1.0d0
        xnew = x + alpha * d
        call evalf(n,xnew,fnew,prob)

        ! Confirm if it's a usefull point
        gtd = dot_product(g,d)
        do while ( fnew > (f + alpha * gamma * gtd ) )
            
            alpha = 0.5d0 * alpha
            xnew = x + alpha * d
            call evalf(n,xnew,fnew,prob)
        
        end do

        ! Update the point
        f = fnew
        x = xnew
        call evalg(n,x,g,prob)
        gnorm = maxval(abs(g))
        
        ! Output
        iter = iter + 1

    end do

    ! Print final conditions
    write(*,*) 'Iterations:', iter, 'f:', f
    write(*,*) 'The solution is: ', x 
   
end subroutine linearSearch

!*****************************************************************************
!*****************************************************************************
subroutine calculateDMaxDesc(n,x,d,prob)

    use mgh_eval

    ! Arguments
    integer, intent(in) :: n
    integer, intent(in) :: prob
    real*8, intent(in) :: x(n)
    real*8, intent(out) :: d(n)

    call grdfcn(n, x, d, prob)
    d = -d

end subroutine calculateDMaxDesc

!*****************************************************************************
!*****************************************************************************
subroutine calculateDNewton(n,x,d,prob)

    use mgh_eval

    ! Arguments
    integer, intent(in) :: n
    integer, intent(in) :: prob
    real*8, intent(in) :: x(n)
    real*8, intent(out) :: d(n)

    call grdfcn(n, x, d, prob)
    d = -d

end subroutine calculateDNewton

!*****************************************************************************
!*****************************************************************************
subroutine evalf(n,x,f,prob)

    use mgh_eval

    ! Arguments
    integer, intent(in) :: n
    integer, intent(in) :: prob
    real*8, intent(in) :: x(n)
    real*8, intent(out) :: f

    call objfcn(n, x, f, prob)

end subroutine evalf

!*****************************************************************************
!*****************************************************************************
subroutine evalg(n,x,g,prob)

    use mgh_eval
    
    ! Arguments
    integer, intent(in) :: n
    integer, intent(in) :: prob
    real*8,  intent(in) :: x(n)
    real*8,  intent(out) :: g(n)

    call grdfcn(n, x, g, prob)

end subroutine evalg

!*****************************************************************************
!*****************************************************************************
subroutine evalh(n,x,hd,hod,prob)

    use mgh_eval
    
    ! Arguments
    integer, intent(in) :: n
    integer, intent(in) :: prob
    real*8,  intent(in) :: x(n)
    real*8,  intent(out) :: hd(n) 
    real*8,  intent(out) :: hod( (n*n -n)/2 )
    
    call hesfcn(n, x, hd, hod, prob)

end subroutine evalh
