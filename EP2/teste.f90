program teste

  use mgh_eval

  implicit none

  integer,parameter::  prob = 14
  integer,parameter::  n = 2

  real(kind=8) :: x(n), g(n), f
  real(kind=8) :: hd(n), hod( (n*n -n)/2 )

  call initpt(n, x, prob, 1.0d0)
  call objfcn(n, x, f, prob)
  call grdfcn(n, x, g, prob)
  call hesfcn(n, x, hd, hod, prob)

  print *, "x = ", x
  print *, "f = ", f
  print *, "g = ", g
  print *, "H = ", hd, hod

end program teste
