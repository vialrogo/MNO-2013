program teste

  use mgh_eval

  implicit none

  real(kind=8) :: x(2), g(2), f
  real(kind=8) :: hd(2), hod(1)

  ! Aqui eu testo com o problema 16 - Beale, n = 2

  call initpt(2, x, 16, 1.0d0)

  print *, "x = ", x

  call objfcn(2, x, f, 16)
  call grdfcn(2, x, g, 16)
  call hesfcn(2, x, hd, hod, 16)

  print *, "f = ", f
  print *, "g = ", g
  print *, "H = ", hd, hod

end program teste
