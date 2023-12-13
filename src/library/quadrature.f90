module quadrature
  use math
  use quicksort
  use precision
  implicit none
  
  ! Dummy module
  
end module quadrature


! ======================================================= !
! Recurrence coefficients for monic Jacobi polynomials    !
! Example:                                                !
!  - a= 0.0 and b= 0.0 => Legendre polynomial             !
!  - a=-0.5 and b=-0.5 => Chebyshev polynomial (1st kind) !
!  - a= 0.5 and b= 0.5 => Chebyshev polynomial (2nd kind) !
! The resulting polynomials are orthogonal on [-1,+1]     !
! with respect to w(x)=(1-x)^a*(1+x)^b. The n alpha-coeff !
! are stored in the first column, the n beta-coeff in the !
! second column of the (n,2) array ab.                    !
! ======================================================= !
subroutine rjacobi(n,a,b,ab)
  use quadrature
  implicit none
  
  integer, intent(in) :: n
  real(WP), intent(in) :: a
  real(WP), intent(in) :: b
  real(WP), dimension(n,2), intent(out) :: ab
  real(WP) :: nu,mu
  integer :: i
  
  ! Check the parameters are ok
  if (n.le.0 .or. a.le.-1.0_WP .or. b.le.-1.0_WP) STOP 'In rjacobi: parameters out of range.'
  
  ! Compute the polynomials
  nu = (b-a)/(a+b+2.0_WP)
  mu = 2.0_WP**(a+b+1.0_WP)*gamma_function(a+1.0_WP)*gamma_function(b+1.0_WP)/gamma_function(a+b+2.0_WP)
  if (n.eq.1) then
     ab(1,1) = nu
     ab(1,2) = mu
     return
  end if
  ab(1,1) = nu
  do i=2,n
     ab(i,1)=(b**2-a**2)/((2.0_WP*real(i-1,WP)+a+b)*(2.0_WP*real(i-1,WP)+a+b+2.0_WP))
  end do
  ab(1,2) = mu
  ab(2,2) = 4.0_WP*(a+1.0_WP)*(b+1.0_WP)/((a+b+2.0_WP)**2*(a+b+3.0_WP))
  do i=3,n
     ab(i,2)=4.0_WP*(real(i-1,WP)+a)*(real(i-1,WP)+b)*real(i-1,WP)*(real(i-1,WP)+a+b) / &
          ((2.0_WP*real(i-1,WP)+a+b)**2*(2.0_WP*real(i-1,WP)+a+b+1.0_WP)*(2.0_WP*real(i-1,WP)+a+b-1.0_WP))
  end do
  
  return
end subroutine rjacobi


! ============================================================ !
! Gauss quadrature rule                                        !
! Needs the weight function w encoded by the (n,2) array ab of !
! the first n recurrence coefficients of the given polynomials !
! This is obtained by calling rjacobi with the correct params. !
! ============================================================ !
subroutine gauss(n,ab,x,w)
  use quadrature
  implicit none
  
  integer, intent(in) :: n
  real(WP), dimension(n,2) :: ab
  real(WP), dimension(n)   :: x
  real(WP), dimension(n)   :: w
  real(WP), dimension(n,n) :: eigv
  integer,  dimension(n)   :: order
  integer :: i
  
  ! Form the matrix
  x(:)=ab(:,1)
  w(2:n)=sqrt(ab(2:n,2))
  
  ! Find the eigenvalues
  eigv=0.0_WP
  do i=1,n
     eigv(i,i) = 1.0_WP
  end do
  call tqli(x,n,w,eigv)
  
  ! Sort x and eigv
  do i=1,n
     order(i) = i
  end do
  call quick_sort(x,order)
  w(:)=ab(1,2)*eigv(1,order(:))**2
  
  return
end subroutine gauss


! ============================================================== !
! Gauss-Lobatto quadrature rule                                  !
! Needs the weight function w encoded by the (n+2,2) array ab of !
! the first n+2 recurrence coefficients of the given polynomials !
! This is obtained by calling rjacobi with the correct params.   !
! Following this quadrature rule, two points at endl and endr    !
! will be part of the collocation points (typically -1 and 1).   !
! ============================================================== !
subroutine lobatto(n,ab,endl,endr,x,w)
  use quadrature
  implicit none
  
  integer, intent(in) :: n
  real(WP), dimension(n+2,2) :: ab
  real(WP), intent(in) :: endl,endr
  real(WP), dimension(n+2), intent(out) :: x
  real(WP), dimension(n+2), intent(out) :: w
  real(WP) :: p0l,p0r,p1l,p1r,pm1l,pm1r,det
  integer :: i
  
  ! Fix two points
  p0l=0.0_WP
  p0r=0.0_WP
  p1l=1.0_WP
  p1r=1.0_WP
  do i=1,n+1
     pm1l=p0l
     p0l=p1l
     pm1r=p0r
     p0r=p1r
     p1l=(endl-ab(i,1))*p0l-ab(i,2)*pm1l
     p1r=(endr-ab(i,1))*p0r-ab(i,2)*pm1r
  end do
  det=p1l*p0r-p1r*p0l
  ab(n+2,1)=(endl*p1l*p0r-endr*p1r*p0l)/det
  ab(n+2,2)=(endr-endl)*p1l*p1r/det
  
  ! Compute Gauss quadrature
  call gauss(n+2,ab,x,w)
  
  return
end subroutine lobatto
