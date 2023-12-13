module math

  ! External modules
  use precision

  implicit none

  ! Trigonometric parameters
  real(WP), parameter :: Pi    = 3.1415926535897932385_WP
  real(WP), parameter :: twoPi = 6.2831853071795864770_WP

  ! Bessel first zero
  real(WP), parameter :: besselJ1Zero = 3.8317059702075123115_WP

contains

  ! Compute the cross product of a 3x1 vector
  ! -----------------------------------------
  function cross_product(x, y) result(z)

    implicit none

    ! Arguments
    real(WP), dimension(3), intent(in) :: x, y

    ! Result
    real(WP), dimension(3) :: z

    z(1) = x(2) * y(3) - x(3) * y(2)
    z(2) = x(3) * y(1) - x(1) * y(3)
    z(3) = x(1) * y(2) - x(2) * y(1)

  end function cross_product


  ! Normalize a vector of size 3x1
  ! ------------------------------
  function normalize(v) result(w)

    implicit none

    ! Arguments
    real(WP), dimension(3), intent(in) :: v

    ! Result
    real(WP), dimension(3) :: w

    ! Local variables
    real(WP) :: norm

    norm = 1.0_WP / (sqrt(dot_product(v, v)) + epsilon(1.0_WP))
    w = v * norm

  end function normalize

  ! Normalize a vector of size N x nDimensions
  ! -----------------------------------------
   subroutine normalize_vector(vector)

    implicit none

    ! Arguments
    real(WP), intent(inout) :: vector(:,:)

    ! Local variables
    integer :: i
    real(WP) :: norm, invnorm, eps

    eps = 10.0_WP * epsilon(1.0_WP)

    select Case (size(vector,2))
       
    case (1)
       do i = 1, size(vector, 1)
          norm = abs(vector(i,1)) 
          if (norm .gt. eps) Then
             invnorm = 1.0_WP / norm
          else
             invnorm = 0.0_WP
          end if
          vector(i,1) = vector(i,1) * invnorm
       end do
       
    case (2)
       do i = 1, size(vector, 1)
          norm = sqrt(vector(i,1)**2 + vector(i,2)**2)
          if (norm .gt. eps) then
             invnorm = 1.0_WP / norm
          else
             invnorm = 0.0_WP
          end if
          vector(i,1) = vector(i,1) * invnorm
          vector(i,2) = vector(i,2) * invnorm
       end do
       
    case (3)
       do i = 1, size(vector, 1)
          norm = sqrt(vector(i,1)**2 + vector(i,2)**2 + vector(i,3)**2)
          if (norm .gt. eps) Then
             invnorm = 1.0_WP / norm
          else
             invnorm = 0.0_WP
          end if
          vector(i,1) = vector(i,1) * invnorm
          vector(i,2) = vector(i,2) * invnorm
          vector(i,3) = vector(i,3) * invnorm
       end do
    end select

    return
  end subroutine normalize_vector


  ! Returns the gamma function
  ! --------------------------
  function gamma_function(xx)
    
    implicit none

    ! Arguments
    real(WP) :: gamma_function
    real(WP), intent(in) :: xx
    
    gamma_function = exp(gammaln(xx))
    
    return
  end function gamma_function
  
  
  ! Returns the log of the gamma function
  ! ------------------------------------
  function gammaln(xx)
    
    implicit none

    ! Arguments
    real(WP) :: gammaln
    real(WP), intent(in) :: xx

    ! Local variables
    real(WP), parameter :: stp = 2.5066282746310005_WP
    real(WP), dimension(6), parameter :: cof = (/ 76.18009172947146_WP, &
         -86.50532032941677_WP, 24.01409824083091_WP,-1.231739572450155_WP, &
         .1208650973866179E-2_WP, -.5395239384953E-5_WP /)
    
    real(WP) :: ser,tmp,x,y
    integer :: j

    x = xx
    y = x
    tmp = x + 5.5_WP
    tmp = (x + 0.5_WP) * log(tmp) - tmp
    ser = 1.000000000190015_WP
    do j = 1, 6
       y = y + 1.0_WP
       ser = ser+cof(j)/y
    end do
    gammaln = tmp + log(stp*ser / x)
    
    return
  end function gammaln


  ! Compute the distance between the point 'x' and an infinite line
  ! defined by the two points 'x1' and 'x2'
  ! -------------------------------------------------------------
  function distance_from_point_to_line(x1, x2, x) result(d)

    implicit none

    ! Arguments
    real(WP), dimension(3), intent(in) :: x1
    real(WP), dimension(3), intent(in) :: x2
    real(WP), dimension(3), intent(in) :: x

    ! Result
    real(WP) :: d

    ! Local variables
    real(WP), dimension(3) :: x21
    real(WP) :: x21_mag

    x21 = x2 - x1
    x21_mag = norm2(x21)

    if (x21_mag .gt. 0.0_WP) then

        d = norm2( cross_product( x21, x1 - x ) ) / x21_mag

    else

        d = norm2(x1 - x)

    end if

  end function distance_from_point_to_line


  ! Compute the distance between the point 'x' and the line segment
  ! defined by the two points 'x1' and 'x2'
  ! -------------------------------------------------------------
  function distance_from_point_to_line_segment(x1, x2, x) result(d)

    implicit none

    ! Arguments
    real(WP), dimension(:), intent(in) :: x1
    real(WP), dimension(:), intent(in) :: x2
    real(WP), dimension(:), intent(in) :: x

    ! Result
    real(WP) :: d

    ! Local variables
    real(WP), dimension(size(x1)) :: x12
    real(WP) :: s, x12_mag

    x12 = x1 - x2
    x12_mag = norm2(x12)

    if (x12_mag .lt. epsilon(1.0_WP)) then

       d = norm2(x1 - x)

    else

       s = dot_product(x1-x, x12) / dot_product(x12, x12)

       ! If the projection is not on the segment then use the closest end point:
       if (s .lt. 0.0_WP) then
          s = 0.0_WP
       else if (s .gt. 1.0_WP) then
          s = 1.0_WP
       end if

       d = norm2(x - x1 + s * x12)

    end if

  end function distance_from_point_to_line_segment


  ! Order 0 Bessel of the first kind
  ! --------------------------------
  function besselJ0(x)

    implicit none

    ! Arguments
    real(WP) :: besselJ0
    real(WP), intent(in) :: x

    ! Local variables
    real(WP), dimension(6), parameter :: r = (/57568490574.0_WP, -13362590354.0_WP,          &
         651619640.7_WP, -11214424.18_WP,77392.33017_WP,-184.9052456_WP/)
    real(WP), dimension(6), parameter :: s = (/57568490411.0_WP, 1029532985.0_WP,            &
         9494680.718_WP, 59272.64853_WP, 267.8532712_WP, 1.0_WP/)
    real(WP), dimension(5), parameter :: p = (/1.0_WP, -.1098628627E-2_WP,                   &
         0.2734510407E-4_WP, -0.2073370639E-5_WP, 0.2093887211E-6_WP/)
    real(WP), dimension(5), parameter :: q = (/-0.1562499995E-1_WP, 0.1430488765E-3_WP,      &
         -0.6911147651E-5_WP, 0.7621095161E-6_WP, -0.934945152E-7_WP/)
    real(WP) :: ax, xx, z, y

    if(abs(x) < 8.0_WP) then
       y = x ** 2
       besselJ0 = (r(1) + y * (r(2) + y*(r(3) + y * (r(4) + y * (r(5) + y * r(6)))))) /      &
            (s(1) + y * (s(2) + y * (s(3) + y * (s(4) + y * (s(5) + y * s(6))))))
    else
       ax = abs(x)
       z = 8.0_WP / ax
       y = z ** 2
       xx = ax - 0.785398164_WP
       besselJ0 = sqrt(0.636619772_WP / ax) * (cos(xx) * (p(1) + y * (p(2) + y* (p(3) + y *  &
            (p(4) + y * p(5))))) - z * sin(xx) * (q(1) + y * (q(2) + y * (q(3) + y * (q(4) + &
            y * q(5))))))
    end if

  end function besselJ0


  ! Order 1 Bessel of the first kind
  ! --------------------------------
  function besselJ1(x)

    implicit none

    ! Arguments
    real(WP) :: besselJ1
    real(WP), intent(in) :: x

    ! Local variables
    real(WP), dimension(6), parameter :: r = (/72362614232.0_WP, -7895059235.0_WP,           &
         242396853.1_WP, -2972611.439_WP, 15704.48260_WP, -30.16036606_WP/)
    real(WP), dimension(6), parameter :: s = (/144725228442.0_WP, 2300535178.0_WP,           &
         18583304.74_WP, 99447.43394_WP, 376.9991397_WP, 1.0_WP/)
    real(WP), dimension(5), parameter :: p = (/1.0_WP, 0.183105E-2_WP, 0-.3516396496E-4_WP,  &
         0.2457520174E-5_WP, -0.240337019E-6_WP/)
    real(WP), dimension(5), parameter :: q = (/0.04687499995_WP, -0.2002690873E-3_WP,        &
         0.8449199096E-5_WP, -.88228987E-6_WP, 0.105787412E-6_WP/)
    real(WP) :: ax, xx, z, y

    if (abs(x) < 8.0_WP) then
       y = x ** 2
       besselJ1 = x * (r(1) + y * (r(2) + y * (r(3) + y * (r(4) + y * (r(5) + y * r(6)))))) /&
            (s(1) + y * (s(2) + y * (s(3) + y * (s(4) + y * (s(5) + y * s(6))))))
    else
       ax = abs(x)
       z = 8.0_WP / ax
       y = z ** 2
       xx = ax - 2.356194491_WP
       besselJ1 = sqrt(0.636619772_WP / ax) * (cos(xx) * (p(1) + y * (p(2) + y * (p(3) + y * &
            (p(4) + y * p(5))))) - z * sin(xx) * (q(1) + y * (q(2) + y * (q(3) + y * (q(4) + &
            y * q(5)))))) * sign(1.0_WP, x)
    end if

  end function besselJ1


  ! Compute the arctan
  ! ------------------
  function arctan(dx, dy)

    implicit none

    ! Arguments
    real(WP), intent(in) :: dx, dy
    real(WP) :: arctan

    ! Evaluate atan
    if (abs(dx) + abs(dy) < 1.0e-9_WP) then
       arcTan = 0.0_WP
    else
       arctan = atan(dy / dx)
    end if

    ! Quadrant correction
    if (dx <= 0.0_WP) then
       arctan = Pi + arctan
    elseif (dy <= 0.0_WP .and. dx > 0.0_WP) then
       arctan = twoPi + arctan
    end if

  end function arctan


  ! Compute the hyperbolic secant
  ! -----------------------------
  function sech(x)

    implicit none

    ! Arguments
    real(WP), intent(in) :: x
    real(WP) :: sech

    sech = 1.0_WP / cosh(x)

  end function sech


  ! Partition items via the pigeonhole principle
  ! --------------------------------------------
  subroutine pigeon_hole(nPigeons, nHoles, holeIndex, offset, nPigeonsInHole)

    implicit none

    ! Arguments
    integer, intent(in) :: nPigeons, nHoles, holeIndex
    integer, intent(out) :: offset, nPigeonsInHole

    offset = holeIndex * (nPigeons / nHoles) + min(holeIndex, mod(nPigeons, nHoles))
    nPigeonsInHole = nPigeons / nHoles
    if (holeIndex < mod(nPigeons, nHoles)) nPigeonsInHole = nPigeonsInHole + 1

  end subroutine pigeon_hole


  ! Find a root via the bisection method
  ! ------------------------------------
  subroutine bisection(xLoc, iLoc, xArray, iStart, iEnd)

    implicit none

    ! Arguments
    real(WP), intent(in) :: xLoc
    integer, intent(out) :: iLoc
    integer, intent(in) :: iStart
    integer, intent(in) :: iEnd
    real(WP), dimension(iStart:iEnd), intent(in) :: xArray

    ! Local variables
    integer :: il, im, iu

    ! Initialize lower and upper limits
    il = iStart - 1
    iu = iEnd

    ! While not done
    do while (iu - il > 1)
       ! Compute a mid-point
       im = (iu + il) / 2
       ! Replace lower of upper limit as appropriate
       if (xLoc >= xArray(im)) then
          il = im
       else
          iu = im
       end if
    end do

    ! Finalize the output
    if (abs(xLoc - xArray(iStart)) < epsilon(1.0_WP)) then
       iLoc = iStart
    else if (xLoc >= xArray(iEnd)) then
       iLoc = iEnd
    else
       iLoc = il
    end if

  end subroutine bisection


  ! B-spline basis function of degree 2
  ! -----------------------------------
  function bspline2(x) result(B02)

    implicit none

    ! Arguments
    real(WP), intent(in) :: x

    ! Result
    real(WP) :: B02

    if (x .gt. 0.0_WP .and. x .le. 1.0_WP / 3.0_WP) then
       B02 = 4.5_WP * x**2
    else if (x .gt. 1.0_WP / 3.0_WP .and. x .le. 2.0_WP / 3.0_WP) then
       B02 = -9.0_WP * x**2 + 9.0_WP * x - 1.5_WP
    else if (x .gt. 2.0_WP / 3.0_WP .and. x .le. 1.0_WP) then
       B02 = 4.5_WP * (1 - x)**2
    else
       B02 = 0.0_WP
    end if
    
  end function bspline2
  

  ! Regularized a Heaviside function
  ! phi: signed-distance (levelset)
  ! eps: regularization parameters (length)
  ! ---------------------------------------
  function regularize_heaviside(phi, eps) result(H)

    implicit none

    ! Arguments
    real(WP), intent(in) :: phi, eps
    real(WP)             :: H

    if (phi .le. -eps) then
       H = 0.0_WP
    else if (abs(phi) .lt. eps) then
       H = 0.5_WP * (1.0_WP + phi / eps + 1.0_WP / pi * sin(phi * pi / eps))
    else
       H = 1.0_WP
    end if

    return
  end function regularize_heaviside

  
  ! Inverse a matrix
  ! ----------------
  subroutine inverse_matrix(A,B)
    implicit none

    real(WP), dimension(:,:), intent(inout) :: A
    real(WP), dimension(:,:), intent(out) :: B
    integer :: i,l,m,n

    ! Get matrix size
    m=size(A,1)
    n=size(A,2)
    if (m.ne.n) call die('inverse_matrix: matrix must be square!')
    i=size(B,1)
    l=size(B,2)
    if (m*n.ne.i*l) call die('inverse_matrix: A /= B')

    ! Zero out B
    B=0.0_WP

    ! Forward elimination
    do i=1,n
       B(i,i) = 1.0_WP
       B(i,:) = B(i,:) / A(i,i)
       A(i,:) = A(i,:) / A(i,i)
       do l=i+1,n
          B(l,:) = B(l,:) - A(l,i)*B(i,:)
          A(l,:) = A(l,:) - A(l,i)*A(i,:)
       end do
    end do

    ! Backward substitution
    do i=n,1,-1
       do l=i+1,n
          B(i,:) = B(i,:) - A(i,l)*B(l,:)
       end do
    end do

    return
  end subroutine inverse_matrix

  
  ! Determinant of a matrix
  ! -----------------------
  function matdet(A,n) result(d)
    implicit none
    
    integer, intent(in) :: n
    real(WP), dimension(:,:) :: A
    real(WP) :: d
    integer :: i

    ! Perform LU decomposition
    call lu_decomp(A,n,d)

    ! Product of diagonal terms
    do i=1,n
       d=d*A(i,i)
    end do

    return
  end function matdet

  ! LU decomposition of A(n,n)
  ! d tracks parity
  ! index provides permutations
  ! ---------------------------
  subroutine lu_decomp(A,n,d)
    implicit none
    
    integer, intent(in) :: n
    real(WP), dimension(:,:) :: A
    real(WP), intent(out) :: d
    integer :: i,j,k,imax
    real(WP) :: Amax,dum,sum
    real(WP), dimension(n) :: vv
    real(WP), parameter :: tiny=1.0e-20_WP
    
    ! No row interchange yet
    d=1.0_WP
    
    ! Get scaling information in first pass
    do i=1,n
       ! Find maximum along row
       Amax=0.0_WP
       do j=1,n
          if (abs(A(i,j)).gt.Amax) Amax=abs(A(i,j))
       end do
       ! Stop if zero row
       if (abs(Amax).lt.epsilon(1.0_WP)) STOP 'Singular matrix in lu_decomp'
       ! Save the scaling
       vv(i)=1.0_WP/Amax
    end do
    
    ! Crout's method
    do j=1,n
       ! Create sum
       do i=1,j-1
          sum=A(i,j)
          do k=1,i-1
             sum=sum-A(i,k)*A(k,j) 
          end do
          A(i,j)=sum
       end do
       ! Find pivot element
       Amax=0.0_WP
       do i=j,n
          sum=A(i,j)
          do k=1,j-1
             sum=sum-A(i,k)*A(k,j) 
          end do
          A(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.Amax) then
             imax=i
             Amax=dum
          end if
       end do
       ! Interchange rows if necessary
       if (j.ne.imax) then
          do k=1,n
             dum=A(imax,k)
             A(imax,k)=A(j,k)
             A(j,k)=dum
          end do
          ! Track parity if necessary
          d=-d
          ! Update scaling
          vv(imax)=vv(j)
       end if
              
       ! Deal with singular pivot
       if (abs(A(j,j)).lt.tiny) A(j,j)=tiny
       
       ! Divide by the pivot element
       if (j.ne.n) then
          dum=1.0_WP/A(j,j)
          do i=j+1,n
             A(i,j)=A(i,j)*dum
          end do
       end if
    end do
    
    return
  end subroutine lu_decomp

  ! Evaluates a 1D Lagrange basis
  ! n: The number of data points
  ! x: The interpolation nodes
  ! xp: The evluation point
  ! L: Lagrange basis
  ! -------------------------------------------
  subroutine lagrange_basis(n, x, xp, L)
    
    implicit none

    ! Arguments
    integer,  intent(in) :: n
    real(WP), intent(in), dimension(n) :: x
    real(WP), intent(in) :: xp
    real(WP), intent(out), dimension(n) :: L
 
    ! Local variables
    integer :: i

    select case (n)

    case (1)
       call die('lagrange_basis: number of data points must be > 1!')

    case (2)
       ! Second-order
       L(1) = (xp - x(2)) / (x(1) - x(2))
       L(2) = (xp - x(1)) / (x(2) - x(1))

    case (3)
       ! Third-order
       L(1) =                                                                                &
            (xp - x(2)) / (x(1) - x(2)) *                                                    &
            (xp - x(3)) / (x(1) - x(3))
       L(2) =                                                                                &
            (xp - x(1)) / (x(2) - x(1)) *                                                    &
            (xp - x(3)) / (x(2) - x(3))
       L(3) =                                                                                &
            (xp - x(1)) / (x(3) - x(1)) *                                                    &
            (xp - x(2)) / (x(3) - x(2))

    case (4)
       ! Fourth-order
       L(1) =                                                                                &
            (xp - x(2)) / (x(1) - x(2)) *                                                    &
            (xp - x(3)) / (x(1) - x(3)) *                                                    &
            (xp - x(4)) / (x(1) - x(4))
       L(2) =                                                                                &
            (xp - x(1)) / (x(2) - x(1)) *                                                    &
            (xp - x(3)) / (x(2) - x(3)) *                                                    &
            (xp - x(4)) / (x(2) - x(4))
       L(3) =                                                                                &
            (xp - x(1)) / (x(3) - x(1)) *                                                    &
            (xp - x(2)) / (x(3) - x(2)) *                                                    &
            (xp - x(4)) / (x(3) - x(4))
       L(4) =                                                                                &
            (xp - x(1)) / (x(4) - x(1)) *                                                    &
            (xp - x(2)) / (x(4) - x(2)) *                                                    &
            (xp - x(3)) / (x(4) - x(3))

    case (5)
       ! Fifth-order
       L(1) =                                                                                &
            (xp - x(2)) / (x(1) - x(2)) *                                                    &
            (xp - x(3)) / (x(1) - x(3)) *                                                    &
            (xp - x(4)) / (x(1) - x(4)) *                                                    &
            (xp - x(5)) / (x(1) - x(5))
       L(2) =                                                                                &
            (xp - x(1)) / (x(2) - x(1)) *                                                    &
            (xp - x(3)) / (x(2) - x(3)) *                                                    &
            (xp - x(4)) / (x(2) - x(4)) *                                                    &
            (xp - x(5)) / (x(2) - x(5))
       L(3) =                                                                                &
            (xp - x(1)) / (x(3) - x(1)) *                                                    &
            (xp - x(2)) / (x(3) - x(2)) *                                                    &
            (xp - x(4)) / (x(3) - x(4)) *                                                    &
            (xp - x(5)) / (x(3) - x(5))
       L(4) =                                                                                &
            (xp - x(1)) / (x(4) - x(1)) *                                                    &
            (xp - x(2)) / (x(4) - x(2)) *                                                    &
            (xp - x(3)) / (x(4) - x(3)) *                                                    &
            (xp - x(5)) / (x(4) - x(5))

       L(5) =                                                                                &
            (xp - x(1)) / (x(5) - x(1)) *                                                    &
            (xp - x(2)) / (x(5) - x(2)) *                                                    &
            (xp - x(3)) / (x(5) - x(3)) *                                                    &
            (xp - x(4)) / (x(5) - x(4))

    case (6)
       ! Sixth-order
       L(1) =                                                                                &
            (xp - x(2)) / (x(1) - x(2)) *                                                    &
            (xp - x(3)) / (x(1) - x(3)) *                                                    &
            (xp - x(4)) / (x(1) - x(4)) *                                                    &
            (xp - x(5)) / (x(1) - x(5)) *                                                    &
            (xp - x(6)) / (x(1) - x(6))
       L(2) =                                                                                &
            (xp - x(1)) / (x(2) - x(1)) *                                                    &
            (xp - x(3)) / (x(2) - x(3)) *                                                    &
            (xp - x(4)) / (x(2) - x(4)) *                                                    &
            (xp - x(5)) / (x(2) - x(5)) *                                                    &
            (xp - x(6)) / (x(2) - x(6))
       L(3) =                                                                                &
            (xp - x(1)) / (x(3) - x(1)) *                                                    &
            (xp - x(2)) / (x(3) - x(2)) *                                                    &
            (xp - x(4)) / (x(3) - x(4)) *                                                    &
            (xp - x(5)) / (x(3) - x(5)) *                                                    &
            (xp - x(6)) / (x(3) - x(6))
       L(4) =                                                                                &
            (xp - x(1)) / (x(4) - x(1)) *                                                    &
            (xp - x(2)) / (x(4) - x(2)) *                                                    &
            (xp - x(3)) / (x(4) - x(3)) *                                                    &
            (xp - x(5)) / (x(4) - x(5)) *                                                    &
            (xp - x(6)) / (x(4) - x(6))
       L(5) =                                                                                &
            (xp - x(1)) / (x(5) - x(1)) *                                                    &
            (xp - x(2)) / (x(5) - x(2)) *                                                    &
            (xp - x(3)) / (x(5) - x(3)) *                                                    &
            (xp - x(4)) / (x(5) - x(4)) *                                                    &
            (xp - x(6)) / (x(5) - x(6))
       L(6) =                                                                                &
            (xp - x(1)) / (x(6) - x(1)) *                                                    &
            (xp - x(2)) / (x(6) - x(2)) *                                                    &
            (xp - x(3)) / (x(6) - x(3)) *                                                    &
            (xp - x(4)) / (x(6) - x(4)) *                                                    &
            (xp - x(5)) / (x(6) - x(5))

    case default
       ! High-order
       do i = 1, n
          L(i) = product( (xp - x(1:i-1)) / (x(i) - x(1:i-1)) )                              &
               * product( (xp - x(i+1:n)) / (x(i) - x(i+1:n)) )
       end do

    end select

    return
  end subroutine lagrange_basis

  ! Interpolation based on Lagrange polynomials
  ! n: The number of data points
  ! x: The data points
  ! y: The data values
  ! xp: The interpolation point
  ! yp: The interpolation value
  ! -------------------------------------------
  subroutine lagrange_interp(n, x, y, xp, yp)
    
    implicit none

    ! Arguments
    integer,  intent(in) :: n
    real(WP), intent(in), dimension(n) :: x, y
    real(WP), intent(in) :: xp
    real(WP), intent(out) :: yp
 
    ! Local variables
    real(WP), dimension(n) :: L

    ! Get the Lagrange polynomial basis
    call lagrange_basis(n, x, xp, L)
    
    ! Interpolate
    yp = sum(L * y)       
    
    return
  end subroutine lagrange_interp

end module math
