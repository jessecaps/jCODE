subroutine polydiagonal(nd,A,R,n,lot,dir,b,per)

  ! External modules
  use precision
  use parallel

  implicit none

  ! Arguments
  integer, intent(in) :: dir
  integer, intent(in) :: n,lot,nd
  real(WP), intent(inout), dimension(lot,n,-nd:nd) :: A
  real(WP), intent(inout), dimension(lot,n)        :: R
  real(WP), intent(inout), dimension(lot,n,nd)     :: b
  logical, intent(in), dimension(3) :: per

  ! Local variables
  integer :: proc,nper
  
  ! Safety check
  if (dir.lt.1 .or. dir .gt.3) stop 'Unknown direction'
  
  ! Get parallel info
  proc = nProcsDir(dir)
  if (per(dir)) then
     nper = 1
  else
     nper = 0
  end if
  
  ! If serial
  if (proc.eq.1) then
     if (nper.eq.0) then
        call polydiagonal_serial(nd,A,R,n,lot)
     else
        call polydiagonal_periodic_serial(nd,A,R,n,lot,b)
     end if
     return
  else
     stop 'polydiagonal: Implemented only in serial.'
  end if
  
  return
end subroutine polydiagonal


subroutine polydiagonal_serial(nd,A,R,n,lot)

  ! External modules
  use precision

  implicit none

  ! Arguments
  integer,  intent(in) :: n,lot,nd
  real(WP), intent(inout), dimension(lot,n,-nd:nd) :: A
  real(WP), intent(inout), dimension(lot,n)        :: R

  ! Local variables
  real(WP), dimension(lot) :: const,pivot
  integer :: i,j,k
  
  ! Forward elimination
  do i=1,n-1
     pivot(:) = 1.0_WP/A(:,i,0)
     do j=1,min(nd,n-i)
        const(:) = A(:,i+j,-j)*pivot(:)
        do k=1,min(nd,n-i)
           A(:,i+j,-j+k) = A(:,i+j,-j+k) - A(:,i,k)*const(:)
        end do
        R(:,i+j) = R(:,i+j) - R(:,i)*const(:)
     end do
  end do
  
  ! Back-substitution
  do i=n,1,-1
     do j=1,min(nd,n-i)
        R(:,i) = R(:,i) - A(:,i,j)*R(:,i+j)
     end do
     R(:,i) = R(:,i)/A(:,i,0)
  end do
  
  return
end subroutine polydiagonal_serial


subroutine polydiagonal_periodic_serial(nd,A,R,n,lot,b)

  ! External modules
  use precision

  implicit none

  ! Arguments
  integer,  intent(in) :: n,lot,nd
  real(WP), intent(inout), dimension(lot,n,-nd:nd) :: A
  real(WP), intent(inout), dimension(lot,n)        :: R
  real(WP), intent(inout), dimension(lot,n,nd)     :: b

  ! Local variables
  real(WP), dimension(lot) :: const,pivot
  integer :: i,j,k,i0
  
  if (n.eq.1) then
     const(:) = 0.0_WP
     do j=-nd,nd
        const(:) = const(:) + A(:,j,1)
     end do
     R(:,1) = R(:,1) / const(:)
     return
  end if
  
  ! Setup bridge array
  do j=1,nd
     do i=1,n
        b(:,i,j) = 0.0_WP
     end do
  end do
  do i=1,nd
     do j=i,nd
        b(:,i,j) = A(:,i,-nd+j-i)
        b(:,n-nd-i+1,j-i+1) = A(:,n-nd-i+1,j)
     end do
  end do
  do i=n-nd+1,n
     do j=1,nd
        b(:,i,j) = A(:,i,-nd+j-(i-n))
     end do
  end do
  
  ! Forward elimination
  do i=1,n-nd
     pivot(:) = 1.0_WP/A(:,i,0)
     do j=1,nd
        const(:) = A(:,i+j,-j)*pivot(:)
        do k=1,nd
           A(:,i+j,-j+k) = A(:,i+j,-j+k) - A(:,i,k)*const(:)
           b(:,i+j,k) = b(:,i+j,k) - b(:,i,k)*const(:)
        end do
        R(:,i+j) = R(:,i+j) - R(:,i)*const(:)
     end do
  end do
  
  ! Backward elimination
  do i=n-nd,1,-1
     pivot(:) = 1.0_WP/A(:,i,0)
     do j=1,min(nd,i-1)
        const(:) = A(:,i-j,j)*pivot(:)
        do k=1,nd
           b(:,i-j,k) = b(:,i-j,k) - b(:,i,k)*const(:)
        end do
        R(:,i-j) = R(:,i-j) - R(:,i)*const(:)
     end do
  end do
  
  ! Eliminate oddball region
  do i=1,nd
     pivot(:) = 1.0_WP/A(:,i,0)
     do j=i,nd
        const(:) = A(:,n-j+i,j)*pivot(:)
        do k=1,nd
           b(:,n-j+i,k) = b(:,n-j+i,k) - b(:,i,k)*const(:)
        end do
        R(:,n-j+i) = R(:,n-j+i) - R(:,i)*const(:)
     end do
  end do
  
  ! Elimination for corner matrix
  i0 = n-nd
  do i=1,nd-1
     pivot(:) = 1.0_WP/b(:,i+i0,i)
     do j=i+1,nd
        const(:) = b(:,j+i0,i)*pivot(:)
        do k=i+1,nd
           b(:,j+i0,k) = b(:,j+i0,k) - b(:,i+i0,k)*const(:)
        end do
        R(:,j+i0) = R(:,j+i0) - R(:,i+i0)*const(:)
     end do
  end do
  
  ! Back-substitution for corner matrix
  i0 = n-nd
  do i=nd,1,-1
     do j=i+1,nd
        R(:,i+i0) = R(:,i+i0) - b(:,i+i0,j)*R(:,j+i0)
     end do
     R(:,i+i0) = R(:,i+i0)/b(:,i+i0,i)
  end do
  
  ! Back-substitution for bridge
  do i=n-nd,1,-1
     do j=1,nd
        R(:,i) = R(:,i) - b(:,i,j)*R(:,n-nd+j)
     end do
     R(:,i) = R(:,i)/A(:,i,0)
  end do
  
  return
end subroutine polydiagonal_periodic_serial
