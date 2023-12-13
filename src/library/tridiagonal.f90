subroutine tridiagonal(A,B,C,R,n,lot,dir,per)

  ! External modules
  use precision
  use parallel

  implicit none
  
  ! Direction
  integer, intent(in) :: dir
  ! Size of the problems
  integer, intent(in) :: n
  ! Number of problems
  integer, intent(in) :: lot
  ! Matrix
  real(WP), dimension(lot,n) :: A     ! LOWER
  real(WP), dimension(lot,n) :: B     ! DIAGONAL
  real(WP), dimension(lot,n) :: C     ! UPPER
  real(WP), dimension(lot,n) :: R     ! RHS - RESULT
  ! Local
  real(WP), dimension(lot)   :: const
  real(WP), dimension(lot)   :: r1
  real(WP), dimension(lot)   :: r2
  real(WP), dimension(lot,n) :: s1
  real(WP), dimension(lot,n) :: s2
  ! Periodicity
  logical, intent(in), dimension(3) :: per
  ! Communication
  integer :: proc,rank,ncom,nper
  integer :: nremain,nlot
  real(WP), dimension(:),     allocatable :: sendbuf
  real(WP), dimension(:,:,:), allocatable :: recvbuf1
  real(WP), dimension(:,:),   allocatable :: recvbuf2
  integer,  dimension(:),     allocatable :: ngroup
  ! Stuff
  integer :: i,igroup
  integer :: k1,L,k2,nk
  integer :: ierr
  
  ! Safety check
  if (dir.lt.1 .or. dir .gt.3) stop 'Unknown direction'
  
  ! Get parallel info
  proc = nProcsDir(dir)
  rank = iRankDir(dir)
  ncom = commDir(dir)
  if (per(dir)) then
     nper = 1
  else
     nper = 0
  end if
  
  ! If serial
  if (proc .eq. 1) then
     if (nper.eq.0) then
        call tridiagonal_serial(A,B,C,R,n,lot)
     else
        call tridiagonal_periodic_serial(A,B,C,R,n,lot)
     end if
     return
  end if
  
  ! Partition the lot
  if (lot .lt. proc) stop 'Tridiagonal solver cannot handle so many proc for such a small problem.'
  allocate(ngroup(proc))
  ngroup(:) = lot/proc
  nremain = mod(lot,proc)
  ngroup(1:nremain) = ngroup(1:nremain) + 1
  nlot = ngroup(1)
  allocate(sendbuf(nlot*12*proc))
  allocate(recvbuf1(nlot,6,2*proc))
  allocate(recvbuf2(nlot,2*proc))
  
  ! Initialize boundary values
  s1(:,1) = a(:,1)
  s2(:,n) = c(:,n)
  if (nper .eq. 0) then
     if (rank .eq. 0)      s1(:,1) = 0.0_WP
     if (rank .eq. proc-1) s2(:,n) = 0.0_WP
  end if
  
  ! Forward elimination
  ! Upper boundary in s1(i)
  do i=2,n
     const(:) = a(:,i)/b(:,i-1)
     b(:,i)   = b(:,i) - c(:,i-1)*const(:)
     r(:,i)   = r(:,i) - r(:,i-1)*const(:)
     s1(:,i)  = -s1(:,i-1)*const(:)
  end do

  ! Backward elimination
  ! Lower boundary in s2(i)
  do i=n-1,1,-1
     const(:) = c(:,i)/b(:,i+1)
     r(:,i)   = r(:,i) - r(:,i+1)*const(:)
     s1(:,i)  = s1(:,i) - s1(:,i+1)*const(:)
     s2(:,i)  = -s2(:,i+1)*const(:)
  end do
  
  ! All dependence has been shifted to the boundary elements
  ! Communicate boundary values to root process
  ! and solve reduced pentadiagonal system
  ! Use of pentadiagonal system is more robust than the
  ! reordered (removes zeros) tridiagonal system
  
  ! Send rows of pentadiagonal system
  ! (0, s1, b, 0, s2; r)
  !    (s1, 0, b, s2, 0; r)
  
  L = 1
  k1 = 1
  do igroup=1,proc
     k2 = k1+ngroup(igroup)-1
     nk = ngroup(igroup)
     
     sendbuf(L:L+nk-1) = 0.0_WP      ; L = L + nlot
     sendbuf(L:L+nk-1) = s1(k1:k2,1) ; L = L + nlot
     sendbuf(L:L+nk-1) = b(k1:k2,1)  ; L = L + nlot
     sendbuf(L:L+nk-1) = 0.0_WP      ; L = L + nlot
     sendbuf(L:L+nk-1) = s2(k1:k2,1) ; L = L + nlot
     sendbuf(L:L+nk-1) = r(k1:k2,1)  ; L = L + nlot
     sendbuf(L:L+nk-1) = s1(k1:k2,n) ; L = L + nlot
     sendbuf(L:L+nk-1) = 0.0_WP      ; L = L + nlot
     sendbuf(L:L+nk-1) = b(k1:k2,n)  ; L = L + nlot
     sendbuf(L:L+nk-1) = s2(k1:k2,n) ; L = L + nlot
     sendbuf(L:L+nk-1) = 0.0_WP      ; L = L + nlot
     sendbuf(L:L+nk-1) = r(k1:k2,n)  ; L = L + nlot
     
     k1 = k2 + 1
  end do
  
  ! Gather the boundary data
  call MPI_ALLTOALL(sendbuf,nlot*12,MPI_REAL_WP,recvbuf1,nlot*12,MPI_REAL_WP,ncom,ierr)
  
  ! Clear unused values
  nk = ngroup(rank+1)
  recvbuf1(nk+1:nlot, :, :) = 0.0_WP
  recvbuf1(nk+1:nlot, 3, :) = 1.0_WP
  
  ! Solve reduced systems
  if (nper .eq. 0) then
     call pentadiagonal_serial(&
          recvbuf1(:,1,2:2*proc-1),recvbuf1(:,2,2:2*proc-1),recvbuf1(:,3,2:2*proc-1),&
          recvbuf1(:,4,2:2*proc-1),recvbuf1(:,5,2:2*proc-1),recvbuf1(:,6,2:2*proc-1),&
          2*proc-2,nlot)
  else
     call pentadiagonal_periodic_serial(&
          recvbuf1(:,1,:),recvbuf1(:,2,:),recvbuf1(:,3,:),&
          recvbuf1(:,4,:),recvbuf1(:,5,:),recvbuf1(:,6,:),&
          2*proc,nlot,sendbuf(1),sendbuf(2*proc*nlot+1))
  end if
  
  ! Move solution to first slot
  do i=1,2*proc
     recvbuf2(:,i) = recvbuf1(:,6,i)
  end do
  
  ! Permute the order
  do i=1,proc-1
     const(1:nlot) = recvbuf2(:,2*i)
     recvbuf2(:,2*i) = recvbuf2(:,2*i+1)
     recvbuf2(:,2*i+1) = const(1:nlot)
  end do
  
  ! If periodic, don't forget the end points
  if (nper .eq. 1) then
     const(1:nlot) = recvbuf2(:,1)
     recvbuf2(:,1) = recvbuf2(:,2*proc)
     recvbuf2(:,2*proc) = const(1:nlot)
  end if
  
  ! Scatter back the solution
  call MPI_ALLTOALL(recvbuf2,nlot*2,MPI_REAL_WP,sendbuf,nlot*2,MPI_REAL_WP,ncom,ierr)
  
  L = 1
  k1 = 1
  do igroup=1,proc
     k2 = k1+ngroup(igroup)-1
     nk = k2-k1+1
     
     r1(k1:k2) = sendbuf(L:L+nk-1) ; L = L + nlot
     r2(k1:k2) = sendbuf(L:L+nk-1) ; L = L + nlot
     
     k1 = k2 + 1
  end do
  
  ! Only if not periodic
  if (nper .eq. 0) then
     if (rank .eq. 0)      r1 = 0.0_WP
     if (rank .eq. proc-1) r2 = 0.0_WP
  end if
  
  do i=1,n
     r(:,i) = (r(:,i) - s1(:,i)*r1(:) - s2(:,i)*r2(:))/b(:,i)
  end do
  
  deallocate(sendbuf)
  deallocate(recvbuf2)
  deallocate(recvbuf1)
  deallocate(ngroup)
  
  return
end subroutine tridiagonal


! =========================================== !
! TriDiagonal Solver - Serial Case - Periodic !
! =========================================== !
subroutine tridiagonal_periodic_serial(a,b,c,r,n,lot)

  ! External modules
  use precision

  implicit none

  ! Arguments
  integer, intent(in) :: n,lot
  real(WP), intent(inout), dimension(lot,n) :: a,b,c,r

  ! Local variables
  real(WP), dimension(lot) :: const
  integer :: i

  if (n .eq. 1) then
     r = r/(a + b + c)
     return
  else if (n .eq. 2) then
     ! Solve 2x2 system
     c(:,1) = c(:,1) + a(:,1)
     a(:,2) = a(:,2) + c(:,2)
     const(:) = a(:,2)/b(:,1)
     b(:,2) = b(:,2) - c(:,1)*const(:)
     r(:,2) = r(:,2) - r(:,1)*const(:)
     r(:,2) = r(:,2)/b(:,2)
     r(:,1) = (r(:,1) - c(:,1)*r(:,2))/b(:,1)
     return
  end if

  ! Forward elimination
  do i=2,n-2
     const(:) = a(:,i)/b(:,i-1)
     b(:,i) = b(:,i) - c(:,i-1)*const(:)
     r(:,i) = r(:,i) - r(:,i-1)*const(:)
     ! Boundary is stored in a(i)
     a(:,i) = -a(:,i-1)*const(:)
  end do
  i=n-1
  const(:) = a(:,i)/b(:,i-1)
  b(:,i) = b(:,i) - c(:,i-1)*const(:)
  r(:,i) = r(:,i) - r(:,i-1)*const(:)
  a(:,i) = c(:,i) - a(:,i-1)*const(:)
  i=n
  const(:) = a(:,i)/b(:,i-1)
  r(:,i) = r(:,i) - r(:,i-1)*const(:)
  a(:,i) = b(:,i) - a(:,i-1)*const(:)

  ! Backward elimination
  do i=n-2,1,-1
     const(:) = c(:,i)/b(:,i+1)
     r(:,i) = r(:,i) - r(:,i+1)*const(:)
     a(:,i) = a(:,i) - a(:,i+1)*const(:)
  end do

  ! Eliminate oddball
  const(:) = c(:,n)/b(:,1)
  r(:,n) = r(:,n) - r(:,1)*const(:)
  a(:,n) = a(:,n) - a(:,1)*const(:)

  ! Backward substitution
  r(:,n) = r(:,n)/a(:,n)
  do i=n-1,1,-1
     r(:,i) = (r(:,i) - a(:,i)*r(:,n))/b(:,i)
  end do

  return
end subroutine tridiagonal_periodic_serial


! =============================================== !
! TriDiagonal Solver - Serial Case - Not Periodic !
! =============================================== !
subroutine tridiagonal_serial(a,b,c,r,n,lot)

  ! External modules
  use precision

  implicit none

  ! Arguments
  integer, intent(in) :: n,lot
  real(WP), intent(inout), dimension(lot,n) :: a,b,c,r

  ! Local variables
  real(WP), dimension(lot) :: const
  integer :: i
  
  ! Forward elimination
  do i=2,n
     const(:) = a(:,i)/(b(:,i-1)+tiny(1.0_WP))
     b(:,i) = b(:,i) - c(:,i-1)*const(:)
     r(:,i) = r(:,i) - r(:,i-1)*const(:)
  end do

  ! Back-substitution
  r(:,n) = r(:,n)/b(:,n)
  do i=n-1,1,-1
     r(:,i) = (r(:,i) - c(:,i)*r(:,i+1))/(b(:,i)+tiny(1.0_WP))
  end do
  
  return
end subroutine tridiagonal_serial
