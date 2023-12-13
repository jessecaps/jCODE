module stat_1d

  ! External modules
  use stat
  use string
  
  implicit none

  ! Global variables
  integer :: stat_nvar, statDirection, nStatPoints
  character(len=str_medium), dimension(:), pointer :: stat_name
  real(WP), dimension(:), allocatable :: grid1D, area1D
  real(WP), dimension(:,:), allocatable :: stat1D
  real(WP), dimension(:,:), allocatable :: buf1D
  real(WP) :: Delta_t
  
end module stat_1d


! ======================== !
! Initialize 1D statistics !
! ======================== !
subroutine stat_1d_setup

  ! Internal modules
  use stat_1d

  ! External modules
  use parser
  use parallel
  use simulation_flags
  use geometry
  use grid_functions

  implicit none

  ! Local variables
  integer :: i, j, k, l, ijk(3), gridIndex
  logical :: fileIsThere    
  
  ! Return if init already done
  if (associated(stat_name)) return
  
  ! Allocate
  stat_nvar = 10
  if (twoWayCoupling) stat_nvar = stat_nvar + 1
  allocate(stat_name(stat_nvar))

  ! Choose direction to average
  call parser_read('stat direction', statDirection)
  if (statDirection.lt.1 .or. statDirection.gt.nDimensions)                                  &
       call die('stat_1d_setup: incorrect stat direction')
 
  ! Velocity statistics
  i = 0
  if (twoWayCoupling) then
     i = 1
     stat_name(1) = 'alpha'
  end if
  stat_name(i+1 ) = 'rho'
  stat_name(i+2 ) = 'rhoU'
  stat_name(i+3 ) = 'rhoV'
  stat_name(i+4 ) = 'rhoW'
  stat_name(i+5 ) = 'rhoUU'
  stat_name(i+6 ) = 'rhoVV'
  stat_name(i+7 ) = 'rhoWW'
  stat_name(i+8 ) = 'rhoUV'
  stat_name(i+9 ) = 'rhoVW'
  stat_name(i+10) = 'rhoWU'
  
  ! Simplifty
  nStatPoints = globalGridSize(statDirection)

  ! Allocate stat arrays
  allocate(stat1D(1:nStatPoints, stat_nvar)); stat1D = 0.0_WP
  allocate(buf1D(1:nStatPoints, stat_nvar)); buf1D = 0.0_WP

  ! Store the grid and local surface area before hand
  allocate(grid1D(1:nStatPoints)); grid1D = 0.0_WP
  allocate(area1D(1:nStatPoints)); area1D = 0.0_WP
  do k = iStart(3), iEnd(3)
     do j = iStart(2), iEnd(2)
        do i = iStart(1), iEnd(1)   
           gridIndex = grid_index(i, j, k)
           ijk(1) = i; ijk(2) = j; ijk(3) = k
           l = ijk(statDirection)
           area1D(l) = area1D(l) + sum(metrics(gridIndex,1+nDimensions*(statDirection-1):    &
                nDimensions*statDirection))
           grid1D(l) = coordinates(gridIndex, statDirection)
        end do
     end do
  end do
  call parallel_sum_dir(grid1D, statDirection)
  call parallel_sum(area1D)
  
  ! Open file if it exists
  inquire(file='stat/stat-1D', exist = fileIsThere)
  if (fileIsThere) then
     call stat_1d_read
  else
     Delta_t = 0.0_WP
  end if
 
  return
end subroutine stat_1d_setup


! ===================== !
! Cleanup 1D statistics !
! ===================== !
subroutine stat_1d_cleanup

  ! Internal modules
  use stat_1d
  
  implicit none

  if (allocated(area1D)) deallocate(area1D)
  if (allocated(stat1D)) deallocate(stat1D)
  if (allocated(buf1D)) deallocate(buf1D)
  if (allocated(grid1D)) deallocate(grid1D)

  return
end subroutine stat_1d_cleanup


! ===================== !
! Compute 1D statistics !
! ===================== !
subroutine stat_1d_compute

  ! Internal modules
  use stat_1d

  ! External modules
  use parallel
  use geometry
  use grid_functions
  use state
  use time_info

  implicit none

  ! Local variables
  integer  :: i, j, k, l, n, ijk(3), gridindex
  real(WP) :: norm

  Delta_t = Delta_t + timeStepSize

  ! Gather the 1D stats
  do k = iStart(3), iEnd(3)
     do j = iStart(2), iEnd(2)
        do i = iStart(1), iEnd(1)   
           gridIndex = grid_index(i, j, k)
           ijk(1) = i; ijk(2) = j; ijk(3) = k
           l = ijk(statDirection)
           n = 0

           ! Get norm
           norm = timeStepSize * sum(metrics(gridIndex,1+nDimensions*(statDirection-1):      &
                nDimensions*statDirection))
           
           ! alpha
           if (twoWayCoupling) then
              stat1D(l, 1) = stat1D(l, 1) + norm * volumeFraction(gridIndex,1)
              n = n + 1
           end if
           ! rho
           stat1D(l, n+1) = stat1D(l, n+1) + norm * conservedVariables(gridIndex,1)
           n = n + 1
           ! rhoU, rhoV, rhoW
           stat1D(l,n+1) = stat1D(l,n+1) + norm * conservedVariables(gridIndex,2)
           if (nDimensions.gt.1) stat1D(l,n+2) = stat1D(l,n+2) + norm *                      &
                conservedVariables(gridIndex,3)
           if (nDimensions.gt.2) stat1D(l,n+3) = stat1D(l,n+3) + norm *                      &
                conservedVariables(gridIndex,4)
           n = n + 3
           ! rhoUU, rhoVV, rhoWW
           stat1D(l,n+1) = stat1D(l,n+1) + norm *                                            &
                conservedVariables(gridIndex,2) * velocity(gridIndex, 1)
           if (nDimensions.gt.1) stat1D(l,n+2) = stat1D(l,n+2) + norm *                      &
                conservedVariables(gridIndex,3) * velocity(gridIndex, 2)
           if (nDimensions.gt.2) stat1D(l,n+3) = stat1D(l,n+3) + norm *                      &
                conservedVariables(gridIndex,4) * velocity(gridIndex, 3)
           n = n  + 3
           ! rhoUV, rhoVW, rhoWV
           if (nDimensions.gt.1) stat1D(l,n+1) = stat1D(l,n+1) + norm *                      &
                conservedVariables(gridIndex,2) * velocity(gridIndex, 2)
           if (nDimensions.gt.2) stat1D(l,n+2) = stat1D(l,n+2) + norm *                      &
                conservedVariables(gridIndex,3) * velocity(gridIndex, 3)
           if (nDimensions.gt.2) stat1D(l,n+3) = stat1D(l,n+3) + norm *                      &
                conservedVariables(gridIndex,4) * velocity(gridIndex, 2)
           n = n  + 3
        end do
     end do
  end do

  return
end subroutine stat_1d_compute


! ================== !
! Read 1D statistics !
! ================== !
subroutine stat_1d_read

  ! Internal modules
  use stat_1d

  ! External modules
  use parallel
  use geometry

  implicit none

  ! Local variables
  real(WP) :: time_
  integer :: i
  integer, dimension(2) :: dims
  integer, dimension(MPI_STATUS_SIZE) :: status
  character(len=str_medium)  :: name
  character(len=str_medium) :: filename
  integer :: var,ierror,ifile
  
  ! Open the file to write
  filename = trim(mpiiofs) // "stat/stat-1D"
  call MPI_FILE_OPEN(comm,filename,MPI_MODE_RDONLY,mpiInfo,ifile,ierror)
  
  ! Read dimensions from header
  call MPI_FILE_READ_ALL(ifile,dims,2,MPI_INTEGER,status,ierror)
  if (dims(1).ne.nStatPoints) then
     print*, 'expected = ', nStatPoints
     print*, 'stat = ', dims(1)
     call die('stat_1d_read: The size of the stat file is incorrect')
  end if
  if (dims(2).ne.stat_nvar) call die('stat_1d_read: Wrong number of variables in stat file')
  
  ! Read some headers
  call MPI_FILE_READ_ALL(ifile,Delta_t,1,MPI_REAL_WP,status,ierror)
  call MPI_FILE_READ_ALL(ifile,time_,1,MPI_REAL_WP,status,ierror)
  
  ! Read variable names
  do var=1,stat_nvar
     call MPI_FILE_READ_ALL(ifile,name,str_medium,MPI_CHARACTER,status,ierror)
     if (name.ne.stat_name(var)) then
        call die('stat_1d_read: Variables names in stat and data files do not correspond')
     end if
  end do
  
  ! Read
  call MPI_FILE_READ_ALL(ifile, buf1D, nStatPoints*stat_nvar, MPI_REAL_WP, status, ierror)

  ! Close the file
  call MPI_FILE_CLOSE(ifile,ierror)
  
  ! Restart the stats
  do i = 1, stat_nvar
     stat1D(iStart(statDirection):iEnd(statDirection),i) =                                   &
          buf1D(iStart(statDirection):iEnd(statDirection),i) *                               &
          area1D(iStart(statDirection):iEnd(statDirection)) * Delta_t
  end do

  return
end subroutine stat_1d_read


! ================== !
! Writ 1D statistics !
! ================== !
subroutine stat_1d_write

  ! Internal modules
  use stat_1d

  ! External modules
  use time_info
  use parallel
  use fileio
  use grid_functions  

  implicit none

  ! Local variables
  character(len=str_short) :: dir
  character(len=str_medium) :: filename
  integer  :: i, iunit, ierror, var
  
  ! Nothing to write
  if (Delta_t .le. epsilon(1.0_WP)) return

  ! Identify direction
  select case (statDirection)
  case (1)
     dir = 'x'
  case (2)
     dir = 'y'
  case (3)
     dir = 'z'
  end select

  ! Gather the data
  call parallel_sum(stat1D, buf1D)

  ! Root process writes
  if (irank.eq.iroot) then

     do i = 1, stat_nvar
        buf1D(:,i) = buf1D(:,i) /area1D(:) / Delta_t
     end do
     
     ! ------- WRITE THE BINARY FILE ------
     ! Open the file
     filename = "stat/stat-1D"
     call BINARY_FILE_OPEN(iunit,trim(filename),"w",ierror)

     ! Write dimensions
     call BINARY_FILE_WRITE(iunit, nStatPoints, 1, kind(nStatPoints), ierror)
     call BINARY_FILE_WRITE(iunit, stat_nvar, 1, kind(stat_nvar), ierror)
     call BINARY_FILE_WRITE(iunit, Delta_t, 1, kind(Delta_t), ierror)
     call BINARY_FILE_WRITE(iunit, time, 1, kind(time), ierror)

     ! Write variable names
     do var=1,stat_nvar
        call BINARY_FILE_WRITE(iunit, stat_name(var), str_medium, kind(stat_name), ierror)
     end do

     ! Write the stats
     call BINARY_FILE_WRITE(iunit, buf1D, nStatPoints*stat_nvar, kind(buf1D), ierror)

     ! Close the file
     call BINARY_FILE_CLOSE(iunit,ierror)

     ! ------- WRITE THE ASCII FILE ------
     iunit = iopen()
     write(filename,'(a28)') "stat/stat-1D.txt"
     open (iunit, file=adjustl(trim(filename)), form="formatted",iostat=ierror)
     write(iunit,'(10000a20)') trim(dir),(trim(adjustl(stat_name(var))),var=1,stat_nvar)
     do i = 1, nStatPoints
        write(iunit,'(10000ES20.12)') grid1D(i), buf1D(i,:)
     end do
     close(iclose(iunit)) 
  end if

  ! Log
  call monitor_log("1D STATISTICS FILE WRITTEN")
  
  return
end subroutine stat_1d_write
