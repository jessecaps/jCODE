module grid_io

  ! External modules
  use string
  use parallel
  use grid

  implicit none

end module grid_io


! ========================== !
! Read grid data in parallel !
! ========================== !
subroutine grid_read_parallel(filename_)

  ! Internal modules
  use grid_io

  ! External modules
  use geometry
  use solver_options

  implicit none

  ! Arguments
  character(len = *), intent(in) :: filename_

  ! Local variables
  integer :: i, ifile, ierror
  integer, dimension(MPI_STATUS_SIZE) :: status
  integer, dimension(8) :: dims
  integer(kind=MPI_Offset_kind) :: disp
  integer(kind=MPI_Offset_kind) :: nx_MOK, ny_MOK, nz_MOK
  integer(kind=MPI_Offset_kind) :: WP_MOK
  real(WP) :: buf
  character(len=str_medium) :: filename

  ! Prefix the filename to prevent potential file locking
  filename = trim(mpiiofs) // trim(filename_)

  ! Resize some integers so MPI can read even the biggest files
  nx_MOK    = int(globalGridSize(1), MPI_Offset_kind)
  ny_MOK    = int(globalGridSize(2), MPI_Offset_kind)
  nz_MOK    = int(globalGridSize(3), MPI_Offset_kind)
  WP_MOK    = int(WP,                MPI_Offset_kind)

  ! Open the file
  call MPI_FILE_OPEN(comm, filename, MPI_MODE_RDONLY, mpiInfo, ifile, ierror)

  ! Check if there was an error opening the file
  if (ierror.ne.0) call die('Error reading ' // trim(filename_))

  ! Read dimensions and periodicity from header
  ! dims(1:3) => globalGridSize
  ! dims(4)   => isDomainCurvilinear
  ! dims(5)   => useIblank
  ! dims(6:8) => periodicityType (NONE=0, PLANE=1, OVERLAP=2, POLAR=3)
  call MPI_FILE_READ_ALL(ifile, dims, 8, MPI_INTEGER, status, ierror)
  disp = 8*4
  do i = 1, nDimensions
     if (dims(i+5).eq.PLANE) then
        call MPI_FILE_READ_ALL(ifile, buf, 1, MPI_REAL_WP, status, ierror)
        disp = disp + WP_MOK
     end if
  end do

  ! Read the coordinates for each dimension
  do i = 1, nDimensions
     call MPI_FILE_SET_VIEW(ifile, disp, MPI_REAL_WP, mpiDerivedTypeRealSubarray,            &
          "native", mpiInfo, ierror)
     call MPI_FILE_READ_ALL(ifile, coordinates(:,i), nGridPoints, MPI_REAL_WP,               &
          status, ierror)
     disp = disp + nx_MOK*ny_MOK*nz_MOK*WP_MOK
  end do

  ! Read the Iblank array
  if (useIblank) then
     call MPI_FILE_SET_VIEW(ifile, disp, MPI_INTEGER, mpiDerivedTypeIntegerSubarray,         &
          "native", mpiInfo, ierror)
     call MPI_FILE_READ_ALL(ifile, iblank(:), nGridPoints, MPI_INTEGER, status, ierror)
  else
     iblank = 1
  end if

  ! Close the file
  call MPI_FILE_CLOSE(ifile, ierror)

  return
end subroutine grid_read_parallel


! ======================== !
! Read grid data in serial !
! ======================== !
subroutine grid_read_serial(filename_)

  ! Internal modules
  use grid_io

  ! External modules
  use fileio
  use geometry
  use solver_options

  implicit none

  ! Arguments
  character(len = *), intent(in) :: filename_

  ! Local variables
  integer :: i, ifile, ierror
  integer, dimension(9) :: dims
  character(len=str_medium) :: filename
  logical :: fileIsThere

  ! Get the name of the header file to read in
  filename = trim(filename_) // '/grid.header'

  ! If file is not there revert to parallel read
  inquire(file = filename, exist = fileIsThere)
  if (.not. fileIsThere) then
     call grid_read_parallel(trim(filename_))
     return
  end if

  ! Root reads the header file
  ! dims(1)   => nProcs
  ! dims(2:4) => globalGridSize
  ! dims(5)   => isDomainCurvilinear
  ! dims(6)   => useIblank
  ! dims(7:9) => periodicityType (NONE=0, PLANE=1, OVERLAP=2, POLAR=3)
  if (irank.eq.iroot) then
     ! Open the file
     call BINARY_FILE_OPEN(ifile, trim(filename), "r", ierror)
     call BINARY_FILE_READ(ifile, dims, 9, kind(dims), ierror)
     if (dims(1).ne.nProcs) then
        print*, 'Expected ', dims(1), ' processor(s)'
        print*, 'Currently ', nProcs, ' processor(s)'
        call die('Number of processors incompatible with serial grid file')
     end if
     ! Close the file
     call BINARY_FILE_CLOSE(ifile, ierror)
  end if

  ! Build actual filename to read for each processor
  filename = trim(filename_) // '/grid.'
  write(filename(len_trim(filename)+1:len_trim(filename)+6),'(i6.6)') irank

  ! Open each grid file and read the grid data
  call BINARY_FILE_OPEN(ifile, trim(filename), "r", ierror)
  call BINARY_FILE_READ(ifile, dims(1:3), 3, kind(dims), ierror)
  if (localGridSize(1).ne.dims(1)) call die('wrong size in grid read serial')
  if (localGridSize(2).ne.dims(2)) call die('wrong size in grid read serial')
  if (localGridSize(3).ne.dims(3)) call die('wrong size in grid read serial')
  do i = 1, nDimensions
     call BINARY_FILE_READ(ifile, coordinates(:,i), nGridPoints, kind(coordinates), ierror)
  end do
  if (useIblank) call BINARY_FILE_READ(ifile, iblank, nGridPoints, kind(iblank), ierror)

  ! Close the file
  call BINARY_FILE_CLOSE(ifile, ierror)

  return
end subroutine grid_read_serial


! =========================== !
! Write grid data in parallel !
! =========================== !
subroutine grid_write_parallel(filename_)

  ! Internal modules
  use grid_io

  ! External modules
  use geometry
  use solver_options

  implicit none

  ! Arguments
  character(len = *), intent(in) :: filename_

  ! Local variables
  integer :: i, useIb, ifile, ierror
  integer, dimension(MPI_STATUS_SIZE) :: status
  integer, dimension(8) :: dims
  integer(kind=MPI_Offset_kind) :: disp
  integer(kind=MPI_Offset_kind) :: nx_MOK, ny_MOK, nz_MOK
  integer(kind=MPI_Offset_kind) :: WP_MOK
  character(len=str_medium) :: filename
  logical :: fileIsThere

  ! Prefix the filename to prevent potential file locking
  filename = trim(mpiiofs) // trim(filename_)

  ! Resize some integers so MPI can read even the biggest files
  nx_MOK    = int(globalGridSize(1), MPI_Offset_kind)
  ny_MOK    = int(globalGridSize(2), MPI_Offset_kind)
  nz_MOK    = int(globalGridSize(3), MPI_Offset_kind)
  WP_MOK    = int(WP,                MPI_Offset_kind)

  ! Determine if iblank should be written
  useIb = 0
  if (allocated(iblank)) then
     if (minval(iblank).ne.1 .or. maxval(iblank).ne.1) useIb = 1
  end if
  call parallel_max(useIb)

  ! Open the file to write
  inquire(file = filename, exist = fileIsThere)
  if (fileIsThere .and. irank.eq.iroot) call MPI_FILE_DELETE(filename, mpiInfo, ierror)
  call MPI_FILE_OPEN(comm, filename, IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), mpiInfo,         &
       ifile, ierror)

  ! Root process writes the header
  if (irank.eq.iRoot) then
     ! dims(1:3) => globalGridSize
     ! dims(4)   => isDomainCurvilinear
     ! dims(5)   => useIblank
     ! dims(6:8) => periodicityType (NONE=0, PLANE=1, OVERLAP=2, POLAR=3)
     dims(1:3) = globalGridSize
     if (isDomainCurvilinear) then
        dims(4) = 1
     else
        dims(4) = 0
     end if
     dims(5)   = useIb
     dims(6:8) = periodicityType
     call MPI_FILE_WRITE(ifile, dims, 8, MPI_INTEGER, status, ierror)
     disp = 8*4
     do i = 1, nDimensions
        if (dims(i+5).eq.PLANE) then
           call MPI_FILE_WRITE(ifile, periodicLength(i), 1, MPI_REAL_WP, status, ierror)
           disp = disp + WP_MOK
        end if
     end do
  end if

  ! Communicate displacement
  i = int(disp, kind=4)
  call parallel_bc(i)
  disp = int(i, MPI_Offset_kind)

  ! Write the coordinates for each dimension
  do i = 1, nDimensions
     call MPI_FILE_SET_VIEW(ifile, disp, MPI_REAL_WP, mpiDerivedTypeRealSubarray, "native",  &
          mpiInfo, ierror)
     call MPI_FILE_WRITE_ALL(ifile, coordinates(:,i), nGridPoints, MPI_REAL_WP,              &
          status, ierror)
     disp = disp + nx_MOK*ny_MOK*nz_MOK*WP_MOK
  end do

  ! Write the iblank array
  if (useIb.eq.1) then
     call MPI_FILE_SET_VIEW(ifile, disp, MPI_INTEGER, mpiDerivedTypeIntegerSubarray,         &
          "native", mpiInfo, ierror)
     call MPI_FILE_WRITE_ALL(ifile, iblank(:), nGridPoints, MPI_INTEGER, status, ierror)
  end if

  ! Close the file
  call MPI_FILE_CLOSE(ifile, ierror)

  return
end subroutine grid_write_parallel


! ========================= !
! Write grid data in serial !
! ========================= !
subroutine grid_write_serial(filename_)

  ! Internal modules
  use grid_io

  ! External modules
  use fileio
  use geometry
  use solver_options

  implicit none

  ! Arguments
  character(len = *), intent(in) :: filename_

  ! Local variables
  integer :: i, useIb, ifile, ierror
  integer, dimension(9) :: dims
  character(len=str_medium) :: filename, dirname

  ! Determine if iblank should be written
  useIb = 0
  if (allocated(iblank)) then
     if (minval(iblank).ne.1 .or. maxval(iblank).ne.1) useIb = 1
  end if
  call parallel_max(useIb)

  ! Set the name of the directory
  dirname = trim(adjustl(filename_))

  ! Create directory
  if (irank.eq.iroot) call CREATE_FOLDER(trim(dirname))
  call MPI_BARRIER(comm, ierror)

  ! Root dumps header file
  if (irank.eq.iroot) then
     ! Create header name
     filename = trim(adjustl(dirname)) // "/" // "grid.header"
     ! Open file to write
     call BINARY_FILE_OPEN(ifile, trim(adjustl(filename)), "w", ierror)
     ! Write common integer info
     dims(1) = nProcs
     dims(2) = globalGridSize(1)
     dims(3) = globalGridSize(2)
     dims(4) = globalGridSize(3)
     if (isDomainCurvilinear) then
        dims(5) = 1
     else
        dims(5) = 0
     end if
     dims(6) = useIb
     dims(7) = periodicityType(1)
     dims(8) = periodicityType(2)
     dims(9) = periodicityType(3)
     call BINARY_FILE_WRITE(ifile, dims, 9, kind(dims), ierror)
     ! Add periodic length
     do i = 1, nDimensions
        if (periodicityType(i).eq.PLANE) call BINARY_FILE_WRITE(ifile, periodicLength(i), 1, &
             kind(periodicLength(i)), ierror)
     end do
     ! Close the file
     call BINARY_FILE_CLOSE(ifile, ierror)
  end if

  ! Create filename
  filename = trim(adjustl(dirname)) // "/grid."
  write(filename(len_trim(filename)+1:len_trim(filename)+6),'(i6.6)') irank

  ! Open the file to write
  call BINARY_FILE_OPEN(ifile, trim(adjustl(filename)), "w", ierror)

  ! Write header
  dims(1) = localGridSize(1)
  dims(2) = localGridSize(2)
  dims(3) = localGridSize(3)
  call BINARY_FILE_WRITE(ifile, dims(1:3), 3, kind(dims), ierror)

  ! Write the coordinates for each dimension
  do i = 1, nDimensions
     call BINARY_FILE_WRITE(ifile, coordinates(:,i), nGridPoints, kind(coordinates), ierror)
  end do

  ! Write the iblank array
  if (useIb.eq.1) call BINARY_FILE_WRITE(ifile, iblank, nGridPoints, kind(iblank), ierror)

  ! Close the file
  call BINARY_FILE_CLOSE(ifile, ierror)

  return
end subroutine grid_write_serial
