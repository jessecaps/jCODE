module levelset_io

  ! External modules
  use string
  use parallel
  use grid_levelset

  implicit none

end module levelset_io


! ============================= !
! Read the levelset in parallel !
! ============================= !
subroutine levelset_read_parallel(filename_)

  ! Internal modules
  use levelset_io

  ! External modules
  use geometry
  use solver_options

  implicit none

  ! Arguments
  character(len = *), intent(in) :: filename_

  ! Local variables
  integer :: ifile, ierror
  integer, dimension(MPI_STATUS_SIZE) :: status
  integer, dimension(4) :: dims
  integer(kind=MPI_Offset_kind) :: disp
  character(len=str_medium) :: filename

  ! Prefix the filename to prevent potential file locking
  filename = trim(mpiiofs) // trim(filename_)

  ! Open the file
  call MPI_FILE_OPEN(comm, filename, MPI_MODE_RDONLY, mpiInfo, ifile, ierror)

  ! Check if there was an error opening the file
  if (ierror.ne.0) call die('Error reading ' // trim(filename_))

  ! Read dimensions from header
  call MPI_FILE_READ_ALL(ifile, dims, 4, MPI_INTEGER, status, ierror)
  if ((dims(1).ne.globalGridSize(1)) .or. (dims(2).ne.globalGridSize(2)) .or.                &
       (dims(3).ne.globalGridSize(3))) then
     print*, 'grid = ', globalGridSize(1), globalGridSize(2), globalGridSize(3)
     print*, trim(filename_), ' = ',dims(1), dims(2), dims(3)
     call die('The size of the levelset file does not correspond to the grid file')
  end if

  ! Read the levelset
  disp = 4*4
  call MPI_FILE_SET_VIEW(ifile, disp, MPI_REAL_WP, mpiDerivedTypeRealSubarray,               &
       "native", mpiInfo, ierror)
  call MPI_FILE_READ_ALL(ifile, levelset(:,1), nGridPoints, MPI_REAL_WP, status, ierror)

  ! Close the file
  call MPI_FILE_CLOSE(ifile, ierror)  

  return
end subroutine levelset_read_parallel


! =========================== !
! Read the levelset in serial !
! =========================== !
subroutine levelset_read_serial(filename_)

  ! Internal modules
  use levelset_io

  ! External modules
  use fileio
  use geometry
  use solver_options

  implicit none

  ! Arguments
  character(len = *), intent(in) :: filename_

  ! Local variables
  integer :: ifile, ierror, ibuf
  integer, dimension(5) :: dims
  character(len=str_medium) :: filename
  real(WP) :: rbuf
  logical :: fileIsThere

  ! Get the name of the header file to read in
  filename = trim(filename_) // '/levelset.header'

  ! If file is not there revert to parallel read
  inquire(file = filename, exist = fileIsThere)
  if (.not. fileIsThere) then
     call levelset_read_parallel(trim(filename_))
     return
  end if

  ! Root reads the header file
  if (irank.eq.iroot) then
     ! Open the file
     call BINARY_FILE_OPEN(ifile, trim(filename), "r", ierror)
     ! Read dimensions from header
     call BINARY_FILE_READ(ifile, dims, 5, kind(dims), ierror)
     if (dims(1).ne.nProcs) then
        print*, 'Expected ', dims(1), ' processor(s)'
        print*, 'Currently ', nProcs, ' processor(s)'
        call die('Number of processors incompatible with serial solution file')
     end if
     if ((dims(2).ne.globalGridSize(1)) .or. (dims(3).ne.globalGridSize(2)) .or.             &
          (dims(4).ne.globalGridSize(3))) then
        print*, 'levelset = ',dims(2), dims(3), dims(4)
        print*, 'grid = ',globalGridSize(1), globalGridSize(2), globalGridSize(3)
        call die('The size of the levelset file does not correspond to the grid file')
     end if
     ! Read additional stuff
     call BINARY_FILE_READ(ifile, rbuf, 1, kind(rbuf), ierror)
     call BINARY_FILE_READ(ifile, ibuf, 1, kind(ibuf), ierror)
     ! Close the file
     call BINARY_FILE_CLOSE(ifile, ierror)
  end if

  ! Build actual filename to read for each processor
  filename = trim(filename_) // '/levelset.'
  write(filename(len_trim(filename)+1:len_trim(filename)+6),'(i6.6)') irank

  ! Open each file and read the levelset
  call BINARY_FILE_OPEN(ifile, trim(filename), "r", ierror)
  call BINARY_FILE_READ(ifile, dims(1:3), 3, kind(dims), ierror)
  if (localGridSize(1).ne.dims(1)) call die('wrong size in solution read serial')
  if (localGridSize(2).ne.dims(2)) call die('wrong size in solution read serial')
  if (localGridSize(3).ne.dims(3)) call die('wrong size in solution read serial')
  call BINARY_FILE_READ(ifile, levelset(:,1), nGridPoints, kind(levelset), ierror)

  ! Close the file
  call BINARY_FILE_CLOSE(ifile, ierror)

  return
end subroutine levelset_read_serial


! ========================== !
! Write levelset in parallel !
! ========================== !
subroutine levelset_write_parallel(filename_)

  ! Internal modules
  use levelset_io

  ! External modules
  use geometry
  use solver_options

  implicit none

  ! Arguments
  character(len = *), intent(in) :: filename_

  ! Local variables
  integer :: ifile, ierror
  integer, dimension(MPI_STATUS_SIZE) :: status
  integer, dimension(4) :: dims
  integer(kind=MPI_Offset_kind) :: disp
  character(len=str_medium) :: filename
  logical :: fileIsThere

  ! Prefix the filename to prevent potential file locking
  filename = trim(mpiiofs) // trim(filename_)

  ! Open the file to write
  inquire(file = filename, exist = fileIsThere)
  if (fileIsThere .and. irank.eq.iroot) call MPI_FILE_DELETE(filename, mpiInfo, ierror)
  call MPI_FILE_OPEN(comm, filename, IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), mpiInfo,         &
       ifile, ierror)

  ! Write header
  if (irank.eq.iroot) then
     ! Write dimensions
     dims(1:3) = globalGridSize
     dims(4) = 1
     call MPI_FILE_WRITE(ifile, dims, 4, MPI_INTEGER, status, ierror)
  end if

  ! Write the levelset
  disp = 4*4
  call MPI_FILE_SET_VIEW(ifile, disp, MPI_REAL_WP, mpiDerivedTypeRealSubarray, "native",     &
       mpiInfo, ierror)
  call MPI_FILE_WRITE_ALL(ifile, levelset(:,1), nGridPoints, MPI_REAL_WP, status, ierror)

  ! Close the file
  call MPI_FILE_CLOSE(ifile, ierror)

  return
end subroutine levelset_write_parallel


! ============================ !
! Write the levelset in serial !
! ============================ !
subroutine levelset_write_serial(filename_)

  ! Internal modules
  use levelset_io

  ! External modules
  use fileio
  use geometry
  use solver_options

  implicit none

  ! Arguments
  character(len = *), intent(in) :: filename_

  ! Local variables
  integer :: ifile, ierror
  integer, dimension(5) :: dims
  character(len=str_medium) :: filename, dirname

  ! Set the name of the directory
  dirname = trim(adjustl(filename_))

  ! Create directory
  if (irank.eq.iroot) call CREATE_FOLDER(trim(dirname))
  call MPI_BARRIER(comm, ierror)

  ! Root dumps header file
  if (irank.eq.iroot) then
     ! Create header name
     filename = trim(adjustl(dirname)) // "/" // "levelset.header"
     ! Open file to write
     call BINARY_FILE_OPEN(ifile, trim(adjustl(filename)), "w", ierror)
     ! Write common integer info
     dims(1) = nProcs
     dims(2) = globalGridSize(1)
     dims(3) = globalGridSize(2)
     dims(4) = globalGridSize(3)
     dims(5) = 1
     call BINARY_FILE_WRITE(ifile, dims, 5, kind(dims), ierror)
     ! Close the file
     call BINARY_FILE_CLOSE(ifile, ierror)
  end if

  ! Create filename
  filename = trim(adjustl(dirname)) // "/levelset."
  write(filename(len_trim(filename)+1:len_trim(filename)+6),'(i6.6)') irank

  ! Open the file to write
  call BINARY_FILE_OPEN(ifile, trim(adjustl(filename)), "w", ierror)

  ! Write header
  dims(1) = localGridSize(1)
  dims(2) = localGridSize(2)
  dims(3) = localGridSize(3)
  call BINARY_FILE_WRITE(ifile, dims(1:3), 3, kind(dims), ierror)

  ! Write the levelset
  call BINARY_FILE_WRITE(ifile, levelset(:,1), nGridPoints, kind(levelset), ierror)

  ! Close the file
  call BINARY_FILE_CLOSE(ifile, ierror)

  return
end subroutine levelset_write_serial
