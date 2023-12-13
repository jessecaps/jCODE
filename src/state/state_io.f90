module state_io

  ! External modules
  use string
  use parallel
  use state

  implicit none

end module state_io


! =========================== !
! Read state data in parallel !
! =========================== !
subroutine state_read_parallel(stateVector, filename_)

  ! Internal modules
  use state_io

  ! External modules
  use geometry
  use solver_options
  use time_info

  implicit none

  ! Arguments
  real(WP), dimension(nGridPoints, nUnknowns), intent(out) :: stateVector
  character(len = *), intent(in) :: filename_

  ! Local variables
  integer :: ifile, ierror, nvars, var
  integer, dimension(MPI_STATUS_SIZE) :: status
  integer, dimension(4) :: dims
  integer(kind=MPI_Offset_kind) :: disp
  integer(kind=MPI_Offset_kind) :: nx_MOK, ny_MOK, nz_MOK
  integer(kind=MPI_Offset_kind) :: WP_MOK, var_MOK
  integer(kind=MPI_Offset_kind) :: NVARS_MOK
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
     call die('The size of the solution file does not correspond to the grid file')
  end if

  ! Check that number of variables is consistent
  nvars = dims(4)
  if (nvars .ne. nUnknowns) then
     print *, 'Expected ', nvars, ' unknowns'
     print *, 'Currently ', nUnknowns, ' unknowns'
     call die('The size of the solution file is not consistent with the simulation')
  end if

  ! Read the current time
  call MPI_FILE_READ_ALL(ifile, time, 1, MPI_REAL_WP, status, ierror)
  call MPI_FILE_READ_ALL(ifile, timestep, 1, MPI_INTEGER, status, ierror)

  ! Resize some integers so MPI can read even the biggest files
  nx_MOK    = int(globalGridSize(1), MPI_Offset_kind)
  ny_MOK    = int(globalGridSize(2), MPI_Offset_kind)
  nz_MOK    = int(globalGridSize(3), MPI_Offset_kind)
  WP_MOK    = int(WP,                MPI_Offset_kind)
  NVARS_MOK = int(nUnknowns,         MPI_Offset_kind)

  ! Read the data for each variable
  do var = 1, nUnknowns
     var_MOK = int(var, MPI_Offset_kind)
     disp = 4*4 + WP_MOK + 4 + nx_MOK*ny_MOK*nz_MOK*WP_MOK*(var_MOK-1)
     call MPI_FILE_SET_VIEW(ifile, disp, MPI_REAL_WP, mpiDerivedTypeRealSubarray,            &
          "native", mpiInfo, ierror)
     call MPI_FILE_READ_ALL(ifile, stateVector(:,var), nGridPoints, MPI_REAL_WP,             &
          status, ierror)
  end do

  ! Close the file
  call MPI_FILE_CLOSE(ifile, ierror)

  return
end subroutine state_read_parallel


! ========================= !
! Read state data in serial !
! ========================= !
subroutine state_read_serial(stateVector, filename_)

  ! Internal modules
  use state_io

  ! External modules
  use fileio
  use geometry
  use solver_options
  use time_info

  implicit none

  ! Arguments
  real(WP), dimension(nGridPoints, nUnknowns), intent(out) :: stateVector
  character(len = *), intent(in) :: filename_

  ! Local variables
  integer :: ifile, ierror, nvars, var
  integer, dimension(5) :: dims
  character(len=str_medium) :: filename
  logical :: fileIsThere

  ! Get the name of the header file to read in
  filename = trim(filename_) // '/data.header'

  ! If file is not there revert to parallel read
  inquire(file = filename, exist = fileIsThere)
  if (.not. fileIsThere) then
     call state_read_parallel(stateVector, trim(filename_))
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
        print*, 'data = ',dims(2), dims(3), dims(4)
        print*, 'grid = ',globalGridSize(1), globalGridSize(2), globalGridSize(3)
        call die('The size of the solution file does not correspond to the grid file')
     end if
     nvars = dims(5)
     if (nvars .ne. nUnknowns) then
        print *, 'Expected ', nvars, ' unknowns'
        print *, 'Currently ', nUnknowns, ' unknowns'
        call die('The size of the solution file is not consistent with the simulation')
     end if
     ! Read additional stuff
     call BINARY_FILE_READ(ifile, timestep, 1, kind(timestep), ierror)
     call BINARY_FILE_READ(ifile, time, 1, kind(time), ierror)
     ! Close the file
     call BINARY_FILE_CLOSE(ifile, ierror)
  end if

  ! Broadcast information
  call parallel_bc(timestep)
  call parallel_bc(time)
  call parallel_bc(nvars)

  ! Build actual filename to read for each processor
  filename = trim(filename_) // '/data.'
  write(filename(len_trim(filename)+1:len_trim(filename)+6),'(i6.6)') irank

  ! Open each file and read the data for each variable
  call BINARY_FILE_OPEN(ifile, trim(filename), "r", ierror)
  call BINARY_FILE_READ(ifile, dims(1:3), 3, kind(dims), ierror)
  if (localGridSize(1).ne.dims(1)) call die('wrong size in solution read serial')
  if (localGridSize(2).ne.dims(2)) call die('wrong size in solution read serial')
  if (localGridSize(3).ne.dims(3)) call die('wrong size in solution read serial')
  do var = 1, nUnknowns
     call BINARY_FILE_READ(ifile, stateVector(:,var), nGridPoints, kind(stateVector), ierror)
  end do

  ! Close the file
  call BINARY_FILE_CLOSE(ifile, ierror)

  return
end subroutine state_read_serial


! ============================ !
! Write state data in parallel !
! ============================ !
subroutine state_write_parallel(stateVector, filename_)

  ! Internal modules
  use state_io

  ! External modules
  use geometry
  use solver_options
  use time_info

  implicit none

  ! Arguments
  real(WP), dimension(nGridPoints, nUnknowns), intent(in) :: stateVector
  character(len = *), intent(in) :: filename_

  ! Local variables
  integer :: ifile, ierror, var
  integer, dimension(MPI_STATUS_SIZE) :: status
  integer, dimension(4) :: dims
  integer(kind=MPI_Offset_kind) :: disp
  integer(kind=MPI_Offset_kind) :: nx_MOK, ny_MOK, nz_MOK
  integer(kind=MPI_Offset_kind) :: WP_MOK, var_MOK
  integer(kind=MPI_Offset_kind) :: NVARS_MOK
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
     dims(4) = nUnknowns
     call MPI_FILE_WRITE(ifile, dims, 4, MPI_INTEGER, status, ierror)
     ! Write the current time
     call MPI_FILE_WRITE(ifile, time, 1, MPI_REAL_WP, status, ierror)
     call MPI_FILE_WRITE(ifile, timestep, 1, MPI_INTEGER, status, ierror)
  end if

  ! Resize some integers so MPI can read even the biggest files
  nx_MOK    = int(globalGridSize(1), MPI_Offset_kind)
  ny_MOK    = int(globalGridSize(2), MPI_Offset_kind)
  nz_MOK    = int(globalGridSize(3), MPI_Offset_kind)
  WP_MOK    = int(WP,                MPI_Offset_kind)
  NVARS_MOK = int(nUnknowns,         MPI_Offset_kind)

  ! Write the data for each variable
  do var = 1, nUnknowns
     var_MOK = int(var, MPI_Offset_kind)
     disp = 4*4 + WP_MOK + 4 + nx_MOK*ny_MOK*nz_MOK*WP_MOK*(var_MOK-1)
     call MPI_FILE_SET_VIEW(ifile, disp, MPI_REAL_WP, mpiDerivedTypeRealSubarray, "native",  &
          mpiInfo, ierror)
     call MPI_FILE_WRITE_ALL(ifile, stateVector(:,var), nGridPoints, MPI_REAL_WP,            &
          status, ierror)
  end do

  ! Close the file
  call MPI_FILE_CLOSE(ifile, ierror)

  return
end subroutine state_write_parallel


! ========================== !
! Write state data in serial !
! ========================== !
subroutine state_write_serial(stateVector, filename_)

  ! Internal modules
  use state_io

  ! External modules
  use fileio
  use geometry
  use solver_options
  use time_info

  implicit none

  ! Arguments
  real(WP), dimension(nGridPoints, nUnknowns), intent(in) :: stateVector
  character(len = *), intent(in) :: filename_

  ! Local variables
  integer :: ifile, ierror, var
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
     filename = trim(adjustl(dirname)) // "/" // "data.header"
     ! Open file to write
     call BINARY_FILE_OPEN(ifile, trim(adjustl(filename)), "w", ierror)
     ! Write common integer info
     dims(1) = nProcs
     dims(2) = globalGridSize(1)
     dims(3) = globalGridSize(2)
     dims(4) = globalGridSize(3)
     dims(5) = nUnknowns
     call BINARY_FILE_WRITE(ifile, dims, 5, kind(dims), ierror)
     ! Add time info
     call BINARY_FILE_WRITE(ifile, timestep, 1, kind(timestep), ierror)
     call BINARY_FILE_WRITE(ifile, time, 1, kind(time), ierror)
     ! Close the file
     call BINARY_FILE_CLOSE(ifile, ierror)
  end if

  ! Create filename
  filename = trim(adjustl(dirname)) // "/data."
  write(filename(len_trim(filename)+1:len_trim(filename)+6),'(i6.6)') irank

  ! Open the file to write
  call BINARY_FILE_OPEN(ifile, trim(adjustl(filename)), "w", ierror)

  ! Write header
  dims(1) = localGridSize(1)
  dims(2) = localGridSize(2)
  dims(3) = localGridSize(3)
  call BINARY_FILE_WRITE(ifile, dims(1:3), 3, kind(dims), ierror)

  ! Write the data for each variable
  do var = 1, nUnknowns
     call BINARY_FILE_WRITE(ifile, stateVector(:,var), nGridPoints, kind(stateVector), ierror)
  end do

  ! Close the file
  call BINARY_FILE_CLOSE(ifile, ierror)

  return
end subroutine state_write_serial
