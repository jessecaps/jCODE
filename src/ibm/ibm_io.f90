module ibm_io

  ! External modules
  use ibm
  use parallel

  implicit none

end module ibm_io


! ========================= !
! Read IBM data in parallel !
! ========================= !
subroutine ibm_read_parallel(filename_)

  ! Internal modules
  use ibm_io

  ! External modules
  use math
  use geometry
  use solver_options

  implicit none

  ! Arguments
  character(len = *), intent(in) :: filename_

  ! Local variables
  integer :: objectSize, ifile, ierror
  integer (kind=MPI_OFFSET_KIND) :: offset
  real(WP) :: buf
  character(len=str_medium) :: filename, message

  ! Prefix the filename to prevent potential file locking
  filename = trim(mpiiofs) // trim(filename_)

  ! Open the file
  call MPI_File_open(comm, trim(filename), MPI_MODE_RDONLY, mpiInfo, ifile, ierror)

  ! Check if there was an error opening the file
  if (ierror.ne.0) call die('Error reading ' // trim(filename_))

  ! Number of objects
  call MPI_File_read_all(ifile, nObjects, 1, MPI_INTEGER, MPI_STATUS_IGNORE, ierror)

  ! Object size
  call MPI_File_read_all(ifile, objectSize, 1, MPI_INTEGER, MPI_STATUS_IGNORE, ierror)

  if (objectSize .ne. SIZE_MPI_OBJECT) then
     write(message, '(3A,I0.0,A,I0.0)') "Problem reading '", trim(filename),                 &
          "': Object size = ", objectSize, ", expected size = ", SIZE_MPI_OBJECT
     call die(trim(message))
  end if

  ! Time
  call MPI_File_read_all(ifile, buf, 1, MPI_REAL_WP, MPI_STATUS_IGNORE, ierror)

  if (allocated(object)) deallocate(object)
  allocate(object(nObjects))

  offset = 4*2 + int(WP, MPI_OFFSET_KIND)

  if (nObjects .gt. 0) call MPI_File_read_at(ifile, offset, object, nObjects, MPI_OBJECT,         &
       MPI_STATUS_IGNORE, ierror)

  ! Close the file
  call MPI_File_close(ifile, ierror)
  
  return
end subroutine ibm_read_parallel


! ======================= !
! Read IBM data in serial !
! ======================= !
subroutine ibm_read_serial(filename_)

  ! Internal modules
  use ibm_io

  ! External modules
  use fileio
  use geometry
  use solver_options
  use time_info

  implicit none

  ! Arguments
  character(len = *), intent(in) :: filename_

  ! Local variables
  integer :: ifile, ierror
  integer, dimension(3) :: dims
  character(len=str_medium) :: filename, message
  logical :: fileIsThere

  ! Get the name of the header file to read in
  filename = trim(filename_) // '/ibm.header'

  ! If file is not there revert to parallel read
  inquire(file = filename, exist = fileIsThere)
  if (.not. fileIsThere) then
     call ibm_read_parallel(trim(filename_))
     return
  end if

  ! Root reads the header file
  if (irank.eq.iroot) then
     ! Open the file
     call BINARY_FILE_OPEN(ifile, trim(filename), "r", ierror)
     ! Read the IBM header
     call BINARY_FILE_READ(ifile, dims, 3, kind(dims), ierror)

     ! Set number of objects
     nObjects = dims(2)

     ! Check object size
     if (dims(3) .ne. SIZE_MPI_OBJECT) then
        write(message, '(3A,I0.0,A,I0.0)') "Problem reading '", trim(filename),              &
             "': Object size = ", dims(3), ", expected size = ", SIZE_MPI_OBJECT
        call die(trim(message))
     end if
     ! Close the file
     call BINARY_FILE_CLOSE(ifile, ierror)
  end if

  ! Broadcast information
  call parallel_bc(nObjects)

  ! Build actual filename to read for each processor
  filename = trim(filename_) // '/ibm'

  ! Open each file and the number of objects for this processor
  call BINARY_FILE_OPEN(ifile, trim(filename), "r", ierror)
  call BINARY_FILE_READ(ifile, nObjects, 1, kind(nObjects), ierror)

  ! Allocate object vector and read
  if (allocated(object)) deallocate(object)
  allocate(object(nObjects))
  if (nObjects.gt.0) call binary_file_read(ifile, object(1:nObjects), nObjects,              &
       SIZE_MPI_OBJECT, ierror)

  ! Close the file
  call BINARY_FILE_CLOSE(ifile, ierror)

  return
end subroutine ibm_read_serial


! ========================== !
! Write IBM data in parallel !
! ========================== !
subroutine ibm_write_parallel(filename_)

  ! Internal modules
  use ibm_io

  ! External modules
  use math, only : pigeon_hole
  use geometry
  use solver_options
  use time_info

  implicit none

  ! Arguments
  character(len = *), intent(in) :: filename_

  ! Local variables
  integer :: i, ifile, nObjects_, objectOffset, ierror
  integer, allocatable :: nObjectsProc(:)
  integer (kind=MPI_OFFSET_KIND) :: offset
  character(len=str_medium) :: filename
  logical :: fileIsThere

  ! Prefix the filename to prevent potentials file locking
  filename = trim(mpiiofs) // trim(filename_)

  ! Distribute the objects evenly among the processors
  call pigeon_hole(nObjects, nProcs, iRank, objectOffset, nObjects_)

  ! Get number of objects per processor in `comm`
  allocate(nObjectsProc(nProcs))
  call MPI_Allgather(nObjects_, 1, MPI_INTEGER, nObjectsProc, 1, MPI_INTEGER, comm, ierror)

  ! Open the file to write
  inquire(file = filename, exist = fileIsThere)
  if (fileIsThere .and. irank.eq.iroot) call MPI_FILE_DELETE(filename, mpiInfo, ierror)
  call MPI_FILE_OPEN(comm, filename, IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), mpiInfo,         &
       ifile, ierror)

  ! Root process writes the header
  if (iRank .eq. iRoot) then
     call MPI_FILE_WRITE(ifile, nObjects, 1, MPI_INTEGER, MPI_STATUS_IGNORE, ierror)
     call MPI_FILE_WRITE(ifile, SIZE_MPI_OBJECT, 1, MPI_INTEGER, MPI_STATUS_IGNORE, ierror)
     call MPI_FILE_WRITE(ifile, time, 1, MPI_REAL_WP, MPI_STATUS_IGNORE, ierror)
  end if

  ! Get the local object offset
  offset = 4*2 + int(WP, MPI_OFFSET_KIND)
  do i = 1, iRank
     offset = offset + int(nObjectsProc(i), MPI_OFFSET_KIND) * int(SIZE_MPI_OBJECT,          &
          MPI_OFFSET_KIND)
  end do

  ! Write local object data
  if (nObjects_ .gt. 0) call MPI_file_write_at(ifile, offset,                                &
       object(objectOffset+1:objectOffset+nObjects_), nObjects_, MPI_OBJECT,                 &
       MPI_STATUS_IGNORE, ierror)

  ! Close the file
  call MPI_File_close(ifile, ierror)

  ! Cleanup
  deallocate(nObjectsProc)

  return
end subroutine ibm_write_parallel


! ======================== !
! Write IBM data in serial !
! ======================== !
subroutine ibm_write_serial(filename_)

  ! Internal modules
  use ibm_io

  ! External modules
  use fileio
  use geometry
  use solver_options
  use time_info

  implicit none

  ! Arguments
  character(len = *), intent(in) :: filename_

  ! Local variables
  integer :: ifile, ierror
  integer, dimension(3) :: dims
  character(len=str_medium) :: filename, dirname

  ! Only root processor writes
  if (iRank .ne. iRoot) return
  
  ! Set the name of the directory
  dirname = trim(adjustl(filename_))

  ! Create directory
  call CREATE_FOLDER(trim(dirname))

  ! Create header name
  filename = trim(adjustl(dirname)) // "/" // "ibm.header"
  ! Open file to write
  call BINARY_FILE_OPEN(ifile, trim(adjustl(filename)), "w", ierror)
  ! Write common integer info
  dims(1) = 1
  dims(2) = nObjects
  dims(3) = SIZE_MPI_OBJECT
  call BINARY_FILE_WRITE(ifile, dims, 3, kind(dims), ierror)
  ! Add time info
  call BINARY_FILE_WRITE(ifile, time, 1, kind(time), ierror)
  ! Close the file
  call BINARY_FILE_CLOSE(ifile, ierror)

  ! Create filename
  filename = trim(adjustl(dirname)) // "/ibm"

  ! Open the file to write
  call BINARY_FILE_OPEN(ifile, trim(adjustl(filename)), "w", ierror)

  ! Write header
  call BINARY_FILE_WRITE(ifile, nObjects, 1, kind(nObjects), ierror)

  ! Write the IBMobject data
  if (nObjects.gt.0) call BINARY_FILE_WRITE(ifile, object(1:nObjects), nObjects,             &
       SIZE_MPI_OBJECT, ierror)

  ! Close the file
  call BINARY_FILE_CLOSE(ifile, ierror)

  return
end subroutine ibm_write_serial
