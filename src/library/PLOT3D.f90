module PLOT3D

  use string

  implicit none
  public

  integer, parameter, public ::                                                              &
       PLOT3D_GRID_FILE     = 0,                                                             &
       PLOT3D_SOLUTION_FILE = 1,                                                             &
       PLOT3D_FUNCTION_FILE = 2

  type, public :: t_PLOT3DDescriptor

     integer :: fileType, nGrids, nDimensions, nScalars
     logical :: hasIblank = .true., isEndiannessNative = .true., isCurvilinear = .true.

  end type t_PLOT3DDescriptor

  interface swap_endianness
     module procedure swap_integer_endianness_, swap_scalar_endianness_
  end interface swap_endianness

  interface plot3d_get_offset
     module procedure plot3d_get_offset_from_grid_sizes_, plot3d_get_offset_from_file_
  end interface plot3d_get_offset

  character(len = str_medium) :: plot3dErrorMessage

  private :: swap_endianness, swap_integer_endianness_, swap_scalar_endianness_,             &
       plot3d_get_offset_from_grid_sizes_, plot3d_get_offset_from_file_

contains

  subroutine plot3d_detect_format(filename, success, descriptor, globalGridSizes,            &
       includeFunctionFiles)

    ! Detects the format of a PLOT3D file `filename` and reads the grid size

    ! External modules
    use string
    use parallel
    use, intrinsic :: iso_c_binding

    implicit none

    ! Arguments
    character(len = *), intent(in) :: filename
    logical, intent(out) :: success
    type(t_PLOT3DDescriptor), intent(out), optional :: descriptor
    integer, allocatable, intent(out), optional :: globalGridSizes(:,:)
    logical, intent(in), optional :: includeFunctionFiles

    interface

       function detect_format(filename, includeFunctionFiles, nGrids, nDimensions, nScalars, &
            hasIblank, isEndiannessNative, fileType, globalGridSizes)                        &
            bind(C, name = "detect_format")
         use, intrinsic :: iso_c_binding
         use string
         integer(kind = C_INT) :: detect_format
         character(kind = C_CHAR) :: filename(str_medium + 1)
         integer(kind = C_INT), value :: includeFunctionFiles
         integer(kind = C_INT32_T) :: nGrids
         integer(kind = C_INT) :: nDimensions, nScalars,                                     &
              isEndiannessNative, hasIblank, fileType
         type(C_PTR) :: globalGridSizes
       end function detect_format

       subroutine free_memory(globalGridSizes) bind(C, name = "free_memory")
         use, intrinsic :: iso_c_binding
         type(C_PTR) :: globalGridSizes
       end subroutine free_memory

    end interface

    ! Local variables
    integer :: i, j, ierror
    integer(C_INT) :: nDimensions, nScalars, errorCode, includeFunctionFiles_,               &
         isEndiannessNative, hasIblank, file_type
    character(C_CHAR) :: filename_(str_medium + 1)
    type(C_PTR) :: globalGridSizes_ptr
    integer(kind = C_INT32_T) :: nGrids
    integer(C_INT), pointer :: globalGridSizes_(:)

    ! By default, exclude function files
    includeFunctionFiles_ = 0
    if (present(includeFunctionFiles)) then
       if (includeFunctionFiles) includeFunctionFiles_ = 1
    end if

    ! Only the ``master'' process in `comm` calls the function that implements the detection
    ! algorithm and broadcasts an error code to other processes
    if (irank .eq. iroot) then
       do i = 1, len_trim(filename)
          filename_(i) = filename(i:i)
       end do
       filename_(len_trim(filename)+1) = char(0)
       errorCode = detect_format(filename_, includeFunctionFiles_, nGrids, nDimensions,      &
            nScalars, isEndiannessNative, hasIblank, file_type, globalGridSizes_ptr)
       if (.not. present(globalGridSizes)) call free_memory(globalGridSizes_ptr)
    end if
    call MPI_Bcast(errorCode, 1, MPI_INTEGER, 0, comm, ierror)

    ! If an error occured, update the static `plot3dErrorMessage` variable and return
    if (errorCode /= 0) then
       select case (errorCode)
       case (-1)
          write(plot3dErrorMessage, '(2A)') trim(filename),                                  &
               ": File not found or permission denied."
       case (-2)
          write(plot3dErrorMessage, '(2A)') trim(filename), ": Unexpected end of file."
       case (-3)
          write(plot3dErrorMessage, '(2A)') trim(filename),                                  &
               ": Not a valid multi-block whole-format PLOT3D file."
       case (-4)
          write(plot3dErrorMessage, '(2A)') trim(filename), ": Inconsistent record markers."
       end select
       success = .false.
       return
    end if

    ! If not, broadcast the number of dimensions and number of grids to all processes
    call MPI_Bcast(nDimensions, 1, MPI_INTEGER, 0, comm, ierror)
    call MPI_Bcast(nGrids, 1, MPI_INTEGER, 0, comm, ierror)

    if (present(descriptor)) then

       ! Fill the descriptor array with useful information
       descriptor%nDimensions = nDimensions
       descriptor%nGrids = int(nGrids, kind = kind(descriptor%nGrids))
       call MPI_Bcast(isEndiannessNative, 1, MPI_INTEGER, 0, comm, ierror)
       descriptor%isEndiannessNative = (isEndiannessNative /= 0)
       call MPI_Bcast(hasIblank, 1, MPI_INTEGER, 0, comm, ierror)
       descriptor%hasIblank = (hasIblank /= 0)
       call MPI_Bcast(file_type, 1, MPI_INTEGER, 0, comm, ierror)

       ! Character identifier for the file type
       select case (file_type)
       case (0)
          descriptor%fileType = PLOT3D_GRID_FILE
       case (1)
          descriptor%fileType = PLOT3D_SOLUTION_FILE
       case (2)
          descriptor%fileType = PLOT3D_FUNCTION_FILE
       case default
       end select

       ! If this is a function file, broadcast the number of components to all processes
       if (descriptor%fileType .eq. PLOT3D_FUNCTION_FILE) then
          call MPI_Bcast(nScalars, 1, MPI_INTEGER, 0, comm, ierror)
          descriptor%nScalars = nScalars
       else
          descriptor%nScalars = 0
       end if

    end if

    ! If requested, read and broadcast the grid dimensions:

    if (present(globalGridSizes)) then

       if (allocated(globalGridSizes)) deallocate(globalGridSizes)
       allocate(globalGridSizes(nDimensions, nGrids))

       if (irank .eq. iroot) then
          call C_F_POINTER(globalGridSizes_ptr, globalGridSizes_, [nDimensions * nGrids])
          do j = 1, size(globalGridSizes, 2)
             do i = 1, size(globalGridSizes, 1)
                globalGridSizes(i,j) = globalGridSizes_(i + nDimensions * (j - 1))
             end do
          end do
          call free_memory(globalGridSizes_ptr)
       end if

       call MPI_Bcast(globalGridSizes, size(globalGridSizes), MPI_INTEGER, 0, comm, ierror)

    end if

    ! Successful return code
    success = .true.

    return
  end subroutine plot3d_detect_format

  function plot3d_get_offset_from_grid_sizes_(fileType, globalGridSizes, gridIndex,          &
       hasIblank, success, nScalars) result(offset)

    ! Computes the offset, in bytes, to the beginning of the data corresponding to block
    ! `gridIndex` in a 3D multi-block whole-format PLOT3D file. The offset counts past the
    ! leading record size and points directly to the beginning of the actual data

    ! External modules
    use precision
    use parallel

    implicit none

    ! Arguments
    integer, intent(in) :: fileType, globalGridSizes(:,:), gridIndex
    logical, intent(in) :: hasIblank
    logical, intent(out) :: success
    integer, intent(in), optional :: nScalars

    ! Local variables
    integer, parameter :: SIZEOF_SCALAR = WP
    integer, parameter :: SIZEOF_PLOT3D_OFF = WP

    ! Result
    integer(kind = MPI_OFFSET_KIND) :: offset

    ! Local variables
    integer :: nScalars_

    offset = 0

    nScalars_ = 1
    if (present(nScalars)) nScalars_ = nScalars

    select case (fileType)

    case (PLOT3D_GRID_FILE)
       if (hasIblank) then
          offset = 4 + 5 * SIZEOF_PLOT3D_OFF +                                               &
               12 * int(size(globalGridSizes, 2), MPI_OFFSET_KIND) +                         &
               2 * SIZEOF_PLOT3D_OFF * (int(gridIndex, MPI_OFFSET_KIND) - 1) +               &
               (3 * SIZEOF_SCALAR + 4) *                                                     &
               sum(product(int(globalGridSizes(:,:gridIndex-1), MPI_OFFSET_KIND), dim = 1))
       else
          offset = 4 + 5 * SIZEOF_PLOT3D_OFF +                                               &
               12 * int(size(globalGridSizes, 2), MPI_OFFSET_KIND) +                         &
               2 * SIZEOF_PLOT3D_OFF * (int(gridIndex, MPI_OFFSET_KIND) - 1) +               &
               3 * SIZEOF_SCALAR *                                                           &
               sum(product(int(globalGridSizes(:,:gridIndex-1), MPI_OFFSET_KIND), dim = 1))
       end if

    case (PLOT3D_SOLUTION_FILE)
       offset = 4 + 5 * SIZEOF_PLOT3D_OFF +                                                  &
            12 * int(size(globalGridSizes, 2), MPI_OFFSET_KIND) +                            &
            4 * SIZEOF_PLOT3D_OFF * (int(gridIndex, MPI_OFFSET_KIND) - 1) +                  &
            (5 * SIZEOF_SCALAR) *                                                            &
            sum(product(int(globalGridSizes(:,:gridIndex-1), MPI_OFFSET_KIND), dim = 1)) +   &
            4 * SIZEOF_SCALAR * (int(gridIndex, MPI_OFFSET_KIND) - 1)

    case (PLOT3D_FUNCTION_FILE)
       offset = 4 + 5 * SIZEOF_PLOT3D_OFF +                                                  &
            16 * int(size(globalGridSizes, 2), MPI_OFFSET_KIND) +                            &
            2 * SIZEOF_PLOT3D_OFF * (int(gridIndex, MPI_OFFSET_KIND) - 1) +                  &
            (int(nScalars_, MPI_OFFSET_KIND) * SIZEOF_SCALAR) *                              &
            sum(product(int(globalGridSizes(:,:gridIndex-1), MPI_OFFSET_KIND), dim = 1))

    case default

    end select

    ! Successful return code
    success = .true.

  end function plot3d_get_offset_from_grid_sizes_

  function plot3d_get_offset_from_file_(filename, gridIndex, success) result(offset)

    ! Computes the offset, in bytes, to the beginning of the data corresponding to block
    ! `gridIndex` in a 3D multi-block whole-format PLOT3D file. The offset counts past the
    ! leading record size and points directly to the beginning of the actual data

    ! External modules
    use parallel

    ! Arguments
    character(len = *), intent(in) :: filename
    integer, intent(in) :: gridIndex
    logical, intent(out) :: success

    ! Result
    integer(kind = MPI_OFFSET_KIND) :: offset

    ! Local variables
    integer, allocatable :: globalGridSizes(:,:)
    type(t_PLOT3DDescriptor) :: descriptor

    offset = 0

    call plot3d_detect_format(filename, success, descriptor, globalGridSizes,                &
         includeFunctionFiles = .true.)
    if (.not. success) return

    offset = plot3d_get_offset(descriptor%fileType, globalGridSizes, gridIndex,              &
         descriptor%hasIblank, success, descriptor%nScalars)
    if (allocated(globalGridSizes)) deallocate(globalGridSizes)

    ! Successful return code
    success = .true.

  end function plot3d_get_offset_from_file_

  subroutine plot3d_write_skeleton(filename, fileType, globalGridSizes, success, nScalars)

    ! Writes the skeleton of a PLOT3D file to `filename`, including the header, and the
    ! leading and trailing record sizes. This subroutine must be called collectively from all
    ! processes in the MPI communicator `comm`

    ! External modules
    use string
    use parallel
    use, intrinsic :: iso_c_binding

    implicit none

    ! Arguments
    character(len = *), intent(in) :: filename
    integer, intent(in) :: fileType, globalGridSizes(:,:)
    logical, intent(out) :: success
    integer, intent(in), optional :: nScalars

    interface
       function write_skeleton(filename, nGrids, nComp, fileType, globalGridSizes)           &
            bind(C, name = "write_skeleton")
         use, intrinsic :: iso_c_binding
         use string
         integer(kind = C_INT) :: write_skeleton
         character(kind = C_CHAR) :: filename(str_medium + 1)
         integer(kind = C_INT32_T), value :: nGrids, nComp
         integer(kind = C_INT), value :: fileType
         integer(kind = C_INT32_T) :: globalGridSizes(*)
       end function write_skeleton
    end interface

    ! Local variables
    integer :: i, idim, iGrid, ierror
    integer(C_INT) :: errorCode, fileType_
    character(C_CHAR) :: filename_(str_medium + 1)
    integer(C_INT32_T) :: nGrids, nScalars_
    integer(C_INT32_T), pointer :: globalGridSizes_(:)

    select case (fileType)
    case (PLOT3D_GRID_FILE)
       fileType_ = 0
    case (PLOT3D_SOLUTION_FILE)
       fileType_ = 1
    case (PLOT3D_FUNCTION_FILE)
       fileType_ = 2
    case default
    end select

    ! Only the ``master'' process in `comm` calls the function that writes the skeleton of the
    ! PLOT3D file and broadcasts an error code to other processes
    if (irank .eq. iroot) then

       nGrids = size(globalGridSizes, 2)
       nScalars_ = 0
       if (fileType .eq. PLOT3D_FUNCTION_FILE) nScalars_ = nScalars
       do i = 1, len_trim(filename)
          filename_(i) = filename(i:i)
       end do
       filename_(len_trim(filename)+1) = char(0)

       allocate(globalGridSizes_(3 * nGrids))
       globalGridSizes_ = 1

       do iGrid = 1, size(globalGridSizes, 2)
          do idim = 1, size(globalGridSizes, 1)
             globalGridSizes_(idim + 3 * (iGrid - 1)) = globalGridSizes(idim, iGrid)
          end do
       end do

       errorCode = write_skeleton(filename_, nGrids, nScalars_, fileType_, globalGridSizes_)
       deallocate(globalGridSizes_)

    end if
    call MPI_Bcast(errorCode, 1, MPI_INTEGER, 0, comm, ierror)

    ! If an error occured, print an error message and abort
    if (errorCode /= 0) then
       select case (errorCode)
       case (-1)
          write(plot3dErrorMessage, '(2A)') trim(filename),                                  &
               ": Could not open file for writing."
       case default
          write(plot3dErrorMessage, '(2A)') trim(filename), ": Failed to write to file."
       end select
       success = .false.
       return
    end if

    call MPI_Barrier(comm, ierror)

    ! Successful return code
    success = .true.

    return
  end subroutine plot3d_write_skeleton

  subroutine plot3d_write_single_grid(filename, offset, globalGridSize, coordinates,         &
       iblank, success)

    ! Writes the coordinates and IBLANK values corresponding to a single block at offset
    ! `offset` from the beginning of file `filename`

    ! External modules
    use precision
    use parallel

    implicit none

    ! Arguments
    character(len = *), intent(in) :: filename
    integer(kind = MPI_OFFSET_KIND), intent(inout) :: offset
    integer, intent(in) :: globalGridSize(3)
    real(WP), intent(in) :: coordinates(:,:)
    integer, intent(in) :: iblank(:)
    logical, intent(out) :: success

    ! Local variables
    integer :: i, mpiFileHandle, ierror
    integer, parameter :: SIZEOF_SCALAR = WP
    integer, parameter :: SIZEOF_PLOT3D_OFF = WP

    call MPI_File_open(comm, trim(mpiiofs)//trim(filename)//char(0), MPI_MODE_WRONLY,       &
         mpiInfo, mpiFileHandle, ierror)

    do i = 1, size(coordinates, 2)
       call MPI_File_set_view(mpiFileHandle, offset, MPI_REAL_WP,                            &
            mpiDerivedTypeRealSubarray, "native", mpiInfo, ierror)
       call MPI_File_write_all(mpiFileHandle, coordinates(:,i), size(coordinates, 1),        &
            MPI_REAL_WP, MPI_STATUS_IGNORE, ierror)
       offset = offset + SIZEOF_SCALAR * product(int(globalGridSize, MPI_OFFSET_KIND))
    end do

    do while (i .le. 3)
       offset = offset + SIZEOF_SCALAR * product(int(globalGridSize, MPI_OFFSET_KIND))
       i = i + 1
    end do

    call MPI_File_set_view(mpiFileHandle, offset, MPI_INTEGER,                               &
         mpiDerivedTypeIntegerSubarray, "native", mpiInfo, ierror)
    call MPI_File_write_all(mpiFileHandle, iblank, size(iblank),                             &
         MPI_INTEGER, MPI_STATUS_IGNORE, ierror)
    offset = offset + 4 * product(int(globalGridSize, MPI_OFFSET_KIND)) +                    &
         2 * SIZEOF_PLOT3D_OFF

    call MPI_File_close(mpiFileHandle, ierror)

    ! Successful return code
    success = .true.

    return
  end subroutine plot3d_write_single_grid

  subroutine plot3d_write_single_solution_header(filename, offset, solutionHeader, success)

    ! Writes the auxiliary solution data (consisting of 4 real(WP) values)
    ! `solutionHeader` corresponding to a single block at offset `offset` from the
    ! beginning of a file `filename`

    ! External modules
    use precision
    use parallel

    implicit none

    ! Arguments
    character(len = *), intent(in) :: filename
    integer(kind = MPI_OFFSET_KIND), intent(inout) :: offset
    real(WP), intent(in) :: solutionHeader(4)
    logical, intent(out) :: success

    ! Local variables
    integer :: mpiFileHandle, ierror
    integer, parameter :: SIZEOF_SCALAR = WP
    integer, parameter :: SIZEOF_PLOT3D_OFF = WP

    call MPI_File_open(comm, trim(mpiiofs)//trim(filename)//char(0), MPI_MODE_WRONLY,        &
         mpiInfo, mpiFileHandle, ierror)

    call MPI_File_set_view(mpiFileHandle, offset, MPI_REAL_WP,                               &
         MPI_REAL_WP, "native", mpiInfo, ierror)
    call MPI_File_write_all(mpiFileHandle, solutionHeader, 4,                                &
         MPI_REAL_WP, MPI_STATUS_IGNORE, ierror)
    offset = offset + 4 * SIZEOF_SCALAR + 2 * SIZEOF_PLOT3D_OFF

    call MPI_File_close(mpiFileHandle, ierror)

    ! Successful return code
    success = .true.

    return
  end subroutine plot3d_write_single_solution_header

  subroutine plot3d_write_single_solution(filename, offset, globalGridSize, solutionVector,  &
       success)

    ! Writes the solution `solutionVector` corresponding to a single block at offset `offset`
    ! from the beginning of a file `filename`

    ! External modules
    use precision
    use parallel

    implicit none

    ! Arguments
    character(len = *), intent(in) :: filename
    integer(kind = MPI_OFFSET_KIND), intent(inout) :: offset
    integer, intent(in) :: globalGridSize(3)
    real(WP), intent(in) :: solutionVector(:,:)
    logical, intent(out) :: success

    ! Local variables
    integer :: i, mpiFileHandle, ierror
    integer, parameter :: SIZEOF_SCALAR = WP
    integer, parameter :: SIZEOF_PLOT3D_OFF = WP

    call MPI_File_open(comm, trim(mpiiofs)//trim(filename)//char(0), MPI_MODE_WRONLY,        &
         mpiInfo, mpiFileHandle, ierror)

    do i = 1, size(solutionVector, 2) - 1
       call MPI_File_set_view(mpiFileHandle, offset, MPI_REAL_WP,                            &
            mpiDerivedTypeRealSubarray, "native", mpiInfo, ierror)
       call MPI_File_write_all(mpiFileHandle, solutionVector(:,i), size(solutionVector, 1),  &
            MPI_REAL_WP, MPI_STATUS_IGNORE, ierror)
       offset = offset + SIZEOF_SCALAR * product(int(globalGridSize, MPI_OFFSET_KIND))
    end do

    do while (i .le. 4)
       offset = offset + SIZEOF_SCALAR * product(int(globalGridSize, MPI_OFFSET_KIND))
       i = i + 1
    end do

    i = size(solutionVector, 2)
    call MPI_File_set_view(mpiFileHandle, offset, MPI_REAL_WP,                               &
         mpiDerivedTypeRealSubarray, "native", mpiInfo, ierror)
    call MPI_File_write_all(mpiFileHandle, solutionVector(:,i), size(solutionVector, 1),     &
         MPI_REAL_WP, MPI_STATUS_IGNORE, ierror)
    offset = offset + SIZEOF_SCALAR *                                                        &
         product(int(globalGridSize, MPI_OFFSET_KIND)) + 2 * SIZEOF_PLOT3D_OFF

    call MPI_File_close(mpiFileHandle, ierror)

    ! Successful return code
    success = .true.

    return
  end subroutine plot3d_write_single_solution

  subroutine plot3d_write_single_function(filename, offset, globalGridSize, functionVector,  &
       success)

    ! Writes the function `functionVector` corresponding to a single block at offset `offset`
    ! from the beginning of a file `filename`

    ! External modules
    use precision
    use parallel

    implicit none

    ! Arguments
    character(len = *), intent(in) :: filename
    integer(kind = MPI_OFFSET_KIND), intent(inout) :: offset
    integer, intent(in) :: globalGridSize(3)
    real(WP), intent(in) :: functionVector(:,:)
    logical, intent(out) :: success

    ! Local variables
    integer :: i, mpiFileHandle, ierror
    integer, parameter :: SIZEOF_SCALAR = WP
    integer, parameter :: SIZEOF_PLOT3D_OFF = WP

    call MPI_File_open(comm, trim(mpiiofs)//trim(filename)//char(0), MPI_MODE_WRONLY,        &
         mpiInfo, mpiFileHandle, ierror)

    do i = 1, size(functionVector, 2)
       call MPI_File_set_view(mpiFileHandle, offset, MPI_REAL_WP,                            &
            mpiDerivedTypeRealSubarray, "native", mpiInfo, ierror)
       call MPI_File_write_all(mpiFileHandle, functionVector(:,i), size(functionVector, 1),  &
            MPI_REAL_WP, MPI_STATUS_IGNORE, ierror)
       offset = offset + SIZEOF_SCALAR * product(int(globalGridSize, MPI_OFFSET_KIND))
    end do
    offset = offset + 2 * SIZEOF_PLOT3D_OFF

    call MPI_File_close(mpiFileHandle, ierror)

    ! Successful return code
    success = .true.

    return
  end subroutine plot3d_write_single_function

  subroutine plot3d_read_single_grid(filename, offset, globalGridSize, coordinates, iblank,  &
       success)

    ! Reads the coordinates and IBLANK values corresponding to a single block at offset
    ! `offset` from the beginning of a file `filename`. The offset may be obtained by calling
    ! `plot3dGetOffset`. If IBLANK values are not present in `filename`, `iblank` is set to 1

    ! External modules
    use precision
    use parallel

    implicit none

    ! Arguments
    character(len = *), intent(in) :: filename
    integer(kind = MPI_OFFSET_KIND), intent(inout) :: offset
    integer, intent(in) :: globalGridSize(3)
    real(WP), intent(out) :: coordinates(:,:)
    integer, intent(out) :: iblank(:)
    logical, intent(out) :: success

    ! Local variables
    integer :: i, idim, mpiFileHandle, ierror
    integer, parameter :: SIZEOF_SCALAR = WP
    integer, parameter :: SIZEOF_PLOT3D_OFF = WP
    type(t_PLOT3DDescriptor) :: descriptor

    call plot3d_detect_format(filename, success, descriptor)
    if (.not. success) return

    if (descriptor%fileType /= PLOT3D_GRID_FILE) then
       write(plot3dErrorMessage, '(3A)') "'", trim(filename),                                &
            "': Not a valid PLOT3D grid file."
       success = .false.
       return
    end if

    call MPI_File_open(comm, trim(mpiiofs)//trim(filename)//char(0), MPI_MODE_RDONLY,        &
         mpiInfo, mpiFileHandle, ierror)

    do idim = 1, min(size(coordinates, 2), 3)
       call MPI_File_set_view(mpiFileHandle, offset, MPI_REAL_WP,                            &
            mpiDerivedTypeRealSubarray, "native", mpiInfo, ierror)
       call MPI_File_read_all(mpiFileHandle, coordinates(:,idim), size(coordinates, 1),      &
            MPI_REAL_WP, MPI_STATUS_IGNORE, ierror)
       if (.not. descriptor%isEndiannessNative) then
          do i = 1, size(coordinates, 1)
             call swap_endianness(coordinates(i, idim))
          end do
       end if
       offset = offset + SIZEOF_SCALAR * product(int(globalGridSize, MPI_OFFSET_KIND))
    end do

    do while (idim .le. 3)
       offset = offset + SIZEOF_SCALAR * product(int(globalGridSize, MPI_OFFSET_KIND))
       idim = idim + 1
    end do

    if (descriptor%hasIblank) then
       call MPI_File_set_view(mpiFileHandle, offset, MPI_INTEGER,                            &
            mpiDerivedTypeIntegerSubarray, "native", mpiInfo, ierror)
       call MPI_File_read_all(mpiFileHandle, iblank, size(iblank),                           &
            MPI_INTEGER, MPI_STATUS_IGNORE, ierror)
       if (.not. descriptor%isEndiannessNative) then
          do i = 1, size(iblank)
             call swap_endianness(iblank(i))
          end do
       end if
       offset = offset + product(int(globalGridSize, MPI_OFFSET_KIND)) + 2 * SIZEOF_PLOT3D_OFF
    else
       iblank = 1
    end if

    call MPI_File_close(mpiFileHandle, ierror)

    ! Successful return code
    success = .true.

    return
  end subroutine plot3d_read_single_grid

  subroutine plot3d_read_single_solution_header(filename, offset, solutionHeader, success)

    ! Reads the auxiliary solution data (consisting of 4 real(WP) values)
    ! `solutionHeader` corresponding to a single block at offset `offset` from the
    ! beginning of a file `filename`. The offset may be obtained by calling
    ! `plot3dGetOffset`

    ! External modules
    use precision
    use parallel

    implicit none

    ! Arguments
    character(len = *), intent(in) :: filename
    integer(kind = MPI_OFFSET_KIND), intent(inout) :: offset
    real(WP), intent(out) :: solutionHeader(4)
    logical, intent(out) :: success

    ! Local variables
    integer :: i, mpiFileHandle, ierror
    integer, parameter :: SIZEOF_SCALAR = WP
    integer, parameter :: SIZEOF_PLOT3D_OFF = WP
    type(t_PLOT3DDescriptor) :: descriptor

    call plot3d_detect_format(filename, success, descriptor)
    if (.not. success) return

    if (descriptor%fileType .ne. PLOT3D_SOLUTION_FILE) then
       write(plot3dErrorMessage, '(3A)') "'", trim(filename),                                &
            "': Not a valid PLOT3D solution file."
       success = .false.
       return
    end if

    call MPI_File_open(comm, trim(mpiiofs)//trim(filename)//char(0), MPI_MODE_RDONLY,        &
         mpiInfo, mpiFileHandle, ierror)

    call MPI_File_set_view(mpiFileHandle, offset, MPI_REAL_WP,                               &
         MPI_REAL_WP, "native", mpiInfo, ierror)
    call MPI_File_read_all(mpiFileHandle, solutionHeader, 4,                                 &
         MPI_REAL_WP, MPI_STATUS_IGNORE, ierror)
    if (.not. descriptor%isEndiannessNative) then
       do i = 1, 4
          call swap_endianness(solutionHeader(i))
       end do
    end if
    offset = offset + 4 * SIZEOF_SCALAR + 2 * SIZEOF_PLOT3D_OFF

    call MPI_File_close(mpiFileHandle, ierror)

    ! Successful return code
    success = .true.

    return
  end subroutine plot3d_read_single_solution_header

  subroutine plot3d_read_single_solution(filename, offset, globalGridSize, solutionVector,   &
       success)

    ! Reads the solution `solutionVector` corresponding to a single block at offset `offset`
    ! from the beginning of a file `filename`. The offset may be obtained by calling
    ! `PLOT3D_get_offset`

    ! External modules
    use precision
    use parallel

    implicit none

    ! Arguments
    character(len = *), intent(in) :: filename
    integer(kind = MPI_OFFSET_KIND), intent(inout) :: offset
    integer, intent(in) :: globalGridSize(3)
    real(WP), intent(out) :: solutionVector(:,:)
    logical, intent(out) :: success

    ! Local variables
    integer :: i, j, mpiFileHandle, ierror
    integer, parameter :: SIZEOF_SCALAR = WP
    integer, parameter :: SIZEOF_PLOT3D_OFF = WP
    type(t_PLOT3DDescriptor) :: descriptor

    call plot3d_detect_format(filename, success, descriptor)
    if (.not. success) return

    if (descriptor%fileType /= PLOT3D_SOLUTION_FILE) then
       write(plot3dErrorMessage, '(3A)') "'", trim(filename),                                &
            "': Not a valid PLOT3D solution file."
       success = .false.
       return
    end if

    call MPI_File_open(comm, trim(mpiiofs)//trim(filename)//char(0), MPI_MODE_RDONLY,        &
         mpiInfo, mpiFileHandle, ierror)

    do i = 1, min(size(solutionVector, 2) - 1, 4)
       call MPI_File_set_view(mpiFileHandle, offset, MPI_REAL_WP,                            &
            mpiDerivedTypeRealSubarray, "native", mpiInfo, ierror)
       call MPI_File_read_all(mpiFileHandle, solutionVector(:,i), size(solutionVector, 1),   &
            MPI_REAL_WP, MPI_STATUS_IGNORE, ierror)
       if (.not. descriptor%isEndiannessNative) then
          do j = 1, size(solutionVector, 1)
             call swap_endianness(solutionVector(j,i))
          end do
       end if
       offset = offset + SIZEOF_SCALAR * product(int(globalGridSize, MPI_OFFSET_KIND))
    end do

    do while (i .le. 4)
       offset = offset + SIZEOF_SCALAR * product(int(globalGridSize, MPI_OFFSET_KIND))
       i = i + 1
    end do

    i = size(solutionVector, 2)
    call MPI_File_set_view(mpiFileHandle, offset, MPI_REAL_WP,                               &
         mpiDerivedTypeRealSubarray, "native", mpiInfo, ierror)
    call MPI_File_read_all(mpiFileHandle, solutionVector(:,i), size(solutionVector, 1),      &
         MPI_REAL_WP, MPI_STATUS_IGNORE, ierror)
    if (.not. descriptor%isEndiannessNative) then
       do j = 1, size(solutionVector, 1)
          call swap_endianness(solutionVector(j,i))
       end do
    end if
    offset = offset + SIZEOF_SCALAR * product(int(globalGridSize, MPI_OFFSET_KIND)) +        &
         2 * SIZEOF_PLOT3D_OFF

    call MPI_File_close(mpiFileHandle, ierror)

    ! Successful return code
    success = .true.

    return
  end subroutine plot3d_read_single_solution

  subroutine plot3d_read_single_function(filename, offset, globalGridSize, functionVector,   &
       success)

    ! Reads the function `functionVector` corresponding to a single block at offset `offset`
    ! from the beginning of a file `filename`. The offset may be obtained by calling
    ! `PLOT3D_get_offset`

    ! External modules
    use precision
    use parallel

    implicit none

    ! Arguments
    character(len = *), intent(in) :: filename
    integer(kind = MPI_OFFSET_KIND), intent(inout) :: offset
    integer, intent(in) :: globalGridSize(3)
    real(WP), intent(out) :: functionVector(:,:)
    logical, intent(out) :: success

    ! Local variables
    integer :: i, j, mpiFileHandle, ierror
    integer, parameter :: SIZEOF_SCALAR = WP
    integer, parameter :: SIZEOF_PLOT3D_OFF = WP
    type(t_PLOT3DDescriptor) :: descriptor

    call plot3d_detect_format(filename, success, descriptor,                                 &
         includeFunctionFiles = .true.)
    if (.not. success) return

    if (descriptor%fileType /= PLOT3D_FUNCTION_FILE) then
       write(plot3dErrorMessage, '(3A)') "'", trim(filename),                                &
            "': Not a valid PLOT3D function file."
       success = .true.
       return
    end if

    call MPI_File_open(comm, trim(mpiiofs)//trim(filename)//char(0), MPI_MODE_RDONLY,        &
         mpiInfo, mpiFileHandle, ierror)

    do i = 1, min(descriptor%nScalars, size(functionVector, 2))
       call MPI_File_set_view(mpiFileHandle, offset, MPI_REAL_WP,                            &
            mpiDerivedTypeRealSubarray, "native", mpiInfo, ierror)
       call MPI_File_read_all(mpiFileHandle, functionVector(:,i), size(functionVector, 1),   &
            MPI_REAL_WP, MPI_STATUS_IGNORE, ierror)
       if (.not. descriptor%isEndiannessNative) then
          do j = 1, size(functionVector, 1)
             call swap_endianness(functionVector(j,i))
          end do
       end if
       offset = offset + SIZEOF_SCALAR * product(int(globalGridSize, MPI_OFFSET_KIND))
    end do
    offset = offset + 2 * SIZEOF_PLOT3D_OFF

    call MPI_File_close(mpiFileHandle, ierror)

    ! Successful return code
    success = .true.

    return
  end subroutine plot3d_read_single_function

  pure subroutine swap_integer_endianness_(a)

    ! Arguments
    integer, intent(inout) :: a

    ! Local variables
    integer :: b

    call mvbits(a,  0, 8, b, 24)
    call mvbits(a,  8, 8, b, 16)
    call mvbits(a, 16, 8, b,  8)
    call mvbits(a, 24, 8, b,  0)

    a = b

  end subroutine swap_integer_endianness_

  subroutine swap_scalar_endianness_(a)

    ! External modules
    use precision

    ! Arguments
    real(WP), intent(inout) :: a

    ! Local variables
    integer :: i
    integer, parameter :: SIZEOF_SCALAR = WP
    integer(kind = 1) :: b(SIZEOF_SCALAR), c(SIZEOF_SCALAR)

    b = transfer(a, b)
    do i = 1, SIZEOF_SCALAR
       c(i) = b(SIZEOF_SCALAR+1-i)
    end do
    a = transfer(c, a)

  end subroutine swap_scalar_endianness_

end module PLOT3D
