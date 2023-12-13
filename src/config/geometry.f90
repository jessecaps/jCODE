module geometry

  ! External modules
  use precision

  implicit none

  ! Periodicity options
  integer, parameter ::                                                                      &
       NONE    = 0,                                                                          &
       PLANE   = 1,                                                                          &
       OVERLAP = 2,                                                                          &
       POLAR   = 3

  ! Grid-specific parameters
  integer :: nDimensions, globalGridSize(3), localGridSize(3), gridOffset(3),                &
       nGridPoints = 0, periodicityType(3), iStart(3), iEnd(3)
  real(WP) :: periodicLength(3), minGridSpacing, globalGridVolume, domainExtent(3,2)
  logical :: isPeriodic(3), isDomainCurvilinear, useIblank

contains

  ! Get the number of dimensions from the grid file
  ! -----------------------------------------------
  subroutine get_dimensions(reportFilename)

    ! External modules
    use string
    use parser
    use fileio
    use parallel
    use simulation_flags

    implicit none

    ! Arguments
    character(len = *), intent(out), optional :: reportFilename

    ! Local variables
    integer :: i, ifile, ierror
    integer, dimension(3) :: dims
    character(len = str_medium) :: filename, filename_
    logical :: fileIsThere

    ! Get the grid file name
    call parser_read('grid file', filename_)

    ! Get the name of the header file in case serial i/o is used
    fileIsThere = .false.
    if (useSerialIO) then
       filename = trim(filename_) // '/grid.header'
       inquire(file = filename, exist = fileIsThere)
    end if
    if (.not. fileIsThere) filename = trim(filename_)

    ! Root reads the header file
    if (irank.eq.iroot) then

       ! Open the file
       call BINARY_FILE_OPEN(ifile, trim(filename), "r", ierror)

       ! Check number of processors are correct if serial file is used
       if (fileIsThere) then
          call BINARY_FILE_READ(ifile, dims(1), 1, kind(dims), ierror)
          if (dims(1) .ne. nProcs) then
             print*, 'Expected ', dims(1), ' processor(s)'
             print*, 'Currently ', nProcs, ' processor(s)'
             call die('Number of processors incompatible with serial grid file')
          end if
       end if

       ! Get the dimensions from header
       call BINARY_FILE_READ(ifile, dims, 3, kind(dims), ierror)
       globalGridSize = dims
       nDimensions = 0
       do i = 1, 3
          if (globalGridSize(i) .gt. 1) nDimensions = nDimensions + 1
       end do

       ! Is the domain curvilinear?
       call BINARY_FILE_READ(ifile, dims(1), 1, kind(dims), ierror)
       isDomainCurvilinear = (dims(1).ne.0)

       ! Use iblank?
       call BINARY_FILE_READ(ifile, dims(1), 1, kind(dims), ierror)
       useIblank = (dims(1) .ne. 0)

       ! Get the periodicity from header
       call BINARY_FILE_READ(ifile, dims, 3, kind(dims), ierror)
       periodicityType = NONE
       periodicityType(1:nDimensions) = dims(1:nDimensions)
       do i = 1, nDimensions
          if (periodicityType(i).eq.PLANE) call BINARY_FILE_READ(ifile, periodicLength(i),&
               1, kind(periodicLength(i)), ierror)
       end do

       ! Close the file
       call BINARY_FILE_CLOSE(ifile, ierror)
    end if

    ! Broadcast information
    call parallel_bc(nDimensions)
    call parallel_bc(globalGridSize)
    call parallel_bc(isDomainCurvilinear)
    call parallel_bc(periodicityType)
    call parallel_bc(periodicLength)
    call parallel_bc(useIblank)

    ! Set periodicity logic
    isPeriodic = (periodicityType .ne. NONE)

    ! Report the filename if present
    if (present(reportFilename)) reportFilename = trim(filename_)

    return
  end subroutine get_dimensions

end module geometry


! ================================== !
! Initialize the geometry parameters !
! ================================== !
subroutine geometry_setup

  ! Internal modules
  use geometry

  ! External modules
  use string
  use parser
  use parallel
  use math, only : pigeon_hole
  use simulation_flags

  implicit none

  ! Local variables
  integer :: i

  ! Initialize parallel i/o info and topology
  call parallel_init_io
  call parallel_init_topology(nDimensions, isPeriodic)

  ! Find the size of the part of the grid distributed to the current process
  gridOffset = 0
  localGridSize = 1
  do i = 1, nDimensions
     call pigeon_hole(globalGridSize(i), nProcsDir(i), procCoords(i), gridOffset(i),         &
          localGridSize(i))
  end do
  nGridPoints = product(localGridSize)

  ! Store the starting and ending indices
  iStart = gridOffset + 1
  iEnd = gridOffset + localGridSize

  ! Derived types for MPI I/O
  call parallel_subarray(globalGridSize, localGridSize, gridOffset)
  if (useParticles) call prepare_mpi_particle
  if (useIBM) call prepare_mpi_object

  return
end subroutine geometry_setup



