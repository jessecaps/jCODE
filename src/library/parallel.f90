module parallel

  ! External modules
  !use MPI
  use precision
  use string

  implicit none

  include 'mpif.h' ! Not compatible with Fortran versions >= 2003

  integer :: nProcs, nProcsDir(3)
  integer :: iRank, iRoot, iRankDir(3)
  integer :: comm, commDir(3)
  integer :: procCoords(3)
  integer :: MPI_REAL_WP, MPI_REAL_SP, MPI_COMPLEX_WP
  integer :: mpiDerivedTypeRealSubarray = MPI_DATATYPE_NULL,                                 &
       mpiDerivedTypeIntegerSubarray = MPI_DATATYPE_NULL
  integer :: globalComm, globalCommRank, globalCommSize, globalCommRoot
  integer :: externalComm, externalCommRank, externalCommSize, externalCommRoot
  logical :: disableManualDecomp

  ! MPI info for parallel file system
  character(len=str_medium) :: mpiiofs
  integer :: mpiInfo = MPI_INFO_NULL

  ! Compute the global maximum of anything
  interface parallel_max
     module procedure parallel_max_real_3d
     module procedure parallel_max_real_0d
     module procedure parallel_max_int_0d
     module procedure parallel_max_real_0d_in_place
     module procedure parallel_max_int_0d_in_place
  end interface parallel_max

  ! Compute the directional maximum of anything
  interface parallel_max_dir
     module procedure parallel_max_real_1d
     module procedure parallel_max_real_1d_in_place
     module procedure parallel_max_real_2d
  end interface parallel_max_dir

  ! Compute the global minimum of anything
  interface parallel_min
     module procedure parallel_min_real_3d
     module procedure parallel_min_real_0d
     module procedure parallel_min_int_0d
     module procedure parallel_min_real_0d_in_place
     module procedure parallel_min_int_0d_in_place
  end interface parallel_min

  ! Compute the global sum of anything
  interface parallel_sum
     module procedure parallel_sum_int_0d
     module procedure parallel_sum_int_0d_in_place
     module procedure parallel_sum_int_1d
     module procedure parallel_sum_int_1d_in_place
     module procedure parallel_sum_real_0d
     module procedure parallel_sum_real_0d_in_place
     module procedure parallel_sum_real_1d
     module procedure parallel_sum_real_1d_in_place
     module procedure parallel_sum_real_2d
     module procedure parallel_sum_real_2d_in_place
     module procedure parallel_sum_real_3d
  end interface parallel_sum

  ! Perform broadcast
  interface parallel_bc
     module procedure parallel_bc_char
     module procedure parallel_bc_log
     module procedure parallel_bc_int_0d
     module procedure parallel_bc_int_1d
     module procedure parallel_bc_int_2d
     module procedure parallel_bc_int_3d
     module procedure parallel_bc_real_0d
     module procedure parallel_bc_real_1d
     module procedure parallel_bc_real_2d
     module procedure parallel_bc_real_3d
  end interface parallel_bc

  ! Perform directional summation
  interface parallel_sum_dir
     module procedure parallel_sum_dir_real_0d
     module procedure parallel_sum_dir_real_0d_in_place
     module procedure parallel_sum_dir_real_1d
     module procedure parallel_sum_dir_real_1d_in_place
     module procedure parallel_sum_dir_real_2d
     module procedure parallel_sum_dir_real_2d_in_place
     module procedure parallel_sum_dir_real_3d
     module procedure parallel_sum_dir_real_3d_in_place
  end interface parallel_sum_dir

  ! Perform directional gathering
  interface parallel_gather_dir
     module procedure parallel_gather_dir_real_1d
  end interface parallel_gather_dir

  ! Compute the global logical or of anything
  interface parallel_lor
     module procedure parallel_lor_0d
     module procedure parallel_lor_0d_in_place
  end interface parallel_lor

  ! Compute the directional logical or of anything
  interface parallel_lor_dir
     module procedure parallel_lor_dir_0d
  end interface parallel_lor_dir

  ! Fill ghost points
  interface fill_ghost_points
     module procedure fill_ghost_points_int
     module procedure fill_ghost_points_real
  end interface fill_ghost_points

contains

  ! MPI ALLGATHER
  subroutine parallel_gather_dir_real_1d(A, B, direction)

    implicit none

    real(WP), dimension(:), intent(in)  :: A
    real(WP), dimension(:), intent(out) :: B
    integer, intent(in) :: direction
    integer :: ierror, n

    n = size(A)

    call MPI_ALLGATHER(A, n, MPI_REAL_WP, B, n, MPI_REAL_WP, commDir(direction), ierror)

    return
  end subroutine parallel_gather_dir_real_1d

  ! MPI SUM-DIR
  subroutine parallel_sum_dir_real_0d(A, B, direction)

    implicit none

    real(WP), intent(in)  :: A
    real(WP), intent(out) :: B
    integer, intent(in) :: direction
    integer :: ierror, n

    n = 1

    call MPI_ALLREDUCE(A, B, n, MPI_REAL_WP, MPI_SUM, commDir(direction), ierror)

    return
  end subroutine parallel_sum_dir_real_0d

  subroutine parallel_sum_dir_real_0d_in_place(A, direction)

    implicit none

    real(WP), intent(inout)  :: A
    integer, intent(in) :: direction
    integer :: ierror, n

    n = 1

    call MPI_ALLREDUCE(MPI_IN_PLACE, A, n, MPI_REAL_WP, MPI_SUM, commDir(direction), ierror)

    return
  end subroutine parallel_sum_dir_real_0d_in_place

  ! MPI SUM-DIR
  subroutine parallel_sum_dir_real_1d(A, B, direction)

    implicit none

    real(WP), dimension(:), intent(in)  :: A
    real(WP), dimension(:), intent(out) :: B
    integer, intent(in) :: direction
    integer :: ierror, n

    n = size(A)

    call MPI_ALLREDUCE(A, B, n, MPI_REAL_WP, MPI_SUM, commDir(direction), ierror)
    
    return
  end subroutine parallel_sum_dir_real_1d

  subroutine parallel_sum_dir_real_1d_in_place(A, direction)

    implicit none

    real(WP), dimension(:), intent(inout) :: A
    integer, intent(in) :: direction
    integer :: ierror, n
    
    n = size(A)

    call MPI_ALLREDUCE(MPI_IN_PLACE, A, n, MPI_REAL_WP, MPI_SUM, commDir(direction), ierror)
    
    return
  end subroutine parallel_sum_dir_real_1d_in_place
  
  subroutine parallel_sum_dir_real_2d(A, B, direction)

    implicit none

    real(WP), dimension(:,:), intent(in)  :: A
    real(WP), dimension(:,:), intent(out) :: B
    integer, intent(in) :: direction
    integer :: ierror, n

    n = size(A)

    call MPI_ALLREDUCE(A, B, n, MPI_REAL_WP, MPI_SUM, commDir(direction), ierror)

    return
  end subroutine parallel_sum_dir_real_2d

  subroutine parallel_sum_dir_real_2d_in_place(A, direction)

    implicit none

    real(WP), dimension(:,:), intent(inout) :: A
    integer, intent(in) :: direction
    integer :: ierror, n
    
    n = size(A)

    call MPI_ALLREDUCE(MPI_IN_PLACE, A, n, MPI_REAL_WP, MPI_SUM, commDir(direction), ierror)
    
    return
  end subroutine parallel_sum_dir_real_2d_in_place
  
  subroutine parallel_sum_dir_real_3d(A, B, direction)

    implicit none

    real(WP), dimension(:,:,:), intent(in)  :: A
    real(WP), dimension(:,:,:), intent(out) :: B
    integer, intent(in) :: direction
    integer :: ierror, n
    
    n = size(A)

    call MPI_ALLREDUCE(A, B, n, MPI_REAL_WP, MPI_SUM, commDir(direction), ierror)

    return
  end subroutine parallel_sum_dir_real_3d

  subroutine parallel_sum_dir_real_3d_in_place(A, direction)

    implicit none

    real(WP), dimension(:,:,:), intent(inout) :: A
    integer, intent(in) :: direction
    integer :: ierror, n
    
    n = size(A)

    call MPI_ALLREDUCE(MPI_IN_PLACE, A, n, MPI_REAL_WP, MPI_SUM, commDir(direction), ierror)
    
    return
  end subroutine parallel_sum_dir_real_3d_in_place
  
  ! MPI MAX
  subroutine parallel_max_real_3d(A, B)

    implicit none

    real(WP), dimension(:,:,:), intent(in) :: A
    real(WP), intent(out) :: B
    real(WP) :: C
    integer :: ierror

    C = maxval(A)

    call MPI_ALLREDUCE(C, B, 1, MPI_REAL_WP, MPI_MAX, comm, ierror)

    return
  end subroutine parallel_max_real_3d
  
  subroutine parallel_max_real_0d(A, B)

    implicit none

    real(WP), intent(in)  :: A
    real(WP), intent(out) :: B
    integer :: ierror

    call MPI_ALLREDUCE(A, B, 1, MPI_REAL_WP, MPI_MAX, comm, ierror)

    return
  end subroutine parallel_max_real_0d
  
  subroutine parallel_max_int_0d(A, B)

    implicit none

    integer, intent(in)  :: A
    integer, intent(out) :: B
    integer :: ierror

    call MPI_ALLREDUCE(A, B, 1, MPI_INTEGER, MPI_MAX, comm, ierror)

    return
  end subroutine parallel_max_int_0d

subroutine parallel_max_real_0d_in_place(A)

    implicit none

    real(WP), intent(inout) :: A
    integer :: ierror

    call MPI_ALLREDUCE(MPI_IN_PLACE, A, 1, MPI_REAL_WP, MPI_MAX, comm, ierror)

    return
  end subroutine parallel_max_real_0d_in_place
  
  subroutine parallel_max_int_0d_in_place(A)

    implicit none

    integer, intent(inout) :: A
    integer :: ierror

    call MPI_ALLREDUCE(MPI_IN_PLACE, A, 1, MPI_INTEGER, MPI_MAX, comm, ierror)

    return
  end subroutine parallel_max_int_0d_in_place
  
  ! MPI DIRECTIONAL MAX
  subroutine parallel_max_real_1d(A,B,direction)

    implicit none
    real(WP), dimension(:), intent(in)  :: A
    real(WP), dimension(:), intent(out) :: B
    integer, intent(in) :: direction
    integer :: ierror, n

    n = size(A)

    call MPI_ALLREDUCE(A, B, n, MPI_REAL_WP, MPI_MAX, commDir(direction), ierror)
    
    return
  end subroutine parallel_max_real_1d

  subroutine parallel_max_real_1d_in_place(A,direction)

    implicit none
    real(WP), dimension(:), intent(inout)  :: A
    integer, intent(in) :: direction
    integer :: ierror, n

    n = size(A)

    call MPI_ALLREDUCE(MPI_IN_PLACE, A, n, MPI_REAL_WP, MPI_MAX, commDir(direction), ierror)
    
    return
  end subroutine parallel_max_real_1d_in_place
  
  subroutine parallel_max_real_2d(A, B, direction)

    implicit none

    real(WP), dimension(:,:), intent(in)  :: A
    real(WP), dimension(:,:), intent(out) :: B
    integer, intent(in) :: direction
    integer :: ierror, n
    
    n = size(A)

    call MPI_ALLREDUCE(A, B, n, MPI_REAL_WP, MPI_MAX, commDir(direction), ierror)

    return
  end subroutine parallel_max_real_2d
  
  ! MPI MIN
  subroutine parallel_min_real_3d(A, B)

    implicit none

    real(WP), dimension(:,:,:), intent(in) :: A
    real(WP), intent(out) :: B
    real(WP) :: C
    integer :: ierror

    C = minval(A)

    call MPI_ALLREDUCE(C, B, 1, MPI_REAL_WP, MPI_MIN, comm, ierror)

    return
  end subroutine parallel_min_real_3d
  
  subroutine parallel_min_real_0d(A, B)

    implicit none
    real(WP), intent(in)  :: A
    real(WP), intent(out) :: B
    integer :: ierror

    call MPI_ALLREDUCE(A, B, 1, MPI_REAL_WP, MPI_MIN, comm, ierror)

    return
  end subroutine parallel_min_real_0d
  
  subroutine parallel_min_int_0d(A, B)
    implicit none

    integer, intent(in)  :: A
    integer, intent(out) :: B
    integer :: ierror

    call MPI_ALLREDUCE(A, B, 1, MPI_INTEGER, MPI_MIN, comm, ierror)

    return
  end subroutine parallel_min_int_0d

  subroutine parallel_min_real_0d_in_place(A)

    implicit none
    real(WP), intent(inout) :: A
    integer :: ierror

    call MPI_ALLREDUCE(MPI_IN_PLACE, A, 1, MPI_REAL_WP, MPI_MIN, comm, ierror)

    return
  end subroutine parallel_min_real_0d_in_place
  
  subroutine parallel_min_int_0d_in_place(A)
    implicit none

    integer, intent(inout) :: A
    integer :: ierror

    call MPI_ALLREDUCE(MPI_IN_PLACE, A, 1, MPI_INTEGER, MPI_MIN, comm, ierror)

    return
  end subroutine parallel_min_int_0d_in_place
  
  ! MPI SUM
  subroutine parallel_sum_real_0d(A, B)

    implicit none

    real(WP), intent(in)  :: A
    real(WP), intent(out) :: B
    integer :: ierror

    call MPI_ALLREDUCE(A, B, 1, MPI_REAL_WP, MPI_SUM, comm, ierror)

    return
  end subroutine parallel_sum_real_0d
  
  subroutine parallel_sum_int_0d(A,B)

    implicit none

    integer, intent(in)  :: A
    integer, intent(out) :: B
    integer :: ierror

    call MPI_ALLREDUCE(A, B, 1, MPI_INTEGER, MPI_SUM, comm, ierror)

    return
  end subroutine parallel_sum_int_0d

  subroutine parallel_sum_real_0d_in_place(A)

    implicit none

    real(WP), intent(inout) :: A
    integer :: ierror

    call MPI_ALLREDUCE(MPI_IN_PLACE, A, 1, MPI_REAL_WP, MPI_SUM, comm, ierror)

    return
  end subroutine parallel_sum_real_0d_in_place
  
  subroutine parallel_sum_int_0d_in_place(A)

    implicit none

    integer, intent(inout) :: A
    integer :: ierror

    call MPI_ALLREDUCE(MPI_IN_PLACE, A, 1, MPI_INTEGER, MPI_SUM, comm, ierror)

    return
  end subroutine parallel_sum_int_0d_in_place

  subroutine parallel_sum_int_1d(A, B)
    implicit none

    integer, dimension(:), intent(in)  :: A
    integer, dimension(:), intent(out) :: B
    integer :: ierror

    call MPI_ALLREDUCE(A, B, size(A), MPI_INTEGER, MPI_SUM, comm, ierror)

    return
  end subroutine parallel_sum_int_1d

  subroutine parallel_sum_int_1d_in_place(A)
    implicit none

    integer, dimension(:), intent(inout)  :: A
    integer :: ierror

    call MPI_ALLREDUCE(MPI_IN_PLACE, A, size(A), MPI_INTEGER, MPI_SUM, comm, ierror)

    return
  end subroutine parallel_sum_int_1d_in_place
  
  subroutine parallel_sum_real_1d(A, B)

    implicit none

    real(WP), dimension(:), intent(in)  :: A
    real(WP), dimension(:), intent(out) :: B
    integer :: ierror

    call MPI_ALLREDUCE(A, B, size(A), MPI_REAL_WP, MPI_SUM, comm, ierror)

    return
  end subroutine parallel_sum_real_1d

  subroutine parallel_sum_real_1d_in_place(A)

    implicit none

    real(WP), dimension(:), intent(inout)  :: A
    integer :: ierror

    call MPI_ALLREDUCE(MPI_IN_PLACE, A, size(A), MPI_REAL_WP, MPI_SUM, comm, ierror)

    return
  end subroutine parallel_sum_real_1d_in_place
  
  subroutine parallel_sum_real_2d(A, B)

    implicit none

    real(WP), dimension(:,:), intent(in)  :: A
    real(WP), dimension(:,:), intent(out) :: B
    integer :: ierror

    call MPI_ALLREDUCE(A, B, size(A), MPI_REAL_WP, MPI_SUM, comm, ierror)

    return
  end subroutine parallel_sum_real_2d

  subroutine parallel_sum_real_2d_in_place(A)

    implicit none

    real(WP), dimension(:,:), intent(inout)  :: A
    integer :: ierror

    call MPI_ALLREDUCE(MPI_IN_PLACE, A, size(A), MPI_REAL_WP, MPI_SUM, comm, ierror)

    return
  end subroutine parallel_sum_real_2d_in_place
  
  subroutine parallel_sum_real_3d(A, B)

    implicit none

    real(WP), dimension(:,:,:), intent(in)  :: A
    real(WP), dimension(:,:,:), intent(out) :: B
    integer :: ierror

    call MPI_ALLREDUCE(A, B, size(A), MPI_REAL_WP, MPI_SUM, comm, ierror)

    return
  end subroutine parallel_sum_real_3d
  
  ! MPI_BCAST
  subroutine parallel_bc_char(A)

    implicit none

    integer :: ierror

    character(len=*) :: A

    call MPI_BCAST(A, len(A), MPI_CHARACTER, iRoot, comm, ierror)

    return
  end subroutine parallel_bc_char

  subroutine parallel_bc_log(A)

    implicit none

    integer :: ierror
    logical :: A

    call MPI_BCAST(A, 1, MPI_LOGICAL, iRoot, comm, ierror)

    return
  end subroutine parallel_bc_log
  
  subroutine parallel_bc_int_0d(A)

    implicit none

    integer :: ierror
    integer :: A

    call MPI_BCAST(A, 1, MPI_INTEGER, iRoot, comm, ierror)

    return
  end subroutine parallel_bc_int_0d
  
  subroutine parallel_bc_int_1d(A)

    implicit none

    integer :: ierror
    integer, dimension(:) :: A

    call MPI_BCAST(A, size(A), MPI_INTEGER, iRoot, comm, ierror)

    return
  end subroutine parallel_bc_int_1d
  
  subroutine parallel_bc_int_2d(A)

    implicit none

    integer :: ierror
    integer, dimension(:,:) :: A

    call MPI_BCAST(A, size(A), MPI_INTEGER, iRoot, comm, ierror)

    return
  end subroutine parallel_bc_int_2d
  
  subroutine parallel_bc_int_3d(A)

    implicit none

    integer :: ierror
    integer, dimension(:,:,:) :: A

    call MPI_BCAST(A, size(A), MPI_INTEGER, iRoot, comm, ierror)

    return
  end subroutine parallel_bc_int_3d
  
  subroutine parallel_bc_real_0d(A)

    implicit none

    integer :: ierror
    real(WP) :: A

    call MPI_BCAST(A, 1, MPI_REAL_WP, iRoot, comm, ierror)

    return
  end subroutine parallel_bc_real_0d
  
  subroutine parallel_bc_real_1d(A)

    implicit none

    integer :: ierror
    real(WP), dimension(:) :: A

    call MPI_BCAST(A, size(A), MPI_REAL_WP, iRoot, comm, ierror)

    return
  end subroutine parallel_bc_real_1d
  
  subroutine parallel_bc_real_2d(A)

    implicit none

    integer :: ierror
    real(WP), dimension(:,:) :: A

    call MPI_BCAST(A, size(A), MPI_REAL_WP, iRoot, comm, ierror)

    return
  end subroutine parallel_bc_real_2d
  
  subroutine parallel_bc_real_3d(A)

    implicit none

    integer :: ierror
    real(WP), dimension(:,:,:) :: A

    call MPI_BCAST(A, size(A), MPI_REAL_WP, iRoot, comm, ierror)

    return
  end subroutine parallel_bc_real_3d

  ! MPI LOGICAL OR
  subroutine parallel_lor_0d(A, B)

    implicit none

    logical, intent(in)  :: A
    logical, intent(out) :: B
    integer :: ierror

    call MPI_ALLREDUCE(A, B, 1, MPI_LOGICAL, MPI_LOR, comm, ierror)

    return
  end subroutine parallel_lor_0d

  ! MPI LOGICAL OR IN PLACE
  subroutine parallel_lor_0d_in_place(A)

    implicit none

    logical, intent(inout)  :: A
    integer :: ierror

    call MPI_ALLREDUCE(MPI_IN_PLACE, A, 1, MPI_LOGICAL, MPI_LOR, comm, ierror)

    return
  end subroutine parallel_lor_0d_in_place

  ! MPI DIRECTIONAL LOGICAL OR
  subroutine parallel_lor_dir_0d(A,B,direction)

    implicit none
    logical, intent(in)  :: A
    logical, intent(out) :: B
    integer, intent(in) :: direction
    integer :: ierror

    call MPI_ALLREDUCE(A, B, 1, MPI_LOGICAL, MPI_LOR, commDir(direction), ierror)
    
    return
  end subroutine parallel_lor_dir_0d

  subroutine fill_ghost_points_int(arrayWithGhostPoints, direction, nGhostPoints,            &
       periodicOffset)

    implicit none

    ! Arguments
    integer, intent(inout) :: arrayWithGhostPoints(:,:,:,:)
    integer, intent(in) :: direction, nGhostPoints(2)
    integer, intent(in), optional :: periodicOffset(2)

    ! Local variables
    integer :: i, j, k, l, j2, nScalars, ny, gridSizeWithGhostPoints(3), gridSize(3),        &
         normalPlaneSize, periodicOffset_(2), request(4), tagPrev, tagNext,                  &
         rankOfPreviousProcess, rankOfNextProcess, ierror
    integer, dimension(MPI_STATUS_SIZE, 4) :: statuses
    integer, allocatable :: receivedFromPreviousProcess(:,:,:),                              &
         receivedFromNextProcess(:,:,:),                                                     &
         sentToPreviousProcess(:,:,:),                                                       &
         sentToNextProcess(:,:,:)
    logical :: polarCommunication(2)

    if (all(nGhostPoints .le. 0)) return
    if (nProcsDir(direction) .eq. 1 .and. nGhostPoints(1) .ne. nGhostPoints(2)) return

    ! Get array dimensions
    nScalars = size(arrayWithGhostPoints, 4)
    gridSizeWithGhostPoints = shape(arrayWithGhostPoints(:,:,:,1))
    normalPlaneSize = product(gridSizeWithGhostPoints) / size(arrayWithGhostPoints, direction)
    gridSize = gridSizeWithGhostPoints - sum(nGhostPoints)
    periodicOffset_ = 0
    if (present(periodicOffset)) periodicOffset_ = periodicOffset

    ! Initialize MPI variables
    request = MPI_REQUEST_NULL
    tagPrev = 1
    tagNext = 2
    call MPI_Cart_shift(comm, direction - 1, 1, rankOfPreviousProcess, rankOfNextProcess,    &
         ierror)

    ! Overwrite send/recv ranks for spherical geometry if first or last proc
    ! Assume: direction 1: radial, direction 2: azimuthal, direction 3: polar
    ! Assume: not decomposed in azimuthal direction
    polarCommunication = .false.
    if (periodicOffset_(1).eq.2 .or. periodicOffset_(2).eq.2) then
       if (direction .ne. 3) then
          call die('fill_ghost_points: periodic type POLAR must be direction 3!')
       end if
       if (nProcsDir(2) .ne. 1) then
          call die('fill_ghost_points: periodic type POLAR requires 1 proc in direction 2!')
       end if
       if (periodicOffset_(1) .eq. 2) then
          polarCommunication(1) = .true.
          ! Adjust rank of previous process
          rankOfPreviousProcess = iRank
       end if
       if (periodicOffset_(2) .eq. 2) then
          ! Adjust rank of next process
          polarCommunication(2) = .true.
          rankOfNextProcess = iRank
       end if
       ny = size(arrayWithGhostPoints, 2)
       periodicOffset_ = 0
    end if

    if (nGhostPoints(1) .gt. 0) then

       ! Allocate temporary buffers which hold send/receive data
       allocate(receivedFromPreviousProcess(nGhostPoints(1), normalPlaneSize, nScalars))
       allocate(sentToPreviousProcess(nGhostPoints(1), normalPlaneSize, nScalars))

       ! Non-blocking receive
       if (polarCommunication(1)) then
          call MPI_Irecv(receivedFromPreviousProcess, size(receivedFromPreviousProcess),     &
               MPI_INTEGER, rankOfPreviousProcess, tagPrev, comm, request(1), ierror)
       else
          call MPI_Irecv(receivedFromPreviousProcess, size(receivedFromPreviousProcess),     &
               MPI_INTEGER, rankOfPreviousProcess, tagNext, comm, request(1), ierror)
       end if

    end if

    if (nGhostPoints(2) .gt. 0) then

       ! Allocate temporary buffers which hold send/receive data
       allocate(receivedFromNextProcess(nGhostPoints(2), normalPlaneSize, nScalars))
       allocate(sentToNextProcess(nGhostPoints(2), normalPlaneSize, nScalars))

       ! Non-blocking receive
       if (polarCommunication(2)) then
          call MPI_Irecv(receivedFromNextProcess, size(receivedFromNextProcess),             &
               MPI_INTEGER, rankOfNextProcess, tagNext, comm, request(2), ierror)
       else
          call MPI_Irecv(receivedFromNextProcess, size(receivedFromNextProcess),             &
               MPI_INTEGER, rankOfNextProcess, tagPrev, comm, request(2), ierror)
       end if

    end if

    ! Copy ghost point data to send buffers:

    select case (direction)

    case (1)

       do l = 1, nScalars
          do k = 1, gridSizeWithGhostPoints(3)
             do j = 1, gridSizeWithGhostPoints(2)
                do i = 1, nGhostPoints(1)
                   sentToPreviousProcess(i, j + gridSizeWithGhostPoints(2) * (k - 1), l) =   &
                        arrayWithGhostPoints(i + nGhostPoints(1) +                           &
                        periodicOffset_(2), j, k, l)
                end do
             end do
          end do
       end do

       do l = 1, nScalars
          do k = 1, gridSizeWithGhostPoints(3)
             do j = 1, gridSizeWithGhostPoints(2)
                do i = 1, nGhostPoints(2)
                   sentToNextProcess(nGhostPoints(2) + 1 - i, j +                            &
                        gridSizeWithGhostPoints(2) * (k - 1), l) =                           &
                        arrayWithGhostPoints(gridSize(1) + 1 - i + nGhostPoints(1) -         &
                        periodicOffset_(1), j, k, l)
                end do
             end do
          end do
       end do

    case (2)

       do l = 1, nScalars
          do k = 1, gridSizeWithGhostPoints(3)
             do i = 1, gridSizeWithGhostPoints(1)
                do j = 1, nGhostPoints(1)
                   sentToPreviousProcess(j, i + gridSizeWithGhostPoints(1) * (k - 1), l) =   &
                        arrayWithGhostPoints(i, j + nGhostPoints(1) +                        &
                        periodicOffset_(2), k, l)
                end do
             end do
          end do
       end do

       do l = 1, nScalars
          do k = 1, gridSizeWithGhostPoints(3)
             do i = 1, gridSizeWithGhostPoints(1)
                do j = 1, nGhostPoints(2)
                   sentToNextProcess(nGhostPoints(2) + 1 - j, i +                            &
                        gridSizeWithGhostPoints(1) * (k - 1), l) =                           &
                        arrayWithGhostPoints(i, gridSize(2) + 1 - j + nGhostPoints(1) -      &
                        periodicOffset_(1), k, l)
                end do
             end do
          end do
       end do

    case (3)

       if (polarCommunication(1)) then

          do l = 1, nScalars
             do j = 1, gridSizeWithGhostPoints(2)
                j2 = mod(floor(0.5_WP * real(ny, WP)) + j - 1, ny - 1) + 1
                do i = 1, gridSizeWithGhostPoints(1)
                   do k = 1, nGhostPoints(1)
                      sentToPreviousProcess(nGhostPoints(1) + 1 - k, i +                     &
                           gridSizeWithGhostPoints(1) * (j - 1), l) =                        &
                           arrayWithGhostPoints(i, j2, k + nGhostPoints(1), l)
                   end do
                end do
             end do
          end do

       else

          do l = 1, nScalars
             do j = 1, gridSizeWithGhostPoints(2)
                do i = 1, gridSizeWithGhostPoints(1)
                   do k = 1, nGhostPoints(1)
                      sentToPreviousProcess(k, i + gridSizeWithGhostPoints(1) * (j - 1), l) =&
                           arrayWithGhostPoints(i, j, k + nGhostPoints(1) +                  &
                           periodicOffset_(2), l)
                   end do
                end do
             end do
          end do

       end if

       if (polarCommunication(2)) then

          do l = 1, nScalars
             do j = 1, gridSizeWithGhostPoints(2)
                j2 = mod(floor(0.5_WP * real(ny, WP)) + j - 1, ny - 1) + 1
                do i = 1, gridSizeWithGhostPoints(1)
                   do k = 1, nGhostPoints(2)
                      sentToNextProcess(k, i + gridSizeWithGhostPoints(1) * (j - 1), l) =    &
                           arrayWithGhostPoints(i, j2, gridSize(3) + 1 - k +                 &
                           nGhostPoints(1), l)
                   end do
                end do
             end do
          end do

       else

          do l = 1, nScalars
             do j = 1, gridSizeWithGhostPoints(2)
                do i = 1, gridSizeWithGhostPoints(1)
                   do k = 1, nGhostPoints(2)
                      sentToNextProcess(nGhostPoints(2) + 1 - k, i +                         &
                           gridSizeWithGhostPoints(1) * (j - 1), l) =                        &
                           arrayWithGhostPoints(i, j, gridSize(3) + 1 - k + nGhostPoints(1)  &
                           - periodicOffset_(1), l)
                   end do
                end do
             end do
          end do

       end if

    end select

    ! Non-blocking send followed by `MPI_Waitall`
    if (nGhostPoints(2) .gt. 0) then
       call MPI_Isend(sentToNextProcess, size(sentToNextProcess), MPI_INTEGER,               &
            rankOfNextProcess, tagNext, comm, request(3), ierror)
    end if
    if (nGhostPoints(1) .gt. 0) then
       call MPI_Isend(sentToPreviousProcess, size(sentToPreviousProcess), MPI_INTEGER,       &
            rankOfPreviousProcess, tagPrev, comm, request(4), ierror)
    end if
    call MPI_Waitall(4, request, statuses, ierror)

    ! Copy receive buffer data to ghost points:

    select case (direction)

    case (1)

       do l = 1, nScalars
          do k = 1, gridSizeWithGhostPoints(3)
             do j = 1, gridSizeWithGhostPoints(2)
                do i = 1, nGhostPoints(1)
                   arrayWithGhostPoints(i,j,k,l) =                                           &
                        receivedFromPreviousProcess(i, j +                                   &
                        gridSizeWithGhostPoints(2) * (k - 1), l)
                end do
             end do
          end do
       end do

       do l = 1, nScalars
          do k = 1, gridSizeWithGhostPoints(3)
             do j = 1, gridSizeWithGhostPoints(2)
                do i = 1, nGhostPoints(2)
                   arrayWithGhostPoints(gridSizeWithGhostPoints(1) + 1 - i, j, k, l) =       &
                        receivedFromNextProcess(nGhostPoints(2) + 1 - i,                     &
                        j + gridSizeWithGhostPoints(2) * (k - 1), l)
                end do
             end do
          end do
       end do

    case (2)

       do l = 1, nScalars
          do k = 1, gridSizeWithGhostPoints(3)
             do i = 1, gridSizeWithGhostPoints(1)
                do j = 1, nGhostPoints(1)
                   arrayWithGhostPoints(i,j,k,l) =                                           &
                        receivedFromPreviousProcess(j, i +                                   &
                        gridSizeWithGhostPoints(1) * (k - 1), l)
                end do
             end do
          end do
       end do

       do l = 1, nScalars
          do k = 1, gridSizeWithGhostPoints(3)
             do i = 1, gridSizeWithGhostPoints(1)
                do j = 1, nGhostPoints(2)
                   arrayWithGhostPoints(i, gridSizeWithGhostPoints(2) + 1 - j, k, l) =       &
                        receivedFromNextProcess(nGhostPoints(2) + 1 - j,                     &
                        i + gridSizeWithGhostPoints(1) * (k - 1), l)
                end do
             end do
          end do
       end do

    case (3)

       do l = 1, nScalars
          do j = 1, gridSizeWithGhostPoints(2)
             do i = 1, gridSizeWithGhostPoints(1)
                do k = 1, nGhostPoints(1)
                   arrayWithGhostPoints(i,j,k,l) =                                           &
                        receivedFromPreviousProcess(k, i +                                   &
                        gridSizeWithGhostPoints(1) * (j - 1), l)
                end do
             end do
          end do
       end do

       do l = 1, nScalars
          do j = 1, gridSizeWithGhostPoints(2)
             do i = 1, gridSizeWithGhostPoints(1)
                do k = 1, nGhostPoints(2)
                   arrayWithGhostPoints(i, j, gridSizeWithGhostPoints(3) + 1 - k, l) =       &
                        receivedFromNextProcess(nGhostPoints(2) + 1 - k,                     &
                        i + gridSizeWithGhostPoints(1) * (j - 1), l)
                end do
             end do
          end do
       end do

    end select

    if (allocated(sentToPreviousProcess)) deallocate(sentToPreviousProcess)
    if (allocated(sentToNextProcess)) deallocate(sentToNextProcess)
    if (allocated(receivedFromPreviousProcess)) deallocate(receivedFromPreviousProcess)
    if (allocated(receivedFromNextProcess)) deallocate(receivedFromNextProcess)

  end subroutine fill_ghost_points_int

  subroutine fill_ghost_points_real(arrayWithGhostPoints, direction, nGhostPoints,           &
       periodicOffset)

    implicit none

    ! Arguments
    real(WP), intent(inout) :: arrayWithGhostPoints(:,:,:,:)
    integer, intent(in) :: direction, nGhostPoints(2)
    integer, intent(in), optional :: periodicOffset(2)

    ! Local variables
    integer :: i, j, k, l, j2, nScalars, ny, gridSizeWithGhostPoints(3), gridSize(3),        &
         normalPlaneSize, periodicOffset_(2), request(4), tagPrev, tagNext,                  &
         rankOfPreviousProcess, rankOfNextProcess, ierror
    integer, dimension(MPI_STATUS_SIZE, 4) :: statuses
    real(WP), allocatable :: receivedFromPreviousProcess(:,:,:),                             &
         receivedFromNextProcess(:,:,:),                                                     &
         sentToPreviousProcess(:,:,:),                                                       &
         sentToNextProcess(:,:,:)
    logical :: polarCommunication(2)

    if (all(nGhostPoints .le. 0)) return
    if (nProcsDir(direction) .eq. 1 .and. nGhostPoints(1) .ne. nGhostPoints(2)) return

    ! Get array dimensions
    nScalars = size(arrayWithGhostPoints, 4)
    gridSizeWithGhostPoints = shape(arrayWithGhostPoints(:,:,:,1))
    normalPlaneSize = product(gridSizeWithGhostPoints) / size(arrayWithGhostPoints, direction)
    gridSize = gridSizeWithGhostPoints - sum(nGhostPoints)
    periodicOffset_ = 0
    if (present(periodicOffset)) periodicOffset_ = periodicOffset

    ! Initialize MPI variables
    request = MPI_REQUEST_NULL
    tagPrev = 1
    tagNext = 2
    call MPI_Cart_shift(comm, direction - 1, 1, rankOfPreviousProcess, rankOfNextProcess,    &
         ierror)

    ! Overwrite send/recv ranks for spherical geometry if first or last proc
    ! Assume: direction 1: radial, direction 2: azimuthal, direction 3: polar
    ! Assume: not decomposed in azimuthal direction
    polarCommunication = .false.
    if (periodicOffset_(1).eq.2 .or. periodicOffset_(2).eq.2) then
       if (direction .ne. 3) then
          call die('fill_ghost_points: periodic type POLAR must be direction 3!')
       end if
       if (nProcsDir(2) .ne. 1) then
          call die('fill_ghost_points: periodic type POLAR requires 1 proc in direction 2!')
       end if
       if (periodicOffset_(1) .eq. 2) then
          polarCommunication(1) = .true.
          ! Adjust rank of previous process
          rankOfPreviousProcess = iRank
       end if
       if (periodicOffset_(2) .eq. 2) then
          ! Adjust rank of next process
          polarCommunication(2) = .true.
          rankOfNextProcess = iRank
       end if
       ny = size(arrayWithGhostPoints, 2)
       periodicOffset_ = 0
    end if

    if (nGhostPoints(1) .gt. 0) then

       ! Allocate temporary buffers which hold send/receive data
       allocate(receivedFromPreviousProcess(nGhostPoints(1), normalPlaneSize, nScalars))
       allocate(sentToPreviousProcess(nGhostPoints(1), normalPlaneSize, nScalars))

       ! Non-blocking receive
       if (polarCommunication(1)) then
          call MPI_Irecv(receivedFromPreviousProcess, size(receivedFromPreviousProcess),     &
               MPI_REAL_WP, rankOfPreviousProcess, tagPrev, comm, request(1), ierror)
       else
          call MPI_Irecv(receivedFromPreviousProcess, size(receivedFromPreviousProcess),     &
               MPI_REAL_WP, rankOfPreviousProcess, tagNext, comm, request(1), ierror)
       end if

    end if

    if (nGhostPoints(2) .gt. 0) then

       ! Allocate temporary buffers which hold send/receive data
       allocate(receivedFromNextProcess(nGhostPoints(2), normalPlaneSize, nScalars))
       allocate(sentToNextProcess(nGhostPoints(2), normalPlaneSize, nScalars))

       ! Non-blocking receive
       if (polarCommunication(2)) then
          call MPI_Irecv(receivedFromNextProcess, size(receivedFromNextProcess),             &
               MPI_REAL_WP, rankOfNextProcess, tagNext, comm, request(2), ierror)
       else
          call MPI_Irecv(receivedFromNextProcess, size(receivedFromNextProcess),             &
               MPI_REAL_WP, rankOfNextProcess, tagPrev, comm, request(2), ierror)
       end if

    end if

    ! Copy ghost point data to send buffers:

    select case (direction)

    case (1)

       do l = 1, nScalars
          do k = 1, gridSizeWithGhostPoints(3)
             do j = 1, gridSizeWithGhostPoints(2)
                do i = 1, nGhostPoints(1)
                   sentToPreviousProcess(i, j + gridSizeWithGhostPoints(2) * (k - 1), l) =   &
                        arrayWithGhostPoints(i + nGhostPoints(1) +                           &
                        periodicOffset_(2), j, k, l)
                end do
             end do
          end do
       end do

       do l = 1, nScalars
          do k = 1, gridSizeWithGhostPoints(3)
             do j = 1, gridSizeWithGhostPoints(2)
                do i = 1, nGhostPoints(2)
                   sentToNextProcess(nGhostPoints(2) + 1 - i, j +                            &
                        gridSizeWithGhostPoints(2) * (k - 1), l) =                           &
                        arrayWithGhostPoints(gridSize(1) + 1 - i + nGhostPoints(1) -         &
                        periodicOffset_(1), j, k, l)
                end do
             end do
          end do
       end do

    case (2)

       do l = 1, nScalars
          do k = 1, gridSizeWithGhostPoints(3)
             do i = 1, gridSizeWithGhostPoints(1)
                do j = 1, nGhostPoints(1)
                   sentToPreviousProcess(j, i + gridSizeWithGhostPoints(1) * (k - 1), l) =   &
                        arrayWithGhostPoints(i, j + nGhostPoints(1) +                        &
                        periodicOffset_(2), k, l)
                end do
             end do
          end do
       end do

       do l = 1, nScalars
          do k = 1, gridSizeWithGhostPoints(3)
             do i = 1, gridSizeWithGhostPoints(1)
                do j = 1, nGhostPoints(2)
                   sentToNextProcess(nGhostPoints(2) + 1 - j, i +                            &
                        gridSizeWithGhostPoints(1) * (k - 1), l) =                           &
                        arrayWithGhostPoints(i, gridSize(2) + 1 - j + nGhostPoints(1) -      &
                        periodicOffset_(1), k, l)
                end do
             end do
          end do
       end do
       
    case (3)

       if (polarCommunication(1)) then

          do l = 1, nScalars
             do j = 1, gridSizeWithGhostPoints(2)
                j2 = mod(floor(0.5_WP * real(ny, WP)) + j - 1, ny - 1) + 1
                do i = 1, gridSizeWithGhostPoints(1)
                   do k = 1, nGhostPoints(1)
                      sentToPreviousProcess(nGhostPoints(1) + 1 - k, i +                     &
                           gridSizeWithGhostPoints(1) * (j - 1), l) =                        &
                           arrayWithGhostPoints(i, j2, k + nGhostPoints(1), l)
                   end do
                end do
             end do
          end do

       else

          do l = 1, nScalars
             do j = 1, gridSizeWithGhostPoints(2)
                do i = 1, gridSizeWithGhostPoints(1)
                   do k = 1, nGhostPoints(1)
                      sentToPreviousProcess(k, i + gridSizeWithGhostPoints(1) * (j - 1), l) =&
                           arrayWithGhostPoints(i, j, k + nGhostPoints(1) +                  &
                           periodicOffset_(2), l)
                   end do
                end do
             end do
          end do

       end if

       if (polarCommunication(2)) then

          do l = 1, nScalars
             do j = 1, gridSizeWithGhostPoints(2)
                j2 = mod(floor(0.5_WP * real(ny, WP)) + j - 1, ny - 1) + 1
                do i = 1, gridSizeWithGhostPoints(1)
                   do k = 1, nGhostPoints(2)
                      sentToNextProcess(k, i + gridSizeWithGhostPoints(1) * (j - 1), l) =    &
                           arrayWithGhostPoints(i, j2, gridSize(3) + 1 - k +                 &
                           nGhostPoints(1), l)
                   end do
                end do
             end do
          end do

       else

          do l = 1, nScalars
             do j = 1, gridSizeWithGhostPoints(2)
                do i = 1, gridSizeWithGhostPoints(1)
                   do k = 1, nGhostPoints(2)
                      sentToNextProcess(nGhostPoints(2) + 1 - k, i +                         &
                           gridSizeWithGhostPoints(1) * (j - 1), l) =                        &
                           arrayWithGhostPoints(i, j, gridSize(3) + 1 - k + nGhostPoints(1)  &
                           - periodicOffset_(1), l)
                   end do
                end do
             end do
          end do

       end if

    end select

    ! Non-blocking send followed by `MPI_Waitall`
    if (nGhostPoints(2) .gt. 0) then
       call MPI_Isend(sentToNextProcess, size(sentToNextProcess), MPI_REAL_WP,               &
            rankOfNextProcess, tagNext, comm, request(3), ierror)
    end if
    if (nGhostPoints(1) .gt. 0) then
       call MPI_Isend(sentToPreviousProcess, size(sentToPreviousProcess), MPI_REAL_WP,       &
            rankOfPreviousProcess, tagPrev, comm, request(4), ierror)
    end if
    call MPI_Waitall(4, request, statuses, ierror)

    ! Copy receive buffer data to ghost points:

    select case (direction)

    case (1)

       do l = 1, nScalars
          do k = 1, gridSizeWithGhostPoints(3)
             do j = 1, gridSizeWithGhostPoints(2)
                do i = 1, nGhostPoints(1)
                   arrayWithGhostPoints(i,j,k,l) =                                           &
                        receivedFromPreviousProcess(i, j +                                   &
                        gridSizeWithGhostPoints(2) * (k - 1), l)
                end do
             end do
          end do
       end do

       do l = 1, nScalars
          do k = 1, gridSizeWithGhostPoints(3)
             do j = 1, gridSizeWithGhostPoints(2)
                do i = 1, nGhostPoints(2)
                   arrayWithGhostPoints(gridSizeWithGhostPoints(1) + 1 - i, j, k, l) =       &
                        receivedFromNextProcess(nGhostPoints(2) + 1 - i,                     &
                        j + gridSizeWithGhostPoints(2) * (k - 1), l)
                end do
             end do
          end do
       end do

    case (2)

       do l = 1, nScalars
          do k = 1, gridSizeWithGhostPoints(3)
             do i = 1, gridSizeWithGhostPoints(1)
                do j = 1, nGhostPoints(1)
                   arrayWithGhostPoints(i,j,k,l) =                                           &
                        receivedFromPreviousProcess(j, i +                                   &
                        gridSizeWithGhostPoints(1) * (k - 1), l)
                end do
             end do
          end do
       end do

       do l = 1, nScalars
          do k = 1, gridSizeWithGhostPoints(3)
             do i = 1, gridSizeWithGhostPoints(1)
                do j = 1, nGhostPoints(2)
                   arrayWithGhostPoints(i, gridSizeWithGhostPoints(2) + 1 - j, k, l) =       &
                        receivedFromNextProcess(nGhostPoints(2) + 1 - j,                     &
                        i + gridSizeWithGhostPoints(1) * (k - 1), l)
                end do
             end do
          end do
       end do

    case (3)

       do l = 1, nScalars
          do j = 1, gridSizeWithGhostPoints(2)
             do i = 1, gridSizeWithGhostPoints(1)
                do k = 1, nGhostPoints(1)
                   arrayWithGhostPoints(i,j,k,l) =                                           &
                        receivedFromPreviousProcess(k, i +                                   &
                        gridSizeWithGhostPoints(1) * (j - 1), l)
                end do
             end do
          end do
       end do

       do l = 1, nScalars
          do j = 1, gridSizeWithGhostPoints(2)
             do i = 1, gridSizeWithGhostPoints(1)
                do k = 1, nGhostPoints(2)
                   arrayWithGhostPoints(i, j, gridSizeWithGhostPoints(3) + 1 - k, l) =       &
                        receivedFromNextProcess(nGhostPoints(2) + 1 - k,                     &
                        i + gridSizeWithGhostPoints(1) * (j - 1), l)
                end do
             end do
          end do
       end do

    end select

    if (allocated(sentToPreviousProcess)) deallocate(sentToPreviousProcess)
    if (allocated(sentToNextProcess)) deallocate(sentToNextProcess)
    if (allocated(receivedFromPreviousProcess)) deallocate(receivedFromPreviousProcess)
    if (allocated(receivedFromNextProcess)) deallocate(receivedFromNextProcess)

  end subroutine fill_ghost_points_real

  subroutine border_summation(arrayWithGhostPoints, direction, nGhostPoints, periodicOffset)

    implicit none

    ! Arguments
    real(WP), intent(inout) :: arrayWithGhostPoints(:,:,:,:)
    integer, intent(in) :: direction, nGhostPoints(2)
    integer, intent(in), optional :: periodicOffset(2)

    ! Local variables
    integer :: i, j, k, l, nScalars, nDimensions, gridSizeWithGhostPoints(3), gridSize(3),   &
         normalPlaneSize, periodicOffset_(2), request(4), tagPrev, tagNext,                  &
         rankOfPreviousProcess, rankOfNextProcess, ierror
    integer, dimension(MPI_STATUS_SIZE, 4) :: statuses
    real(WP), allocatable :: receivedFromPreviousProcess(:,:,:),                             &
         receivedFromNextProcess(:,:,:),                                                     &
         sentToPreviousProcess(:,:,:),                                                       &
         sentToNextProcess(:,:,:)

    ! Get number of dimensions from communicator
    call MPI_Cartdim_get(comm, nDimensions, ierror)

    if (all(nGhostPoints .le. 0)) return
    if (nProcsDir(direction) .eq. 1 .and. nGhostPoints(1) .ne. nGhostPoints(2)) return

    ! Get array dimensions
    nScalars = size(arrayWithGhostPoints, 4)
    gridSizeWithGhostPoints = shape(arrayWithGhostPoints(:,:,:,1))
    normalPlaneSize = product(gridSizeWithGhostPoints) / size(arrayWithGhostPoints, direction)
    gridSize = gridSizeWithGhostPoints - sum(nGhostPoints)
    periodicOffset_ = 0
    if (present(periodicOffset)) periodicOffset_ = periodicOffset

    ! Initialize MPI variables
    request = MPI_REQUEST_NULL
    tagPrev = 1
    tagNext = 2
    call MPI_Cart_shift(comm, direction - 1, 1, rankOfPreviousProcess, rankOfNextProcess,    &
         ierror)

    if (nGhostPoints(1) .gt. 0) then

       ! Allocate temporary buffers which hold send/receive data
       allocate(receivedFromPreviousProcess(nGhostPoints(1), normalPlaneSize, nScalars))
       allocate(sentToPreviousProcess(nGhostPoints(1), normalPlaneSize, nScalars))

       ! Non-blocking receive
       call MPI_Irecv(receivedFromPreviousProcess, size(receivedFromPreviousProcess),        &
            MPI_REAL_WP, rankOfPreviousProcess, tagNext, comm, request(1), ierror)

    end if

    if (nGhostPoints(2) .gt. 0) then

       ! Allocate temporary buffers which hold send/receive data
       allocate(receivedFromNextProcess(nGhostPoints(2), normalPlaneSize, nScalars))
       allocate(sentToNextProcess(nGhostPoints(2), normalPlaneSize, nScalars))

       ! Non-blocking receive
       call MPI_Irecv(receivedFromNextProcess, size(receivedFromNextProcess), MPI_REAL_WP,   &
            rankOfNextProcess, tagPrev, comm, request(2), ierror)

    end if

    ! Copy ghost point data to send buffers:

    select case (direction)

    case (1)

       do l = 1, nScalars
          do k = 1, gridSizeWithGhostPoints(3)
             do j = 1, gridSizeWithGhostPoints(2)
                do i = 1, nGhostPoints(1)
                   sentToPreviousProcess(i, j + gridSizeWithGhostPoints(2) * (k - 1), l) =   &
                        arrayWithGhostPoints(i + periodicOffset_(2), j, k, l)
                end do
             end do
          end do
       end do

       do l = 1, nScalars
          do k = 1, gridSizeWithGhostPoints(3)
             do j = 1, gridSizeWithGhostPoints(2)
                do i = 1, nGhostPoints(2)
                   sentToNextProcess(nGhostPoints(2) + 1 - i, j +                          &
                        gridSizeWithGhostPoints(2) * (k - 1), l) =                           &
                        arrayWithGhostPoints(gridSize(1) + 1 - i + sum(nGhostPoints) -     &
                        periodicOffset_(1), j, k, l)
                end do
             end do
          end do
       end do

    case (2)

       do l = 1, nScalars
          do k = 1, gridSizeWithGhostPoints(3)
             do i = 1, gridSizeWithGhostPoints(1)
                do j = 1, nGhostPoints(1)
                   sentToPreviousProcess(j, i + gridSizeWithGhostPoints(1) * (k - 1), l) =   &
                        arrayWithGhostPoints(i, j + periodicOffset_(2), k, l)
                end do
             end do
          end do
       end do

       do l = 1, nScalars
          do k = 1, gridSizeWithGhostPoints(3)
             do i = 1, gridSizeWithGhostPoints(1)
                do j = 1, nGhostPoints(2)
                   sentToNextProcess(nGhostPoints(2) + 1 - j, i +                            &
                        gridSizeWithGhostPoints(1) * (k - 1), l) =                           &
                        arrayWithGhostPoints(i, gridSize(2) + 1 - j + sum(nGhostPoints) -    &
                        periodicOffset_(1), k, l)
                end do
             end do
          end do
       end do

    case (3)

       do l = 1, nScalars
          do j = 1, gridSizeWithGhostPoints(2)
             do i = 1, gridSizeWithGhostPoints(1)
                do k = 1, nGhostPoints(1)
                   sentToPreviousProcess(k, i + gridSizeWithGhostPoints(1) * (j - 1), l) =   &
                        arrayWithGhostPoints(i, j, k + periodicOffset_(2), l)
                end do
             end do
          end do
       end do

       do l = 1, nScalars
          do j = 1, gridSizeWithGhostPoints(2)
             do i = 1, gridSizeWithGhostPoints(1)
                do k = 1, nGhostPoints(2)
                   sentToNextProcess(nGhostPoints(2) + 1 - k, i +                            &
                        gridSizeWithGhostPoints(1) * (j - 1), l) =                           &
                        arrayWithGhostPoints(i, j, gridSize(3) + 1 - k + sum(nGhostPoints)   &
                        - periodicOffset_(1), l)
                end do
             end do
          end do
       end do

    end select

    ! Non-blocking send followed by `MPI_Waitall`
    if (nGhostPoints(2) .gt. 0) then
       call MPI_Isend(sentToNextProcess, size(sentToNextProcess), MPI_REAL_WP,               &
            rankOfNextProcess, tagNext, comm, request(3), ierror)
    end if
    if (nGhostPoints(1) .gt. 0) then
       call MPI_Isend(sentToPreviousProcess, size(sentToPreviousProcess), MPI_REAL_WP,       &
            rankOfPreviousProcess, tagPrev, comm, request(4), ierror)
    end if
    call MPI_Waitall(4, request, statuses, ierror)

    ! Copy receive buffer data to ghost points:

    select case (direction)

    case (1)

       do l = 1, nScalars
          do k = 1, gridSizeWithGhostPoints(3)
             do j = 1, gridSizeWithGhostPoints(2)
                do i = 1, nGhostPoints(1)
                   arrayWithGhostPoints(i+nGhostPoints(1),j,k,l) =                           &
                        arrayWithGhostPoints(i+nGhostPoints(1),j,k,l) +                      &
                        receivedFromPreviousProcess(i, j +                                   &
                        gridSizeWithGhostPoints(2) * (k - 1), l)
                end do
             end do
          end do
       end do

       do l = 1, nScalars
          do k = 1, gridSizeWithGhostPoints(3)
             do j = 1, gridSizeWithGhostPoints(2)
                do i = 1, nGhostPoints(2)
                   arrayWithGhostPoints(gridSize(1) + nGhostPoints(1) + 1 - i, j, k, l) =    &
                        arrayWithGhostPoints(gridSize(1) + &
                        nGhostPoints(1) + 1 - i, j, k, l)                                    &
                        + receivedFromNextProcess(nGhostPoints(2) + 1 - i,                   &
                        j + gridSizeWithGhostPoints(2) * (k - 1), l)
                end do
             end do
          end do
       end do

    case (2)

       do l = 1, nScalars
          do k = 1, gridSizeWithGhostPoints(3)
             do i = 1, gridSizeWithGhostPoints(1)
                do j = 1, nGhostPoints(1)
                   arrayWithGhostPoints(i,j+nGhostPoints(1),k,l) =                           &
                        arrayWithGhostPoints(i,j+nGhostPoints(1),k,l) +                      &
                        receivedFromPreviousProcess(j, i +                                   &
                        gridSizeWithGhostPoints(1) * (k - 1), l)
                end do
             end do
          end do
       end do

       do l = 1, nScalars
          do k = 1, gridSizeWithGhostPoints(3)
             do i = 1, gridSizeWithGhostPoints(1)
                do j = 1, nGhostPoints(2)
                   arrayWithGhostPoints(i, gridSize(2) + nGhostPoints(1) + 1 - j, k, l) =    &
                        arrayWithGhostPoints(i, gridSize(2) +                                &
                        nGhostPoints(1) + 1 - j, k, l)                                       &
                        + receivedFromNextProcess(nGhostPoints(2) + 1 - j,                   &
                        i + gridSizeWithGhostPoints(1) * (k - 1), l)
                end do
             end do
          end do
       end do

    case (3)

       do l = 1, nScalars
          do j = 1, gridSizeWithGhostPoints(2)
             do i = 1, gridSizeWithGhostPoints(1)
                do k = 1, nGhostPoints(1)
                   arrayWithGhostPoints(i,j,k+nGhostPoints(1),l) =                           &
                        arrayWithGhostPoints(i,j,k+nGhostPoints(1),l) +                      &
                        receivedFromPreviousProcess(k, i +                                   &
                        gridSizeWithGhostPoints(1) * (j - 1), l)
                end do
             end do
          end do
       end do

       do l = 1, nScalars
          do j = 1, gridSizeWithGhostPoints(2)
             do i = 1, gridSizeWithGhostPoints(1)
                do k = 1, nGhostPoints(2)
                   arrayWithGhostPoints(i, j, gridSize(3) + nGhostPoints(1) + 1 - k, l) =    &
                        arrayWithGhostPoints(i, j, gridSize(3) +                             &
                        nGhostPoints(1) + 1 - k, l)                                          &
                        + receivedFromNextProcess(nGhostPoints(2) + 1 - k,                   &
                        i + gridSizeWithGhostPoints(1) * (j - 1), l)
                end do
             end do
          end do
       end do

    end select

    if (allocated(sentToPreviousProcess)) deallocate(sentToPreviousProcess)
    if (allocated(sentToNextProcess)) deallocate(sentToNextProcess)
    if (allocated(receivedFromPreviousProcess)) deallocate(receivedFromPreviousProcess)
    if (allocated(receivedFromNextProcess)) deallocate(receivedFromNextProcess)

  end subroutine border_summation

  subroutine gather_along_direction(localArray, localSize, direction, offsetAlongDirection,  &
       gatheredArray)

    implicit none

    ! Arguments
    real(WP), intent(in) :: localArray(:,:)
    integer, intent(in) :: localSize(3), direction, offsetAlongDirection
    real(WP), intent(out) :: gatheredArray(:,:)

    ! Local variables
    integer :: i, j, k, l, nScalars, globalSizeAlongDirection, nProcsInPencil, ierror
    integer, allocatable :: localSizes(:), offsets(:)
    real(WP), allocatable :: sendBuffer(:,:), receiveBuffer(:,:)

    ! Get the dimensions of the local and gathered arrays.
    globalSizeAlongDirection = size(gatheredArray, 1) /                                      &
         (size(localArray, 1) / localSize(direction))
    nScalars = size(localArray, 2)

    ! Gather the offsets and sizes of `localArray` along dimension `direction` of all
    ! processes in the pencil communicator.
    call MPI_Comm_size(commDir(direction), nProcsInPencil, ierror)
    allocate(localSizes(nProcsInPencil), offsets(nProcsInPencil))
    call MPI_Allgather(localSize(direction), 1, MPI_INTEGER,                                 &
         localSizes, 1, MPI_INTEGER, commDir(direction), ierror)
    call MPI_Allgather(offsetAlongDirection, 1, MPI_INTEGER, offsets,                        &
         1, MPI_INTEGER, commDir(direction), ierror)

    ! Find the sizes and displacements of data to be gathered
    localSizes = localSizes * size(localArray) / localSize(direction)
    offsets = offsets * size(localArray) / localSize(direction)

    ! Allocate send and receive buffers with shapes that can be gathered easily using
    ! `MPI_Allgatherv`.
    allocate(sendBuffer(nScalars * size(localArray) / localSize(direction),                  &
         localSize(direction)))
    allocate(receiveBuffer(nScalars * size(localArray) / localSize(direction),               &
         globalSizeAlongDirection))

    ! Pack the data to be gathered:

    select case (direction)

    case (1)
       do l = 1, nScalars
          do k = 1, localSize(3)
             do j = 1, localSize(2)
                do i = 1, localSize(1)
                   sendBuffer(j + localSize(2) * (k - 1 + localSize(3) * (l - 1)), i) =      &
                        localArray(i + localSize(1) * (j - 1 + localSize(2) * (k - 1)), l)
                end do
             end do
          end do
       end do

    case (2)
       do l = 1, nScalars
          do k = 1, localSize(3)
             do j = 1, localSize(2)
                do i = 1, localSize(1)
                   sendBuffer(i + localSize(1) * (k - 1 + localSize(3) * (l - 1)), j) =      &
                        localArray(i + localSize(1) * (j - 1 + localSize(2) * (k - 1)), l)
                end do
             end do
          end do
       end do

    case (3)
       do l = 1, nScalars
          do k = 1, localSize(3)
             do j = 1, localSize(2)
                do i = 1, localSize(1)
                   sendBuffer(i + localSize(1) * (j - 1 + localSize(2) * (l - 1)), k) =      &
                        localArray(i + localSize(1) * (j - 1 + localSize(2) * (k - 1)), l)
                end do
             end do
          end do
       end do

    end select

    call MPI_Allgatherv(sendBuffer, size(sendBuffer), MPI_REAL_WP, receiveBuffer,            &
         localSizes, offsets, MPI_REAL_WP, commDir(direction), ierror)

    ! Unpack the data that was gathered:

    select case (direction)

    case (1)
       do l = 1, nScalars
          do k = 1, localSize(3)
             do j = 1, localSize(2)
                do i = 1, globalSizeAlongDirection
                   gatheredArray(i + globalSizeAlongDirection * (j - 1 +                     &
                        localSize(2) * (k - 1)), l) =                                        &
                        receiveBuffer(j + localSize(2) * (k - 1 + localSize(3) * (l - 1)), i)
                end do
             end do
          end do
       end do

    case (2)
       do l = 1, nScalars
          do k = 1, localSize(3)
             do j = 1, globalSizeAlongDirection
                do i = 1, localSize(1)
                   gatheredArray(i + localSize(1) * (j - 1 +                                 &
                        globalSizeAlongDirection * (k - 1)), l) =                            &
                        receiveBuffer(i + localSize(1) * (k - 1 + localSize(3) * (l - 1)), j)
                end do
             end do
          end do
       end do

    case (3)
       do l = 1, nScalars
          do k = 1, globalSizeAlongDirection
             do j = 1, localSize(2)
                do i = 1, localSize(1)
                   gatheredArray(i + localSize(1) * (j - 1 + localSize(2) * (k - 1)), l) =   &
                        receiveBuffer(i + localSize(1) * (j - 1 + localSize(2) * (l - 1)), k)
                end do
             end do
          end do
       end do

    end select
    
    if (allocated(receiveBuffer)) deallocate(receiveBuffer)
    if (allocated(sendBuffer)) deallocate(sendBuffer)
    if (allocated(offsets)) deallocate(offsets)
    if (allocated(localSizes)) deallocate(localSizes)

  end subroutine gather_along_direction

  subroutine generate_comm

    implicit none

    ! Local variables
    integer :: ierror, color, key

    color = 0 ! Procs with the same color will be in comm
    key = globalCommRank ! Determines the ordering (rank) of the procs in comm from
    ! the smallest to the greatest values of key

    call MPI_Comm_Split(globalComm, color, key, comm, ierror)

    ! Store processor information
    call MPI_Comm_Rank(comm, iRank, ierror)
    call MPI_Comm_Size(comm, nProcs, ierror) 
    iRoot = 0

    return
  end subroutine generate_comm

  subroutine generate_external_communicator

    implicit none

    ! Local variables
    integer :: ierror, color, key

    color = 0 ! Procs with the same color will be in externalComm
    key = 0 ! Determines the ordering (rank) of the procs in externalcomm from
    ! the smallest to the greatest values of key

    ! Only iRoot from the Fortran code will communicate with externalComm
    if (iRank.ne.iRoot) color = MPI_UNDEFINED
    call MPI_Comm_Split(globalComm, color, key, externalComm, ierror)

    ! only iRoot can use externalComm
    if (iRank.eq.iRoot) then
       call MPI_Comm_Size(externalComm, externalCommSize, ierror)
       if (externalCommSize.ne.2) call parallel_kill('generate_external_communicator: &
            &Failed to split into 2-rank Comm!')

       call MPI_Comm_Rank(externalComm, externalCommRank, ierror)
       externalCommRoot = 0
       if (externalCommRank.ne.externalCommRoot)                                               &
            call parallel_kill('generate_external_communicator: Failed to get externalCommRank&
            &=0 in the externalComm')
    else
       externalComm = MPI_COMM_NULL
       externalCommSize = 0
       externalCommRank = -1
    end if

    print *,'jCode: rank/size of globalComm: ', globalCommRank, '/', globalCommSize,           &
         'rank/size of comm: ', iRank, '/', nProcs,                                            &
         'rank/size of externalComm: ', externalCommRank, '/', externalCommSize

    return
  end subroutine generate_external_communicator

  subroutine parallel_kill(errorMessage)

    implicit none

    ! Local variables
    integer :: ierror
    character(len = *), intent(in), optional :: errorMessage

    ! Specify who sends the abort signal and what it means
    if (present(errorMessage)) then
       write(*,*) char(27) // '[1;31mERROR' // char(27) // '[0m: ' // trim(errorMessage)
    else
       write(*,*) char(27) // '[1;31mERROR' // char(27) //                                   &
            '[0m: initiated general abort signal due to an unknown error'
    end if

    ! Call general abort
    call MPI_ABORT(MPI_COMM_WORLD, iRoot, ierror)

    return
  end subroutine parallel_kill

end module parallel

subroutine parallel_get_inputname(filename)

  ! Internal modules
  use parallel

  ! External modules
  use parser
  use string

  implicit none

  ! Arguments
  character(len = str_medium), intent(out) :: filename

  ! Local variables
  integer :: ierror

  if (irank .eq. iroot) then
     if (command_argument_count() .ge. 1) then
        call get_command_argument(1, filename)
        if (filename(1:1) .eq. '-' .or. len_trim(filename) .eq. 0) then
           print *, 'No input file name was detected, using "input".'
           filename = 'input'
        end if
     else
        print *, 'No input file name was detected, using "input".'
        filename = 'input'
     end if
  end if

  call MPI_BCAST(filename, str_medium, MPI_CHARACTER, 0, comm, ierror)

end subroutine parallel_get_inputname

subroutine parallel_init

  ! Internal modules
  use parallel

  implicit none

  ! Local variables
  integer :: ierror
  integer :: sizeReal, sizeDP, sizeComplex, sizeComplexDP

  ! Initialize a basic MPI environment
  call MPI_Init(ierror)

  ! MPI_COMM_WORLD contains the total procs called by jcode and external programs (if any)
  globalComm = MPI_COMM_WORLD

  ! Store processor information
  call MPI_Comm_Rank(globalComm, globalCommRank, ierror)
  call MPI_Comm_Size(globalComm, globalCommSize, ierror)
  globalCommRoot = 0

  ! Generate jcode communication (comm) and external communication (if any)
  call generate_comm
  externalComm = MPI_COMM_NULL
  externalCommSize = 0
  externalCommRank = -1
  if (globalCommSize - nProcs .gt. 0) call generate_external_communicator

  ! Set MPI working precision - WP
  call MPI_TYPE_SIZE(MPI_REAL, sizeReal, ierror)
  call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, sizeDP, ierror)
  if (WP .eq. sizeReal) then
     MPI_REAL_WP = MPI_REAL
  else if (WP .eq. sizeDP) then
     MPI_REAL_WP = MPI_DOUBLE_PRECISION
  else
     call parallel_kill('Error in parallel_init: no WP equivalent in MPI')
  end if

  ! Set MPI single precision
  call MPI_TYPE_SIZE(MPI_REAL, sizeReal, ierror)
  call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, sizeDP, ierror)
  if (SP .eq. sizeReal) then
     MPI_REAL_SP = MPI_REAL
  else if (SP .eq. sizeDP) then
     MPI_REAL_SP = MPI_DOUBLE_PRECISION
  else
     call parallel_kill('Error in parallel_init: no SP equivalent in MPI')
  end if

  ! Set MPI working precision (for complex types)
  call MPI_TYPE_SIZE(MPI_COMPLEX,sizeComplex,ierror)
  call MPI_TYPE_SIZE(MPI_DOUBLE_COMPLEX,sizeComplexDP,ierror)
  if (2*WP .eq. sizeComplex) then
     MPI_COMPLEX_WP = MPI_COMPLEX
  else if (2*WP .eq. sizeComplexDP) then
     MPI_COMPLEX_WP = MPI_DOUBLE_COMPLEX
  else
     call parallel_kill('Error in parallel_init: no complex WP equivalent in MPI')
  end if

  ! By default allow for manual processor decomposition (init_flow will turn this off)
  disableManualDecomp = .false.

  return
end subroutine parallel_init

subroutine parallel_init_topology(nDimensions, isPeriodic)

  ! Internal modules
  use parallel

  ! External modules
  use string
  use parser

  implicit none

  ! Arguments
  integer, intent(in) :: nDimensions
  logical, intent(in) :: isPeriodic(3)

  ! Local variables
  integer :: i, oldComm, ierror
  logical :: manualProcDecomp, keepThisDirection(3)
  character(len = str_long) :: message

  ! Get the process distribution
  nProcsDir = 1
  call parser_read('use manual domain decomposition', manualProcDecomp, .false.)
  if (manualProcDecomp .and. .not.disableManualDecomp) then
     call parser_read('processor decomposition', nProcsDir)
  else
     nProcsDir(1:nDimensions) = 0
     call MPI_Dims_create(nProcs, nDimensions, nProcsDir(1:nDimensions), ierror)
  end if

  ! Validate the process distribution
  if (product(nProcsDir) .ne. nProcs) then
     write(message, '(2(A,I0.0),A)') 'Invalid process distribution: expected a total of ',   &
          nProcs, ' procs, got ', product(nProcsDir), '!'
     call parallel_kill(trim(message))
  end if

  ! Create a Cartesian communicator
  oldComm = comm
  call MPI_Cart_create(oldComm, nDimensions, nProcsDir(1:nDimensions),                       &
       isPeriodic(1:nDimensions), .true., comm, ierror)

  ! Find process coordinates in Cartesian topology
  procCoords = 0
  call MPI_Comm_rank(comm, iRank, ierror)
  call MPI_Cart_coords(comm, iRank, nDimensions, procCoords(1:nDimensions), ierror)

  ! Create a pencil communicator and find process coordinates in each direction
  commDir = MPI_COMM_NULL
  iRankDir = -1
  do i = 1, nDimensions
     keepThisDirection = .false.; keepThisDirection(i) = .true.
     call MPI_Cart_sub(comm, keepThisDirection(1:nDimensions), commDir(i), ierror)
     call MPI_Comm_Rank(commDir(i), iRankDir(i), ierror)
  end do

  return
end subroutine parallel_init_topology

subroutine parallel_init_io

  ! Internal modules
  use parallel

  ! External modules
  use parser

  implicit none

  ! Local variables
  integer :: ierror
  
  ! Create info for parallel i/o
  call parser_read('MPIIO fs', mpiiofs, 'win')
  select case(trim(mpiiofs))
  case('win')
     mpiiofs = ""
     mpiInfo = MPI_INFO_NULL
  case('ufs')
     mpiiofs = "ufs:"
     mpiInfo = MPI_INFO_NULL
  case('lustre')
     mpiiofs = "lustre:"
     call MPI_INFO_CREATE(mpiInfo, ierror)
     call MPI_INFO_SET(mpiInfo, "romio_ds_write", "disable", ierror)     
  case('marvin')
     mpiiofs = "ufs:"
     call MPI_INFO_CREATE(mpiInfo, ierror)
     call MPI_INFO_SET(mpiInfo, "romio_ds_write", "disable", ierror)
  case('panfs')
     mpiiofs = "panfs:"
     call MPI_INFO_CREATE(mpiInfo, ierror)
     call MPI_INFO_SET(mpiInfo, "panfs_concurrent_write", "1", ierror)
  case('mira')
     mpiiofs = "bglockless:"
     mpiInfo = MPI_INFO_NULL
  case default
     call parallel_kill("parallel_init: unknown MPIIO fs: '"//trim(mpiiofs)//"'")
  end select 

  return
end subroutine parallel_init_io

subroutine parallel_subarray(globalGridSize, localGridSize, gridOffset)

  ! Internal modules
  use parallel

  implicit none

  ! Arguments
  integer :: globalGridSize(3), localGridSize(3), gridOffset(3)

  ! Local variables
  integer :: ierror

  call MPI_Type_create_subarray(3, globalGridSize, localGridSize, gridOffset,                &
       MPI_ORDER_FORTRAN, MPI_REAL_WP, mpiDerivedTypeRealSubarray, ierror)
  call MPI_Type_commit(mpiDerivedTypeRealSubarray, ierror)
  call MPI_Type_create_subarray(3, globalGridSize, localGridSize, gridOffset,                &
       MPI_ORDER_FORTRAN, MPI_INTEGER, mpiDerivedTypeIntegerSubarray, ierror)
  call MPI_Type_commit(mpiDerivedTypeIntegerSubarray, ierror)

  return
end subroutine parallel_subarray

subroutine parallel_time(wtime)

  ! Internal modules
  use parallel

  implicit none

  ! Arguments
  real(WP), intent(out) :: wtime
  
  ! Get time
  wtime = MPI_WTIME()
  
  return
end subroutine parallel_time

subroutine parallel_finalize

  ! Internal modules
  use parallel

  implicit none

  ! Local variables
  integer :: ierror

  if (mpiDerivedTypeIntegerSubarray .ne. MPI_DATATYPE_NULL)                                  &
       call MPI_Type_free(mpiDerivedTypeIntegerSubarray, ierror)
  if (mpiDerivedTypeRealSubarray .ne. MPI_DATATYPE_NULL)                                     &
       call MPI_Type_free(mpiDerivedTypeRealSubarray, ierror)

  mpiDerivedTypeRealSubarray = MPI_DATATYPE_NULL
  mpiDerivedTypeIntegerSubarray = MPI_DATATYPE_NULL

  if (comm .ne. MPI_COMM_NULL) call MPI_Comm_free(comm, ierror)
  comm = MPI_COMM_NULL

  ! Finalize MPI
  call MPI_FINALIZE(ierror)

  return
end subroutine parallel_finalize
