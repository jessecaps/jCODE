module ibm_marker

  ! External modules
  use precision

  implicit none

  ! Define the Lagrangian marker type
  type :: t_Marker
     integer(KIND=8) :: id
     real(WP) :: area
     real(WP) :: temperature
     real(WP), dimension(3) :: position
     real(WP), dimension(3) :: velocity
     real(WP), dimension(3) :: normalVector
     integer :: gridIndex(3)
     integer :: stop
  end type t_Marker
  type(t_Marker), dimension(:), allocatable :: marker, ghostMarker

  ! Lagrangian marker parameters
  integer :: nMarkers = 0, nMarkersGlobal = 0, nGhostMarkers = 0
  integer :: MPI_MARKER, SIZE_MPI_MARKER, markerOffset

contains

  ! Send `stopped` Lagrangian markers to the end and resize
  ! -------------------------------------------------------
  subroutine recycle_markers(markerVector)

    implicit none

    ! Arguments
    type(t_Marker), allocatable, intent(inout) :: markerVector(:)

    ! Local variables
    integer :: i, newSize

    ! Compact real particles at the beginning of the array
    newSize = 0
    if (allocated(markerVector)) then
       do i = 1, size(markerVector)
          if (markerVector(i)%stop .eq. 0) then
             newSize = newSize + 1
             if (i .ne. newSize) then
                markerVector(newSize) = markerVector(i)
                markerVector(i)%stop = 1
             end if
          end if
       end do
    end if

    ! Resize Lagrangian markers
    call resize_markers(markerVector, newSize)

    return
  end subroutine recycle_markers


  ! Resize the Lagrangian markers based on new size `n`
  ! ---------------------------------------------------
  subroutine resize_markers(markerVector, n)

    implicit none

    ! Arguments
    type(t_Marker), allocatable, intent(inout) :: markerVector(:)
    integer, intent(in) :: n

    ! Local variables
    integer :: i, nMarkersOld
    type(t_Marker), allocatable :: tempMarker(:)

    ! Resize marker array to size n
    if (.not. allocated(markerVector)) then
       ! Marker is of size 0
       if (n .eq. 0) then
          ! Nothing to do, that's what we want
       else
          ! Allocate directly of size n
          allocate(markerVector(n))
          markerVector(1:n)%stop = 1
       end if
    else if (n .eq. 0) then
       ! Empty the marker vector
       deallocate(markerVector); allocate(markerVector(0))
    else
       ! Update non zero size to another non zero size
       nMarkersOld = size(markerVector)
       if (n .gt. nMarkersOld) then
          ! Increase from `nMarkersOld` to `n`
          allocate(tempMarker(n))
          do i = 1, nMarkersOld
             tempMarker(i) = markerVector(i)
          end do
          deallocate(markerVector)
          allocate(markerVector(size(tempMarker)))
          markerVector = tempMarker
          markerVector(nMarkersOld + 1:n)%stop = 1
       else if (n .lt. nMarkersOld) then
          ! Decrease from `nMarkersOld` to `n`
          allocate(tempMarker(n))
          do i = 1, n
             tempMarker(i) = markerVector(i)
          end do
          deallocate(markerVector)
          allocate(markerVector(size(tempMarker)))
          markerVector = tempMarker
       end if
    end if

    return
  end subroutine resize_markers

end module ibm_marker


! ========================================== !
! Setup MPI marker for communication and i/o !
! ========================================== !
subroutine prepare_mpi_marker

  ! Internal modules
  use ibm_marker

  ! External modules
  use parallel

  implicit none

  ! Local variables
  integer, dimension(17) :: types, lengths, displacement ! ... size based on t_Marker
  integer :: ierror

  ! Create the MPI structure to send Lagrangian markers
  types( 1: 1) = MPI_INTEGER8
  types( 2:12) = MPI_REAL_WP
  types(13:16) = MPI_INTEGER4
  types(17)    = MPI_UB
  lengths(:)   = 1

  ! Hard-code marker type displacement here
  displacement( 1) = 0                     ! ... ID
  displacement( 2) = displacement( 1) + 8  ! ... Area
  displacement( 3) = displacement( 2) + WP ! ... Temperature
  displacement( 4) = displacement( 3) + WP ! ... Position(1)
  displacement( 5) = displacement( 4) + WP ! ... Position(2)
  displacement( 6) = displacement( 5) + WP ! ... Position(3)
  displacement( 7) = displacement( 6) + WP ! ... Velocity(1)
  displacement( 8) = displacement( 7) + WP ! ... Velocity(2)
  displacement( 9) = displacement( 8) + WP ! ... Velocity(3)
  displacement(10) = displacement( 9) + WP ! ... NormalVector(1)
  displacement(11) = displacement(10) + WP ! ... NormalVector(2)
  displacement(12) = displacement(11) + WP ! ... NormalVector(3)
  displacement(13) = displacement(12) + WP ! ... GridIndex(1)
  displacement(14) = displacement(13) + 4  ! ... GridIndex(2)
  displacement(15) = displacement(14) + 4  ! ... GridIndex(3)
  displacement(16) = displacement(15) + 4  ! ... Stop
  displacement(17) = displacement(16) + 4  ! ... UpperBound

  ! Finalize by creating and commiting the new type
  call MPI_Type_struct(17, lengths, displacement, types, MPI_MARKER, ierror)
  call MPI_Type_commit(MPI_MARKER, ierror)

  ! If problem, say it
  if (ierror .ne. 0) call die('Problem with MPI_MARKER!')

  ! Get the size of this type
  call MPI_type_size(MPI_MARKER, SIZE_MPI_MARKER, ierror)

  return
end subroutine prepare_mpi_marker
