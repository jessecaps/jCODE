module ibm_comm

  ! External modules
  use ibm_marker
  use ibm_exchange

  implicit none

contains

  ! Communicate ibm markers in a given direction
  ! --------------------------------------------
  subroutine communicate_direction(dir)

    ! External modules
    use parallel

    implicit none

    ! Arguments
    integer, intent(in) :: dir

    ! Local variables
    integer :: i, i1, i2, rankLeft, rankRight, nSendLeft, nSendRight, nRecv, ierror
    type(t_Marker), dimension(:), allocatable :: markersToRecv, markersToSendLeft,           &
         markersToSendRight

    ! Initialize MPI variables
    call MPI_Cart_shift(commDir(dir), 0, 1, rankLeft, rankRight, ierror)

    ! Get processor extents
    i1 = iStart(dir)
    i2 = iEnd(dir)

    ! Determine the rank the marker should belong to
    nSendLeft = 0; nSendRight = 0
    do i = 1, nMarkers
       if (marker(i)%gridIndex(dir) .lt. i1) then
          nSendLeft = nSendLeft + 1
       else if (marker(i)%gridIndex(dir) .gt. i2) then
          nSendRight = nSendRight + 1
       end if
    end do

    ! Allocate the send buffers
    allocate(markersToSendLeft(nSendLeft))
    allocate(markersToSendRight(nSendRight))

    ! Copy markers to buffer
    nSendLeft = 0; nSendRight = 0
    do i = 1, nMarkers
       if (marker(i)%gridIndex(dir) .lt. i1) then
          nSendLeft = nSendLeft + 1
          markersToSendLeft(nSendLeft) = marker(i)
          marker(i)%stop = 1 ! ... remove the marker
       else if (marker(i)%gridIndex(dir) .gt. i2) then
          nSendRight = nSendRight + 1
          markersToSendRight(nSendRight) = marker(i)
          marker(i)%stop = 1 ! ... remove the marker
       end if
    end do

     ! <<<<<<<< Communicate LEFT <<<<<<<<<<<<

     ! Send to left, receive from right
     nRecv = 0
     call MPI_SendRecv(nSendLeft, 1, MPI_INTEGER, rankLeft, 0, nRecv, 1, MPI_INTEGER,        &
          rankRight, 0, commDir(dir), MPI_STATUS_IGNORE, ierror)
     allocate(markersToRecv(nRecv))
     call MPI_SendRecv(markersToSendLeft, nSendLeft, MPI_MARKER, rankLeft, 1,                &
          markersToRecv, nRecv, MPI_MARKER, rankRight,1, commDir(dir),                       &
          MPI_STATUS_IGNORE, ierror)
     deallocate(markersToSendLeft)

     ! Add to the marker vector
     if (nRecv .gt. 0) then
        call resize_markers(marker, nMarkers + nRecv)
        marker(nMarkers + 1 : nMarkers + nRecv) = markersToRecv(1:nRecv)
        nMarkers = nMarkers + nRecv
     end if
     deallocate(markersToRecv)

     ! >>>>>>>> Communicate RIGHT >>>>>>>>>

     ! Send to right, receive from left
     nRecv = 0
     call MPI_SendRecv(nSendRight, 1, MPI_INTEGER, rankRight, 2, nRecv, 1, MPI_INTEGER,      &
          rankLeft, 2, commDir(dir), MPI_STATUS_IGNORE, ierror)
     allocate(markersToRecv(nRecv))
     call MPI_SendRecv(markersToSendRight, nSendRight, MPI_MARKER, rankRight, 3,             &
          markersToRecv, nRecv, MPI_MARKER, rankLeft, 3, commDir(dir),                       &
          MPI_STATUS_IGNORE, ierror)
     deallocate(markersToSendRight)

     ! Add to the marker vector
     if (nRecv .gt. 0) then
        call resize_markers(marker, nMarkers + nRecv)
        marker(nMarkers + 1:nMarkers + nRecv) = markersToRecv(1:nRecv)
        nMarkers = nMarkers + nRecv
     end if
     deallocate(markersToRecv)

     ! Recycle
     call recycle_markers(marker)

     ! Update number of markers
     call update_marker_size

    return
  end subroutine communicate_direction


  ! Communicate ghost markers in a given direction
  ! ----------------------------------------------
  subroutine communicate_border(dir)

    ! External modules
    use parallel

    implicit none

    ! Arguments
    integer, intent(in) :: dir

    ! Local variables
    integer :: n, nSend, nRecv, rankLeft, rankRight, tag(3,2), ierror
    integer :: iStatus(MPI_STATUS_SIZE)
    type(t_Marker), dimension(:), allocatable :: markersToSend, markersToRecv

    ! Define unitque tags for send/ecv
    tag(1,1) = 0
    tag(1,2) = 2
    tag(2,1) = 4
    tag(2,2) = 6
    tag(3,1) = 8
    tag(3,2) = 10

    ! <<<<<<<< Communicate LEFT <<<<<<<<<<<<

    ! Count number of markers to send left
    nSend = 0
    do n = 1, nMarkers
       if (marker(n)%gridIndex(dir) .lt. gridOffset(dir) + 1 + nOverlap(dir,1))              &
            nSend = nSend + 1
    end do
    do n = 1, nGhostMarkers
       if (ghostMarker(n)%gridIndex(dir) .lt. gridOffset(dir) + 1 + nOverlap(dir,1))         &
            nSend = nSend + 1
    end do

    ! Allocate send buffer
    allocate(markersToSend(nSend))

    ! Copy markers in buffer
    nSend = 0
    do n = 1, nMarkers
       if (marker(n)%gridIndex(dir) .lt. gridOffset(dir) + 1 + nOverlap(dir,1)) then
          nSend = nSend + 1
          markersToSend(nSend) = marker(n)
          if (isPeriodic(dir) .and. procCoords(dir) .eq. 0) then
             markersToSend(nSend)%position(dir) = markersToSend(nSend)%position(dir) +       &
                  periodicLength(dir)
             markersToSend(nSend)%gridIndex(dir) = markersToSend(nSend)%gridIndex(dir) +     &
                  globalGridSize(dir)
          end if
       end if
    end do
    do n = 1, nGhostMarkers
       if (ghostMarker(n)%gridIndex(dir) .lt. gridOffset(dir) +1 + nOverlap(dir,1)) then
          nSend = nSend + 1
          markersToSend(nSend) = ghostMarker(n)
          if (isPeriodic(dir) .and. procCoords(dir) .eq. 0) then
             markersToSend(nSend)%position(dir) = markersToSend(nSend)%position(dir) +       &
                  periodicLength(dir)
             markersToSend(nSend)%gridIndex(dir) = markersToSend(nSend)%gridIndex(dir) +     &
                  globalGridSize(dir)
          end if
       end if
    end do

    ! Communicate sizes
    nRecv = 0
    call MPI_Cart_shift(commDir(dir), 0, -1, rankRight, rankLeft, ierror)
    call MPI_SendRecv(nSend, 1, MPI_INTEGER, rankLeft, tag(dir,1), nRecv, 1, MPI_INTEGER,    &
         rankRight, tag(dir,1), commDir(dir), MPI_STATUS_IGNORE, ierror)

    ! Allocate recv buffer
    allocate(markersToRecv(nRecv))

    ! Communicate
    call MPI_SendRecv(markersToSend, nSend, MPI_MARKER, rankLeft, 0, markersToRecv,          &
         nRecv, MPI_MARKER, rankRight, 0, commDir(dir), MPI_STATUS_IGNORE, ierror)

    ! Add to ghost markers
    if (nRecv .gt. 0) then
       call resize_markers(ghostMarker, nGhostMarkers+nRecv)
       ghostMarker(nGhostMarkers+1:nGhostMarkers + nRecv) = markersToRecv(1:nRecv)
       nGhostMarkers = nGhostMarkers + nRecv
    end if

    ! Clean up
    deallocate(markersToSend)
    deallocate(markersToRecv)

    ! >>>>>>>> Communicate RIGHT >>>>>>>>>

    ! Count number of markers to send right
    nSend = 0
    do n = 1, nMarkers
       if (marker(n)%gridIndex(dir) .gt. localGridSize(dir) + gridOffset(dir) -              &
            nOverlap(dir,2)) nSend = nSend + 1
    end do
    do n = 1, nGhostMarkers - nRecv
       if (ghostMarker(n)%gridIndex(dir) .gt. localGridSize(dir) + gridOffset(dir) -         &
            nOverlap(dir,2)) nSend = nSend + 1
    end do

    ! Allocate send buffer
    allocate(markersToSend(nSend))

    ! Copy markers in buffer
    nSend = 0
    do n = 1, nMarkers
       if (marker(n)%gridIndex(dir) .gt. localGridSize(dir) + gridOffset(dir) -              &
            nOverlap(dir,2)) then
          nSend = nSend + 1
          markersToSend(nSend) = marker(n)
          if (isPeriodic(dir) .and. procCoords(dir) .eq. nProcsDir(dir)-1) then
             markersToSend(nSend)%position(dir) = markersToSend(nSend)%position(dir) -       &
                  periodicLength(dir)
             markersToSend(nSend)%gridIndex(dir) = markersToSend(nSend)%gridIndex(dir) -     &
                  globalGridSize(dir)
          end if
       end if
    end do
    do n = 1, nGhostMarkers - nRecv
       if (ghostMarker(n)%gridIndex(dir) .gt. localGridSize(dir)  + gridOffset(dir) -        &
            nOverlap(dir,2)) then
          nSend = nSend + 1
          markersToSend(nSend) = ghostMarker(n)
          if (isPeriodic(dir) .and. procCoords(dir) .eq. nProcsDir(dir)-1) then
             markersToSend(nSend)%position(dir) = markersToSend(nSend)%position(dir) -       &
                  periodicLength(dir)
             markersToSend(nSend)%gridIndex(dir) = markersToSend(nSend)%gridIndex(dir) -     &
                  globalGridSize(dir)
          end if
       end if
    end do

    ! Communicate sizes
    nRecv = 0
    call MPI_Cart_shift(commDir(dir), 0, +1, rankLeft, rankRight, ierror)
    call MPI_SendRecv(nSend, 1, MPI_INTEGER, rankRight, tag(dir,2), nRecv, 1, MPI_INTEGER,   &
         rankLeft, tag(dir,2), commDir(dir), iStatus, ierror)

    ! Allocate recv buffer
    allocate(markersToRecv(nRecv))

    ! Communicate
    call MPI_SendRecv(markersToSend, nSend, MPI_MARKER, rankRight, 0, markersToRecv,         &
         nRecv, MPI_MARKER, rankLeft, 0, commDir(dir), iStatus, ierror)

    ! Add to ghost markers
    if (nRecv .gt. 0) then
       call resize_markers(ghostMarker, nGhostMarkers+nRecv)
       ghostMarker(nGhostMarkers+1:nGhostMarkers + nRecv) = markersToRecv(1:nRecv)
       nGhostMarkers = nGhostMarkers + nRecv
    end if

    ! Clean up
    deallocate(markersToSend)
    deallocate(markersToRecv)

    return
  end subroutine communicate_border

  
  ! Check if a marker left the domain
  ! -----------------------------------
  subroutine check_marker_bounds(marker_)

    implicit none

    ! Arguments
    type(t_Marker), intent(inout) :: marker_

    if (marker_%stop .eq. 1) return

    select case (nDimensions)

    case(1)

       if (.not. isPeriodic(1)) then
          if (marker_%position(1).lt.cartesianCoordinates(1,1) .or. marker_%position(1).gt.  &
               cartesianCoordinates(globalGridSize(1),1)) marker_%stop = 1
       end if

    case(2)

       if (.not. isPeriodic(1)) then
          if (marker_%position(1).lt.cartesianCoordinates(1,1) .or. marker_%position(1).gt.  &
               cartesianCoordinates(globalGridSize(1),1)) marker_%stop = 1
       end if

       if (.not. isPeriodic(2)) then
          if (marker_%position(2).lt.cartesianCoordinates(1,2) .or. marker_%position(2).gt.  &
               cartesianCoordinates(globalGridSize(2),2)) marker_%stop = 1
       end if

    case(3)

       if (.not. isPeriodic(1)) then
          if (marker_%position(1).lt.cartesianCoordinates(1,1) .or. marker_%position(1).gt.  &
               cartesianCoordinates(globalGridSize(1),1)) marker_%stop = 1
       end if

       if (.not. isPeriodic(2)) then
          if (marker_%position(2).lt.cartesianCoordinates(1,2) .or. marker_%position(2).gt.  &
               cartesianCoordinates(globalGridSize(2),2)) marker_%stop = 1
       end if

       if (.not. isPeriodic(3)) then
          if (marker_%position(3).lt.cartesianCoordinates(1,3) .or. marker_%position(3).gt.  &
               cartesianCoordinates(globalGridSize(3),3)) marker_%stop = 1
       end if

    end select

    return
  end subroutine check_marker_bounds


  ! Localize a marker to the grid
  ! -----------------------------
  subroutine localize_marker(marker_)

    implicit none

    ! Arguments
    type(t_Marker), intent(inout) :: marker_

    ! Local variables
    integer :: i, ip

    if (marker_%stop .eq. 1) return

    do i = 1, nDimensions
       ip = marker_%gridIndex(i)
       do while (marker_%position(i) - cartesianCoordinates(ip, i) .lt. 0.0_WP)
          ip = ip - 1
       end do
       do while (cartesianCoordinates(ip + 1, i) - marker_%position(i) .le. 0.0_WP)
          ip = ip + 1
       end do
       marker_%gridIndex(i) = ip
    end do

    return
  end subroutine localize_marker


  ! Correct a marker position for periodicity
  ! -----------------------------------------
  subroutine correct_marker_position(marker_)

    implicit none

    ! Arguments
    type(t_Marker), intent(inout) :: marker_

    if (marker_%stop .eq. 1) return

    select case (nDimensions)

    case(1)

       if (isPeriodic(1)) marker_%position(1) = cartesianCoordinates(1,1) +                  &
            modulo(marker_%position(1) - cartesianCoordinates(1,1), periodicLength(1))

    case(2)

       if (isPeriodic(1)) marker_%position(1) = cartesianCoordinates(1,1) +                  &
            modulo(marker_%position(1) - cartesianCoordinates(1,1), periodicLength(1))

       if (isPeriodic(2)) marker_%position(2) = cartesianCoordinates(1,2) +                  &
            modulo(marker_%position(2) - cartesianCoordinates(1,2), periodicLength(2))

    case(3)

       if (isPeriodic(1)) marker_%position(1) = cartesianCoordinates(1,1) +                  &
            modulo(marker_%position(1) - cartesianCoordinates(1,1), periodicLength(1))

       if (isPeriodic(2)) marker_%position(2) = cartesianCoordinates(1,2) +                  &
            modulo(marker_%position(2) - cartesianCoordinates(1,2), periodicLength(2))

       if (isPeriodic(3)) marker_%position(3) = cartesianCoordinates(1,3) +                  &
            modulo(marker_%position(3) - cartesianCoordinates(1,3), periodicLength(3))

    end select

    call localize_marker(marker_)

    return
  end subroutine correct_marker_position

  
  ! Update the marker size
  ! -----------------------
  subroutine update_marker_size

    ! External modules
    use parallel

    implicit none

    ! Update number of Lagrangian markers
    nMarkers = size(marker)
    call parallel_sum(nMarkers, nMarkersGlobal)

    return
  end subroutine update_marker_size

end module ibm_comm


! ================================================= !
! Communicate IBM markers to neighboring processors !
! ================================================= !
subroutine communicate_ibm_markers

  ! Internal modules
  use ibm_comm

  ! External modules
  use parallel

  implicit none

  ! Local variables
  integer :: i

  ! Return if not used
  if (nMarkersGlobal .eq. 0) return

  ! Recycle
  call recycle_markers(marker)

  ! Update number of markers
  call update_marker_size

  select case (nDimensions)

  case (1)

     call communicate_direction(1)

  case (2)

     call communicate_direction(1)
     call communicate_direction(2)

  case (3)

     call communicate_direction(1)
     call communicate_direction(2)
     call communicate_direction(3)

  end select

  ! Update number of markers
  nMarkers = size(marker)
  call parallel_sum(nMarkers, nMarkersGlobal)

  ! Correct marker position
  do i = 1, nMarkers
     call correct_marker_position(marker(i))
  end do

  return
end subroutine communicate_ibm_markers


! ================================================ !
! Communicate ghost markers for computing levelset !
! ================================================ !
subroutine communicate_ghost_markers

  ! Internal modules
  use ibm_comm

  ! External modules
  use parallel

  implicit none

  ! Clean ghost marker structure
  call resize_markers(ghostMarker, 0)
  nGhostMarkers = 0

  select case (nDimensions)

  case (1)

     call communicate_border(1)

  case (2)

     call communicate_border(1)
     call communicate_border(2)

  case (3)

     call communicate_border(1)
     call communicate_border(2)
     call communicate_border(3)

  end select

  ! Update number of markers
  if (allocated(ghostMarker)) nGhostMarkers = size(ghostMarker)

  return
end subroutine communicate_ghost_markers
