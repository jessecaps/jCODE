module fileio

  implicit none
  
  ! File index
  integer :: fileIndex
  integer, dimension(128) :: iunits
  
contains
  
  ! ====================== !
  ! File index management: !
  !   - open a file        !
  !   - add it to the list !
  ! ====================== !
  integer function iopen()
    implicit none

    ! Local variables
    integer, save :: icall = 1
    integer :: i

    if (icall .eq. 1) then
       fileIndex = 1
       icall = 0
    end if
    iunits(fileIndex) = 0
    do i=1,fileIndex
       if (iunits(i) .eq. 0) exit
    end do
    if (i .eq. fileIndex) then
       fileIndex = fileIndex + 1
       if (fileIndex .ge. 128) stop "iopen: maximum units number exceeded"
    end if
    iunits(i) = 1
    iopen = i + 10
    return
  end function iopen
  
  ! ======================= !
  ! File index management:  !
  !   - close a file        !
  !   - remove it from list !
  ! ======================= !
  integer function iclose(iu)
    implicit none

    ! Arguments
    integer :: iu

    iu = iu - 10
    if (iu .gt. 0 .and. iu .lt. fileIndex) then
       iunits(iu) = 0
       iclose = iu + 10
    else
       iclose = -1
    end if
    return
  end function iclose
  
end module fileio
