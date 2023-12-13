module timing

  ! External modules
  use precision
  use string
  use parallel

  implicit none

  ! Number of timer
  integer :: nTimers

  ! Maximum of timers to be used
  integer, parameter :: nTimersMax = 32
  
  ! Definition of a timer
  type :: t_Timer
     character(len=str_medium) :: name
     real(WP) :: time
     real(WP) :: timeIn
     logical  :: started
  end type t_Timer
  
  ! Array of timers
  type(t_Timer), dimension(:), allocatable, public :: timers

  ! Number and List of timers started
  integer :: nStarted
  integer, dimension(:), allocatable :: started

  ! Timing for full timestep
  real(WP) :: fullTimeIn

  ! Array of values to transfer to monitor
  real(WP), dimension(:), allocatable :: monitorValue

  ! Name of the timing monitor file
  character(len = str_medium) :: monitorName

contains
  
  ! Find a timer of a given name
  ! ----------------------------
  subroutine timing_find_timer(name, iTimer)

    implicit none

    ! Arguments
    character(len = *), intent(in) :: name
    integer, intent(out) :: iTimer

    loop: do iTimer = 1, nTimers
       if (trim(timers(iTimer)%name) .eq. trim(name)) exit loop
    end do loop
    if (iTimer .gt. nTimers) then
       print *, 'Timer : ', name
       call die('timing_find_timer: unknown timer')
    end if
    
    return
  end subroutine timing_find_timer

end module timing


! ============================ !
! Initialize the timing module !
! ============================ !
subroutine timing_setup

  ! Internal modules
  use timing

  implicit none
  
  ! Start with no timers
  nTimers  = 0
  nStarted = 0

  ! Get time for full time step
  fullTimeIn = MPI_WTIME()

  ! Safely deallocate the timers
  if (allocated(timers )) deallocate(timers)
  if (allocated(started)) deallocate(started)
  if (allocated(monitorValue)) deallocate(monitorValue)

  ! Allocate the timers
  allocate(timers (nTimersMax))
  allocate(started(nTimersMax))

  return
end subroutine timing_setup


! ========================= !
! Cleanup the timing module !
! ========================= !
subroutine timing_cleanup

  ! Internal modules
  use timing

  implicit none

  ! Start with no timers
  nTimers  = 0
  nStarted = 0

  if (allocated(timers )) deallocate(timers)
  if (allocated(started)) deallocate(started)
  if (allocated(monitorValue)) deallocate(monitorValue)

  return
end subroutine timing_cleanup


! ================================ !
! Setup the timing monitor routine !
! ================================ !
subroutine timing_setup_monitor(mode, controlIteration)

  ! Internal modules
  use timing

  ! External modules
  use solver_options

  implicit none

  ! Arguments
  integer, intent(in) :: mode, controlIteration

  ! Local variables
  integer :: iTimer

  ! Allocate the array of data to transfer to monitor
  allocate(monitorValue(2 * nTimers + 4))

  ! Output timing based on run status
  select case (mode)

  case (FORWARD)
     monitorName = 'timing'
     if (controlIteration .gt. 0) write(monitorName, '(A,I2.2)')                             &
          'timing_', controlIteration

  case (ADJOINT)
     monitorName = 'adjoint_timing'
     if (controlIteration .gt. 0) write(monitorName, '(A,I2.2)')                             &
          'adjoint_timing_', controlIteration

  end select

  ! Create the file to monitor
  call monitor_create(trim(monitorName), 2 * ntimers + 4)
  call monitor_set_header(1,'Total [s]','r')
  call monitor_set_header(2,'Time/points [us]','r')
  do iTimer = 1, nTimers
     call monitor_set_header(iTimer*2+1,trim(timers(iTimer)%name)//' [s]','r')
     call monitor_set_header(iTimer*2+2,trim(timers(iTimer)%name)//' [%]','r')
  end do
  call monitor_set_header(2*nTimers+3,'Rest [s]','r')
  call monitor_set_header(2*nTimers+4,'Rest [%]','r')

  return
end subroutine timing_setup_monitor


! ========================= !
! Create a new timer object !
! ========================= !
subroutine timing_create(name)

  ! Internal modules
  use timing

  implicit none

  ! Arguments
  character(len=*), intent(in) :: name

  if (.not. allocated(timers)) return

  ! Add a timer
  ntimers = ntimers + 1
  if (nTimers .gt. nTimersMax) call die ('timing_create: too many timers!')

  ! Preset the values
  timers(ntimers)%name = name
  timers(ntimers)%time = 0.0_WP
  timers(ntimers)%timeIn = 0.0_WP
  timers(ntimers)%started = .false.
  
  return
end subroutine timing_create


! =============== !
! Start the timer !
! =============== !
subroutine timing_start(name)

  ! Internal modules
  use timing

  implicit none

  ! Arguments
  character(len=*), intent(in) :: name

  ! Local variables
  integer :: iTimer
  real(WP) :: timeIn

  if (.not. allocated(timers)) return

  call timing_find_timer(name, iTimer)

  if (timers(itimer)%started) then
     print *, 'Timer : ', name
     call die('timing_start: timer already started')
  end if
  
  timeIn = MPI_WTIME()

  timers(iTimer)%timeIn = timeIn
  timers(iTimer)%started = .true.

  if (nStarted .ne. 0) then
     timers(started(nStarted))%time = timers(started(nStarted))%time                         &
          + timeIn - timers(started(nStarted))%timeIn
     timers(started(nStarted))%timeIn = 0.0_WP
  end if
  nStarted = nStarted + 1
  started(nStarted) = iTimer
  
  return
end subroutine timing_start


! ============== !
! Stop the timer !
! ============== !
subroutine timing_stop(name)

  ! Internal modules
  use timing

  implicit none

  ! Arguments
  character(len=*), intent(in) :: name

  ! Local variables
  integer :: iTimer
  real(WP) :: timeOut

  if (.not. allocated(timers)) return

  call timing_find_timer(name, iTimer)

  if (.not.timers(iTimer)%started) then
     print *, 'Timer : ', name
     call die('timing_stop: timer already stopped')
  end if
  
  timeOut = MPI_WTIME()

  timers(iTimer)%time = timers(iTimer)%time + timeOut - timers(iTimer)%timeIn
  timers(iTimer)%timeIn = 0.0_WP
  timers(iTimer)%started = .false.

  nStarted = nStarted - 1  
  if (nStarted .ne. 0) then
     timers(started(nStarted))%timeIn = timeOut
  end if
  
  return
end subroutine timing_stop


! ================== !
! Monitor the timers !
! ================== !
subroutine timing_monitor

  ! Internal modules
  use timing

  ! External modules
  use geometry, only : nGridPoints

  implicit none

  ! Local variables
  integer :: iTimer
  real(WP) :: fullTimeOut, restTime

  ! Get time for full time step
  fullTimeOut = MPI_WTIME()

  ! Create the array of values
  restTime = 0.0_WP
  monitorValue(1) = fullTimeOut - fullTimeIn + epsilon(1.0_WP)
  monitorValue(2) = 1e6_WP * monitorValue(1) * real(nProcs, WP) / real(nGridPoints, WP)
  do iTimer = 1, nTimers
     restTime = restTime + timers(iTimer)%time
     monitorValue(2 * iTimer + 1) = timers(iTimer)%time
     monitorValue(2 * iTimer + 2) = 100.0_WP * timers(iTimer)%time / monitorValue(1)
  end do
  monitorValue(2 * nTimers + 3) = monitorValue(1) - restTime
  monitorValue(2 * nTimers + 4) = 100.0_WP * (monitorValue(1) - restTime) / monitorValue(1)
  
  ! Transfer data to monitor
  call monitor_select(trim(monitorName))
  call monitor_set_array_values(monitorValue)
  
  ! Reset the timers
  fullTimeIn = fullTimeOut
  do iTimer = 1, nTimers
     timers(iTimer)%time = 0.0_WP
     timers(iTimer)%timeIn = 0.0_WP
  end do
  
  return
end subroutine timing_monitor

