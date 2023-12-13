module monitor

  use precision
  use string

  implicit none
  public

  type :: t_Monitor
     integer :: nColumns, iunit
     character(len = str_medium) :: filename
     character(len = str_medium), dimension(:), allocatable :: header
     character(len = 1), dimension(:), allocatable :: columnType
     real(WP), dimension(:), allocatable :: val
  end type t_Monitor

  type(t_Monitor), allocatable, public :: monitorFiles(:)
  integer, public :: fileLog = -1, iFile, nFiles
  integer, parameter, public :: maxFiles = 32
  integer, parameter, public :: columnLength = 14
  integer :: reportInterval
  character(len = *), parameter, public :: f1 = '(a12)'
  character(len = *), parameter, public :: f2 = '(i12)'
  character(len = *), parameter, public :: f3 = '(ES12.5)'

contains

  ! Print some values to the screen
  ! -------------------------------
  subroutine monitor_to_screen(mode, controlIteration)

    ! External modules
    use parallel
    use string
    use simulation_flags
    use solver_options
    use state
    use time_info
    use functional, only : instantaneousCostFunctional
    use controller, only : instantaneousSensitivity

    implicit none

    ! Arguments
    integer, intent(in) :: mode, controlIteration

    ! Local variables
    character(len = str_long) :: str1, str2

    ! Only the root process does something
    if (iRank .ne. iRoot) return

    if (useConstantCfl) then
       write(str1, '(I8,2(A,ES12.5))') timestep, ', time = ', abs(time),                     &
            ', dt = ', timeStepSize
    else
       write(str1, '(I8,2(A,ES12.5))') timestep, ', time = ', abs(time),                     &
            ', CFL = ', cfl
    end if

    if (.not. predictionOnly) then

       select case (mode)
       case (FORWARD)
          if (controlIteration .gt. 0) then
             write(str2, '(A,1ES12.5,A,I2)') ', cost = ', instantaneousCostFunctional,       &
                  ', iteration = ',controlIteration
          else
             write(str2, '(A,1ES12.5)') ', cost = ', instantaneousCostFunctional
          end if
       case (ADJOINT)
          write(str2, '(A,1ES12.5)') ', gradient = ', sum(instantaneousSensitivity**2)
       end select

       str1 = trim(str1) // trim(str2)

    end if

    write(*, '(A)') trim(str1)

    return
  end subroutine monitor_to_screen


  ! Monitor new solver run
  ! ----------------------
  subroutine monitor_new_run(mode)

    ! External modules
    use parallel
    use solver_options 

    implicit none

    ! Arguments
    integer, intent(in) :: mode

    if (iRank .eq. iRoot) then

       write(*, *)

       select case (mode)

       case (FORWARD)

          write(*, '(A)') "New forward run:"
          write(*, '(A)') repeat('-', 32)

       case (ADJOINT)

          write(*, '(A)') "New adjoint run:"
          write(*, '(A)') repeat('-', 32)

       end select

       write(*, *)

    end if

    return
  end subroutine monitor_new_run

end module monitor


! ======================== !
! Setup the monitor module !
! ======================== !
subroutine monitor_setup

  ! Internal modules
  use monitor

  ! External modules
  use parser
  use parallel
  use fileio

  implicit none

  ! Local variables
  integer :: ierror

  ! Set number of files to zero
  nFiles = 0

  ! But preallocate the number of files to some value
  allocate(monitorFiles(maxFiles))

  ! Only the root process does the following
  if (irank .eq. iroot) then

     ! Create the directory
     call CREATE_FOLDER("monitor")

     ! Open log file
     fileLog = iOpen()
     open(fileLog, file = 'monitor/log', form = 'formatted', iostat = ierror,                &
          status = 'REPLACE')

  end if

  ! Read in the report interval
  call parser_read('report interval', reportInterval, 1)
  if (reportInterval .eq. 0) reportInterval = -1

  return
end subroutine monitor_setup


! ======================== !
! Reset the monitor module !
! ======================== !
subroutine monitor_reset

  ! Internal modules
  use monitor

  ! External modules
  use parallel
  use fileio

  implicit none

  ! Local variables
  integer :: i

  if (irank .eq. iroot .and. allocated(monitorFiles)) then
     ! Close the files
     do i = 1, nfiles
       close(iClose(monitorFiles(i)%iunit))
     end do
  end if

  nFiles = 0

  if (allocated(monitorFiles)) deallocate(monitorFiles)
  allocate(monitorFiles(maxFiles))

  return
end subroutine monitor_reset


! ============================================== !
! Create a new file to monitor at each iteration !
! ============================================== !
subroutine monitor_create(filename, nColumns)

  ! Internal modules
  use monitor

  implicit none

  ! Arguments
  integer, intent(in) :: nColumns
  character(len = *), intent(in) :: filename

  ! Add a file
  nFiles = nFiles + 1
  if (nFiles .gt. maxFiles) call die('monitor_create: too many files to monitor')

  iFile  = nFiles

  ! Preset the values
  monitorFiles(iFile)%filename  = trim(adjustl(filename))
  monitorFiles(iFile)%nColumns  = nColumns

  ! Allocate the arrays
  allocate(monitorFiles(ifile)%header(nColumns))
  allocate(monitorFiles(ifile)%columnType(nColumns))
  allocate(monitorFiles(ifile)%val(nColumns))

  return
end subroutine monitor_create

subroutine monitor_select(filename)

  ! Internal modules
  use monitor

  implicit none

  ! Arguments
  character(len = *), intent(in) :: filename

  ! Find the file.
  loop: do iFile = 1, nFiles
     if (trim(monitorFiles(ifile)%filename) .eq. trim(filename)) exit loop
  end do loop
  if (iFile .gt. nFiles) call die('monitor_select: unknown file to monitor: ' // trim(filename))

  return
end subroutine monitor_select


! ========================== !
! Set the header of the file !
! ========================== !
subroutine monitor_set_header(icol, header, columnType)

  ! Internal modules
  use monitor

  implicit none

  ! Arguments
  integer, intent(in) :: icol
  character(len = *), intent(in) :: header
  character(len = 1), intent(in) :: columnType

  ! Test if there are enough columns
  if (icol .gt. monitorFiles(iFile)%nColumns)                                                &
       call die('monitor_set_header: column index too large for file: ' //                   &
       trim(monitorFiles(iFile)%filename))

  ! Set the header
  monitorFiles(iFile)%header(icol)     = trim(header)
  monitorFiles(iFile)%columnType(icol) = columnType

  return
end subroutine monitor_set_header


! ================================== !
! Set the values to dump in the file !
! ================================== !
subroutine monitor_set_single_value(int, val)

  ! Internal modules
  use monitor

  ! External modules
  use precision

  implicit none

  ! Arguments
  integer, intent(in) :: int
  real(WP), intent(in) :: val

  ! Local variables
  real(WP) :: trimVal

  if (int .gt. monitorFiles(iFile)%nColumns) then
     call die('set_monitor_single_value: too  many values')
  else
     if (abs(val) .lt. 1E-99_WP) then
        trimVal = 0.0_WP
     else
        trimVal = val
     end if
     monitorFiles(iFile)%val(int) = trimVal
  end if

  return
end subroutine monitor_set_single_value


! ================================== !
! Set the values to dump in the file !
! ================================== !
subroutine monitor_set_array_values(val)

  ! Internal modules
  use monitor

  ! External modules
  use precision

  implicit none

  ! Arguments
  real(WP), dimension(monitorFiles(iFile)%nColumns) :: val

  monitorFiles(iFile)%val = val

  return
end subroutine monitor_set_array_values


! ============================= !
! Monitor some additional stuff !
! ============================= !
subroutine monitor_log(text)

  ! Internal modules
  use monitor

  ! External modules
  use parallel

  implicit none

  ! Arguments
  character(len = *), intent(in) :: text

  ! Local variables
  integer, dimension(8) :: date

  if (irank .eq. iroot .and. fileLog .ge. 0) then
     call date_and_time(VALUES = date)
     write(fileLog, 100) date(2), date(3), date(1), date(5), date(6), date(7), text
     flush(fileLog)
  end if

100 format(i2.2,'/',i2.2,'/',i4.4,' @ ',i2.2,':',i2.2,':',i2.2,a59)

  return
end subroutine monitor_log


! ============================= !
! Dump the headers of the files !
! ============================= !
subroutine monitor_dump_header

  ! Internal modules
  use monitor

  ! External modules
  use parallel
  use string
  use fileio

  implicit none

  ! Local variables
  integer :: i, j, offset, ierror
  character(len = columnLength) :: col
  character(len = str_long) :: header
  character(len = str_medium) :: filename, buffer
  logical :: twoLines
  integer :: index1

  ! Only the root process does something
  if (irank .ne. iroot) return

  do i = 1, nFiles

     twoLines = .false.
     filename = 'monitor/'//trim(monitorFiles(i)%filename)

     ! Open the file
     monitorFiles(i)%iunit = iOpen()
     open(monitorFiles(i)%iunit, file = filename, form = "formatted",                        &
          iostat = ierror, status = "REPLACE")

     ! Create the header
     write(header(1 + 0 * columnLength:), f1) 'Step'
     write(header(1 + 1 * columnLength:), f1) 'Time'
     offset = 2 * columnLength

     ! Extract the first line
     do j = 1, monitorFiles(i)%nColumns
        read (monitorFiles(i)%header(j), '(a)') buffer
        index1 = index(trim(buffer), ' ')
        if (index1 .ne. 0 .and. index1 .lt. columnLength - 1) then
           twoLines = .true.
           read(buffer(1 : index1), f1) col
        else
           read(buffer, f1) col
        end if
        write(header(offset + 1 + (j - 1) * columnLength:), f1) trim(col)
     end do

     ! Dump the header
     write(monitorFiles(i)%iunit, '(a)') trim(header)

     ! Dump second line if necessary
     if (twoLines) then
        header = ''
        do j = 1, monitorFiles(i)%nColumns
           read (monitorFiles(i)%header(j), '(a)') buffer
           index1 = index(trim(buffer), ' ')
           if (index1 .ne. 0 .and. index1 .lt. columnLength - 1) then
              read(buffer(index1:), f1) col
           else
              col = ''
           end if
           write(header(offset + 1 + (j - 1) * columnLength:), f1) trim(col)
        end do
        write(monitorFiles(i)%iunit,'(a)') trim(header)
     end if

  end do

  return
end subroutine monitor_dump_header


! ============================================= !
! Dump the values to the files at each timestep !
! ============================================= !
subroutine monitor_dump_values

  ! Internal modules
  use monitor

  ! External modules
  use parallel
  use string
  use time_info

  implicit none

  ! Local variables
  integer :: i, j, offset
  character(len = str_long) :: line

  ! Only the root process writes
  if (irank .ne. iroot) return

  do i = 1, nFiles

     ! Create the line to dump
     write (line(1 + 0 * columnLength:), f2) timestep
     write (line(1 + 1 * columnLength:), f3) time
     offset = 2 * columnLength

     do j = 1, monitorFiles(i)%nColumns
        select case(monitorFiles(i)%columnType(j))
        case ('i')
           write (line(offset + 1 + (j - 1) * columnLength:), f2) int(monitorFiles(i)%val(j))
        case ('r')
           write (line(offset + 1 + (j - 1) * columnLength:), f3) monitorFiles(i)%val(j) 
        end select
     end do

     ! Dump the header
     write(monitorFiles(i)%iunit,'(a)') trim(line)

     ! Flush
     flush(monitorFiles(i)%iunit)

  end do

  return
end subroutine monitor_dump_values


! ======================= !
! Finalize the monitoring !
! ======================= !
subroutine monitor_finalize

  ! Internal modules
  use monitor

  ! External modules
  use parallel
  use fileio

  implicit none

  ! Local variables
  integer :: i

  if (.not. allocated(monitorFiles)) return

  if (irank .eq. iroot) then
     ! Close the files
     do i = 1, nfiles
        close(iClose(monitorFiles(i)%iunit))
     end do
     ! Close the log file
     close(iClose(fileLog))
  end if

  if (allocated(monitorFiles)) deallocate(monitorFiles)

  return
end subroutine monitor_finalize


! ============ !
! Welcome text !
! ============ !
subroutine monitor_welcome

  ! Internal modules
  use monitor

  ! External modules
  use parallel

  implicit none

  if (iRank .eq. iRoot) then
     write(*, '(A)') "jCODE"
     write(*, '(A)') repeat('-', 32)
     if (nProcs .gt. 1) then
        write(*, '(I4, A)') nProcs, " processes reporting for duty"
     else
        write(*, '(I4, A)') nProcs, " process reporting for duty"
     end if
     write(*, *)
  end if

  return
end subroutine monitor_welcome


! ============================================ !
! Output useful grid information to the screen !
! ============================================ !
subroutine monitor_grid_diagnostics

  ! Internal modules
  use monitor

  ! External modules
  use parallel
  use geometry
  use grid_functions

  implicit none

  ! Local variables
  integer :: iGlobal, jGlobal, kGlobal
  real(WP) :: buf

  if (iRank .eq. iRoot) then
     write(*, '(A,I0.0,A)') "Simulation is ", nDimensions, "D"
     write(*, '(A)') repeat('-', 32)
     write(*, '(2X,A,3(A,I4),A)') "Grid size ", ": ",                                       &
          globalGridSize(1), " x ", globalGridSize(2), " x ", globalGridSize(3), " points"
  end if

  call find_minimum(minval(gridSpacing, dim = 2), buf, iGlobal, jGlobal, kGlobal)
  if (iRank .eq. iRoot) write(*, '(4X,A,(SS,ES9.2E2),3(A,I4),A)') "min grid spacing = ",     &
       buf, " at (", iGlobal, ", ", jGlobal, ", ", kGlobal, ")"

  call find_maximum(maxval(gridSpacing, dim = 2), buf, iGlobal, jGlobal, kGlobal)
  if (iRank .eq. iRoot) then
     write(*, '(4X,A,(SS,ES9.2E2),3(A,I4),A)') "max grid spacing = ",                        &
          buf, " at (", iGlobal, ", ", jGlobal, ", ", kGlobal, ")"
     write(*, *)
  end if

  ! Report processor decomposition
  call parallel_min(localGridSize(1) / nProcsDir(1), iGlobal)
  call parallel_min(localGridSize(2) / nProcsDir(2), jGlobal)
  call parallel_min(localGridSize(3) / nProcsDir(3), kGlobal)
  if (iRank .eq. iRoot) then
     write(*, '(4X,3(A,I4),A)') "number of procs = (",                                       &
          nProcsDir(1), ", ", nProcsDir(2), ", ", nProcsDir(3), ")"
          write(*, '(4X,3(A,I4),A)') "min points/proc = (",                                  &
               iGlobal, ", ", jGlobal, ", ", kGlobal, ")"
     write(*, *)
  end if

  return
end subroutine monitor_grid_diagnostics


! ======================================= !
! Reset / setup the monitor for a new run !
! ======================================= !
subroutine monitor_setup_new_run(mode, controlIteration)

  ! Internal modules
  use monitor

  ! External modules
  use simulation_flags
  use solver_options
  use dissipation
  use filter
  use combustion

  implicit none

  ! Arguments
  integer, intent(in) :: mode, controlIteration

  ! Write to the screen
  call monitor_new_run(mode)

  ! Reset the monitor in case it was previouly used
  call monitor_reset

  ! Reset the timing module
  call timing_setup

  ! Create the timers
  call timing_create('i/o')
  call timing_create('fluxes')
  call timing_create('boundary')
  call timing_create('sources')
  call timing_create('update')
  call timing_create('monitor')
  if (useDissipation) call timing_create('dissipation')
  if (useFilter)  call timing_create('filter')
  if (nSpecies.gt.0 .and. nReactions.gt.0) call timing_create('combustion')
  if (useParticles) call timing_create('particles')
  if (useIBM) call timing_create('ibm')
  if (.not. predictionOnly) then
     call timing_create('functional')
     call timing_create('controller')
  end if
  if (mode .eq. ADJOINT) call timing_create('checkpoint')

  ! Setup the sub modules
  call timing_setup_monitor(mode, controlIteration)
  call monitor_state_setup(mode, controlIteration)
  call monitor_time_info_setup(mode, controlIteration)
  call monitor_particles_setup(mode, controlIteration)
  call monitor_ibm_setup(mode, controlIteration)
  call monitor_adjoint_setup(mode, controlIteration)
  call monitor_cases_setup(mode, controlIteration)

  ! Dump the file headers
  call monitor_dump_header

  return
end subroutine monitor_setup_new_run


! ================== !
! Monitor a timestep !
! ================== !
subroutine monitor_timestep(mode, controlIteration)

  ! Internal modules
  use monitor

  ! External modules
  use solver_options
  use time_info

  implicit none

  ! Arguments
  integer, intent(in) :: mode, controlIteration

  if (reportInterval .le. 0 .or. mod(timestep, max(1, reportInterval)) .ne. 0) return

  ! Monitor the timing per timestep
  call timing_monitor

  ! Start the monitor timing
  call timing_start('monitor')

  ! Monitor the sub-modules
  select case (mode)

  case (FORWARD)

     call monitor_state_forward(controlIteration)
     call monitor_particles_forward(controlIteration)
     call monitor_ibm_forward(controlIteration)
     call monitor_time_info_forward(controlIteration)
     call monitor_adjoint_functional(controlIteration)
     call monitor_cases_timestep(controlIteration)

  case (ADJOINT)

     call monitor_adjoint_sensitivity(controlIteration)

  end select

  ! Dump the data
  call monitor_dump_values

  ! Print data to screen
  call monitor_to_screen(mode, controlIteration)

  ! Stop the monitor timing
  call timing_stop('monitor')

  return
end subroutine monitor_timestep
