module reverse_migrator

  ! External modules
  use precision
  use state, only : update_state, filter_solution, conservedVariables
  use time_integrator, only : nStages
  use solver, only : saveInterval, adjointIterations, solutionFile

  implicit none

  integer :: startTimestep, endTimestep, numIntermediateStates, loadedTimestep
  real(WP), allocatable :: checkpointTimes(:), checkpointBuffer(:,:,:)

end module reverse_migrator


! ========================== !
! Setup the reverse migrator !
! ========================== !
subroutine reverse_migrator_setup(start, end, n)

  ! Internal modules
  use reverse_migrator

  ! External modules
  use string
  use solver_options
  use geometry

  implicit none

  ! Arguments
  integer, intent(in) :: start, end, n

  ! Local variables
  character(len = str_long) :: message

  ! Set the reverse migrator options for uniform checkpointing
  startTimestep = start
  endTimestep = end
  numIntermediateStates = min(n, saveInterval * nStages)
  loadedTimestep = -1

  if (numIntermediateStates .ne. saveInterval * nStages) then
     write(message, '(2(A,I0.0),A)') 'Uniform checkpointing requires storage for ',          &
          saveInterval * nStages, ' intermediate states, but only ',                         &
          numIntermediateStates, ' are available!'
     call die(trim(message))
  end if

  allocate(checkpointBuffer(nGridPoints, nUnknowns, numIntermediateStates))
  allocate(checkpointTimes(numIntermediateStates))

  return
end subroutine reverse_migrator_setup


! ============================ !
! Cleanup the reverse migrator !
! ============================ !
subroutine reverse_migrator_cleanup

  ! Internal modules
  use reverse_migrator

  implicit none

  if (allocated(checkpointTimes)) deallocate(checkpointTimes)
  if (allocated(checkpointBuffer)) deallocate(checkpointBuffer)

  loadedTimestep = -1

  return
end subroutine reverse_migrator_cleanup


! ======================== !
! Run the reverse migrator !
! ======================== !
subroutine reverse_migrator_run(currentTimestep, stage)

  ! Internal modules
  use reverse_migrator

  ! External modules
  use string
  use simulation_flags
  use time_info
  use state_io

  implicit none

  ! Arguments
  integer, intent(in) :: currentTimestep, stage

  ! Local variables
  integer :: i, savedTimestep, timestep_, stage_
  character(len = str_medium) :: filename
  real(WP) :: savedTime

  ! Start the timer
  call timing_start('checkpoint')

  savedTime = time
  savedTimestep = timestep

  ! Check if requested sub-step is in memory
  if (loadedTimestep .eq. -1 .or.                                                            &
       currentTimestep .lt. loadedTimestep .or.                                              &
       currentTimestep .gt. loadedTimestep + saveInterval .or.                               &
       (currentTimestep .eq. loadedTimestep .and. stage .lt. nStages) .or.                   &
       (currentTimestep .eq. loadedTimestep + saveInterval .and. stage .eq. nStages)) then

     ! Find the last timestep at which solution is available on disk
     if (mod(currentTimestep, saveInterval) .eq. 0 .and. stage .eq. nStages) then
        timestep_ = currentTimestep
     else if (mod(currentTimestep, saveInterval) .eq. 0) then
        timestep_ = currentTimestep - saveInterval
     else
        timestep_ = currentTimestep - mod(currentTimestep, saveInterval)
     end if

     ! Load the saved solution
     write(filename, '(2A,I8.8)') trim(solutionFile), ".", timestep_
     call simulation_read(IO_FORWARD_STATE, trim(filename))
     loadedTimestep = timestep_ !... register timestep at which solution was loaded

     ! Update the first entry in memory
     i = 1
     checkpointBuffer(:,:,i) = conservedVariables
     checkpointTimes(i) = time

     if (loadedTimestep .ne. endTimestep) then

        ! Update the state
        call update_state

        do timestep_ = loadedTimestep + 1, loadedTimestep + saveInterval

           timestep = timestep_
           call get_timestep_size

           do stage_ = 1, nStages

              ! Take a single sub-step using the time integrator
              call substep_forward(stage_)

              ! Update the state
              call update_state

              if (timestep_ .eq. loadedTimestep + saveInterval .and. stage_ .eq. nStages) exit

              i = i + 1
              checkpointBuffer(:,:,i) = conservedVariables
              checkpointTimes(i) = time

           end do

           ! Filter solution if required
           call filter_solution(conservedVariables, timestep_)

        end do

     end if

  end if

  ! Sub-step requested is available in memory
  i = (currentTimestep - 1 - loadedTimestep) * nStages + stage + 1
  conservedVariables = checkpointBuffer(:,:,i)
  time = checkpointTimes(i)
  call update_state

  time = savedTime
  timestep = savedTimestep

  ! Stop the timer
  call timing_stop('checkpoint')

  return
end subroutine reverse_migrator_run
