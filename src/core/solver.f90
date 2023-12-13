module solver

  ! External modules
  use precision
  use string
  use simulation_flags
  use solver_options
  use grid
  use state
  use time_info
  use time_integrator

  implicit none

  integer :: startTimestep, nTimesteps, saveInterval, adjointIterations
  real(WP) :: finalTime, maxWallTime, initWallTime
  real(WP), dimension(2) :: densityRange, temperatureRange, massFractionRange
  logical :: outputOn = .true., dataOverwrite, enableSolutionLimits, interpFluidToParticle
  character(len = str_medium) :: solutionFile, adjointFile, partFile, ibmFile

contains

  ! Forward (predictive) run - returns a cost functional
  ! -----------------------------------------------------
  function run_forward(actuationAmount, controlIteration) result(costFunctional)

    ! External modules
    use state_io
    use functional

    implicit none

    ! Arguments
    integer, intent(in), optional :: controlIteration
    real(WP), intent(in), optional :: actuationAmount

    ! Result
    real(WP) :: costFunctional

    ! Local variables
    character(len = str_medium) :: filename
    integer :: controlIteration_
    real(WP) :: actuationAmount_

    ! Initialize the cost functionals
    costFunctional = 0.0_WP
    currentCostFunctional = 0.0_WP
    auxilaryCostFunctional = 0.0_WP
    instantaneousCostFunctional = 0.0_WP

    ! Set the control iteration
    controlIteration_ = 0
    if (present(controlIteration)) controlIteration_ = controlIteration

    ! Set the actuation amount
    actuationAmount_ = 0.0_WP
    if (present(actuationAmount) .and. .not. predictionOnly)                                 &
         actuationAmount_ = actuationAmount

    ! Setup the monitor and timing for the new run
    call monitor_setup_new_run(FORWARD, controlIteration_)

    ! Load the initial condition
    call load_initial_condition(FORWARD)

    ! Overwrite the initial conditions using the adjoint gradient if t=0
    if (.not. predictionOnly .and. timestep.eq.0)                                            &
         call controller_update_ic(actuationAmount_, controlIteration_)

    ! Set the time
    startTimestep = timestep

    ! Save the initial condition with the appropriate filename
    if (.not. predictionOnly) then
       write(filename, '(2A,I8.8)') trim(solutionFile), ".", timestep
       if (trim(filename) .ne. trim(solutionFile)) then
          call simulation_write(IO_FORWARD_STATE, trim(filename))
       end if
    end if

    ! Call controller hooks before time marching starts
    if (.not. predictionOnly .and. abs(actuationAmount_) .gt. 0.0_WP)                        &
         call controller_hook_before_timemarch(FORWARD)

    ! Update the state
    call update_state

    ! Initialize the CFL & timestep
    call get_timestep_size
    call get_cfl

    ! Setup the data dumping for postprocessing
    if (outputOn) then
       call dump_cleanup
       call dump_setup(FORWARD)
       call dump_result(FORWARD)
       call stat_cleanup
       call stat_setup
    end if
    
    ! Monitor the initial data
    call monitor_timestep(FORWARD, controlIteration_)

    ! Begin time marching
    do while (.not. done())

       ! Update the current time step
       timestep = timestep + 1
       call get_timestep_size(timeStepSize)
       timeStepSize = min(timeStepSize, finalTime - time + epsilon(1.0_WP))

       do timeStage = 1, nStages

          ! Check if physical quantities are within allowed limits
          call check_solution_limits(FORWARD)

          ! Update control parameters
          if (.not. predictionOnly .and. abs(actuationAmount_) .gt. 0.0_WP)                  &
               call controller_update_parameters(actuationAmount_)

          ! Take a single sub-step using the time integrator
          call substep_forward(timeStage)

          ! Update the state
          call update_state

          ! Compute the cost functional
          if (.not. predictionOnly) then
             call functional_compute
             currentCostFunctional = currentCostFunctional + timeNorm(timeStage) *           &
                  timeStepSize * instantaneousCostFunctional
          end if

       end do

       ! Compute cost functional at terminal condition
       if (.not. predictionOnly .and. timestep.eq.startTimestep+nTimesteps) then
          instantaneousCostFunctional = 0.0_WP
          call functional_compute_tc
          currentCostFunctional = currentCostFunctional + instantaneousCostFunctional *      &
               timeStepSize / 6.0_WP ! ... RK4 only
       end if

       ! Filter solution if required
       call filter_solution(conservedVariables, timestep)

       ! Write restart files
       call write_restart_files(FORWARD, controlIteration_)

       ! Dump results for visualization
       if (outputOn) then
          call dump_result(FORWARD)
          call dump_statistics
       end if

       ! Monitor the current timestep
       call monitor_timestep(FORWARD, controlIteration_)

    end do

    ! Call controller hooks after time marching ends
    if (predictionOnly .and. abs(actuationAmount_) .gt. 0.0_WP)                              &
         call controller_hook_after_timemarch(FORWARD)

    ! Cleanup the timer
    call timing_cleanup

    ! Set the cost functional
    if (.not. predictionOnly) costFunctional = currentCostFunctional

    return
  end function run_forward


  ! Adjoint (inverse) run - returns a cost sensitivity gradient
  ! -----------------------------------------------------------
  function run_adjoint(controlIteration) result(costSensitivity)

    ! External modules
    use state_io
    use functional
    use controller

    implicit none

    ! Arguments
    integer, intent(in), optional :: controlIteration

    ! Result
    real(WP), dimension(:), allocatable :: costSensitivity

    ! Local variables
    character(len = str_medium) :: filename
    integer :: i, j, startTimestep, timemarchDirection, controlIteration_

    ! Set the control iteration and adjoint correction
    controlIteration_ = 0
    if (present(controlIteration)) controlIteration_ = controlIteration
    adjointCorrection = 1.0_WP

    ! Setup the monitor and timing for the new run
    call monitor_setup_new_run(ADJOINT, controlIteration_)

    ! Load the initial condition
    call load_initial_condition(FORWARD) !... for control horizon end timestep

    ! Load the adjoint coefficients corresponding to the end of the control time horizon
    write(filename, '(2A,I8.8)') trim(solutionFile), ".", timestep + adjointIterations
    call simulation_read(IO_FORWARD_STATE, trim(filename))

    ! Update the state
    call update_state

    ! Set the time
    startTimestep = timestep

    ! Setup the revese-time migrator
    call reverse_migrator_setup(startTimestep - adjointIterations, startTimestep,            &
         saveInterval * nStages)
    timemarchDirection = -1

    ! Adjoint initial condition (if specified)
    call load_initial_condition(ADJOINT)

    write(filename, '(2A,I8.8)') trim(adjointFile), ".", timestep
    call simulation_write(IO_ADJOINT_STATE, trim(filename))

    ! Call controller hooks before time marching starts
    call controller_hook_before_timemarch(ADJOINT)

    ! Initialize the cost sensitivity
    allocate(costSensitivity(nControlParameters))
    costSensitivity = 0.0_WP
    currentSensitivity = 0.0_WP
    instantaneousSensitivity = 0.0_WP

    ! Special case: check if the cost functional is a function of the control AND
    ! if the cost functional is measured only at the final time only
    if (costFunctionalType .eq. DATA_ASSIMILATION_COST .and.                                 &
         controllerType .eq. CHEMICAL_ACTUATOR) then
       call data_assimilation_compute_sensitivity
       currentSensitivity = instantaneousSensitivity * timeStepSize / 6.0_WP
    end if

    ! Setup the data dumping for postprocessing
    if (outputOn) then
       call dump_cleanup
       call dump_setup(ADJOINT)
       call dump_result(ADJOINT)
    end if

    ! Monitor the initial data
    call monitor_timestep(ADJOINT, controlIteration_)

    do i = startTimestep + sign(1, timemarchDirection),                                      &
         startTimestep + sign(adjointIterations, timemarchDirection), timemarchDirection

       ! Update the timestep
       timestep = i
       call get_timestep_size

       do j = nStages, 1, -1

          ! Load adjoint coefficients
          if (j .eq. 1) then
             call reverse_migrator_run(i, nStages)
          else
             call reverse_migrator_run(i + 1, j - 1)
          end if

          ! Update the sensitivity gradient and save if necessary
          call controller_compute_sensitivity
          currentSensitivity = currentSensitivity + timeNorm(j) * timeStepSize *             &
               instantaneousSensitivity

          ! Do not force adjoint at the initial time (needed for IC sensitivity)
          if (i.eq.startTimestep + sign(adjointIterations, timemarchDirection) .and. j.eq.1) &
               adjointCorrection = 0.0_WP

          ! Take a single adjoint sub-step using the time integrator
          call substep_adjoint(j)

          ! TODO: how to enforce limits on adjoint variables... check for NaN?
          call check_solution_limits(ADJOINT)

       end do

       ! Compute the sensitivity of the initial condition
       if (i .eq. startTimestep + sign(adjointIterations, timemarchDirection)) then
          call controller_compute_sensitivity_ic
          currentSensitivity = currentSensitivity + instantaneousSensitivity
       end if

       ! Filter solution if required
       call filter_solution(adjointVariables, i)

       ! Write restart files
       call write_restart_files(ADJOINT, controlIteration_)

       ! Dump results for visualization
       if (outputOn) call dump_result(ADJOINT)

       ! Monitor the current timestep
       call monitor_timestep(ADJOINT, controlIteration_)

    end do

    ! Call controller hooks after time marching ends
    call controller_hook_after_timemarch(ADJOINT)

    ! Cleanup the reverse migrator
    call reverse_migrator_cleanup

    ! Cleanup the timer
    call timing_cleanup

    ! Set the cost sensitivity
    costSensitivity = currentSensitivity

    return
  end function run_adjoint


  ! Load the initial conditions for the forward/adjoint run
  ! -------------------------------------------------------
  subroutine load_initial_condition(mode)

    ! External modules
    use parser
    use state_io
    use functional, only : adjointForcingFactor

    implicit none

    ! Arguments
    integer, intent(in) :: mode

    ! Local variables
    integer :: i
    character(len = str_medium) :: filename
    logical :: fileIsThere

    select case (mode)

    case (FORWARD) !... initialize conserved variables

       ! Initialize from file
       call parser_read('solution file to read', filename)
       call simulation_read(IO_FORWARD_STATE, trim(adjustl(filename)))

       ! Read in particles
       if (useParticles) then
          call parser_read('particle file to read', filename)
          call simulation_read(IO_PARTICLE, trim(adjustl(filename)))
          if (interpFluidToParticle) call overwrite_particle_velocity
       end if

       ! Read IBM and assemble objects
       if (useIBM) then
          call parser_is_defined('ibm file to read', fileIsThere)
          if (fileIsThere) then
             call parser_read('ibm file to read', filename)
             call simulation_read(IO_IBM, trim(adjustl(filename)))
          end if
          call ibm_setup_objects
       end if

       ! Adjust for particles if using two-way coupling
       if (twoWayCoupling) then
          call compute_volume_fraction
          do i = 1, nUnknowns
             conservedVariables(:,i) = conservedVariables(:,i) * volumeFraction(:,1)
          end do
       end if

       ! Store initial momentum
       call momentum_source_init

       ! Correct the state variables
       call correct_state(conservedVariables)

    case (ADJOINT) !... initialize adjoint variables

       if (.not. predictionOnly) then
          call parser_read('adjoint file to read', filename, '')
          if (len_trim(filename) .eq. 0) then
             adjointVariables = 0.0_WP
          else
             call simulation_read(IO_ADJOINT_STATE, trim(filename))
             return
          end if
       end if

       if (.not. useContinuousAdjoint) then

          call get_timestep_size

          adjointForcingFactor = - timeStepSize / 6.0_WP !... RK4 only

          ! Update the adjoint coefficient time
          adjointCoefficientTime = time

          ! Account for adjoint forcing from cost functional
          rightHandSide = 0.0_WP
          call functional_adjoint_source(rightHandSide)
          call functional_adjoint_source_tc(rightHandSide)
          
          adjointVariables = rightHandSide
          adjointForcingFactor = 1.0_WP !... restore

       end if

    end select

    return
  end subroutine load_initial_condition


  ! Write restart files for forward/adjoint runs
  ! --------------------------------------------
  subroutine write_restart_files(mode, controlIteration)

    ! External modules
    use state_io
    use ibm

    implicit none

    ! Arguments
    integer, intent(in) :: mode, controlIteration

    ! Local variables
    integer :: i
    character(len = str_medium) :: filename

    ! Return if output is turned off
    if (.not. outputOn) return

    if (saveInterval .gt. 0 .and. mod(timestep, max(1, saveInterval)) .eq. 0) then

       select case (mode)

       case (FORWARD)

          ! Remove volume fraction before writing
          if (twoWayCoupling) then
             call compute_volume_fraction
             do i = 1, nUnknowns
                conservedVariables(:,i) = conservedVariables(:,i) / volumeFraction(:,1)
             end do
          end if

          ! Write the state variables
          if (dataOverwrite) then
             filename = trim(solutionFile)
          else
             !if (controlIteration .gt. 0) then
             !   write(filename, '(2A,I2.2,A,I8.8)') trim(solutionFile), "_",                 &
             !        controlIteration, ".", timestep
             !else
                write(filename, '(2A,I8.8)') trim(solutionFile), ".", timestep
             !end if
          end if
          call simulation_write(IO_FORWARD_STATE, trim(filename))

          ! Multiply back the volume fraction
          if (twoWayCoupling) then
             do i = 1, nUnknowns
                conservedVariables(:,i) = conservedVariables(:,i) * volumeFraction(:,1)
             end do
          end if

          ! Write the particles
          if (useParticles) then
             if (dataOverwrite) then
                filename = trim(partFile)
             else
                !if (controlIteration .gt. 0) then
                !   write(filename, '(2A,I2.2,A,I8.8)') trim(partFile), "_",                  &
                !        controlIteration, ".", timestep
                !else
                   write(filename, '(2A,I8.8)') trim(partFile), ".", timestep
                !end if
             end if
             call simulation_write(IO_PARTICLE, trim(filename))
          end if

          ! Write IBM objects
          if (useIBM .and. nObjects.gt.0) then
             if (dataOverwrite) then
                filename = trim(ibmFile)
             else
                write(filename, '(2A,I8.8)') trim(ibmFile), ".", timestep
             end if
             call simulation_write(IO_IBM, trim(filename))
          end if

       case (ADJOINT)
          !if (controlIteration .gt. 0) then
          !   write(filename, '(2A,I2.2,A,I8.8)') trim(adjointFile), "_",                     &
          !        controlIteration, ".", timestep
          !else
             write(filename, '(2A,I8.8)') trim(adjointFile), ".", timestep
          !end if
          call simulation_write(IO_ADJOINT_STATE, trim(filename))

       end select
    end if

    return
  end subroutine write_restart_files

  ! Check if the solution is within a user-defined range
  ! ----------------------------------------------------
  subroutine check_solution_limits(mode)

    ! External modules
    use parallel
    use grid_functions, only : isVariableWithinRange
    use state_io

    implicit none

    ! Arguments
    integer, intent(in) :: mode

    ! Local variables
    integer :: i, iGlobal, jGlobal, kGlobal, rankReportingError, ierror
    character(len = str_long) :: message
    real(WP) :: fOutsideRange

    if (.not. enableSolutionLimits) return

    rankReportingError = -1

    ! Adjust for volume fraction
    if (twoWayCoupling) then
       do i = 1, nUnknowns
          conservedVariables(:, i) = conservedVariables(:, i) / volumeFraction(:,1)
       end do
    end if

    do

       ! Check density
       if (.not. isVariableWithinRange(conservedVariables(:,1), fOutsideRange,               &
            iGlobal, jGlobal, kGlobal,                                                       &
            minValue = densityRange(1), maxValue = densityRange(2))) then
          write(message, '(3(A,I0.0),3(A,(SS,ES9.2E2)),A)') "Density at (",                  &
               iGlobal, ", ", jGlobal, ", ", kGlobal, "): ", fOutsideRange,                  &
               " out of range (", densityRange(1), ", ", densityRange(2), ")!"
          rankReportingError = iRank
          exit
       end if

       ! Check temperature
       if (.not. isVariableWithinRange(temperature(:,1), fOutsideRange,                      &
            iGlobal, jGlobal, kGlobal,                                                       &
            minValue = temperatureRange(1), maxValue = temperatureRange(2))) then
          write(message, '(3(A,I0.0),3(A,(SS,ES9.2E2)),A)') "Temperature at (",              &
               iGlobal, ", ", jGlobal, ", ", kGlobal, "): ", fOutsideRange,                  &
               " out of range (", temperatureRange(1), ", ", temperatureRange(2), ")!"
          rankReportingError = iRank
          exit
       end if

       ! Check mass fraction
       do i = 1, nSpecies
          if (.not. isVariableWithinRange(massFraction(:,i), fOutsideRange,                  &
               iGlobal, jGlobal, kGlobal,                                                    &
               minValue = massFractionRange(1), maxValue = massFractionRange(2))) then
             write(message, '(3(A,I0.0),3(A,(SS,ES9.2E2)),A)') "Mass fraction at (",         &
                  iGlobal, ", ", jGlobal, ", ", kGlobal, "): ", fOutsideRange,               &
                  " out of range (", massFractionRange(1), ", ", massFractionRange(2), ")!"
             rankReportingError = iRank
             exit
          end if
       end do

       exit

    end do

    call parallel_max(rankReportingError)

    if (rankReportingError .ne. -1) then

       if (iRank .eq. iRoot .and. rankReportingError .ne. iRoot)                             &
            call MPI_Recv(message, str_long, MPI_CHARACTER, rankReportingError,              &
            rankReportingError, comm, MPI_STATUS_IGNORE, ierror)
       if (iRank .eq. rankReportingError .and. rankReportingError .ne. iRoot)                &
            call MPI_Send(message, str_long, MPI_CHARACTER, iRoot, iRank, comm, ierror)

       select case (mode)
       case (FORWARD)
          call simulation_write(IO_FORWARD_STATE, trim(solutionFile) // '.crashed')
       end select

       if (iRank .eq. iRoot) call die(trim(message))

    end if

    ! Adjust for volume fraction
    if (twoWayCoupling) then
       do i = 1, nUnknowns
          conservedVariables(:, i) = conservedVariables(:, i) * volumeFraction(:,1)
       end do
    end if

    return
  end subroutine check_solution_limits

  ! Check whether the simulation is done or not
  ! -------------------------------------------
  function done()

    ! External modules
    use parallel
    use time_info

    implicit none

    ! Arguments
    logical :: done

    ! Local arguments
    real(WP) :: wallTime
    
    done = .false.
    
    ! Check if number of iterations have been reached
    if (timestep .ge. startTimestep + nTimesteps) then
       done = .true.
       call monitor_log("MAXIMUM NUMBER OF ITERATION REACHED")
    end if

    if (predictionOnly) then

       ! Check if final time has been reached
       if (time .ge. finalTime) then
          done = .true.
          call monitor_log("MAXIMUM SIMULATION TIME REACHED")
       end if

       ! Check if wall clock time is reached
       call parallel_time(wallTime)
       wallTime = wallTime - initWallTime
       if (wallTime .ge. maxWallTime) then
          done = .true.
          call monitor_log("MAXIMUM WALL TIME REACHED")
       end if
       
    end if
    
    return
  end function done

end module solver


! ================ !
! Setup the solver !
! ================ !
subroutine solver_setup

  ! Internal modules
  use solver

  ! External modules
  use parser
  use grid_io
  use state_io

  implicit none

  ! Local variables
  character(len = str_medium) :: filename

  ! Clean slate
  call solver_cleanup

  ! Overwrite data output?
  call parser_read('data overwrite', dataOverwrite, .true.)
  if (.not. predictionOnly) dataOverwrite = .false.

  ! Read filenames from the input
  call parser_read('solution file to write', solutionFile)
  if (useParticles) call parser_read('particle file to write', partFile)
  if (useIBM) call parser_read('ibm file to write', ibmFile, '')
  if (.not. predictionOnly) call parser_read('adjoint file to write', adjointFile)

  ! Enable solution limits?
  call parser_read('enable solution limits', enableSolutionLimits, .false.)
  if (enableSolutionLimits) then
     call parser_read('minimum density', densityRange(1))
     call parser_read('maximum density', densityRange(2))
     call parser_read('minimum temperature', temperatureRange(1))
     call parser_read('maximum temperature', temperatureRange(2))
     if (nSpecies .gt. 0) then
        call parser_read('minimum mass fraction', massFractionRange(1))
        call parser_read('maximum mass fraction', massFractionRange(2))
     end if
  end if

  ! Overwrite particle velocity with local fluid velocity
  call parser_read('overwrite particle velocity', interpFluidToParticle, .false.)

  ! Read in time limits
  call parser_read('number of timesteps', nTimesteps, huge(1))
  nTimesteps = max(0, nTimesteps)
  if (predictionOnly) then
     call parser_read('final time', finalTime, huge(1.0_WP))
     call parser_read('maximum wall time', maxWallTime, huge(1.0_WP))
     maxWallTime = max(0.0_WP, maxWallTime)
     maxWallTime = maxWallTime * 3600.0_WP
     call parallel_time(initWallTime)
  end if

  call parser_read('save interval', saveInterval, -1)
  if (saveInterval .eq. 0) saveInterval = -1

  call parser_read('adjoint iterations', adjointIterations, nTimesteps)
  adjointIterations = max(0, adjointIterations)
  adjointIterations = min(nTimesteps, adjointIterations)

  ! Setup the time integrator
  call time_integrator_setup

  ! If a target state file was specified, read the target state. Otherwise, initialize the
  ! target state to a quiescent state by default
  if (useTargetState) then
     call parser_read('target state file', filename, '')
     if (len_trim(filename) .eq. 0) then
        call make_quiescent(targetState)
     else
        call simulation_read(IO_TARGET_STATE, trim(filename))
     end if
  end if

  ! Setup the boundary conditions
  call boundary_setup

  if (.not. predictionOnly) then

     ! Setup the controller
     call controller_setup

     ! Setup the cost functional
     call functional_setup

  end if

  return
end subroutine solver_setup


! ================== !
! Cleanup the solver !
! ================== !
subroutine solver_cleanup

  ! Internal modules
  use solver

  implicit none

  call timing_cleanup
  call dump_cleanup
  call time_integrator_cleanup
  call functional_cleanup
  call controller_cleanup

  return
end subroutine solver_cleanup
