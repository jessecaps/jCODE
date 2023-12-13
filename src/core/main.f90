program jcode

  call main_setup
  call main_solve
  call main_finalize

end program jcode




subroutine main_setup

  ! External modules
  use string

  implicit none

  ! Local variables
  character(len = str_medium) :: input

  ! Initialize parallel environment
  call parallel_init

  ! Initialize the random number generator
  call random_setup

  ! Parse the command line
  call parallel_get_inputname(input)
  
  ! Parse the input file
  call parser_init
  call parser_parsefile(input)

  ! Setup the monitoring system
  call monitor_setup
  call monitor_log('Begin simulation')
  call monitor_welcome

  ! Setup the simulation
  call simulation_setup

  ! Report some useful information
  call monitor_grid_diagnostics

  return
end subroutine main_setup




subroutine main_solve

  ! External modules
  use parser
  use precision
  use simulation_flags
  use solver_options
  use controller
  use solver, only : run_forward, run_adjoint

  implicit none

  ! Local variables
  real(WP) :: forwardSolution, actuationAmount
  real(WP), dimension(:), pointer :: adjointSolution
  logical :: checkGradientAccuracy, findOptimalForcing, singleControlledPrediction

  call parser_read('check gradient accuracy', checkGradientAccuracy, .false.)
  call parser_read('find optimal forcing', findOptimalForcing, .false.)
  call parser_read('single controlled prediction', singleControlledPrediction, .false.)
  call parser_read('actuation amount', actuationAmount, 1.0_WP)

  ! Main code logic
  if (predictionOnly) then !... just a predictive simulation
     forwardSolution = run_forward()
  else if (checkGradientAccuracy) then !... verify gradient is exact
     call check_gradient_accuracy
  else if (findOptimalForcing) then !... adjoint optimization
     call find_optimal_forcing
  else if (singleControlledPrediction) then
     forwardSolution = run_forward(actuationAmount = actuationAmount)
  else
     if (.not. isBaselineAvailable) forwardSolution = run_forward()
     if (.not. isGradientAvailable) then
        allocate(adjointSolution(nControlParameters))
        adjointSolution = run_adjoint()
        deallocate(adjointSolution)
     end if
  end if

  return
end subroutine main_solve



subroutine main_finalize

  ! External modules
  use parallel

  implicit none

  ! Cleanup the simulation
  call simulation_cleanup

  ! Finalize output
  call monitor_log('Simulation complete')
  call monitor_finalize

  ! Finalize the parallel environment
  call parallel_finalize

  return
end subroutine main_finalize




! -------------------------------------------
subroutine die(errorText)

  ! External modules
  use parallel

  implicit none

  ! Arguments
  character(len = *), intent(in) :: errorText
  
  call monitor_log('KILLED')
  call parallel_kill(errorText)
  
  return
end subroutine die
