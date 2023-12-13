module gradient_accuracy

 ! External modules
  use solver

  implicit none

end module gradient_accuracy


! ============================================ !
! Adjoint vericiation: check gradient accuracy !
! ============================================ !
subroutine check_gradient_accuracy

  ! Internal modules
  use gradient_accuracy

  ! External modules
  use fileio
  use parallel
  use parser
  use controller

  implicit none

  ! Local variables
  integer :: i, j, nIterations, restartIteration, fileUnit, iostat
  character(len = str_medium) :: filename
  real(WP) :: actuationAmount, baselineCostFunctional, costFunctional, costSensitivity,      &
       initialActuationAmount, geometricGrowthFactor, gradientError, dummyValue
  real(WP), dimension(:), allocatable :: individualSensitivities
  logical :: minimizeCost

  call parser_read('number of control iterations', nIterations)
  if (nIterations .lt. 0) call die('check_gradient_accuracy:&
       &Number of control iterations must be a non-negative number!')

  call parser_read('restart control iteration', restartIteration, 0)
  restartIteration = max(restartIteration, 0)

  if (restartIteration .gt. 0 .and. .not. isBaselineAvailable)                               &
       call die("check_gradient_accuracy: Can't restart with controlled prediction&
       & without baseline!")

  filename = 'gradient_error.txt'
  if (iRank .eq. iRoot) then
     if (restartIteration .eq. 0 .and. .not. isBaselineAvailable) then
        fileUnit = iOpen()
        open(unit = fileUnit, file = trim(filename), action = 'write', status = 'unknown',   &
             iostat = iostat)
     else
        fileUnit = iOpen()
        open(unit = fileUnit, file = trim(filename), action = 'readwrite', status = 'old',   &
             position = 'rewind', iostat = iostat)
     end if
  end if

  call parallel_bc(iostat)
  if (iostat .ne. 0) call die('check_gradient_accuracy: Failed to open file for writing!')

  if (nIterations .gt. 0) then
     call parser_read('initial actuation amount', initialActuationAmount)
     gradientDirection = -1
     call parser_read('minimize cost functional', minimizeCost, .true.)
     if (.not. minimizeCost) gradientDirection = 1
     call parser_read('actuation amount geometric growth', geometricGrowthFactor,            &
          1.0_WP / (10.0_WP ** 0.25_WP)) ! ... default yields 4 iterations per decade
  end if

  ! Find (or load from file) the cost functional & sensitivity for the baseline prediction
  allocate(individualSensitivities(nControlParameters))
  if (isBaselineAvailable) then
     if (iRank .eq. iRoot) then
        read(fileUnit, *, iostat = iostat)
        read(fileUnit, *, iostat = iostat) i, actuationAmount,                               &
             baselineCostFunctional, costSensitivity, gradientError
     end if
     call parallel_bc(iostat)
     if (iostat .ne. 0) call die(trim(filename) //                                           &
          ': Failed to read baseline cost functional from file!')
     call parallel_bc(baselineCostFunctional)
     if (.not.spaceTimeGradient) call load_cost_sensitivity
  else
     baselineCostFunctional = run_forward()
  end if

  ! Find the sensitivity gradient (this is the only time the adjoint simulation will be run)
  if (restartIteration .eq. 0) then
     individualSensitivities = run_adjoint()
  else if (spaceTimeGradient) then
     individualSensitivities(1) = costSensitivity
  else
     individualSensitivities = currentSensitivity
  end if
  if (spaceTimeGradient) then
     costSensitivity = individualSensitivities(1)
  else
     costSensitivity = sum(individualSensitivities**2)
  end if

  ! Store gradient to be used for control forcing
  controlGradient = individualSensitivities

  if (iRank .eq. iRoot .and. .not. isBaselineAvailable) then
     write(fileUnit, '(A4,4A24)') 'i', 'Actuation amount', 'Cost functional',                &
          'Cost sensitivity', 'Gradient error'
     write(fileUnit, '(I4,4(1X,SP,(SP,ES23.15E3)))') 0, 0.0_WP, baselineCostFunctional,      &
          costSensitivity, 0.0_WP
     flush(fileUnit)
  end if

  if (nIterations .eq. 0) return

  ! Turn off output for controlled predictions
  call parser_read('output control iterations', outputOn, .false.)

  if (restartIteration .eq. 0) restartIteration = restartIteration + 1

  do i = 1, restartIteration - 1
     if (iRank .eq. iRoot)                                                                   &
          read(fileUnit, *, iostat = iostat) j, actuationAmount, costFunctional,             &
          dummyValue, gradientError
     call parallel_bc(iostat)
     if (iostat .ne. 0) call die(trim(filename) //                                           &
          ': Cost functional history is too short for the specified restart iteration!')
  end do

  do i = restartIteration, restartIteration + nIterations - 1
     actuationAmount = initialActuationAmount * geometricGrowthFactor ** real(i - 1, WP)
     costFunctional = run_forward(actuationAmount = actuationAmount, controlIteration = i)
     if (spaceTimeGradient) then
        gradientError = abs((costFunctional - baselineCostFunctional) / actuationAmount +    &
             costSensitivity)
     else
        gradientError = abs( abs(costFunctional - baselineCostFunctional) /                  &
             actuationAmount - costSensitivity )
     end if
     if (iRank .eq. iRoot) then
        write(fileUnit, '(I4,4(1X,SP,(SP,ES23.15E3)))') i, actuationAmount, costFunctional,  &
             (costFunctional - baselineCostFunctional) / actuationAmount, gradientError
        flush(fileUnit)
     end if
  end do

  if (iRank .eq. iRoot) close(iclose(fileUnit))

  return
end subroutine check_gradient_accuracy
