module ignition_boundary

  ! External modules
  use gradient_search

  implicit none

end module ignition_boundary


! ===================================== !
! Find a point on the ignition boundary !
! ===================================== !
subroutine run_ignition_boundary

  ! Internal modules
  use ignition_boundary

  ! External modules
  use parser
  use parallel
  use fileio
  use functional
  use controller

  implicit none

  ! Local variables
  integer :: i, j, nIterations, restartIteration, controlIteration, nForward, nAdjoint,      &
       iunit, iostat
  character(len = str_medium) :: filename
  real(WP) :: baselineCostFunctional, costFunctional, indicatorFunction, costSensitivity,    &
       actuationAmount, minimumTolerance, burnValue, previousSensitivity
  real(WP), dimension(:), allocatable :: individualSensitivities, parameters
  logical:: done, burning

  ! Read in optimization parameters
  call parser_read('minimum actuation tolerance', minimumTolerance, 1.0E-9_WP)

  call parser_read('number of control iterations', nIterations)
  if (nIterations .lt. 0) call die('ignition_boundary: number of control iterations&
       & must be a non-negative number!')

  call parser_read('restart control iteration', restartIteration, 0)
  restartIteration = max(restartIteration, 0)

  if (restartIteration .gt. 0 .and. .not. isBaselineAvailable)                               &
       call die('ignition_boundary: cannot restart with controlled&
       &  prediction without baseline!')

  ! Set the filename to output optimization progress
  filename = 'cost_optimization.txt'
  if (iRank .eq. iRoot) then
     iunit = iopen()
     if (restartIteration .eq. 0 .and. .not. isBaselineAvailable) then
        open(iunit, file = trim(filename), action ='write', status = 'unknown',              &
             iostat = iostat)
     else
        open(iunit, file = trim(filename), action ='readwrite', status = 'old',              &
             position = 'rewind', iostat = iostat)
     end if
  end if

  call parallel_bc(iostat)
  if (iostat .ne. 0) call die('ignition_boundary: failed to open file for writing!')

  ! Find (or load from file) useful data from the baseline prediction
  allocate(parameters(nControlParameters))
  allocate(individualSensitivities(nControlParameters))
  if (isBaselineAvailable) then
     if (iRank .eq. iRoot) then
        read(iunit, *, iostat = iostat)
        read(iunit, *, iostat = iostat) i, actuationAmount, baselineCostFunctional,          &
             indicatorFunction, parameters, individualSensitivities
     end if
     call parallel_bc(iostat)
     if (iostat .ne. 0) call die('ignition_boundary:&
          & failed to read baseline cost function from file!')
     call parallel_bc(baselineCostFunctional)
     call parallel_bc(indicatorFunction)
     call parallel_bc(parameters)
     baselineValue = parameters
  else
     baselineCostFunctional = run_forward()
     indicatorFunction = auxilaryCostFunctional
  end if

  ! Find the initial sensitivity gradient
  if (restartIteration .eq. 0) then
     individualSensitivities = run_adjoint()
  else
     call parallel_bc(individualSensitivities)
  end if
  costSensitivity = sum(individualSensitivities**2)

  if (restartIteration .eq. 0) restartIteration = restartIteration + 1

  ! Number of forward/adjoint runs
  nForward = restartIteration
  nAdjoint = 1

  ! Find the previous data
  if (restartIteration .le. 1) then
     call parser_read('initial actuation amount', actuationAmount)
  else
     if (iRank.eq.iRoot) then
        do i = 1, restartIteration - 1
           previousSensitivity = individualSensitivities(1)
           read(iunit, *, iostat = iostat) j, actuationAmount, costFunctional,                  &
                indicatorFunction, parameters, individualSensitivities
           if (abs(individualSensitivities(1)-previousSensitivity) .lt. epsilon(1.0_WP)) then
              baselineValue = parameters
              nAdjoint = nAdjoint + 1
           end if
        end do
     end if
     call parallel_bc(actuationAmount)
     call parallel_bc(costFunctional)
     call parallel_bc(indicatorFunction)
     call parallel_bc(parameters)
     call parallel_bc(individualSensitivities)
     call parallel_bc(baselineValue)
     call parallel_bc(nAdjoint)
     costSensitivity = sum(individualSensitivities**2)
  end if

  ! Store gradient to be used for control forcing
  controlGradient = individualSensitivities / sqrt(costSensitivity)

  if (iRank .eq. iRoot .and. .not. isBaselineAvailable) then
     write(iunit, '(A4,1000A24)') 'i', 'Actuation amount', 'Cost functional',                &
          'Inidicator functional',                                                           &
          (trim(sensitivityParameter(i)), i = 1, nControlParameters),                        &
          ('dJ/d ' // trim(sensitivityParameter(i)), i = 1, nControlParameters)
     write(iunit, '(I4,1000(1X,SP,(SP,ES23.15E3)))') 0, 0.0_WP,                              &
          baselineCostFunctional, indicatorFunction,                                         &
          (baselineValue(i), i = 1, nControlParameters),                                     &
          (individualSensitivities(i), i = 1, nControlParameters)
     flush(iunit)
  end if

  if (nIterations .eq. 0) return

  ! Turn off output for controlled predictions
  call parser_read('output control iterations', outputOn, .false.)

  ! Determine if the initial run was burning based on the last value of the
  ! instantaneous cost functional
  call parser_read('burn value', burnValue)
  burning = indicatorFunction .gt. burnValue

  ! We have at this point a baseline cost functional and the first gradient with respect
  ! to each parameter
  done = .false.
  controlIteration = restartIteration
  do while (controlIteration .lt. nIterations .and. .not.done)

     ! Perform line search
     do i = controlIteration + 1, restartIteration + nIterations

        ! Choose a direction to march
        gradientDirection = -1
        if (.not. burning) gradientDirection = 1

        ! Compute a new cost functional
        costFunctional = run_forward(actuationAmount = actuationAmount,                      &
             controlIteration = nForward)
        indicatorFunction = auxilaryCostFunctional
        nForward = nForward + 1
        controlIteration = controlIteration + 1

        ! Output progress
        if (iRank .eq. iRoot) then
           write(iunit, '(I4,1000(1X,SP,(SP,ES23.15E3)))') i,                                &
                actuationAmount, costFunctional, indicatorFunction,                          &
                (baselineValue(j) + real(gradientDirection, WP) * actuationAmount *          &
                controlGradient(j), j = 1, nControlParameters),                              &
                (individualSensitivities(j), j = 1, nControlParameters)
           flush(iunit)
        end if

        ! Exit loop and recompute the gradient if we didn't passed the threshold
        if (burning .and. indicatorFunction .gt. burnvalue) then
           exit
        else if (.not.burning .and. indicatorFunction .lt. burnValue) then
           exit
        end if

        ! If we made it this far, reduce the actuation amount and try again
        actuationAmount = 0.5_WP * actuationAmount

        ! Check actuation tolerance
        if (actuationAmount .lt. minimumTolerance) done = .true.

     end do

     ! Update the baseline values and compute a new sensitivity gradient
     if (.not.done .and. controlIteration .lt. nIterations) then
        do i = 1, nControlParameters
           baselineValue(i) = baselineValue(i) + real(gradientDirection, WP) *               &
                actuationAmount * controlGradient(i)
        end do
        individualSensitivities = run_adjoint(controlIteration = nAdjoint)
        nAdjoint = nAdjoint + 1
        costSensitivity = sum(individualSensitivities**2)
        controlGradient = individualSensitivities / sqrt(costSensitivity)! * baselineValue
     end if

  end do

  ! Dump optimizatiom summary and close
  if (iRank .eq. iRoot) then
     write(iunit, *) ''
     write(iunit, '(A28,I4)') 'Number of forward runs:', nForward
     write(iunit, '(A28,I4)') 'Number of adjoint runs:', nAdjoint
     if (minimizeCost) then
        write(iunit, '(A28,L4)') 'Minimization found:', done
     else
        write(iunit, '(A28,L4)') 'Maximization found:', done
     end if
     close(iclose(iunit))
  end if

  ! Clean up
  deallocate(parameters)
  deallocate(individualSensitivities)

  return
end subroutine run_ignition_boundary
