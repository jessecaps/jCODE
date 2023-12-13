module map_ignition

  ! External modules
  use gradient_search

  implicit none

end module map_ignition


! ================================================= !
! Traverse the ignition boundary using the          !
! tangent vector obtained from the adjoint gradient !
! ================================================= !
subroutine run_map_ignition

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
       p1, p2, parameterDirection, iunit, iostat
  character(len = str_medium) :: filename
  real(WP) :: baselineCostFunctional, costFunctional, indicatorFunction, costSensitivity,    &
       actuationAmount, previousActuationAmount, tangentActuationAmount, minTolerance,       &
       burnValue, previousSensitivity
  real(WP), dimension(:), allocatable :: individualSensitivities, parameters, tangentVector
  logical:: done, burning, minimizeParameter, foundNewMinimum

  ! Read in optimization parameters
  call parser_read('minimum actuation tolerance', minTolerance, 1.0E-9_WP)

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

  ! Assume we are at the burn boundary. We have at this point a baseline cost functional
  ! and a gradient with respect to each parameter
  controlIteration = restartIteration
  tangentActuationAmount = actuationAmount
  allocate(tangentVector(nControlParameters))

  ! Determine which parameters to control/adjust while keeping all others constant
  call parser_read('control parameter', p1)
  call parser_read('adjust parameter', p2)
  if (p1.gt.nControlParameters)                                                              &
       call die('map_ignition: control parameter must be <= number of control parameters')
  if (p2.gt.nControlParameters)                                                              &
       call die('map_ignition: adjust parameter must be <= number of control parameters')
  call parser_read('minimize control parameter', minimizeParameter, .true.)
  parameterDirection = 1
  if (minimizeParameter) parameterDirection = -1

  ! Perform line search
  do i = controlIteration, restartIteration + nIterations - 1

     ! Compute the tangent vector
     tangentVector = 0.0_WP
     tangentVector(p1) = real(parameterDirection, WP)
     tangentVector(p2) = - tangentVector(p1) * individualSensitivities(p1) /                 &
          individualSensitivities(p2)
     controlGradient = tangentVector / sqrt(sum(tangentVector**2))

     ! Compute a new cost functional in the tangent space
     gradientDirection = 1
     costFunctional = run_forward(actuationAmount = tangentActuationAmount,                  &
          controlIteration = nForward)
     nForward = nForward + 1
     indicatorFunction = auxilaryCostFunctional
     burning = indicatorFunction .gt. burnValue

       ! Update the baseline values
       do j = 1, nControlParameters
          baselineValue(j) = baselineValue(j) + real(gradientDirection, WP) *                &
               tangentActuationAmount * controlGradient(j)
       end do

       ! Check if we are close to the ignition boundary and correct if needed
       actuationAmount = minTolerance
       previousActuationAmount = 0.0_WP
       foundNewMinimum = .false.
       controlGradient = individualSensitivities / sqrt(costSensitivity)
       if (burning) gradientDirection = -1
       do
          costFunctional = run_forward(actuationAmount = actuationAmount,                    &
               controlIteration = nForward)
          nForward = nForward + 1
          indicatorFunction = auxilaryCostFunctional
          if ((burning .and. indicatorFunction .lt. burnvalue) .or.                          &
               (.not.burning .and. indicatorFunction .gt. burnvalue)) then
             exit
          end if
          foundNewMinimum = .true.
          previousActuationAmount = actuationAmount
          actuationAmount = actuationAmount + minTolerance
          exit
       end do

       if (foundNewMinimum) then
          ! Update the baseline values and compute a new sensitivity gradient
          do j = 1, nControlParameters
             baselineValue(j) = baselineValue(j) + real(gradientDirection, WP) *             &
                  previousActuationAmount * controlGradient(j)
          end do
          individualSensitivities = run_adjoint(controlIteration = nAdjoint)
          nAdjoint = nAdjoint + 1
          costSensitivity = sum(individualSensitivities**2)
       end if

       ! Output progress
       if (iRank .eq. iRoot) then
          write(iunit, '(I4,1000(1X,SP,(SP,ES23.15E3)))') i,                                &
               actuationAmount, costFunctional, indicatorFunction,                          &
               (baselineValue(j) + real(gradientDirection, WP) * actuationAmount *          &
               controlGradient(j), j = 1, nControlParameters),                              &
               (individualSensitivities(j), j = 1, nControlParameters)
          flush(iunit)
       end if

       controlIteration = controlIteration + 1

    end do

    ! Dump optimizatiom summary and close
    if (iRank .eq. iRoot) then
       write(iunit, *) ''
       write(iunit, '(A28,I4)') 'Number of forward runs:', nForward
       write(iunit, '(A28,I4)') 'Number of adjoint runs:', nAdjoint
       close(iunit)
    end if

    ! Clean up
    deallocate(parameters, individualSensitivities, tangentVector)

  return
end subroutine run_map_ignition
