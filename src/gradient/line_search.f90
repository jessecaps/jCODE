module line_search

  ! External modules
  use gradient_search

  implicit none

  ! Line search methods
  integer, parameter :: BRACKET = 0

  ! Global variables
  integer :: nLineSearchIterations, lineSearchMethod
  real(WP) :: initialActuationAmount, defaultInitialActuationAmount, maxActuationAmount,     &
       defaultMaxActuationAmount, mu1, mu2
  real(WP), dimension(:), allocatable :: costFunctional, actuationAmount

end module line_search


! ===================================== !
! Initialize line search optimiziation  !
! ===================================== !
subroutine initialize_line_search

  ! Internal modules
  use line_search

  ! External modules
  use parser

  implicit none

  ! Local variables
  character(len = str_medium) :: key

  ! Read in line search information
  call parser_read('number of line search iterations', nLineSearchIterations, 20)
  nLineSearchIterations = max(1, nLineSearchIterations)
  call parser_read('maximum actuation amount', defaultMaxActuationAmount,huge(1.0_WP))
  defaultMaxActuationAmount = max(0.0_WP, defaultMaxActuationAmount)
  call parser_read('initial actuation amount', defaultInitialActuationAmount, 1.0_WP)
  if (searchDirectionMethod .eq. BFGS) defaultInitialActuationAmount = 1.0_WP
  defaultInitialActuationAmount = min(defaultInitialActuationAmount,defaultMaxActuationAmount)

  call parser_read('line search method', key, 'bracket')
  select case (trim(key))
  case ('bracket')
     lineSearchMethod = BRACKET
  case default
     call die('initialize_line_search: Unknown line search method: "' // trim(key) // '"')
  end select

  ! Some optimization parameters
  mu1 = 1.0E-4_WP ! ... to staisfy a small fraction of the expected decrease
  select case (searchDirectionMethod)
  case (CONJUGATE_GRADIENT)
     mu2 = 0.1_WP ! ... a positive number between 0 and 1 to determine when bracketing happens
  case (BFGS)
     mu2 = 0.9_WP ! ... a positive number between 0 and 1 to determine when bracketing happens
  end select

  return
end subroutine initialize_line_search


! ================================= !
! Perform line search optimiziation !
! ================================= !
subroutine run_line_search(scaledBaselineValue, baselineCostFunctional, baselineGradient,     &
     searchDirection, optActuationAmount)

  ! Internal modules
  use line_search

  implicit none

  ! Arguments
  real(WP), intent(in) :: baselineCostFunctional
  real(WP), dimension(nControlParameters), intent(in) :: scaledBaselineValue,                 &
       baselineGradient, searchDirection
  real(WP), intent(out) :: optActuationAmount

  ! Local variables
  integer :: i

  ! Allocate memory for actuation amount and cost function sub iterations (during line search)
  allocate(actuationAmount(0:nLineSearchIterations), costFunctional(0:nLineSearchIterations))

  ! Compute maximum actuation amount according to control parameter bounds
  maxActuationAmount = defaultMaxActuationAmount
  do i = 1, nControlParameters
     if ((searchDirection(i) .gt. 0.0_WP) .and.                                              &
          (parameterUpperBound(i) - scaledBaselineValue(i) .gt. 0.0_WP .or.                   &
             nControlParameters .eq. 1)) then
        maxActuationAmount = min(maxActuationAmount,                                         &
             (parameterUpperBound(i) - scaledBaselineValue(i)) / searchDirection(i))
     else if ((searchDirection(i) .lt. 0.0_WP) .and.                                         &
          (parameterLowerBound(i) - scaledBaselineValue(i) .lt. 0.0_WP .or.                   &
             nControlParameters .eq. 1)) then
        maxActuationAmount = min(maxActuationAmount,                                         &
             (parameterLowerBound(i) - scaledBaselineValue(i)) / searchDirection(i))
     end if
  end do

  ! Update initial actuation amount
  initialActuationAmount = min(defaultInitialActuationAmount, maxActuationAmount)

  select case (lineSearchMethod)
  case (BRACKET)
     call line_search_bracket(baselineCostFunctional, baselineGradient, searchDirection,     &
          optActuationAmount)
  end select

  ! Clean up
  deallocate(actuationAmount, costFunctional)

  return
end subroutine run_line_search


! An algorithm for doing the line search to satisfy the strong Wolfe conditions
! Based on Algorithm 2 of Multidisciplinary Design Optimization by J.R.R.A. Martins, A. Ning,&
! and J. Hicken
subroutine line_search_bracket(baselineCostFunctional, baselineGradient, searchDirection,    &
     optActuationAmount)

  ! Internal modules
  use line_search

  implicit none

  ! Arguments
  real(WP), intent(in) :: baselineCostFunctional
  real(WP), dimension(nControlParameters), intent(in) :: baselineGradient, searchDirection
  real(WP), intent(out) :: optActuationAmount

  ! Local Variables
  integer :: i, j
  real(WP) :: dphi0, dphiOld, dphi, costFunctionalGradient(nControlParameters),              &
       optCostFunctional

  ! Set/Compute values for zero actuation amount
  i = 0 ! ... line search iteration number
  actuationAmount(0) = 0.0_WP
  costFunctional(0) = baselineCostFunctional 
  dphi0 = sum( baselineGradient * searchDirection ) ! ... directional derivative
  if (i .ge. nLineSearchIterations) then
     optActuationAmount = 0.0_WP
     return
  end if

  ! Initialize 
  i = 1
  actuationAmount(i) = initialActuationAmount ! ... initial step size
  dphiOld = dphi0 ! ... dphi that corrresponds to actuationAmount(i-1)

  do
     ! Evaluate cost functional
     costFunctional(i) = objective_function(actuationAmount = actuationAmount(i),            &
          searchDirection = searchDirection)

     if (costFunctional(i) .gt. costFunctional(0) + mu1 * actuationAmount(i) * dphi0 .or.    &
          (costFunctional(i) .gt. costFunctional(i-1) .and. i .gt. 1 )) then
        call line_search_zoom(i-1, i, dphi0, dphiOld, i, searchDirection, optActuationAmount)
        return
     end if

     ! Evaluate the gradient of cost functional
     costFunctionalGradient = objective_function_gradient(actuationAmount=actuationAmount(i),&
          searchDirection = searchDirection, costFunctional = costFunctional(i))
     dphi = sum( costFunctionalGradient * searchDirection ) ! ... directional derivative
     if (abs(dphi) .le. -mu2*dphi0 ) then ! ... dphi0 is negative!
        optActuationAmount = actuationAmount(i)
        return
     elseif (dphi.ge.0.0_WP) then
        call line_search_zoom(i, i-1, dphi0, dphi, i, searchDirection, optActuationAmount)
        return
     elseif (i .lt. nLineSearchIterations) then
        actuationAmount(i+1) = 2.0_WP * actuationAmount(i)
        if (actuationAmount(i+1).gt.maxActuationAmount)                                      &
             actuationAmount(i+1) = 0.5_WP * (actuationAmount(i) + maxActuationAmount)
        i = i + 1
        dphiOld = dphi
     else ! ... line searching is enough
        ! Determine the actuation amount with smallest cost functional value
        optActuationAmount = actuationAmount(0)
        optCostFunctional = costFunctional(0)
        do j = 1, nLineSearchIterations
           if (costFunctional(j) .le. optCostFunctional) then
              optActuationAmount = actuationAmount(j)
              optCostFunctional = costFunctional(j)
           end if
        end do
        return
     end if

  end do

end subroutine line_search_bracket


! The algorithm of zoom function, which is as a subroutine for 
! the bracketing line search, to satisfy the strong Wolfe conditions
! Based on Algorithm 3 of Multidisciplinary Design Optimization by J.R.R.A. Martins,
! A. Ning, and J. Hicken
subroutine line_search_zoom(iLow, iHigh, dphi0, dphiLow, currentIteration, searchDirection,  &
     optActuationAmount)

  ! Internal modules
  use line_search

  implicit none

  ! Arguments
  integer, intent(in) :: iLow, iHigh
  real(WP), intent(in) :: dphi0
  real(WP), intent(in) :: dphiLow ! ... dphi that corrresponds to iLow
  integer, intent(inout) :: currentIteration
  real(WP), dimension(nControlParameters), intent(in) :: searchDirection
  real(WP), intent(out) :: optActuationAmount
  real(WP), parameter :: smallNumber = 1.0E-9_WP

  ! Local Variables
  integer :: i, j, iLow_, iHigh_
  real(WP) :: optCostFunctional, dphi, dphiLow_, costFunctionalGradient(nControlParameters)

  iLow_ = iLow; iHigh_ = iHigh; dphiLow_ = dphiLow

  j = currentIteration ! ... current line search iteration

  do
     ! Check if line searching is enough
     if (j .ge. nLineSearchIterations) then
        ! Determine actuation amount with smallest cost functional value
        optActuationAmount = actuationAmount(0)
        optCostFunctional = costFunctional(0)
        do i = 1, nLineSearchIterations
           if (costFunctional(i) .le. optCostFunctional) then
              optActuationAmount = actuationAmount(i)
              optCostFunctional = costFunctional(i)
           end if
        end do
        currentIteration = j
        return
     end if

     ! Check if the boundary points are too close
     if (abs(actuationAmount(iHigh_) - actuationAmount(iLow_)).lt.smallNumber) then
        optActuationAmount = actuationAmount(iLow_)
        currentIteration = j
        return
     end if

     j = j + 1

     ! Use the quadratic interpolation to estimate the actuation amount of the minimum,
     ! i.e. Eq. (2.27) of the new MDO notes
     actuationAmount(j) = 2.0_WP * actuationAmount(iLow_) * (costFunctional(iHigh_) -        &
          costFunctional(iLow_)) + dphiLow_ * (actuationAmount(iLow_)**2 -                   &
          actuationAmount(iHigh_)**2)
     actuationAmount(j) = 0.5_WP * actuationAmount(j) / (costFunctional(iHigh_) -            &
          costFunctional(iLow_) + dphiLow_ * (actuationAmount(iLow_) -                       &
          actuationAmount(iHigh_)))

     ! Modify actuationAmount(j) if is too close to the bounds
     if (abs(actuationAmount(j) - actuationAmount(iLow_)) .lt.                               &
          0.05_WP*abs(actuationAmount(iHigh_) - actuationAmount(iLow_)) .or.                 &
          abs(actuationAmount(j) - actuationAmount(iHigh_)) .lt.                             &
          0.05_WP*abs(actuationAmount(iHigh_) - actuationAmount(iLow_)))                     &
          actuationAmount(j) = 0.5_WP * (actuationAmount(iLow_) + actuationAmount(iHigh_))

     ! Evaluate cost functional
     costFunctional(j) = objective_function(actuationAmount = actuationAmount(j),            &
          searchDirection = searchDirection)

     if (costFunctional(j) .gt. costFunctional(0) + mu1 * actuationAmount(j) * dphi0 .or.    &
          costFunctional(j) .gt. costFunctional(iLow_)) then
        iHigh_ = j
     else
        ! Evaluate the gradient of cost functional
        costFunctionalGradient = objective_function_gradient(                                &
             actuationAmount = actuationAmount(j), searchDirection = searchDirection,        &
             costFunctional = costFunctional(j))
        dphi = sum( costFunctionalGradient * searchDirection ) ! ... directional derivative
        if (abs(dphi) .le. -mu2*dphi0 ) then ! ... dphi0 is negative!
           optActuationAmount = actuationAmount(j)
           currentIteration = j
           return
        elseif (dphi * (actuationAmount(iHigh_) - actuationAmount(iLow_)) .ge. 0.0_WP ) then
           iHigh_ = iLow_
        end if
        iLow_ = j
        dphiLow_ = dphi
     end if

  end do
end subroutine line_search_zoom
