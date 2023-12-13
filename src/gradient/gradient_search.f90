module gradient_search

  ! External modules
  use simulation_flags
  use controller
  use fileio

  implicit none

  ! Search direction methods
  integer, parameter ::                                                                      &
       CONJUGATE_GRADIENT = 0,                                                               &
       BFGS               = 1

  ! Beta Developer parameters
  integer, parameter ::                                                                      &
       FLETCHER_REEVES = 1,                                                                  &
       POLAK_RIBIERE   = 2,                                                                  &
       DAI_YUAN        = 3

  ! Global variables
  integer :: iunit, iostat, nIterations, majorIteration, nForward, nAdjoint, nResetGradient, &
       searchDirectionMethod, betaDeveloperCG, restartIteration
  real(WP) :: tolerance, costFunctionalScale, costFunctionalScaleInverse
  real(WP), dimension(:), allocatable :: individualSensitivities, parameterLowerBound,       &
       parameterUpperBound, parameterScale, parameterScaleInverse
  logical :: minimizeCost

contains

  ! Compute the objective function value for optimization
  ! -----------------------------------------------------
  function objective_function(actuationAmount, scaledBaselineValue, searchDirection)

    ! External module
    use solver

    implicit none

    ! Arguments
    real(WP), intent(in), optional :: actuationAmount
    real(WP), intent(in), optional :: scaledBaselineValue(nControlParameters)
    real(WP), intent(in), optional :: searchDirection(nControlParameters)

    ! Result
    real(WP) :: objective_function

    ! Local variables
    integer :: j
    real(WP) :: actuationAmount_, parameters(nControlParameters),                            &
         searchDirection_(nControlParameters)

    ! Set the actuation amount
    actuationAmount_ = 0.0_WP
    if (present(actuationAmount)) actuationAmount_ = actuationAmount

    ! Update baselineValue if a new one has been provided
    if (present(scaledBaselineValue)) baselineValue = scaledBaselineValue * parameterScale

    ! Set search direction
    searchDirection_ = 0.0_WP
    if (present(searchDirection)) searchDirection_ = searchDirection

    ! Set the parameters for dumping into the optimization history file and update 
    ! controlGradient if necessary
    parameters = baselineValue
    if (abs(actuationAmount_).gt.0.0_WP) then
       ! Suppose x as a scaled vector, X is the unscaled one with x = X / xbar
       ! x = x0 + alpha * p
       ! X = xbar * (x0 + alpha * p) = X0 + alpha * (xbar * p) 
       controlGradient = parameterScale * searchDirection_
       parameters = parameters + actuationAmount_ * controlGradient
    end if

    objective_function =                                                                     &
         run_forward(actuationAmount = actuationAmount_, controlIteration = nForward + 1)
    nForward = nForward + 1

    ! Dump the optimization history 
    if (iRank.eq.iRoot) then
       write(iunit,'(I7,I8,I10,1000(1X,SP,(SP,ES23.15E3)))') nForward,nAdjoint,              &
            majorIteration + 1, objective_function, (parameters(j), j = 1, nControlParameters)
       flush(iunit)
    end if

    ! Scale the function
    objective_function = objective_function * costFunctionalScaleInverse
    if (.not.minimizeCost) objective_function = - objective_function

    return
  end function objective_function


  ! Compute the objective function gradient for optimization
  ! --------------------------------------------------------
  function objective_function_gradient(actuationAmount, scaledBaselineValue, searchDirection, &
       costFunctional)

    ! External module
    use solver

    implicit none

    ! Arguments
    real(WP), intent(in), optional :: actuationAmount
    real(WP), intent(in), optional :: scaledBaselineValue(nControlParameters)
    real(WP), intent(in), optional :: searchDirection(nControlParameters)
    real(WP), intent(in), optional :: costFunctional

    ! Result
    real(WP), dimension(nControlParameters) :: objective_function_gradient

    ! Local variables
    integer :: j
    real(WP) :: actuationAmount_, parameters(nControlParameters), costFunctional_,          &
         searchDirection_(nControlParameters)

    ! Set the actuation amount
    actuationAmount_ = 0.0_WP
    if (present(actuationAmount)) actuationAmount_ = actuationAmount

    ! Update baselineValue if a new one has been provided
    if (present(scaledBaselineValue)) baselineValue = scaledBaselineValue * parameterScale

    ! Set search direction
    searchDirection_ = 0.0_WP
    if (present(searchDirection)) searchDirection_ = searchDirection

    ! Set the parameters for dumping into the optimization history file and update 
    ! controlGradient if necessary
    parameters = baselineValue
    if (abs(actuationAmount_).gt.0.0_WP) then
       ! Suppose x as a scaled vector, X is the unscaled one with x = X / xbar
       ! x = x0 + alpha * p
       ! X = xbar * (x0 + alpha * p) = X0 + alpha * (xbar * p) 
       controlGradient = parameterScale * searchDirection_
       parameters = parameters + actuationAmount_ * controlGradient
    end if

    ! Check if a new forward run is required for dumping data for checkpointing
    if (present(costFunctional)) then
       costFunctional_ = costFunctional * costFunctionalScale ! for dumping opt. history
    else
       costFunctional_ = objective_function(actuationAmount_) * costFunctionalScale
    end if

    individualSensitivities = run_adjoint(controlIteration = nAdjoint + 1)
    nAdjoint = nAdjoint + 1

    ! Dump the optimization history 
    if (iRank.eq.iRoot) then
       write(iunit,'(I7,I8,I10,1000(1X,SP,(SP,ES23.15E3)))') nForward,nAdjoint,              &
            majorIteration + 1, costFunctional_, (parameters(j), j = 1, nControlParameters), &
            (individualSensitivities(j), j = 1, nControlParameters)
       flush(iunit)
    end if

    ! Scale the gradient
    objective_function_gradient = individualSensitivities * parameterScale *                 &
         costFunctionalScaleInverse
    if (.not.minimizeCost) objective_function_gradient = - objective_function_gradient

    return  
  end function objective_function_gradient


  ! Read optimization history
  ! -------------------------
  subroutine read_optimization_history(restartIteration)  

    implicit none

    ! Arguments
    integer, intent(in) :: restartIteration

    ! Local variables
    integer :: i, j
    real(WP) :: rBuffer

    if (iRank .eq. iRoot) then
       read(iunit, *, iostat = iostat)
       i = 0
       do while (i.lt.restartIteration)
          read(iunit, *, iostat = iostat) j, rBuffer, rBuffer,         &
               rBuffer, rBuffer
          if (j .ne. -1) i = i + 1
       end do
    end if
    call parallel_bc(iostat)
    if (iostat .ne. 0) call die('read_optimization_hisrory: failed to &
         &read baseline cost function from file!')

    ! TODO : Broadcast read data

    return
  end subroutine read_optimization_history


  ! Compute the beta parameter for conjugate gradient
  ! From Multidisciplinary Design Optimization by J.R.R.A. Martins, A. Ning, and J. Hicken
  ! --------------------------------------------------------------------------------------
  function get_beta(gradX, gradXold, searchDirection)

    implicit none

    ! Arguments
    real(WP), dimension(:), intent(in) :: gradX, gradXold, searchDirection

    ! Result
    real(WP) :: get_beta

    ! Local variables
    integer :: j
    real(WP) :: dotProduct1, dotProduct2

    select case (betaDeveloperCG)

    case (FLETCHER_REEVES)
       dotProduct1 = 0.0_WP
       dotProduct2 = 0.0_WP
       do j = 1, nControlParameters
          dotProduct1 = dotProduct1 + gradX(j) * gradX(j)
          dotProduct2 = dotProduct2 + gradXold(j) * gradXold(j)
       end do
       get_beta = dotProduct1 / dotProduct2
       if (nResetGradient .eq. 0) get_beta = max(0.0_WP, get_beta)

    case (POLAK_RIBIERE)
       dotProduct1 = 0.0_WP
       dotProduct2 = 0.0_WP
       do j = 1, nControlParameters
          dotProduct1 = dotProduct1 + gradX(j) * (gradX(j) - gradXold(j))
          dotProduct2 = dotProduct2 + gradXold(j) * gradXold(j)
       end do
       get_beta = dotProduct1 / dotProduct2
       if (nResetGradient .eq. 0) get_beta = max(0.0_WP, get_beta)

    case (DAI_YUAN)
       dotProduct1 = 0.0_WP
       dotProduct2 = 0.0_WP
       do j = 1, nControlParameters
          dotProduct1 = dotProduct1 + gradX(j) * gradX(j)
          dotProduct2 = dotProduct2 + (gradX(j) - gradXold(j)) * searchDirection(j)
       end do
       get_beta = dotProduct1 / dotProduct2
       if (nResetGradient .eq. 0) get_beta = max(0.0_WP, get_beta)

    end select

    return
  end function get_beta

end module gradient_search


! ================================================ !
! Adjoint-based optimization: find optimal forcing !
! ================================================ !
subroutine find_optimal_forcing

  ! Internal modules
  use gradient_search

  ! External modules
  use parser

  implicit none

  ! Local variables
  integer :: j
  real(WP) :: baselineCostFunctional, previousCostFunctional, actuationAmount, epsilon_g,    &
       epsilon_a, epsilon_r, searchDirectionNorm
  real(WP), dimension(:), allocatable :: scaledBaselineValue, previousParameters,            &
       baselineGradient, previousGradient, searchDirection
  real(WP), dimension(:,:), allocatable :: hessianInverse
  character(len = str_medium) :: fileName, key

  if (spaceTimeGradient)                                                                     &
       call die('find_optimal_forcing: not yet implemented for space-time sensitivities')

  ! For maximizing problems, the cost functional and its gradient are multplied by -1 in 
  ! the objective_function and objective_functional_gradient funciton
  gradientDirection = 1 

  ! Check if it is a minimzing problem
  call parser_read('minimize cost functional', minimizeCost, .true.)

  ! Read in optimization parameters
  call parser_read('optimization tolerance', tolerance, 1.0E-6_WP)
  call parser_read('restart control iteration', restartIteration, 0)
  restartIteration = max(restartIteration, 0)
  call parser_read('restart forward evaluation', nForward, restartIteration)
  nForward = max(nForward, restartIteration)
  call parser_read('restart adjoint evaluation', nAdjoint, 0)
  nAdjoint = max(nAdjoint, 0)

  call parser_read('number of control iterations', nIterations)
  if (nIterations .lt. restartIteration) call die('find_optimal_forcing: number of control &
       &iterations must be greater than or equal to restart control iteration!')

  if (restartIteration.gt.0 .and. .not.isBaselineAvailable) call die('find_optimal_forcing: &
       &cannot restart with controlled prediction without baseline!')

  ! Set the filename to output optimization progress
  filename = 'cost_optimization.txt'
  if (iRank .eq. iRoot) then
     iunit = iopen()
     if (restartIteration .eq. 0 .and. .not. isBaselineAvailable) then
        open(iunit, file = trim(filename), action ='write', status = 'unknown',              &
             iostat = iostat)
        write(iunit, '(A7,A8,A10,1000A24)') 'Forward', 'Adjoint', 'Iteration',               &
             'Cost functional', (trim(sensitivityParameter(j)), j = 1, nControlParameters),  &
             ('dJ/d ' // trim(sensitivityParameter(j)), j = 1, nControlParameters)
     else
        open(iunit, file = trim(filename), action ='readwrite', status = 'old',              &
             position = 'rewind', iostat = iostat)
     end if
  end if

  call parallel_bc(iostat)
  if (iostat .ne. 0) call die('find_optimal_forcing: failed to open file for writing!')

  ! Read initial data for starting optimization, if it is available
  if (isBaselineAvailable) call read_optimization_history(restartIteration)

  ! Allocate sensitivity
  allocate(individualSensitivities(nControlParameters))
  individualSensitivities = 0.0_WP

  ! Allocate, initialize, and set lower and upper bounds and scales of the control parameters
  allocate(parameterLowerBound(nControlParameters), parameterUpperBound(nControlParameters))
  allocate(parameterScale(nControlParameters), parameterScaleInverse(nControlParameters))
  parameterLowerBound = - huge(1.0_WP)
  parameterUpperBound =   huge(1.0_WP)
  parameterScale = 1.0_WP
  call controller_parameter_bound(parameterLowerBound, parameterUpperBound, parameterScale)
  parameterScaleInverse = 1.0_WP / parameterScale

  ! Scale parameter bounds
  parameterLowerBound = parameterLowerBound * parameterScaleInverse
  parameterUpperBound = parameterUpperBound * parameterScaleInverse

  ! Allocate and initialize scaled control parameters, sensitivity, and search direciton vectors
  allocate(scaledBaselineValue(nControlParameters))
  scaledBaselineValue = baselineValue * parameterScaleInverse

  ! Set cost functional scale
  call parser_read('cost functional scale', costFunctionalScale, 1.0_WP)
  costFunctionalScaleInverse = 1.0_WP / costFunctionalScale

  majorIteration = restartIteration
  nForward = nForward - 1
  nAdjoint = nAdjoint - 1

  ! Read in the gradient search method type and call the appropriate routine
  call parser_read('optimization library', key, 'jcode')

  select case (trim(key))

  case ('snopt')
     call run_snopt(scaledBaselineValue)

  case ('python')
     call run_python_optimizer(scaledBaselineValue)

!!$  case('ignition boundary')
!!$     call run_ignition_boundary

  case default!based on Algorithm 4 of Multidisciplinary Design Optimization by Martins et al.

     ! Read in tolerance parameters
     call parser_read('gradient tolerance', epsilon_g, tolerance)
     call parser_read('absolute tolerance', epsilon_a, 1.0E-6_WP)
     call parser_read('relative tolerance', epsilon_r, 0.01_WP)

     ! Read in search direction and line search methods
     call parser_read('search direction method', key, 'conjugate gradient')
     select case (trim(key))
     case ('conjugate gradient')
        searchDirectionMethod = CONJUGATE_GRADIENT
        ! Determine method for computing beta
        call parser_read('conjugate gradient developer', key, 'fletcher reeves')
        select case (trim(key))

        case ('Fletcher Reeves', 'fletcher reeves')
           betaDeveloperCG = FLETCHER_REEVES

        case ('Polak Ribiere', 'polak ribiere')  ! Polak–Ribière
           betaDeveloperCG = POLAK_RIBIERE

        case ('Dai Yuan', 'dai yuan')
           betaDeveloperCG = DAI_YUAN

        case default
           call die("find_optimal_forcing: Unknown conjugate gradient developer: '" //         &
                trim(key) // "'")
        end select

     case ('BFGS', 'bfgs')
        searchDirectionMethod = BFGS

     case default
        call die("find_optimal_forcing: unknown search direction method '" // trim(key) //"'")

     end select

     ! Read in the interval for reset gradient
     ! If it is 1, means it is gradient descent method!
     ! If it is 0 and method is CG, means it is automatically reset whenever is required!
     call parser_read('reset gradient interval',nResetGradient, nControlParameters)
     nResetGradient = max(0,nResetGradient)

     ! Allocate and initialize scaled sensitivity, and search direciton vectors
     allocate(baselineGradient(nControlParameters)); baselineGradient = 0.0_WP
     allocate(searchDirection(nControlParameters)); searchDirection = 0.0_WP
     searchDirectionNorm =0.0_WP

     ! Allocate and initialize the inverse of hessian approximation
     if (searchDirectionMethod .eq. BFGS) then
        allocate(hessianInverse(nControlParameters,nControlParameters))
        hessianInverse = 0.0_WP
        do j = 1, nControlParameters
           hessianInverse(j,j) = 1.0_WP
        end do
     else
        allocate(hessianInverse(0,0))
     end if

     ! Find (or load from file) initial data for starting optimization
     if (.not. isBaselineAvailable) then
        ! Compute the baselineCostFunctional functional and its gradient
        baselineCostFunctional = objective_function()
        baselineGradient = objective_function_gradient(costFunctional=baselineCostFunctional)
     end if

     ! Allocate and initialize vectors for previous scaled control parameters and sensitivity
     allocate(previousParameters(nControlParameters)); previousParameters = 0.0_WP
     allocate(previousGradient(nControlParameters)); previousGradient = 0.0_WP
     ! Initialize previous cost functional with a nonsense value
     if (minimizeCost) then
        previousCostFunctional = huge(1.0_WP)
     else
        previousCostFunctional = -huge(1.0_WP)
     end if

     ! Initialize the line search subroutine
     call initialize_line_search

     if (majorIteration .ge. nIterations) return

     do 
        ! Compute a search direction
        call compute_seach_direction(scaledBaselineValue, previousParameters,                 &
             baselineGradient, previousGradient, searchDirectionNorm, searchDirection,       &
             hessianInverse)

        ! Perform line searching and find a step length
        call run_line_search(scaledBaselineValue, baselineCostFunctional, baselineGradient,   &
             searchDirection, actuationAmount)

        ! Save the current values before updating them
        previousParameters = scaledBaselineValue
        previousCostFunctional = baselineCostFunctional
        previousGradient = baselineGradient

        ! Update the control parameters
        scaledBaselineValue = scaledBaselineValue + actuationAmount * searchDirection

        ! Compute the current cost functional and its gradient
        baselineCostFunctional = objective_function(scaledBaselineValue = scaledBaselineValue)
        baselineGradient = objective_function_gradient(costFunctional=baselineCostFunctional)

        ! Update number of iterations
        majorIteration = majorIteration + 1

        ! Check if enough
        if ( majorIteration .ge. nIterations .or.                                            &
             (abs(baselineCostFunctional - previousCostFunctional) .le.                      &
             epsilon_a + epsilon_r * abs(previousCostFunctional) .and.                       &
             maxval(abs(baselineGradient)) .lt. tolerance) )                                 &
             exit

     end do

     ! Clean up
     deallocate(previousParameters)
     deallocate(baselineGradient, previousGradient, searchDirection)
     if (allocated(hessianInverse)) deallocate(hessianInverse)

  end select

  ! Dump optimizatiom summary and close
  if (iRank .eq. iRoot) then
     write(iunit, *) ''
     write(iunit, '(A28,I8)') 'Number of major iterations:', majorIteration
     write(iunit, '(A28,I8)') 'Number of forward runs:', nForward + 1 ! include the baseline
     write(iunit, '(A28,I8)') 'Number of adjoint runs:', nAdjoint + 1 ! include the baseline
     close(iclose(iunit))
  end if

  ! Clean up
  deallocate(individualSensitivities)
  deallocate(parameterLowerBound, parameterUpperBound)
  deallocate(parameterScale, parameterScaleInverse)
  deallocate(scaledBaselineValue)

  return
end subroutine find_optimal_forcing


! ============================ !
! Compute the search direction !
! ============================ !
! Based on Algorithm 2 of Multidisciplinary Design Optimization by J.R.R.A. Martins, 
! A. Ning, and J. Hicken
subroutine compute_seach_direction(currentParameter, previousParameter, currentGradient,     &
     previousGradient, searchDirectionNorm, searchDirection, hessianInverse)

  ! Internal modules
  use gradient_search

  implicit none

  ! Arguments
  real(WP), dimension(nControlParameters), intent(in) :: currentParameter, previousParameter,&
       currentGradient, previousGradient
  real(WP), intent(inout) :: searchDirectionNorm
  real(WP), dimension(nControlParameters), intent(inout) :: searchDirection
  real(WP), dimension(nControlParameters,nControlParameters), intent(inout) :: hessianInverse

  ! Local variables
  integer :: i, j
  real(WP) :: beta, sTyInv
  real(WP), dimension(:), allocatable :: s, y
  real(WP), dimension(:,:), allocatable :: A, B, C

  ! Initially the search direciton points towards the descent direction
  if (majorIteration .eq. restartIteration) then
     searchDirectionNorm = maxval(abs(currentGradient))
     searchDirection = - currentGradient / searchDirectionNorm ! ... normalized
     return
  end if

  select case (searchDirectionMethod)

  case (CONJUGATE_GRADIENT)
     ! Compute gradient actuator amount and control gradient
     if (nResetGradient.ne.0 .and. mod(majorIteration, nResetGradient).eq.0) then
        ! Reset the search direction
        beta = 0.0_WP
        searchDirectionNorm = maxval(abs(currentGradient))
        searchDirection = - currentGradient / searchDirectionNorm
     else
        beta = get_beta(currentGradient, previousGradient, searchDirection)
        searchDirection = - currentGradient + beta * searchDirectionNorm * searchDirection
        searchDirectionNorm = maxval(abs(searchDirection))
        searchDirection = searchDirection / searchDirectionNorm
     end if

  case (BFGS)     
     if (nResetGradient.ne.0 .and. mod(majorIteration, nResetGradient).eq.0) then
        ! Reset the search direction and hessian inverse approximation
        searchDirectionNorm = maxval(abs(currentGradient))
        searchDirection = - currentGradient / searchDirectionNorm
        hessianInverse = 0.0_WP
        do i = 1, nControlParameters
           hessianInverse(i,i) = 1.0_WP
        end do
     else
        ! Allocate and initialize useful vectors and matrices
        allocate(s(nControlParameters)); s = 0.0_WP
        allocate(y(nControlParameters)); y = 0.0_WP
        allocate(A(nControlParameters, nControlParameters)); A = 0.0_WP
        allocate(B(nControlParameters, nControlParameters)); B = 0.0_WP
        allocate(C(nControlParameters, nControlParameters)); C = 0.0_WP

        s = currentParameter - previousParameter
        y = currentGradient - previousGradient

        sTyInv = sum(s * y) ! ... s^T * y
        sTyInv = 1.0_WP / sTyInv

        do i = 1, nControlParameters
           do j = 1, nControlParameters
              A(i,j) = A(i,j) - s(i) * y(j) * sTyInv ! ... -s * y^T / s^T * y
              B(i,j) = B(i,j) - y(i) * s(j) * sTyInv ! ... -y * s^T / s^T * y
              C(i,j) = C(i,j) + s(i) * s(j) * sTyInv ! ...  s * s^T / s^T * y
              if (i.eq.j) then
                 A(i,j) = A(i,j) + 1.0_WP ! ... I - s * y^T / s^T * y
                 B(i,j) = B(i,j) + 1.0_WP ! ... I - y * s^T / s^T * y
              end if
           end do
        end do

        ! Update hessian inverse approximation
        hessianInverse = matmul(A, hessianInverse)
        hessianInverse = matmul(hessianInverse, B)
        hessianInverse = hessianInverse + C

        ! Compute the search direction and its norm
        searchDirection = - matmul(hessianInverse, currentGradient)
        searchDirectionNorm = maxval(abs(searchDirection))

        deallocate(s, y, A, B, C)
     end if

  end select

  return
end subroutine compute_seach_direction
