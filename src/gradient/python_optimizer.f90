module python_optimizer

  ! External modules
  use gradient_search
  use controller

  implicit none

  ! Communication tags, assigned the same number with the python code
  integer, parameter ::                                                                      &
       TAG_DONE             = 0,                                                             &
       TAG_STOP             = 1,                                                             &
       TAG_ARGSIZE          = 2,                                                             &
       TAG_ARG              = 3,                                                             &
       TAG_CALL_OBJ         = 4,                                                             &
       TAG_RESULT_OBJ       = 5,                                                             &
       TAG_CALL_OBJ_GRAD    = 6,                                                             &
       TAG_RESULT_OBJ_GRAD  = 7,                                                             &
       TAG_CALL_CONS        = 8,                                                             &
       TAG_RESULT_CONS      = 9,                                                             &
       TAG_CALL_CONS_GRAD   = 10,                                                            &
       TAG_RESULT_CONS_GRAD = 11,                                                            &
       TAG_CALLBACK         = 12 

  ! Rank of the python proc in the python code
  integer, parameter :: PYTHON_SERVER_RANK = 1 

  ! Python library, assigned the same numbers with the pyhton code
  integer, parameter ::                                                                      &
       SCIPY = 0

  ! Optimizer method, assigned the same numbers with the pyhton code
  integer, parameter ::                                                                      &
       CG        = 0,                                                                        &
       NEWTON_CG = 1,                                                                        &
       TNC       = 2,                                                                        &
       BFGS_     = 3,                                                                        &
       L_BFGS_B  = 4,                                                                        &     
       SLSQP     = 5

  ! Options index, assigned the same number with the python code
  integer, parameter ::                                                                      &
       iLIBRARY    = 0,                                                                      &
       iOPT_METHOD = 1,                                                                      &
       iMAX_ITR    = 2,                                                                      &
       iTOLERANCE  = 3,                                                                      &
       iIS_MIN     = 4,                                                                      &
       iDISPLAY    = 5

end module python_optimizer


! ====================== !
!  Use Python libraries  !
! ====================== !
subroutine run_python_optimizer(baselineParameters)

  ! Internal modules
  use python_optimizer

  ! External modules
  use parser

  implicit none

  ! Arguments
  real(WP), dimension(nControlParameters), intent(inout) :: baselineParameters

  ! Local parameter
  integer, parameter :: nOptions = 20

  ! Local Variables
  integer :: pythonComm, pythonCommSize, pythonCommRank, pythonCommRoot
  integer :: i, j, index, ierror, pythonLibrary, optimizerMthod
  integer, dimension(MPI_STATUS_SIZE) :: status
  integer :: pythonTemp(1)
  character(len = str_medium) :: key
  real(WP) :: norm, costFunctional
  real(WP), dimension(:), allocatable :: options, parameterVector, costFunctionalGradient,   &
       constraintFunction, constraintGradientVector, previousValue
  real(WP), dimension(:,:), allocatable :: constraintGradient

  ! Point external comm to python comm
  pythonComm = externalComm
  pythonCommSize = externalCommSize
  pythonCommRank = externalCommRank
  pythonCommRoot = externalCommRoot

  if (irank.eq.iroot .and. pythonCommSize.eq.0) call die('run_python_optimizer: no proc was &
       &launched for Python program!')

  ! Initialize cost functional with a nonsense value
  if (minimizeCost) then
     costFunctional = huge(1.0_WP)
  else
     costFunctional = -huge(1.0_WP)
  end if

  ! Read in the optimizer library
  call parser_read('python library', key, 'scipy')
  select case(trim(key))
  case('scipy', 'SciPy')
     pythonLibrary = SCIPY
  case default
     call die('run_python_optimizer: Unknown python library: "' // trim(key) // '"')
  end select

  ! Read in the optimizer method
  call parser_read('optimzer method', key, 'CG')
  select case(trim(key))
  case('cg', 'CG')
     optimizerMthod = CG
  case('newton-cg', 'NEWTON-CG')
     optimizerMthod = NEWTON_CG
  case('tnc', 'TNC')
     optimizerMthod = TNC
  case('bfgs', 'BFGS')
     optimizerMthod = BFGS_
  case('l-bfgs-b', 'L-BFGS-B')
     optimizerMthod = L_BFGS_B
  case('slsqp', 'SLSQP')
     optimizerMthod = SLSQP     
  case default
     call die('run_python_optimizer: Unknown optimzer method: "' // trim(key) // '"')
  end select

  ! Allocate constriant vectors/matrix
  allocate(constraintFunction(nConstraint))
  allocate(constraintGradient(nConstraint,nControlParameters))
  allocate(constraintGradientVector(nConstraint*nControlParameters))

  ! Allocate cost functional gradient
  allocate(costFunctionalGradient(nControlParameters)); costFunctionalGradient = 0.0_WP

  ! Allocate, then initiailize previousValue with a nonsense value
  allocate(previousValue(nControlParameters)); previousValue = huge(1.0_WP)

  allocate(parameterVector(3*nControlParameters)) ! For compressing initial and boundary of parameters
  parameterVector(1:nControlParameters) = baselineParameters
  parameterVector(nControlParameters+1:2*nControlParameters) = parameterLowerBound
  parameterVector(2*nControlParameters+1:3*nControlParameters) = parameterUpperBound

  ! Gather different options for sending to python code
  allocate(options(0:nOptions))
  options = -1
  options(iLIBRARY) = real(pythonLibrary, WP)
  options(iOPT_METHOD) = real(optimizerMthod, WP)
  options(iMAX_ITR) = real(nIterations, WP)
  options(iTOLERANCE) = tolerance
  if (minimizeCost) then
     options(iIS_MIN) = 1.0_WP
  else
     options(iIS_MIN) = 0.0_WP
  end if
  options(iDISPLAY) = 1.0_WP ! Print results

  if (pythonCommRank.eq.pythonCommRoot) then
     ! Sent nOptions and options to python server
     call MPI_Send(nOptions, 1, MPI_INTEGER, PYTHON_SERVER_RANK, TAG_ARGSIZE, pythonComm,    &
          ierror)
     call MPI_SendRecv(options, nOptions, MPI_REAL_WP, PYTHON_SERVER_RANK, TAG_ARG,          &
          pythonTemp, 1, MPI_INTEGER, PYTHON_SERVER_RANK, MPI_ANY_TAG, pythonComm, status,   &
          ierror)
     ! Check received tag
     if (status(MPI_TAG).eq.TAG_STOP)                                                        &
          call die("run_python_optimizer: Unexpected stop was sent by Pyhton server!")
     if (status(MPI_TAG).ne.TAG_DONE) call die("run_python_optimizer: Expected tag DONE; &
          &Unepected tag was sent by Python server: " // char(status(MPI_TAG)) // "!")

     ! Send nControlParameters, initial and bounds of parameters to python server
     call MPI_Send(nControlParameters, 1, MPI_INTEGER, PYTHON_SERVER_RANK, TAG_ARGSIZE,      &
          pythonComm, ierror)
     call MPI_SendRecv(parameterVector, 3*nControlParameters, MPI_REAL_WP,PYTHON_SERVER_RANK,&
          TAG_ARG, pythonTemp, 1, MPI_INTEGER, PYTHON_SERVER_RANK, MPI_ANY_TAG, pythonComm,  &
          status, ierror)
     ! Check received tag
     if (status(MPI_TAG).eq.TAG_STOP)                                                        &
          call die("run_python_optimizer: Unexpected stop was sent by Pyhton server!")
     if (status(MPI_TAG).ne.TAG_DONE) call die("run_python_optimizer: Expected tag DONE; &
          &Unepected tag was sent by Python server: " // char(status(MPI_TAG)) // "!")

     ! Send nConstraint and their type to python server
     call MPI_Send(nConstraint, 1, MPI_INTEGER, PYTHON_SERVER_RANK, TAG_ARGSIZE, pythonComm, &
          ierror)
     call MPI_SendRecv(equalityConstraint, nConstraint, MPI_LOGICAL, PYTHON_SERVER_RANK,     &
          TAG_ARG, pythonTemp, 1, MPI_INTEGER, PYTHON_SERVER_RANK, MPI_ANY_TAG, pythonComm,  &
          status, ierror)
     ! Check received tag
     if (status(MPI_TAG).eq.TAG_STOP)                                                        &
          call die("run_python_optimizer: Unexpected stop was sent by Pyhton server!")
     if (status(MPI_TAG).ne.TAG_DONE) call die("run_python_optimizer: Expected tag DONE; &
          &Unepected tag was sent by Python server: " // char(status(MPI_TAG)) // "!")

     ! Request start optimization from python server
     i = nControlParameters + 1
     call MPI_Send(i, 1, MPI_INTEGER, PYTHON_SERVER_RANK, TAG_ARGSIZE, pythonComm, ierror)

  end if
  deallocate(options)
  deallocate(parameterVector)

  ! Get new parameters from python server, send back cost, gradient, and constrints
  do while (majorIteration .lt. nIterations)
     if (pythonCommRank.eq.pythonCommRoot)                                                   &
          call MPI_Recv(baselineParameters, nControlParameters, MPI_REAL_WP,                 &
          PYTHON_SERVER_RANK, MPI_ANY_TAG, pythonComm, status, ierror)

     ! Make comm slaves aware of the python result
     call parallel_bc(status)
     call parallel_bc(baselineParameters)

     ! Compare with the previous one
     norm = sum(abs(baselineParameters - previousValue))

     ! Check MPI tag
     select case(status(MPI_TAG))
     case (TAG_STOP)
        call die("run_python_optimizer: Unexpected stop was sent by Pyhton server!")

     case (TAG_DONE) ! Optimization has been finished!
        exit

     case(TAG_CALL_OBJ, TAG_CALLBACK) ! Calculate objcative function
        if (norm.gt.0.0_WP) then ! New control parameters, call objective function
           costFunctional = objective_function(scaledBaselineValue = baselineParameters)
           previousValue = baselineParameters
        end if
        if (status(MPI_TAG) .eq. TAG_CALLBACK) majorIteration = majorIteration + 1

        ! Send the cost functional result to Python server
        if (pythonCommRank.eq.pythonCommRoot)                                                &
             call MPI_Send(costFunctional, 1, MPI_REAL_WP, PYTHON_SERVER_RANK,               &
             TAG_RESULT_OBJ, pythonComm, ierror)

     case(TAG_CALL_OBJ_GRAD) ! Calculate gradient of cost
        if (norm.gt.0.0_WP) then
           ! New control parameters, a new forward run is also required for providing 
           ! checkpointing intermediate memory!
           costFunctionalGradient =                                                          &
                objective_function_gradient(scaledBaselineValue = baselineParameters)
           previousValue = baselineParameters
        else
           costFunctionalGradient =                                                          &
                objective_function_gradient(costFunctional = costFunctional)
        end if

        ! Send the gradient result to Python server
        if (pythonCommRank.eq.pythonCommRoot)                                                &
             call MPI_Send(costFunctionalGradient, nControlParameters, MPI_REAL_WP,          &
             PYTHON_SERVER_RANK, TAG_RESULT_OBJ_GRAD, pythonComm, ierror)

     case(TAG_CALL_CONS) ! Calculate constraints
        call controller_constraint(baselineParameters, constraintFunction)
        ! Send the constraint to Python server
        if (pythonCommRank.eq.pythonCommRoot)                                                &
             call MPI_Send(constraintFunction, nConstraint, MPI_REAL_WP, PYTHON_SERVER_RANK, &
             TAG_RESULT_CONS, pythonComm, ierror)

     case(TAG_CALL_CONS_GRAD) ! Calculate constraint gradient
        call controller_constraint_gradient(baselineParameters, constraintGradient)
        ! Substiture constraintGradient into constraintGradientVector
        do i = 1, nConstraint
           do j = 1, nControlParameters
              index = (i-1)*nControlParameters + j
              constraintGradientVector(index) = constraintGradient(i,j)
           end do
        end do
        ! Send the constraint gradient to Python server
        if (pythonCommRank.eq.pythonCommRoot)                                                &
             call MPI_Send(constraintGradientVector, nConstraint*nControlParameters,         &
             MPI_REAL_WP, PYTHON_SERVER_RANK, TAG_RESULT_CONS_GRAD, pythonComm, ierror)

     case default
        call die("run_python_optimizer: Unepected tag was sent by Python server: " //        &
             char(status(MPI_TAG)) // "!")
     end select

  end do

  ! Send a stop request to Python server
  i = nControlParameters + 1
  if (pythonCommRank.eq.pythonCommRoot)                                                      &
       call MPI_Send(i, 1, MPI_INTEGER, PYTHON_SERVER_RANK, TAG_STOP, pythonComm, ierror)

  ! Clean up
  deallocate(costFunctionalGradient)  
  deallocate(previousValue)
  deallocate(constraintFunction, constraintGradient, constraintGradientVector)

  return
end subroutine run_python_optimizer
