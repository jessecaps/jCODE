module snopt

  ! External modules
  use gradient_search
  use controller

  implicit none

  ! Global variables
  real(WP) :: costFunctional
  real(WP), allocatable, dimension(:) :: previousValue

contains

  ! Compute function and its derivatives for snOpt
  ! ----------------------------------------------
  subroutine usrfun(Status, n, x, needf, nF, f, needG, lenG, G, cu, lencu, iu, leniu, ru,    &
       lenru)

    implicit none

    ! Arguments
    integer, intent(in) :: lencu, lenG, leniu, lenru, n, needf, needG, nF, iu(leniu)
    integer, intent(inout) :: Status
    character(len = str_short), intent(in) :: cu(lencu)
    real(WP), intent(in) :: ru(lenru), x(n)
    real(WP), intent(out) :: f(nF), G(lenG)

    ! Local variables
    integer :: Out, i, j, index
    real(WP) :: norm
    real(WP) :: fCons(nConstraint)
    real(WP) :: gradCons(nConstraint, nControlParameters)

    if (majorIteration .ge. nIterations)                                                     &
         call die('snop.f90: maximum number of iterations was exceeded!')

    ! Output unit number
    Out = 15

    if (Status .ge. 2) then
       return
    end if

    ! Compare with the previous one
    norm = sum(abs(x - previousValue))

    if (needF .gt. 0) then
       if (norm.gt.0.0_WP) then ! New control parameters, call objective function
          costFunctional = objective_function(scaledBaselineValue = x)
          previousValue = x
          majorIteration = majorIteration + 1
          norm = 0.0_WP
       end if

       ! f(1) was set as objective functional in run_snopt
       f(1) = costFunctional

       ! Calculate the constraint function
       call controller_constraint(fCons)
       do i = 2, nF
          f(i) = fCons(i-1)
       end do
    end if

    if (needG .gt. 0) then
       if (norm.gt.0.0_WP) then
          ! New baseline values, a new forward run is also required for providing
          ! checkpointing intermediate memory!
          G(1:n) = objective_function_gradient(scaledBaselineValue = x)
          previousValue = x
          majorIteration = majorIteration + 1
          f(1) = costFunctional
          ! Calculate the constraint function
          call controller_constraint(fCons)
          do i = 2, nF
             f(i) = fCons(i-1)
          end do
       else
          G(1:n) = objective_function_gradient(costFunctional = costFunctional)
       end if

       ! Calculate the constraint gradient
       call controller_constraint_gradient(gradCons)
       ! Substiture gradCons into G
       do i = 2, nF
          do j = 1, n
             index = (i-1)*n + j
             ! iGfun(index) = i was set in the run_snopt soubruotine
             ! jGvar(index) = j was set in the run_snopt soubruotine
             G(index) = gradCons(i-1,j)
          end do
       end do
    end if

    return
  end subroutine usrfun

end module snopt


! =========================================================================================!
!  Use SNOPT library                                                                       !
!  [1] Philip, E., Murray, W. and Saunders, M.A., 2008. User’s Guide for SNOPT Version 7:  !
!    Software for Large-Scale Nonlinear Programming.                                       !
! =========================================================================================!
subroutine run_snopt(scaledBaselineValue)

  ! Internal modules
  use snopt

  ! External modules
  use parser

  implicit none

  ! Arguments
  real(WP), dimension(nControlParameters), intent(inout) :: scaledBaselineValue

#ifdef USE_SNOPT
  ! Local Variables
  integer :: i, j, index, INFO, lenA, lencu, lenG, leniu, lenru, lencw, leniw, lenrw, mincw, &
       miniw, minrw, n, neA, neG, nF, nFname, nInf, nS, nxname, ObjRow, Start, iPrint, iSumm,&
       iSpecs, ne
  integer, dimension(:), allocatable :: iw,  iAfun, iGfun, iu, jAvar, jGvar, xstate, Fstate
  character(len = str_short):: Prob
  character(len = str_short), dimension(:), allocatable :: cw, cu, Fnames, xnames
  real(WP) :: ObjAdd, sInf
  real(WP), dimension(:), allocatable :: rw, A, F, Fmul, Flow, Fupp, ru, x, xlow, xmul, xupp

  ! Check if it is a minimzing problem
  !if (.not.minimizeCost)                                                                     &
  !     call die('snopt.f90: maximization problems must be set in the snopt specific file!') 

  ! Initialize cost functional with a nonsense value
  if (minimizeCost) then
     costFunctional = huge(1.0_WP)
  else
     costFunctional = -huge(1.0_WP)
  end if

  ! Allocate, then initiailize previousValue with a nonsense value
  allocate(previousValue(nControlParameters))
  previousValue = huge(1.0_WP)

  ! Set input values for snOptA
  Start = 0 ! "(Cold start) requests that the CRASH procedure be used, unless an  
  ! Old basis, Insert, or Load file is specified [1]."
  ! Start = 1 ! "is the same as 0 but more meaningful when a basis file is given [1]."
  ! Start = 2 ! "(Warm start) means that xstate and Fstate define a valid starting
  ! point (perhaps from an earlier call, though not necessarily) [1]."
  nF = 1 + nConstraint ! number of problem functions
  ! Set number of independent parameters
  n = nControlParameters

  nxname = 1 ! # xnames' components, =1 means names are not used
  allocate(xnames(nxname)) ! names of the variables
  xnames = ''
  nFname = 1 ! # Fnames' components, = 1 means names are not used
  allocate(Fnames(nFname))
  Fnames = '' ! names of the problem functions

  ObjAdd = 0.0_WP ! "is a constant that will be added to the objective row F(Objrow)
  ! for printing puposes. Typically,ObjAdd=0.0d+0 [1]."
  ObjRow = 1 ! "says which row of F(x) is to act as the objective function [1]."
  Prob = '' ! "is an 8-character name for the problem. Prob is used in the printed solution and in
  ! some routines that output basis files. A blank name may be used. [1]"

  ! usrfun  ! "(Section 3.6) is the name of a subroutine that calculates the nonlinear portion f(x)
  ! of the vector of problem functions F (x) = f (x) + Ax, and (optionally) its
  ! Jacobian G(x) for a given vector x. usrfun must be declared external in the routine that
  ! calls snOptA. [1]"

  lenA = 1  ! "is the dimension of the arrays iAfun, jAvar, A that hold (i, j, Aij ). (lenA>=1) [1]"
  neA  = 0  ! "is the number of nonzeros in A such that F(x) = f(x) + Ax. (0<=neA<=lenA) [1]"
  ! "iAfun(lenA), jAvar(lenA), A(lenA) define the coordinates (i, j) and values Aij of the
  ! nonzero elements of the linear part A of the function F (x) = f (x) + Ax. The first
  ! neA entries are valid. (0 <= neA <= lenA and lenA >=  1) [1]"
  allocate(iAfun(lenA), jAvar(lenA), A(lenA))
  A = 0; iAfun = 0; jAvar = 0;

  lenG = n * nF  ! "is the dimension of the coordinate arrays iGfun, jGvar. (lenG >=  1) [1]"
  neG = lenG ! "is the number of nonzero entries in G. (neG >=  0) [1]"
  ! "iGfun(lenG), jGvar(lenG) define the coordinates (i, j) of Gij , the nonzero elements of
  ! the nonlinear part of the derivatives G(x) + A of the function F (x) = f (x) + Ax.
  ! The first neG entries are valid. (0 <= neG <= lenG and lenG >=  1) [1]"
  allocate(iGfun(lenG), jGvar(lenG))
  do i = 1, nF
     do j = 1, n
        index = (i-1)*n + j
        iGfun(index) = i 
        jGvar(index) = j
     end do
  end do

  ! lower and upper bounds of the parameters
  allocate(xlow(nControlParameters), xupp(nControlParameters)) 
  xlow = parameterLowerBound
  xupp = parameterUpperBound

  allocate(Flow(nF), Fupp(nF)) ! lower and upper bounds of F
  Flow(1) = 0.0_WP;   Fupp(1) = huge(1.0_WP)
  do i = 2, nF
     Flow(i) = 0.0_WP;   Fupp(i) = 0.0_WP
  end do

  allocate(x(n)) ! "usually contains a set of initial values for x [1]."
  allocate(xstate(n)) ! "usually contains a set of initial states for each variable x [1]."
  x = scaledBaselineValue
  xstate = 0

  ne = 2*n
  allocate(xmul(ne)) ! On exit: "is the vector of dual variables for the bound constraints ... [1]"

  allocate(F(nF)) ! "sometimes contains a set of initial values for the functions F (x)."
  F(1) = Fupp(1)
  do i = 2, nF
     F(i) = Flow(i)
  end do
  allocate(Fstate(nF)) ! "sometimes contains a set of initial states for the problem functions F [1]."
  Fstate = 0
  allocate(Fmul(nF)) ! "contains an estimate of \lambda, the vector of Lagrange multipliers (shadow prices) for
  ! the constraints lF <= F (x) <= uF [1]."
  Fmul = 0.0_WP

  ! Arrays of workspace
  leniw = 30000 !leniw = 500 + n + nF
  lencw = 500 + n + nF
  lenrw = 40000 !500 + n + nF
  allocate(iw(leniw))
  allocate(cw(lencw))
  allocate(rw(lenrw))

  ! Arrays of user workspace
  leniu = leniw
  lencu = lencw
  lenru = lenrw
  allocate(iu(leniu))
  allocate(cu(lencu))
  allocate(ru(lenru))

  ! Initialize snopt
  iPrint = 9
  iSumm = 6
  call snInit(iPrint, iSumm, cw, lencw, iw, leniw, rw, lenrw)

  iSpecs = 4
  call snSpec(iSpecs, INFO, cw, lencw, iw, leniw, rw, lenrw)

  call snOptA(Start, nF, n, nxname, nFname, ObjAdd, ObjRow, Prob, usrfun, iAfun, jAvar, lenA,&
       neA, A, iGfun, jGvar, lenG, neG, xlow, xupp, xnames, Flow, Fupp, Fnames, x, xstate,   &
       xmul, F, Fstate, Fmul, INFO, mincw, miniw, minrw, nS, nInf, sInf, cu, lencu, iu,      &
       leniu, ru, lenru, cw, lencw, iw, leniw, rw, lenrw)

  ! Clean up
  deallocate(previousValue)
  deallocate(iAfun, jAvar, A, iGfun, jGvar, xlow, xupp, Flow, Fupp, xnames, Fnames, x,       &
       xstate, F, Fstate, Fmul, cu, iu, ru, cw, iw, rw, xmul)
#endif

  return
end subroutine run_snopt
