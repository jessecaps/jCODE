module time_integrator

  ! External modules
  use simulation

  implicit none

  ! Time scheme
  integer :: TIME_SCHEME
  integer, parameter ::                                                                      &
       EULER           = 1,                                                                  &
       RK2             = 2,                                                                  &
       RK3_TVD         = 3,                                                                  &
       RK4             = 4

  ! Time integration parameters
  integer :: nStages
  real(WP), allocatable :: timeNorm(:)

  ! Buffer arrays
  real(WP), allocatable :: rightHandSide(:,:), buffer1(:,:), buffer2(:,:)

end module time_integrator


! ========================= !
! Setup the time integrator !
! ========================= !
subroutine time_integrator_setup

  ! Internal modules
  use time_integrator

  ! External modules
  use string
  use parser
  use simulation_flags

  implicit none

  ! Local parameters
  character(len = str_medium) :: val

  ! Allocate the right-handside array
  allocate(rightHandSide(nGridPoints, nUnknowns))
  
  ! Determine time integration scheme
  call parser_read('time integration scheme', val, 'rk4')
  select case (trim(val))
     
  case ('euler', 'EULER', 'Euler')
     ! First-order Euler scheme
     TIME_SCHEME = EULER
     nStages = 1

  case ('rk2', 'RK2')
     ! Second-order Runge-Kutta scheme
     TIME_SCHEME = RK2
     nStages = 2
     allocate(buffer1(nGridPoints, nUnknowns))

  case ('rk3-tvd', 'RK3-TVD')
     ! Optimal third-order TVD Runge-Kutta scheme
     TIME_SCHEME = RK3_TVD
     nStages = 3
     allocate(buffer1(nGridPoints, nUnknowns))

  case ('rk4', 'RK4')
     ! Fourth-order Runge-Kutta scheme
     TIME_SCHEME = RK4
     nStages = 4
     allocate(timeNorm(nStages))
     timeNorm = (/ 1.0_WP / 6.0_WP, 1.0_WP / 3.0_WP, 1.0_WP / 3.0_WP, 1.0_WP / 6.0_WP /)
     allocate(buffer1(nGridPoints, nUnknowns))
     allocate(buffer2(nGridPoints, nUnknowns))

  case default

     call die("Unknown time integration scheme: '" //  trim(val) // "'")

  end select

  ! Adjoint solver is only implemented for RK4
  if (TIME_SCHEME.ne.RK4 .and. .not.predictionOnly) then
     call die('time_integrator_setup: Adjoint is only implemented with RK4!')
  end if

  return
end subroutine time_integrator_setup


! =========================== !
! Cleanup the time integrator !
! =========================== !
subroutine time_integrator_cleanup

  ! Internal modules
  use time_integrator

  implicit none

  if (allocated(timeNorm)) deallocate(timeNorm)
  if (allocated(rightHandSide)) deallocate(rightHandSide)
  if (allocated(buffer1)) deallocate(buffer1)
  if (allocated(buffer2)) deallocate(buffer2)

  return
end subroutine time_integrator_cleanup


! ==================================== !
! Take a sub-step during a forward run !
! ==================================== !
subroutine substep_forward(stage)

  ! Internal modules
  use time_integrator

  implicit none

  ! Arguments
  integer, intent(in) :: stage

  select case (TIME_SCHEME)

  case (EULER)

     call substep_euler

  case (RK2)

     call substep_rk2(stage)

  case (RK3_TVD)

     call substep_rk3_tvd(stage)
     
  case (RK4)

     call substep_rk4(stage)

  end select

  ! Correct the state variables if necessary
  call correct_state(conservedVariables)

  return
end subroutine substep_forward


! ================= !
! First-order Euler !
! ================= !
subroutine substep_euler

  ! Internal modules
  use time_integrator

  ! External modules
  use solver_options
  use state

  implicit none
  
  ! Integrate the particles
  call particle_substep_euler

  ! Compute the right-hand side
  call state_rhs(FORWARD, rightHandSide)

  ! Integrate IBM objects
  call ibm_substep_euler

  ! Update the conserved variables
  conservedVariables = conservedVariables + timeStepSize * rightHandSide

  ! Update time
  time = time + timeStepSize

  return
end subroutine substep_euler


! ======================== !
! Second-order Runge-Kutta !
! ======================== !
subroutine substep_rk2(stage)

  ! Internal modules
  use time_integrator

  ! External modules
  use solver_options
  use state

  implicit none

  ! Arguments
  integer, intent(in) :: stage

  ! Local variables

  select case (stage)

  case (1)

     ! Store the conserved variables
     buffer1 = conservedVariables

     ! Integrate the particles
     call particle_substep_rk2(stage)

     ! Compute the right-hand side
     call state_rhs(FORWARD, rightHandSide)

     ! Integrate IBM objects
     !call ibm_substep_rk2(stage)

     ! Update the conserved variables
     conservedVariables = buffer1 + 0.5_WP * timeStepSize * rightHandSide

     ! Update time
     time = time + 0.5_WP * timeStepSize

  case (2)

     ! Integrate the particles
     call particle_substep_rk2(stage)

     ! Compute right-hand side
     call state_rhs(FORWARD, rightHandSide)

     ! Integrate IBM objects
     !call ibm_substep_rk2(stage)

     ! Update the conserved variables
     conservedVariables = buffer1 + timeStepSize * rightHandSide

     ! Update time
     time = time + 0.5_WP * timeStepSize

  end select

  return
end subroutine substep_rk2


! ============================ !
! Second-order TVD Runge-Kutta !
! ============================ !
subroutine substep_rk2_tvd(stage)

  ! Internal modules
  use time_integrator

  ! External modules
  use solver_options
  use state

  implicit none

  ! Arguments
  integer, intent(in) :: stage

  ! Local variables

  select case (stage)

  case (1)

     ! Store the conserved variables
     buffer1 = conservedVariables

     ! Integrate the particles
     !call particle_substep_rk2_tvd(stage)

     ! Compute the right-hand side
     call state_rhs(FORWARD, rightHandSide)

     ! Integrate IBM objects
     !call ibm_substep_rk2_tvd(stage)

     ! Update the conserved variables
     conservedVariables = buffer1 + timeStepSize * rightHandSide

     ! Update time
     time = time + timeStepSize

  case (2)

     ! Integrate the particles
     !call particle_substep_rk2_tvd(stage)

     ! Compute right-hand side
     call state_rhs(FORWARD, rightHandSide)

     ! Integrate IBM objects
     !call ibm_substep_rk2_tvd(stage)

     ! Update the conserved variables
     conservedVariables = 0.5_WP * (buffer1 + conservedVariables +                           &
          timeStepSize * rightHandSide)

  end select

  return
end subroutine substep_rk2_tvd


! ================================================================= !
! Optimal third-order total variation diminishing (TVD) Runge-Kutta !
! ================================================================= !
subroutine substep_rk3_tvd(stage)

  ! Internal modules
  use time_integrator

  ! External modules
  use solver_options
  use state

  implicit none

  ! Arguments
  integer, intent(in) :: stage

  select case (stage)

  case (1)

     ! Store the conserved variables
     buffer1 = conservedVariables

     ! Integrate the particles
     call particle_substep_rk3_tvd(stage)

     ! Compute the right-hand side
     call state_rhs(FORWARD, rightHandSide)

     ! Integrate IBM objects
     call ibm_substep_rk3_tvd(stage)

     ! Update the conserved variables
     conservedVariables = buffer1 + timeStepSize * rightHandSide

     ! Update time
     time = time + timeStepSize

  case (2)

     ! Integrate the particles
     call particle_substep_rk3_tvd(stage)

     ! Compute the right-hand side
     call state_rhs(FORWARD, rightHandSide)

     ! Integrate IBM objects
     call ibm_substep_rk3_tvd(stage)

     ! Update the conserved variables
     conservedVariables = 0.75_WP * buffer1 + 0.25_WP * (conservedVariables +                &
          timeStepSize * rightHandSide)

     ! Update time
     time = time - 0.5_WP * timeStepSize
     
  case (3)

     ! Integrate the particles
     call particle_substep_rk3_tvd(stage)

     ! Compute right-hand side
     call state_rhs(FORWARD, rightHandSide)

     ! Integrate IBM objects
     call ibm_substep_rk3_tvd(stage)

     ! Update the conserved variables
     conservedVariables = 1.0_WP / 3.0_WP * buffer1 + 2.0_WP / 3.0_WP *                      &
          (conservedVariables + timeStepSize * rightHandSide)

     ! Update time
     time = time + 0.5_WP * timeStepSize

  end select

  return
end subroutine substep_rk3_tvd


! ======================== !
! Fourth-order Runge-Kutta !
! ======================== !
subroutine substep_rk4(stage)

  ! Internal modules
  use time_integrator

  ! External modules
  use solver_options
  use state

  implicit none

  ! Arguments
  integer, intent(in) :: stage

  ! Local variables
  real(WP), parameter :: oneHalf = 1.0_WP / 2.0_WP
  real(WP), parameter :: oneThird = 1.0_WP / 3.0_WP
  real(WP), parameter :: oneSixth = 1.0_WP / 6.0_WP

  select case (stage)

  case (1)

     ! Store the conserved variables
     buffer1 = conservedVariables

     ! Integrate the particles
     call particle_substep_rk4(stage)

     ! Compute the right-hand side
     call state_rhs(FORWARD, rightHandSide)

     ! Integrate IBM objects
     call ibm_substep_rk4(stage)

     ! Update the conserved variables
     buffer2 = conservedVariables + timeStepSize * rightHandSide * oneSixth
     conservedVariables = buffer1 + timeStepSize * rightHandSide * oneHalf

  case (2)

     ! Update the time
     time = time + timeStepSize * oneHalf

     ! Integrate the particles
     call particle_substep_rk4(stage)

     ! Compute the right-hand side
     call state_rhs(FORWARD, rightHandSide)

     ! Integrate IBM objects
     call ibm_substep_rk4(stage)

     ! Update the conserved variables
     buffer2 = buffer2 + timeStepSize * rightHandSide * oneThird
     conservedVariables = buffer1 + timeStepSize * rightHandSide * oneHalf

  case (3)

     ! Integrate the particles
     call particle_substep_rk4(stage)

     ! Compute the right-hand side
     call state_rhs(FORWARD, rightHandSide)

     ! Integrate IBM objects
     call ibm_substep_rk4(stage)

     ! Update the conserved variables     
     buffer2 = buffer2 + timeStepSize * rightHandSide * oneThird
     conservedVariables = buffer1 + timeStepSize * rightHandSide

  case (4)

     ! Update time
     time = time + timeStepSize * oneHalf

     ! Integrate the particles
     call particle_substep_rk4(stage)

     ! Compute right-hand side
     call state_rhs(FORWARD, rightHandSide)

     ! Integrate IBM objects
     call ibm_substep_rk4(stage)

     ! Update the conserved variables
     conservedVariables = buffer2 + timeStepSize * rightHandSide * oneSixth

  end select

  return
end subroutine substep_rk4


! ===================================== !
! Take a sub-step during an adjoint run !
! Only implemented for RK4              !
! ===================================== !
subroutine substep_adjoint(stage)

  ! Internal modules
  use time_integrator

  ! External modules
  use functional, only : adjointForcingFactor, adjointCorrection

  implicit none

  ! Arguments
  integer, intent(in) :: stage

  ! Local variables
  real(WP), parameter :: oneHalf = 1.0_WP / 2.0_WP
  real(WP), parameter :: oneThird = 1.0_WP / 3.0_WP
  real(WP), parameter :: oneSixth = 1.0_WP / 6.0_WP

  select case (stage)

  case (4)

     ! Current time for adjoint coefficients
     adjointCoefficientTime = time - 0.5_WP * timeStepSize

     ! Store adjoint variables
     buffer1 = adjointVariables

     ! Compute the right-hand side
     adjointForcingFactor = 2.0_WP * adjointCorrection
     call state_rhs(ADJOINT, rightHandSide)

     ! Update the adjoint variables
     buffer2 = adjointVariables - timeStepSize * rightHandSide * oneSixth
     adjointVariables = buffer1 - timeStepSize * rightHandSide * oneHalf

  case (3)

     ! Update the current time
     time = time - timeStepSize * oneHalf
     adjointCoefficientTime = time

     ! Compute the right-hand side
     adjointForcingFactor = 1.0_WP * adjointCorrection
     call state_rhs(ADJOINT, rightHandSide)

     ! Update the adjoint variables
     buffer2 = buffer2 - timeStepSize * rightHandSide * oneThird
     adjointVariables = buffer1 - timeStepSize * rightHandSide * oneHalf

  case (2)

     ! Update the adjoint coefficient time
     adjointCoefficientTime = time - 0.5_WP * timeStepSize

     ! Compute the right-hand side
     adjointForcingFactor = 0.5_WP * adjointCorrection
     call state_rhs(ADJOINT, rightHandSide)

     ! Update the adjoint variables
     buffer2 = buffer2 - timeStepSize * rightHandSide * oneThird
     adjointVariables = buffer1 - timeStepSize * rightHandSide

  case (1)

     ! Update the current time
     time = time - timeStepSize * oneHalf
     adjointCoefficientTime = time

     ! Compute the right-hand side
     adjointForcingFactor = 1.0_WP * adjointCorrection
     call state_rhs(ADJOINT, rightHandSide)

     ! Update the adjoint variables
     adjointVariables = buffer2 - timeStepSize * rightHandSide * oneSixth

  end select

  return
end subroutine substep_adjoint
