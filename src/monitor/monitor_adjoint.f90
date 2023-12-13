module monitor_adjoint

  ! External modules
  use monitor
  use functional
  use controller

  implicit none

end module monitor_adjoint


! ============================================= !
! Setup the routine to monitor the adjoint data !
! ============================================= !
subroutine monitor_adjoint_setup(mode, controlIteration)

  ! Internal modules
  use monitor_adjoint

  ! External modules
  use string
  use simulation_flags
  use solver_options

  implicit none

  ! Arguments
  integer, intent(in) :: mode, controlIteration

  ! Local variables
  integer :: i
  character(len = str_medium) :: filename

  if (predictionOnly) return

  select case (mode)

  case (FORWARD)

     ! Set the name
     filename = 'functional'
     if (controlIteration .gt. 0) write(filename, '(A,I2.2)') 'functional_', controlIteration

     ! Monitor the cost functional
     call monitor_create(trim(filename), 3)
     call monitor_set_header(1, 'integrand', 'r')
     call monitor_set_header(2, 'functional', 'r')
     call monitor_set_header(3, 'auxiliary', 'r')

  case (ADJOINT)

     ! Set the name
     filename = 'adjoint_sensitivity'
     if (controlIteration .gt. 0) write(filename, '(A,I2.2)') 'adjoint_sensitivity_',        &
          controlIteration

     ! Monitor the cost sensitivity
     call monitor_create(trim(filename), 2 * nControlParameters)
     do i = 1, nControlParameters
        call monitor_set_header(i, trim(sensitivityParameter(i)) // ' [t]', 'r')
     end do
     do i = 1, nControlParameters
        call monitor_set_header(nControlParameters + i, trim(sensitivityParameter(i)), 'r')
     end do

  end select

  return
end subroutine monitor_adjoint_setup


! =========================== !
! Monitor the cost functional !
! =========================== !
subroutine monitor_adjoint_functional(controlIteration)

  ! Internal modules
  use monitor_adjoint

  ! External modules
  use simulation_flags

  implicit none

  ! Arguments
  integer, intent(in) :: controlIteration

  ! Local variables
  character(len = str_medium) :: filename

  if (predictionOnly) return

  ! Set the name
  filename = 'functional'
  if (controlIteration .gt. 0) write(filename, '(A,I2.2)') 'functional_', controlIteration

  ! Monitor the cost functional
  call monitor_select(trim(filename))
  call monitor_set_single_value(1, instantaneousCostFunctional)
  call monitor_set_single_value(2, currentCostFunctional)
  call monitor_set_single_value(3, auxilaryCostFunctional)

  return
end subroutine monitor_adjoint_functional



! ============================ !
! Monitor the cost sensitivity !
! ============================ !
subroutine monitor_adjoint_sensitivity(controlIteration)

  ! Internal modules
  use monitor_adjoint

  implicit none

  ! Arguments
  integer, intent(in) :: controlIteration

  ! Local variables
  integer :: i
  character(len = str_medium) :: filename

  ! Set the name
  filename = 'adjoint_sensitivity'
  if (controlIteration .gt. 0) write(filename, '(A,I2.2)') 'adjoint_sensitivity_',           &
       controlIteration

  ! Monitor the cost sensitivity
  call monitor_select(trim(filename))
  do i = 1, nControlParameters
     call monitor_set_single_value(i, instantaneousSensitivity(i))
  end do
  do i = 1, nControlParameters
     call monitor_set_single_value(nControlParameters + i, currentSensitivity(i))
  end do

  return
end subroutine monitor_adjoint_sensitivity
