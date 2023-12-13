module monitor_ibm

  ! External modules
  use monitor
  use ibm

  implicit none

end module monitor_ibm


! ==================================== !
! Setup the routine to monitor the IBM !
! ==================================== !
subroutine monitor_ibm_setup(mode, controlIteration)

  ! Internal modules
  use monitor_ibm

  ! external modules
  use simulation_flags
  use solver_options
  use geometry

  implicit none

  ! Arguments
  integer, intent(in) :: mode, controlIteration

  if (.not. useIBM) return

  select case (mode)

  case (FORWARD)

     if (controlIteration .eq. 0) then

        ! -----------------------------------
        ! Monitor the baseline IBM prediction

        ! Monitor IBM forces
        call monitor_create('ibm_force', 12)
        call monitor_set_header(1, 'Force_x', 'r')
        call monitor_set_header(2, 'Force_y', 'r')
        call monitor_set_header(3, 'Force_z', 'r')
        call monitor_set_header(4, 'Pres_x', 'r')
        call monitor_set_header(5, 'Pres_y', 'r')
        call monitor_set_header(6, 'Pres_z', 'r')
        call monitor_set_header(7, 'Visc_x', 'r')
        call monitor_set_header(8, 'Visc_y', 'r')
        call monitor_set_header(9, 'Visc_z', 'r')
        call monitor_set_header(10, 'Torque_x', 'r')
        call monitor_set_header(11, 'Torque_y', 'r')
        call monitor_set_header(12, 'Torque_z', 'r')


        ! Monitor IBM objects
        if (ibm_move) then
           call monitor_create('ibm_object', 10)
           call monitor_set_header(1 , 'nObjects', 'i')
           call monitor_set_header(2 , 'minU', 'r')
           call monitor_set_header(3 , 'maxU', 'r')
           call monitor_set_header(4 , 'meanU', 'r')
           call monitor_set_header(5 , 'minV', 'r')
           call monitor_set_header(6 , 'maxV', 'r')
           call monitor_set_header(7 , 'meanV', 'r')
           call monitor_set_header(8 , 'minW', 'r')
           call monitor_set_header(9 , 'maxW', 'r')
           call monitor_set_header(10, 'meanW', 'r')
        end if
        
     else

        ! -------------------------------------
        ! Monitor the forward control iteration

     end if

  case (ADJOINT)

     if (controlIteration .eq. 0) then

        ! -------------------------------------
        ! Monitor the baseline adjoint solution

     else

        ! -------------------------------------
        ! Monitor the adjoint control iteration

     end if

  end select

  return
end subroutine monitor_ibm_setup


! ========================================================== !
! Compute various statistics of the IBM during a forward run !
! ========================================================== !
subroutine monitor_ibm_forward(controlIteration)

  ! Internal modules
  use monitor_ibm

  ! External modules
  use simulation_flags
  use parallel
  use solver_options
  use grid

  implicit none

  ! Arguments
  integer, intent(in) :: controlIteration

  ! Local variables
  integer :: i, j
  real(WP) :: pForce(3), vForce(3), meanVelocity(3), minVelocity(3), maxVelocity(3), torque(3)

  if (.not. useIBM) return

  if (controlIteration .eq. 0) then

     ! -------------------------------
     ! Monitor the baseline prediction

     ! Get forces acting on IBM objects computed on the grid
     call ibm_integrate_forces

     pForce = 0.0_WP; vForce = 0.0_WP; torque = 0.0_WP
     do i = 1, nObjects
        pForce(1:nDimensions) = pForce(1:nDimensions) + object(i)%pForce(1:nDimensions)
        vForce(1:nDimensions) = vForce(1:nDimensions) + object(i)%vForce(1:nDimensions)
        torque(1:3) = torque(1:3) + object(i)%hTorque(1:3)
     end do
     
     ! Set IBM force to monitor
     call monitor_select('ibm_force')
     call monitor_set_single_value(1, pForce(1)+vForce(1))
     call monitor_set_single_value(2, pForce(2)+vForce(2))
     call monitor_set_single_value(3, pForce(3)+vForce(3))
     call monitor_set_single_value(4, pForce(1))
     call monitor_set_single_value(5, pForce(2))
     call monitor_set_single_value(6, pForce(3))
     call monitor_set_single_value(7, vForce(1))
     call monitor_set_single_value(8, vForce(2))
     call monitor_set_single_value(9, vForce(3))
     call monitor_set_single_value(10, torque(1))
     call monitor_set_single_value(11, torque(2))
     call monitor_set_single_value(12, torque(3))

     ! Sum IBM object forces and velocities
     if (ibm_move) then
        meanVelocity = 0.0_WP;  minVelocity = 0.0_WP; maxVelocity = 0.0_WP
        minVelocity(1:nDimensions) = huge(1.0_WP)
        maxVelocity(1:nDimensions) = -huge(1.0_WP)
        do i = 1, nObjects
           do j = 1, nDimensions
              meanVelocity(j) = meanVelocity(j) + object(i)%velocity(j)
              minVelocity(j) = min(minVelocity(j),  object(i)%velocity(j))
              maxVelocity(j) = max(maxVelocity(j),  object(i)%velocity(j))
           end do
        end do
        meanVelocity = meanVelocity / real(nObjects, WP)

        ! Set IBM object stats to monitor
        call monitor_select('ibm_object')
        call monitor_set_single_value(1 , real(nObjects, WP))
        call monitor_set_single_value(2 , minVelocity(1))
        call monitor_set_single_value(3 , maxVelocity(1))
        call monitor_set_single_value(4 , meanVelocity(1))
        call monitor_set_single_value(5 , minVelocity(2))
        call monitor_set_single_value(6 , maxVelocity(2))
        call monitor_set_single_value(7 , meanVelocity(2))
        call monitor_set_single_value(8 , minVelocity(3))
        call monitor_set_single_value(9 , maxVelocity(3))
        call monitor_set_single_value(10, meanVelocity(3))
     end if

  else

     ! --------------------------------------
     ! Monitor the forward control iteration


  end if

  return
end subroutine monitor_ibm_forward
