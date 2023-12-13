module ibm_integrator

  ! External modules
  use ibm

  implicit none

  ! RK4 parameters
  type(t_Object), allocatable :: objectBuffer1(:), objectBuffer2(:)

contains

  ! Rigid body motion - Oscillating cylinder
  ! ----------------------------------------
  subroutine rigid_body_motion

    ! External modules
    use math
    use geometry
    use time_info

    implicit none

    ! Local variables
    real(WP) :: y0, f, D, A, U

    ! Get object location and diameter
    U = 0.25_WP
    y0 = object(rigidBodyID)%position(2)
    D = (2.0_WP * real(nDimensions, WP) * object(rigidBodyID)%volume / pi)                   &
         ** (1.0_WP / real(nDimensions, WP))

    ! Update object (assumes Re=185, St=0.195)
    A = 0.2_WP * D
    f = 0.8_WP * 0.191_WP * U / D
    object(rigidBodyID)%position(2) = A * sin(twoPi * f * time)
    object(rigidBodyID)%velocity(2) = A * twoPi * f * cos(twoPi * f * time)
    call correct_object_position(object(rigidBodyID))

    ! Update levelset using IBM info
    call ibm_compute_levelset

    return
  end subroutine rigid_body_motion

end module ibm_integrator


! =============================== !
! Setup the IBM object integrator !
! =============================== !
subroutine ibm_integrator_setup

  ! Internal modules
  use ibm_integrator

  ! External modules
  use simulation_flags

  implicit none

  if (.not. useIBM) return
  if (.not. ibm_move) return

  return
end subroutine ibm_integrator_setup


! ================================= !
! Cleanup the IBM object integrator !
! ================================= !
subroutine ibm_integrator_cleanup

  ! Internal modules
  use ibm_integrator

  implicit none

  if (allocated(objectBuffer1)) deallocate(objectBuffer1)
  if (allocated(objectBuffer2)) deallocate(objectBuffer2)

  return
end subroutine ibm_integrator_cleanup


! ============================================ !
! First-order Euler scheme for the IBM objects !
! ============================================ !
subroutine ibm_substep_euler

  ! Internal modules
  use ibm_integrator

  ! External modules
  use math
  use simulation_flags
  use geometry
  use time_info

  implicit none

  ! Local variables
  integer :: i, p

  if (.not. useIBM) return
  if (.not. ibm_move) return
  
  ! Start the IBM timer
  call timing_start('ibm')

  ! Integrate in time analytically using rigid body motion
  if (ibm_rigid_motion) then
     call rigid_body_motion
     call timing_stop('ibm')
     return
  end if

  ! Compute the right-hand side
  call ibm_solver_rhs

  ! Update objects
  do p = 1, nObjects

     do i = 1, nDimensions

        ! Update object position
        object(p)%position(i) = object(p)%position(i) + timeStepSize * object(p)%velocity(i)

        ! Update object velocity
        object(p)%velocity(i) = object(p)%velocity(i) + timeStepSize * object(p)%dudt(i)

     end do

     do i = 1, 3
        ! Update object angular velocity
        object(p)%angularVelocity(i) = object(p)%angularVelocity(i) +                        &
             timeStepSize * object(p)%dwdt(i)
     end do
     
     ! Correct the position to take into account periodicity
     call correct_object_position(object(p))

     ! Check if object left the domain
     call check_object_bounds(object(p))

  end do

  ! Recycle
  call recycle_objects(object)

  ! Update levelset using IBM info
  call ibm_compute_levelset

  ! Stop the IBM timer
  call timing_stop('ibm')

  return
end subroutine ibm_substep_euler


! ================================================ !
! 2nd-order Runge-Kutta scheme for the IBM objects !
! ================================================ !
subroutine ibm_substep_rk2(stage)

  ! Internal modules
  use ibm_integrator

  ! External modules
  use math
  use simulation_flags
  use geometry
  use time_info

  implicit none

  ! Arguments
  integer, intent(in) :: stage

  ! Local variables
  integer :: i, p

  if (.not. useIBM) return
  if (.not. ibm_move) return

  ! Start the IBM timer
  call timing_start('ibm')

  ! Integrate in time analytically using rigid body motion
  if (ibm_rigid_motion) then
     call rigid_body_motion
     call timing_stop('ibm')
     return
  end if

  select case (stage)

  case (1)

     ! Prepare temporary object vectors
     if (allocated(objectBuffer1)) deallocate(objectBuffer1)
     allocate(objectBuffer1(nObjects))
     objectBuffer1 = object

     ! Compute the right-hand side
     call ibm_solver_rhs

     ! Update objects
     do p = 1, nObjects

        do i = 1, nDimensions

           ! Update object position
           object(p)%position(i) = objectBuffer1(p)%position(i) +                            &
                0.5_WP * timeStepSize * object(p)%velocity(i)

           ! Update object velocity
           object(p)%velocity(i) = objectBuffer1(p)%velocity(i) +                            &
                0.5_WP * timeStepSize * object(p)%dudt(i)

        end do

        do i = 1, 3
           ! Update object angular velocity
           object(p)%angularVelocity(i) = objectBuffer1(p)%angularVelocity(i) +              &
                0.5_WP * timeStepSize * object(p)%dwdt(i)
        end do
        
        ! Correct the position to take into account periodicity
        call correct_object_position(object(p))

        ! Check if object left the domain
        call check_object_bounds(object(p))

     end do

     ! Recycle
     call recycle_objects(object)

     ! Update levelset using IBM info
     call ibm_compute_levelset

  case (2)

     ! Compute the right-hand side
     call ibm_solver_rhs

     ! Update objects
     do p = 1, nObjects

        do i = 1, nDimensions

           ! Update object position
           object(p)%position(i) = objectBuffer1(p)%position(i) + timeStepSize *             &
                object(p)%velocity(i)

           ! Update object velocity
           object(p)%velocity(i) = objectBuffer1(p)%velocity(i) + timeStepSize *             &
                object(p)%dudt(i)

        end do

        do i = 1, 3
           ! Update object angular velocity
           object(p)%angularVelocity(i) = objectBuffer1(p)%angularVelocity(i) +              &
                timeStepSize * object(p)%dwdt(i)
        end do

        ! Correct the position to take into account periodicity
        call correct_object_position(object(p))

        ! Check if object left the domain
        call check_object_bounds(object(p))

     end do

     ! Recycle
     call recycle_objects(object)

     ! Update levelset using IBM info
     call ibm_compute_levelset

  end select

  ! Stop the IBM timer
  call timing_stop('ibm')

  return
end subroutine ibm_substep_rk2


! ============================================================ !
! Optimal 3rd-order TVD Runge-Kutta scheme for the IBM objects !
! ============================================================ !
subroutine ibm_substep_rk3_tvd(stage)

  ! Internal modules
  use ibm_integrator

  ! External modules
  use math
  use simulation_flags
  use geometry
  use time_info

  implicit none

  ! Arguments
  integer, intent(in) :: stage

  ! Local variables
  integer :: i, p

  if (.not. useIBM) return
  if (.not. ibm_move) return

  ! Start the IBM timer
  call timing_start('ibm')

  ! Integrate in time analytically using rigid body motion
  if (ibm_rigid_motion) then
     call rigid_body_motion
     call timing_stop('ibm')
     return
  end if

  select case (stage)

  case (1)

     ! Prepare temporary object vectors
     if (allocated(objectBuffer1)) deallocate(objectBuffer1)
     allocate(objectBuffer1(nObjects))
     objectBuffer1 = object

     ! Compute the right-hand side
     call ibm_solver_rhs

     ! Update objects
     do p = 1, nObjects

        do i = 1, nDimensions

           ! Update object position
           object(p)%position(i) = objectBuffer1(p)%position(i) +                            &
                timeStepSize * object(p)%velocity(i)

           ! Update object velocity
           object(p)%velocity(i) = objectBuffer1(p)%velocity(i) +                            &
                timeStepSize * object(p)%dudt(i)

        end do

        do i = 1, 3
           ! Update object angular velocity
           object(p)%angularVelocity(i) = objectBuffer1(p)%angularVelocity(i) +              &
                timeStepSize * object(p)%dwdt(i)
        end do
        
        ! Correct the position to take into account periodicity
        call correct_object_position(object(p))

        ! Check if object left the domain
        call check_object_bounds(object(p))

     end do

     ! Recycle
     call recycle_objects(object)

     ! Update levelset using IBM info
     call ibm_compute_levelset

  case (2)

     ! Compute the right-hand side
     call ibm_solver_rhs

     ! Update objects
     do p = 1, nObjects

        do i = 1, nDimensions

           ! Update object position
           object(p)%position(i) = 0.75_WP * objectBuffer1(p)%position(i) +                  &
                0.25_WP * (object(p)%position(i) + timeStepSize * object(p)%velocity(i))

           ! Update object velocity
           object(p)%velocity(i) = 0.75_WP * objectBuffer1(p)%velocity(i) +                  &
                0.25_WP * (object(p)%velocity(i) + timeStepSize * object(p)%dudt(i))

        end do

        do i = 1, 3
           ! Update object angular velocity
           object(p)%angularVelocity(i) = 0.75_WP * objectBuffer1(p)%angularVelocity(i) +    &
                0.25_WP * (object(p)%angularVelocity(i) + timeStepSize * object(p)%dwdt(i))
        end do

        ! Correct the position to take into account periodicity
        call correct_object_position(object(p))

        ! Check if object left the domain
        call check_object_bounds(object(p))

     end do

     ! Recycle
     call recycle_objects(object)

     ! Update levelset using IBM info
     call ibm_compute_levelset

  case (3)

     ! Compute the right-hand side
     call ibm_solver_rhs

     ! Update objects
     do p = 1, nObjects

        do i = 1, nDimensions

           ! Update object position
           object(p)%position(i) = 1.0_WP / 3.0_WP * objectBuffer1(p)%position(i) +          &
                2.0_WP / 3.0_WP * (object(p)%position(i) + timeStepSize *                    &
                object(p)%velocity(i))

           ! Update object velocity
           object(p)%velocity(i) = 1.0_WP / 3.0_WP * objectBuffer1(p)%velocity(i) +          &
                2.0_WP / 3.0_WP * (object(p)%velocity(i) + timeStepSize * object(p)%dudt(i))

        end do

        do i = 1, 3
           ! Update object angular velocity
           object(p)%angularVelocity(i) = 1.0_WP / 3.0_WP *                                  &
                objectBuffer1(p)%angularVelocity(i) + 2.0_WP / 3.0_WP *                      &
                (object(p)%angularVelocity(i) + timeStepSize * object(p)%dwdt(i))
        end do

        ! Correct the position to take into account periodicity
        call correct_object_position(object(p))

        ! Check if object left the domain
        call check_object_bounds(object(p))

     end do

     ! Recycle
     call recycle_objects(object)

     ! Update levelset using IBM info
     call ibm_compute_levelset

  end select

  ! Stop the IBM timer
  call timing_stop('ibm')

  return
end subroutine ibm_substep_rk3_tvd


! ================================================ !
! 4th-order Runge-Kutta scheme for the IBM objects !
! ================================================ !
subroutine ibm_substep_rk4(stage)

  ! Internal modules
  use ibm_integrator

  ! External modules
  use math
  use simulation_flags
  use geometry
  use time_info

  implicit none

  ! Arguments
  integer, intent(in) :: stage

  ! Local variables
  integer :: i, p
  real(WP), parameter :: oneHalf = 1.0_WP / 2.0_WP
  real(WP), parameter :: oneThird = 1.0_WP / 3.0_WP
  real(WP), parameter :: oneSixth = 1.0_WP / 6.0_WP

  if (.not. useIBM) return
  if (.not. ibm_move) return

  ! Start the IBM timer
  call timing_start('ibm')

  ! Integrate in time analytically using rigid body motion
  if (ibm_rigid_motion) then
     call rigid_body_motion
     call timing_stop('ibm')
     return
  end if

  select case (stage)

  case (1)

     ! Prepare temporary object vectors
     if (allocated(objectBuffer1)) deallocate(objectBuffer1)
     if (allocated(objectBuffer2)) deallocate(objectBuffer2)
     allocate(objectBuffer1(nObjects))
     allocate(objectBuffer2(nObjects))
     objectBuffer1 = object

     ! Compute the right-hand side
     call ibm_solver_rhs

     ! Update objects
     do p = 1, nObjects

        do i = 1, nDimensions

           ! Update object position
           objectBuffer2(p)%position(i) = object(p)%position(i) +                            &
                timeStepSize * object(p)%velocity(i) * oneSixth
           object(p)%position(i) = objectBuffer1(p)%position(i) +                            &
                timeStepSize * object(p)%velocity(i) * oneHalf

           ! Update object velocity
           objectBuffer2(p)%velocity(i) = object(p)%velocity(i) +                            &
                timeStepSize * object(p)%dudt(i) * oneSixth
           object(p)%velocity(i) = objectBuffer1(p)%velocity(i) +                            &
                timeStepSize * object(p)%dudt(i) * oneHalf

        end do

        do i = 1, 3
           ! Update object angular velocity
           objectBuffer2(p)%angularVelocity(i) = object(p)%angularVelocity(i) +              &
                timeStepSize * object(p)%dwdt(i) * oneSixth
           object(p)%angularVelocity(i) = objectBuffer1(p)%angularVelocity(i) +              &
                timeStepSize * object(p)%dwdt(i) * oneHalf
        end do
        
        ! Correct the position to take into account periodicity
        call correct_object_position(object(p))

        ! Check if object left the domain
        call check_object_bounds(object(p))

     end do

     ! Recycle
     call recycle_objects(object)

     ! Update levelset using IBM info
     call ibm_compute_levelset

  case (2)

     ! Compute the right-hand side
     call ibm_solver_rhs

     ! Update objects
     do p = 1, nObjects

        do i = 1, nDimensions

           ! Update object position
           objectBuffer2(p)%position(i) = objectBuffer2(p)%position(i) +                     &
                timeStepSize * object(p)%velocity(i) * oneThird
           object(p)%position(i) = objectBuffer1(p)%position(i) +                            &
                timeStepSize * object(p)%velocity(i) * oneHalf

           ! Update object velocity
           objectBuffer2(p)%velocity(i) = objectBuffer2(p)%velocity(i) +                     &
                timeStepSize * object(p)%dudt(i) * oneThird
           object(p)%velocity(i) = objectBuffer1(p)%velocity(i) +                            &
                timeStepSize * object(p)%dudt(i) * oneHalf

        end do

        do i = 1, 3
           ! Update object angular velocity
           objectBuffer2(p)%angularVelocity(i) = objectBuffer2(p)%angularVelocity(i) +       &
                timeStepSize * object(p)%dwdt(i) * oneThird
           object(p)%angularVelocity(i) = objectBuffer1(p)%angularVelocity(i) +              &
                timeStepSize * object(p)%dwdt(i) * oneHalf
        end do

        ! Correct the position to take into account periodicity
        call correct_object_position(object(p))

        ! Check if object left the domain
        call check_object_bounds(object(p))

     end do

     ! Recycle
     call recycle_objects(object)

     ! Update levelset using IBM info
     call ibm_compute_levelset

  case (3)

     ! Compute the right-hand side
     call ibm_solver_rhs

     ! Update objects
     do p = 1, nObjects

        do i = 1, nDimensions

           ! Update object position
           objectBuffer2(p)%position(i) = objectBuffer2(p)%position(i) +                     &
                timeStepSize * object(p)%velocity(i) * oneThird
           object(p)%position(i) = objectBuffer1(p)%position(i) +                            &
                timeStepSize * object(p)%velocity(i)

           ! Update object velocity
           objectBuffer2(p)%velocity(i) = objectBuffer2(p)%velocity(i) +                     &
                timeStepSize * object(p)%dudt(i) * oneThird
           object(p)%velocity(i) = objectBuffer1(p)%velocity(i) +                            &
                timeStepSize * object(p)%dudt(i)

        end do

        do i = 1, 3
           ! Update object angular velocity
           objectBuffer2(p)%angularVelocity(i) = objectBuffer2(p)%angularVelocity(i) +       &
                timeStepSize * object(p)%dwdt(i) * oneThird
           object(p)%angularVelocity(i) = objectBuffer1(p)%angularVelocity(i) +              &
                timeStepSize * object(p)%dwdt(i)
        end do

        ! Correct the position to take into account periodicity
        call correct_object_position(object(p))

        ! Check if object left the domain
        call check_object_bounds(object(p))

     end do

     ! Recycle
     call recycle_objects(object)

     ! Update levelset using IBM info
     call ibm_compute_levelset

  case (4)

     ! Compute the right-hand side
     call ibm_solver_rhs

     ! Update objects
     do p = 1, nObjects

        do i = 1, nDimensions

           ! Update object position
           object(p)%position(i) = objectBuffer2(p)%position(i) +                            &
                timeStepSize * object(p)%velocity(i) * oneSixth

           ! Update object velocity
           object(p)%velocity(i) = objectBuffer2(p)%velocity(i) +                            &
                timeStepSize * object(p)%dudt(i) * oneSixth

        end do

        do i = 1, 3
           ! Update object angular velocity
           object(p)%angularVelocity(i) = objectBuffer2(p)%angularVelocity(i) +              &
                timeStepSize * object(p)%dwdt(i) * oneSixth
        end do

        ! Correct the position to take into account periodicity
        call correct_object_position(object(p))

        ! Check if object left the domain
        call check_object_bounds(object(p))

     end do

     ! Recycle
     call recycle_objects(object)

     ! Update levelset using IBM info
     call ibm_compute_levelset

  end select

  ! Stop the IBM timer
  call timing_stop('ibm')

  return
end subroutine ibm_substep_rk4


! ============================================ !
! Compute the CFL based on the object velocity !
! ============================================ !
subroutine get_ibm_cfl(dt)

  ! Internal modules
  use ibm_integrator

  ! External modules
  use grid
  use state
  use time_info
  use ibm

  implicit none

  ! Arguments
  real(WP), intent(in) :: dt

  ! Local variables
  integer :: n

  if (.not.useIBM .or. .not.ibm_move) return

  ibmCfl = 0.0_WP

  do n = 1, nObjects
     ibmCfl = max(ibmCfl, sum(abs(object(n)%velocity(1:nDimensions))) * dt /   &
          (sqrt(real(nDimensions, WP)) * minGridSpacing))
  end do

  return
end subroutine get_ibm_cfl


! ====================================================== !
! Compute the timestep size based on the object velocity !
! ====================================================== !
subroutine get_ibm_timestep_size

  ! Internal modules
  use ibm_integrator

  ! External modules
  use grid
  use state
  use time_info
  use ibm

  implicit none

  ! Local variables
  integer :: n

  if (.not.useConstantCfl .or. .not.useIBM .or. .not.ibm_move) return

  do n = 1, nObjects
     timeStepSize = min(timeStepSize, inputCFL * sqrt(real(nDimensions, WP)) *               &
          minGridSpacing / (sum(abs(object(n)%velocity(1:nDimensions))) + epsilon(1.0_WP)))
  end do

  return
end subroutine get_ibm_timestep_size
