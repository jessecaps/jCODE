module monitor_time_info

  ! External modules
  use monitor

  implicit none

end module monitor_time_info


! ====================================== !
! Setup the routine to monitor time info !
! ====================================== !
subroutine monitor_time_info_setup(mode, controlIteration)

  ! Internal modules
  use monitor_time_info

  ! external modules
  use simulation_flags
  use solver_options
  use ibm

  implicit none

  ! Arguments
  integer, intent(in) :: mode, controlIteration

  ! Local variables
  integer :: n

  select case (mode)

  case (FORWARD)

     if (controlIteration .eq. 0) then

        ! -------------------------------
        ! Monitor the baseline prediction

        n = 4
        if (useParticles) n = n + 1
        if (useIBM .and. ibm_move) n = n + 1

        call monitor_create('timestep', n)
        call monitor_set_header(1, 'timestep', 'r')
        call monitor_set_header(2, 'max_cfl', 'r')
        call monitor_set_header(3, 'acoustic_cfl', 'r')
        call monitor_set_header(4, 'viscous_cfl', 'r')
        n = 4
        if (useParticles) then
           call monitor_set_header(n + 1, 'particle_cfl', 'r')
           n = n + 1
        end if
        if (useIBM .and. ibm_move) then
           call monitor_set_header(n + 1, 'ibm_cfl', 'r')
           n = n + 1
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
end subroutine monitor_time_info_setup


! ============================================= !
! Set the timestep and CFL during a forward run !
! ============================================= !
subroutine monitor_time_info_forward(controlIteration)

  ! Internal modules
  use monitor_particles

  ! External modules
  use parallel
  use simulation_flags
  use solver_options
  use time_info
  use ibm

  implicit none

  ! Arguments
  integer, intent(in) :: controlIteration

  ! Local variables
  integer :: n

  if (controlIteration .eq. 0) then

     ! -------------------------------
     ! Monitor the baseline prediction

     ! Set the values
     call monitor_select('timestep')
     call monitor_set_single_value(1, timeStepSize)
     call get_cfl
     call monitor_set_single_value(2, cfl)
     call monitor_set_single_value(3, acousticCfl)
     call monitor_set_single_value(4, viscousCfl)
     n = 4
     if (useParticles) then
        call monitor_set_single_value(n + 1, particleCfl)
        n = n + 1
     end if
     if (useIBM .and. ibm_move) then
        call monitor_set_single_value(n + 1, ibmCfl)
        n = n + 1
     end if

  else

     ! --------------------------------------
     ! Monitor the forward control iteration


  end if

  return
end subroutine monitor_time_info_forward
