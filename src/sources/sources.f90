module sources

  ! External modules
  use precision

  implicit none

end module sources


! ============================== !
! Setup the various source terms !
! ============================== !
subroutine sources_setup

  ! Internal modules
  use sources

  implicit none

  call gravity_setup
  call acoustic_source_setup
  call ignition_source_setup
  call solenoidal_excitation_setup
  call fuel_source_setup
  call pressure_gradient_setup
  call momentum_source_setup
  call linear_forcing_setup

  return
end subroutine sources_setup


! ================================ !
! Cleanup the various source terms !
! ================================ !
subroutine sources_cleanup

  ! Internal modules
  use sources

  implicit none

  call gravity_cleanup
  call acoustic_source_cleanup
  call ignition_source_cleanup
  call solenoidal_excitation_cleanup
  call fuel_source_cleanup
  call pressure_gradient_cleanup
  call momentum_source_cleanup
  call linear_forcing_cleanup

  return
end subroutine sources_cleanup


! ============================================ !
! Add source terms during forward/adjoint runs !
! ============================================ !
subroutine add_sources(mode, source)

  ! External modules
  use solver_options
  use geometry
  use time_info

  implicit none

  ! Arguments
  integer, intent(in) :: mode
  real(WP), dimension(nGridPoints, nUnknowns), intent(inout) :: source

  ! Start the sources timer
  call timing_start('sources')

  select case (mode)

  case (FORWARD)

     call gravity_forward(source)

     call acoustic_source_forward(time, source)

     call fuel_source_forward(time, source)

     call ignition_source_forward(time, source)

     call solenoidal_excitation_forward(time, source)

     call pressure_gradient_forward(source)

     call momentum_source_forward(source)

     call linear_forcing_forward(source)

  case (ADJOINT)

     call gravity_adjoint(source)

     call momentum_source_adjoint(source)

     call pressure_gradient_adjoint(source)

     call linear_forcing_adjoint(source)

  end select

  ! Stop the sources timer
  call timing_stop('sources')

  return
end subroutine add_sources
