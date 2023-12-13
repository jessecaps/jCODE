module combustion

  ! External modules
  use precision
  use solver_options

  implicit none

  integer, parameter, public ::                                                              &
       NO_CHEMISTRY     = 0,                                                                 &
       ONE_STEP         = 1,                                                                 &
       BOIVIN_SKELETAL  = 2

  integer :: chemistryModel, nReactions = 0
  logical :: wellStirredReactor

end module combustion


! =========================== !
! Setup the combustion module !
! =========================== !
subroutine combustion_setup

  ! Internal modules
  use combustion

  ! External modules
  use parser
  use simulation_flags

  implicit none

  ! Local variables
  character(len = str_medium) :: val

  if (nSpecies .eq. 0) return

  ! Combustion model
  call parser_read('combustion model', val, 'NONE')

  select case (trim(val))

  case ('NONE', 'none')
     ! No combustion
     chemistryModel = NO_CHEMISTRY
     nReactions = 0

  case ('ONE STEP', 'one step', 'ONE-STEP', 'one-step', 'onestep', 'ONESTEP')
     ! One-step irreversible reaction
     chemistryModel = ONE_STEP
     call onestep_setup

  case ('BOIVIN SKELETAL', 'Boivin skeletal', 'boivin skeletal', 'boivin-skeletal')
     ! Skeletal mechanism of Boivin (2011)
     ! `Reduced-Kinetic Mechanisms for Hydrogen and Syngas Combustion Including Autoignition`
     chemistryModel = BOIVIN_SKELETAL
     call boivin_setup

  case default
     call die("combustion_setup: unknown combustion model '" // trim(val) // "' !")

  end select

  ! Well-stirred reactor (homogeneous flows only)
  call parser_read('well stirred reactor', wellStirredReactor, .false.)
  if (wellStirredReactor) call well_stirred_setup

  return
end subroutine combustion_setup


! ============================= !
! Cleanup the combustion module !
! ============================= !
subroutine combustion_cleanup

  ! Internal modules
  use combustion

  implicit none

  call onestep_cleanup
  call boivin_cleanup
  call well_stirred_cleanup

  chemistryModel = NO_CHEMISTRY
  nReactions = 0
  wellStirredReactor = .false.

  return
end subroutine combustion_cleanup


subroutine combustion_source(mode, source)

  ! Internal modules
  use combustion

  ! External modules
  use solver_options
  use geometry

  implicit none

  ! Arguments
  integer, intent(in) :: mode
  real(WP), dimension(nGridPoints, nUnknowns), intent(inout) :: source

  if (nSpecies .eq. 0 .or. nReactions .eq. 0) return

  ! Start the combustion timer
  call timing_start('combustion')

  select case (mode)

  case (FORWARD)

     select case (chemistryModel)

     case (NO_CHEMISTRY)
        ! Nothing to do

     case (ONE_STEP)
        call onestep_forward(source)

     case (BOIVIN_SKELETAL)
        call boivin_forward(source)

     end select

     if (wellStirredReactor) call well_stirred_forward(source)

  case (ADJOINT)

     select case (chemistryModel)

     case (NO_CHEMISTRY)
        ! Nothing to do

     case (ONE_STEP)
        call onestep_adjoint(source)

     case (BOIVIN_SKELETAL)
        call boivin_adjoint(source)

     end select

     if (wellStirredReactor) call well_stirred_adjoint(source)

  end select

  ! Stop the combustion timer
  call timing_stop('combustion')

  return
end subroutine combustion_source
