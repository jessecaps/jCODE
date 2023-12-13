module les

  ! External modules
  use precision
  use simulation_flags

  implicit none

  ! Subgrid-scale models
  integer, parameter ::                                                                      &
       SMAGORINSKY_LILLY = 1,                                                                &
       GERMANO           = 3

  integer :: sgsModel
  real(WP) :: sgsCoefficient

contains

  ! Constant coefficient eddy viscosity model first published by  J. Smagorinsky (1963)
  ! -----------------------------------------------------------------------------------
  subroutine les_smagorinsky_lilly(density, velocityGradient, turbulentViscosity)

    ! External modules
    use geometry
    use grid

    implicit none

    ! Arguments
    real(WP), intent(in) :: density(:), velocityGradient(:,:)
    real(WP), intent(out) :: turbulentViscosity(:)

    ! Local variables
    integer :: i
    real(WP), dimension(nDimensions**2) :: Sij
    real(WP) :: mixingLength

    do i = 1, nGridPoints
       ! Get the mixing length scale
       mixingLength = sgsCoefficient * gridNorm(i, 1) ** (1.0_WP / real(nDimensions, WP))

       ! Get the rate-of-strain tensor
       select case (nDimensions)
       case (1)
          Sij(1) = velocityGradient(i,1)
       case (2)
          Sij(1) = velocityGradient(i,1)
          Sij(2) = 0.5_WP * (velocityGradient(i,2) + velocityGradient(i,3))
          Sij(3) = Sij(2)
          Sij(4) = velocityGradient(i,4)          
       case (3)
          Sij(1) = velocityGradient(i,1)
          Sij(2) = 0.5_WP * (velocityGradient(i,2) + velocityGradient(i,4))
          Sij(3) = 0.5_WP * (velocityGradient(i,3) + velocityGradient(i,7))
          Sij(4) = Sij(2)
          Sij(5) = velocityGradient(i,5)
          Sij(6) = 0.5_WP * (velocityGradient(i,6) + velocityGradient(i,8))
          Sij(7) = Sij(3)
          Sij(8) = Sij(6)
          Sij(9) = velocityGradient(i,9)
       end select

       ! Compute turbulent viscosity
       turbulentViscosity(i) = density(i) * mixingLength**2 * sqrt(2.0_WP * sum(Sij**2))
    end do

    return
  end subroutine les_smagorinsky_lilly

end module les


! =====================  !
!  Setup the LES routine !
! ====================== !
subroutine les_setup

  ! Internal modules
  use les

  ! External modules
  use string
  use parser

  implicit none

  ! Local variables
  character(len = str_medium) :: inputModel

  ! Return if not used
  if (.not. useLES) return

  call parser_read('subgrid-scale model', inputModel, 'smagorinsky')

  ! Choose the model
  select case (trim(inputModel))

  case ('smagorinsky', 'Smagorinsky', 'SMAGORINSKY')
     sgsModel = SMAGORINSKY_LILLY
     call parser_read('viscosity model coefficient', sgsCoefficient, 0.1_WP)

  case ('germano', 'Germano', 'GERMANO')
     sgsModel = GERMANO
     call die('les_setup: Germano SGS model not yet implemented!')

  case default
     call die("les_setup: Unknown model '" //trim(inputModel) // "'!")

  end select

  return
end subroutine les_setup


! =======================  !
!  Cleanup the LES routine !
! ======================== !
subroutine les_cleanup

  ! Internal modules
  use les

  implicit none

  ! Return if not used
  if (.not. useLES) return

  ! Nothing to do

  return
end subroutine les_cleanup


! =============================== !
! Compute the turbulent viscosity !
! =============================== !
subroutine les_compute_viscosity(density, velocityGradient, turbulentViscosity)

  ! Internal modules
  use les

  ! External modules
  use geometry

  implicit none

  ! Arguments
  real(WP), intent(in) :: density(nGridPoints), velocityGradient(nGridPoints, nDimensions**2)
  real(WP), intent(inout) :: turbulentViscosity(nGridPoints)

  ! Return if not used
  if (.not. useLES) return

  select case (sgsModel)
     
  case (SMAGORINSKY_LILLY)
     call les_smagorinsky_lilly(density, velocityGradient, turbulentViscosity)

  case (GERMANO)
     ! Not yet implemented

  end select

  return
end subroutine les_compute_viscosity


! ================ !
! Adjoint solution !
! ================ !
subroutine les_viscosity_adjoint

  ! Internal modules
  use les

  implicit none

  ! Return if not used
  if (.not. useLES) return

  ! Work on this later

  return
end subroutine les_viscosity_adjoint
