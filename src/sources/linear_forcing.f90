! ================================================================================ !
! Linear forcing of homogeneous isotropic turbulence of Lundgren, T. S.,           !
! “Linearly forced isotropic turbulence,” CTR Annual Research Briefs,              !
! Center for Turbulence Research, CTR, Stanford University, 2003, pp. 461–473.     !
!                                                                                  !
! This forcing routine should be modified for compressible flow, see e.g.,         !
!                                                                                  !
! Petersen, M. R., & Livescu, D. (2010). Forcing for statistically                 !
! stationary compressible isotropic turbulence. Physics of Fluids, 22(11), 116101. !
! ================================================================================ !
module linear_forcing

  ! External modules
  use precision

  implicit none

  real(WP) :: forcingCoefficient
  logical :: useLinearForcing

end module linear_forcing


! ======================== !
! Setup the linear forcing !
! ======================== !
subroutine linear_forcing_setup

  ! Internal modules
  use linear_forcing

  ! External modules
  use parser
  use geometry

  implicit none

  ! Only applicable for periodic flows
  useLinearForcing = .false.
  if (.not. all(isPeriodic(1:nDimensions))) return

  call parser_read('use linear forcing', useLinearForcing, .false.)
  if (useLinearForcing) then
     call parser_read('forcing coefficient', forcingCoefficient)
     if (forcingCoefficient .lt. epsilon(1.0_WP)) useLinearForcing = .false.
  end if

  return
end subroutine linear_forcing_setup


! ========================== !
! Cleanup the linear forcing !
! ========================== !
subroutine linear_forcing_cleanup

  ! Internal modules
  use linear_forcing

  implicit none

  return
end subroutine linear_forcing_cleanup


! ================================================= !
! Add the the linear forcing during the forward run !
! ================================================= !
subroutine linear_forcing_forward(source)

  ! Internal modules
  use linear_forcing

  ! External modules
  use parallel
  use geometry
  use grid
  use state

  implicit none

  ! Arguments
  real(WP), dimension(nGridPoints, nUnknowns), intent(inout) :: source

  ! Local variables
  integer :: i
  real(WP), dimension(nDimensions) :: meanMomentum, momentum

  if (.not. useLinearForcing) return

  ! Ensure mean momentum is zero
  if (predictionOnly) then
     do i = 1, nDimensions
        meanMomentum(i) = sum(conservedVariables(:,i+1) * gridNorm(:,1))
        call parallel_sum(meanMomentum(i))
        meanMomentum(i) = meanMomentum(i) / globalGridVolume
     end do
  else
     meanMomentum = 0.0_WP
  end if
     
  ! Apply linear forcing
  do i = 1, nGridPoints
     momentum = conservedVariables(i, 2:nDimensions+1) - meanMomentum
     source(i, 2:nDimensions+1) = source(i, 2:nDimensions+1) + forcingCoefficient *          &
          momentum
     source(i, nDimensions+2) = source(i, nDimensions+2) + forcingCoefficient *              &
          sum(velocity(i,:) * momentum)
  end do
 
  return
end subroutine linear_forcing_forward


! ================================================= !
! Add the the linear forcing during the adjoint run !
! ================================================= !
subroutine linear_forcing_adjoint(source)

  ! Internal modules
  use linear_forcing

  ! External modules
  use geometry
  use state

  implicit none

  ! Arguments
  real(WP), dimension(nGridPoints, nUnknowns), intent(inout) :: source

  if (.not. useLinearForcing) return

  ! TODO
  call die('linear_forcing_adjoint: not yet implemented')
 
  return
end subroutine linear_forcing_adjoint
