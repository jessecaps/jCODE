module reactant_depletion

  ! External modules
  use functional

  implicit none

  integer :: reactant

contains

  subroutine verify_reactant_depletion_patch(patch)
    implicit none

    ! Arguments
    type(t_Patch), intent(in) :: patch

    ! Local variables
    integer :: i, extent(6)

    extent(1) = patch%iMin; extent(2) = patch%iMax
    extent(3) = patch%jMin; extent(4) = patch%jMax
    extent(5) = patch%kMin; extent(6) = patch%kMax

    i = abs(patch%normalDirection)
    if (extent((i-1)*2+1) .ne. extent((i-1)*2+2)) call die('verify_reactant_depletion_patch: &
         &Extends more than 1 grid point along normal direction!')

    return
  end subroutine verify_reactant_depletion_patch

end module reactant_depletion


! ================================== !
! Setup the reactant cost functional !
! ================================== !
subroutine reactant_depletion_setup

  ! Internal modules
  use reactant_depletion

  ! External modules
  use parser
  use solver_options

  implicit none

  ! Verify the patch type
  call verify_reactant_depletion_patch(functionalPatch)

  call parser_read('deplete species number', reactant)
  if (reactant .gt. nSpecies) call die('reactant_depletion_setup: &
       &species number must be <= nSpecies!')
  if (reactant .le. 0) call die('reactant_depletion_setup: &
       &species number must be > 0!')

  return
end subroutine reactant_depletion_setup


! ==================================== !
! Cleanup the reactant cost functional !
! ==================================== !
subroutine reactant_depletion_cleanup

  ! Internal modules
  use reactant_depletion

  implicit none

  ! Nothing to do

  return
end subroutine reactant_depletion_cleanup


! ==================================== !
! Compute the reactant cost functional !
! ==================================== !
subroutine reactant_depletion_compute

  ! Internal modules
  use reactant_depletion

  ! External modules
  use grid_functions, only : inner_product
  use state
  use onestep

  implicit none

  ! Local variables
  real(WP), allocatable :: F(:,:)

  allocate(F(nGridPoints, 1))

  F(:,1) = massFraction(:, reactant)

  instantaneousCostFunctional = inner_product(F, F, targetMollifier(:,1))

  deallocate(F)

  return
end subroutine reactant_depletion_compute


! ==================================== !
! Compute the reactant adjoint forcing !
! ==================================== !
subroutine reactant_depletion_adjoint_source(source)

  ! Internal modules
  use reactant_depletion

  ! External modules
  use simulation_flags
  use solver_options
  use geometry
  use grid
  use state

  implicit none

  ! Arguments
  real(WP), dimension(nGridPoints, nUnknowns), intent(inout) :: source

  ! Local variables
  integer :: i, j, k, gridIndex, patchIndex
  real(WP) :: forcingFactor, F

  if (useContinuousAdjoint) then
     forcingFactor = 1.0_WP
  else
     forcingFactor = adjointForcingFactor
  end if

  do k = functionalPatch%iStart(3), functionalPatch%iEnd(3)
     do j = functionalPatch%iStart(2), functionalPatch%iEnd(2)
        do i = functionalPatch%iStart(1), functionalPatch%iEnd(1)
           gridIndex = i - gridOffset(1) + localGridSize(1) *                                &
                (j - 1 - gridOffset(2) + localGridSize(2) *                                  &
                (k - 1 - gridOffset(3)))
           patchIndex = i - functionalPatch%offset(1) +                                      &
                functionalPatch%localSize(1) * (j - 1 - functionalPatch%offset(2) +          &
                functionalPatch%localSize(2) * (k - 1 - functionalPatch%offset(3)))

           F = - 2.0_WP * forcingFactor * targetMollifier(gridIndex, 1) *                    &
                massFraction(gridIndex, reactant) * specificVolume(gridIndex, 1)

           source(gridIndex, 1) = source(gridIndex, 1) -                                     &
                F * massFraction(gridIndex, reactant)
           source(gridIndex, nDimensions+2+reactant) =                                       &
                source(gridIndex, nDimensions+2+reactant) + F

        end do
     end do
  end do

  return
end subroutine reactant_depletion_adjoint_source
