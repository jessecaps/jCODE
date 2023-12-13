module binary_mixing

  ! External modules
  use functional

  implicit none

contains

  subroutine verify_binary_mixing_patch(patch)

    ! External modules
    use solver_options

    implicit none

    ! Arguments
    type(t_Patch), intent(in) :: patch

    ! Local variables
    integer :: i, n, extent(6)
    character(len = str_long) :: message

    if (nSpecies .ne. 1) call die('verify_binary_mixing_patch: nSpecies must = 1')

    extent(1) = patch%iMin; extent(2) = patch%iMax
    extent(3) = patch%jMin; extent(4) = patch%jMax
    extent(5) = patch%kMin; extent(6) = patch%kMax

    n = nDimensions

    do i = 1, nDimensions
       if (extent((i-1)*2+1) .eq. extent((i-1)*2+2)) n = n - 1
    end do

    n = nDimensions
    if (n .ne. nDimensions) then
       write(message, '(2(A,I0.0),A)') 'verify_binary_mixing_patch: Expected a ',            &
            nDimensions, 'D patch, but extent represents a ', n, 'D patch!'
       call die(trim(message))
    end if

    return
  end subroutine verify_binary_mixing_patch

end module binary_mixing


! ================================ !
! Setup the mixing cost functional !
! ================================ !
subroutine binary_mixing_setup

  ! Internal modules
  use binary_mixing

  ! External modules
  use parser
  use solver_options

  implicit none

  ! Verify the patch type
  call verify_binary_mixing_patch(functionalPatch)

  return
end subroutine binary_mixing_setup


! ================================== !
! Cleanup the mixing cost functional !
! ================================== !
subroutine binary_mixing_cleanup

  ! Internal modules
  use binary_mixing

  implicit none

  ! Nothing to do

  return
end subroutine binary_mixing_cleanup


! ================================== !
! Compute the mixing cost functional !
! ================================== !
subroutine binary_mixing_compute

  ! Internal modules
  use binary_mixing

  ! External modules
  use grid_functions, only : inner_product
  use state
  use onestep

  implicit none

  ! Local variables
  real(WP), allocatable :: F(:)

  allocate(F(nGridPoints))

  ! Define mixing as Y1 * Y2 = Y * (1 - Y) = Y - Y^2
  F = massFraction(:, 1) - massFraction(:, 1)**2

  instantaneousCostFunctional = inner_product(F, targetMollifier(:,1))

  deallocate(F)

  return
end subroutine binary_mixing_compute


! ================================== !
! Compute the mixing adjoint forcing !
! ================================== !
subroutine binary_mixing_adjoint_source(source)

  ! Internal modules
  use binary_mixing

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

           F = - forcingFactor * targetMollifier(gridIndex, 1) *                             &
                specificVolume(gridIndex, 1)

           source(gridIndex, 1) = source(gridIndex, 1) + F *                                 &
                (2.0_WP * massFraction(gridIndex, 1)**2 - massFraction(gridIndex, 1))
           source(gridIndex, nDimensions+2+1) = source(gridIndex, nDimensions+2+1) +         &
                F * (1.0_WP - 2.0_WP * massFraction(gridIndex, 1))

        end do
     end do
  end do

  return
end subroutine binary_mixing_adjoint_source
