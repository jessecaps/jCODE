module acoustic_noise

  ! External modules
  use functional

  implicit none

  real(WP), allocatable :: meanPressure(:)

contains

  subroutine verify_acoustic_noise_patch(patch)

    implicit none

    ! Arguments
    type(t_Patch), intent(in) :: patch

    ! Local variables
    integer :: i, n, extent(6)
    character(len = str_long) :: message

    extent(1) = patch%iMin; extent(2) = patch%iMax
    extent(3) = patch%jMin; extent(4) = patch%jMax
    extent(5) = patch%kMin; extent(6) = patch%kMax

    n = nDimensions
    do i = 1, nDimensions
       if (extent((i-1)*2+1) .eq. extent((i-1)*2+2)) n = n - 1
    end do

    if (n .ne. nDimensions) then
       write(message, '(2(A,I0.0),A)') 'verify_acoustic_noise_patch: Expected a ',           &
            nDimensions, 'D patch, but extent represents a ', n, 'D patch!'
       call die(trim(message))
    end if

    return
  end subroutine verify_acoustic_noise_patch

end module acoustic_noise


! ======================================== !
! Setup the acoustic noise cost functional !
! ======================================== !
subroutine acoustic_noise_setup

  ! Internal modules
  use acoustic_noise

  ! External modules
  use parser
  use state

  implicit none

  ! Local variables
  real(WP) :: constantPressure

  ! Verify the patch type
  call verify_acoustic_noise_patch(functionalPatch)

  ! Get the mean pressure
  allocate(meanPressure(nGridPoints))
  call parser_read('mean pressure', constantPressure, -1.0_WP)
  if (constantPressure .ge. 0.0_WP) then ! ... set mean pressure constant everywhere
     meanPressure = constantPressure
  else
     meanPressure = pressure(:,1) ! ... set mean pressure to quiescent or target solution
  end if

  return
end subroutine acoustic_noise_setup


! ========================================== !
! Cleanup the acoustic noise cost functional !
! ========================================== !
subroutine acoustic_noise_cleanup

  ! Internal modules
  use acoustic_noise

  implicit none

  if (allocated(meanPressure)) deallocate(meanPressure)

  return
end subroutine acoustic_noise_cleanup


! ========================================== !
! Compute the acoustic noise cost functional !
! ========================================== !
subroutine acoustic_noise_compute

  ! Internal modules
  use acoustic_noise

  ! External modules
  use geometry
  use grid
  use grid_functions, only : inner_product
  use state, only : pressure

  implicit none

  ! Local variables
  real(WP), allocatable :: F(:,:)

  allocate(F(nGridPoints, 1))
  F(:,1) = pressure(:,1) - meanPressure
  instantaneousCostFunctional = inner_product(F, F, targetMollifier(:,1))
  deallocate(F)

  return
end subroutine acoustic_noise_compute


! ========================================== !
! Compute the acoustic noise adjoint forcing !
! ========================================== !
subroutine acoustic_noise_adjoint_source(source)

  ! Internal modules
  use acoustic_noise

  ! External modules
  use simulation_flags
  use solver_options
  use geometry
  use grid
  use grid_functions, only : inner_product
  use grid_patch
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
                (ratioOfSpecificHeats - 1.0_WP) *                                            &
                (pressure(gridIndex, 1) - meanPressure(gridIndex))

           source(gridIndex, nDimensions+2) = source(gridIndex, nDimensions+2) + F
           source(gridIndex, 2:nDimensions+1) = source(gridIndex, 2:nDimensions+1) -         &
                velocity(gridIndex,:) * F
           source(gridIndex, 1) = source(gridIndex, 1) +                                     &
                0.5_WP * sum(velocity(gridIndex,:) ** 2) * F

        end do
     end do
  end do

  return
end subroutine acoustic_noise_adjoint_source
