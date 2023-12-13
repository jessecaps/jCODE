module heat_release

  ! External modules
  use functional

  implicit none

  logical :: useTimeRamp
  real(WP) :: rampPeak, rampWidth

contains

  subroutine verify_heat_release_patch(patch)

    ! External modules
    use combustion

    implicit none

    ! Arguments
    type(t_Patch), intent(in) :: patch

    ! Local variables
    integer :: i, n, extent(6)
    character(len = str_long) :: message

    if (chemistryModel .ne. ONE_STEP) call die('verify_heat_release_patch: &
         &flame temperature cost functional requires one-step chemistry!')

    extent(1) = patch%iMin; extent(2) = patch%iMax
    extent(3) = patch%jMin; extent(4) = patch%jMax
    extent(5) = patch%kMin; extent(6) = patch%kMax

    n = nDimensions
    do i = 1, nDimensions
       if (extent((i-1)*2+1) .eq. extent((i-1)*2+2)) n = n - 1
    end do

    if (n .ne. nDimensions) then
       write(message, '(2(A,I0.0),A)') 'verify_heat_release_patch: Expected a ',             &
            nDimensions, 'D patch, but extent represents a ', n, 'D patch!'
       call die(trim(message))
    end if

    return
  end subroutine verify_heat_release_patch

end module heat_release


! ====================================== !
! Setup the heat release cost functional !
! ====================================== !
subroutine heat_release_setup

  ! Internal modules
  use heat_release

  ! External modules
  use parser

  implicit none

  ! Verify the patch type
  call verify_heat_release_patch(functionalPatch)

  call parser_read('heat release use time ramp', useTimeRamp, .false.)
  if (useTimeRamp) then
     call parser_read('heat release ramp peak', rampPeak)
     call parser_read('heat release ramp width', rampWidth)
  end if

  return
end subroutine heat_release_setup


! ======================================== !
! Cleanup the heat release cost functional !
! ======================================== !
subroutine heat_release_cleanup

  ! Internal modules
  use heat_release

  implicit none

  ! Nothing to do

  return
end subroutine heat_release_cleanup


! ======================================== !
! Compute the heat release cost functional !
! ======================================== !
subroutine heat_release_compute

  ! Internal modules
  use heat_release

  ! External modules
  use solver_options
  use geometry
  use grid
  use grid_functions, only : inner_product
  use state
  use time_info
  use onestep

  implicit none

  ! Local variables
  real(WP) :: timeRampFactor, referenceTemperature, flameTemperature, activationTemperature, H
  real(WP), allocatable :: F(:)

  ! Heat release cost functional is based on one-step chemistry
  referenceTemperature = 1.0_WP / (ratioOfSpecificHeats - 1.0_WP)
  flameTemperature = referenceTemperature / (1.0_WP - heatRelease)
  activationTemperature = zelDovich / heatRelease * flameTemperature
  H = heatRelease * flameTemperature / Yfs ! ... heat release coefficient

  timeRampFactor = 1.0_WP
  if (useTimeRamp) timeRampFactor = exp(- 0.5_WP * (time - rampPeak)**2 / rampWidth**2)

  allocate(F(nGridPoints))

  F = H * controlParameter * Damkohler * conservedVariables(:,nDimensions+2+H2) *            &
       conservedVariables(:,nDimensions+2+O2) *                                              &
       exp(- activationTemperature / temperature(:,1))
  instantaneousCostFunctional = inner_product(F, targetMollifier(:,1))

  deallocate(F)

  auxilaryCostFunctional = instantaneousCostFunctional
  instantaneousCostFunctional = instantaneousCostFunctional * timeRampFactor

  return
end subroutine heat_release_compute


! ======================================== !
! Compute the heat release adjoint forcing !
! ======================================== !
subroutine heat_release_adjoint_source(source)

  ! Internal modules
  use heat_release

  ! External modules
  use state_jacobian
  use simulation_flags
  use solver_options
  use geometry
  use state
  use time_info
  use onestep

  implicit none

  ! Arguments
  real(WP), dimension(nGridPoints, nUnknowns), intent(inout) :: source

  ! Local variables
  integer :: i, j, k, gridIndex, patchIndex
  real(WP) :: forcingFactor, flameTemperature, referenceTemperature, activationTemperature,  &
       F, H, timeRampFactor, temp
  real(WP), dimension(:), allocatable :: deltaTemperature

  if (useContinuousAdjoint) then
     forcingFactor = 1.0_WP
  else
     forcingFactor = adjointForcingFactor
  end if

  referenceTemperature = 1.0_WP / (ratioOfSpecificHeats - 1.0_WP)
  flameTemperature = referenceTemperature / (1.0_WP - heatRelease)
  activationTemperature = zelDovich / heatRelease * flameTemperature
  H = heatRelease * flameTemperature / Yfs ! ... heat release coefficient

  timeRampFactor = 1.0_WP
  if (useTimeRamp) timeRampFactor = exp(-0.5_WP * (adjointCoefficienTtime - rampPeak)**2 /   &
       rampWidth **2)

  allocate(deltaTemperature(nDimensions + nSpecies + 2))

  do k = functionalPatch%iStart(3), functionalPatch%iEnd(3)
     do j = functionalPatch%iStart(2), functionalPatch%iEnd(2)
        do i = functionalPatch%iStart(1), functionalPatch%iEnd(1)
           gridIndex = i - gridOffset(1) + localGridSize(1) *                                &
                (j - 1 - gridOffset(2) + localGridSize(2) *                                  &
                (k - 1 - gridOffset(3)))
           patchIndex = i - functionalPatch%offset(1) +                                      &
                functionalPatch%localSize(1) * (j - 1 - functionalPatch%offset(2) +          &
                functionalPatch%localSize(2) * (k - 1 - functionalPatch%offset(3)))

           call compute_delta_variables(conservedVariables(gridIndex,:),                     &
                deltaTemperature = deltaTemperature)

           temp = H * controlParameter * Damkohler *                                         &
                conservedVariables(gridIndex,nDimensions+2+H2) *                             &
                conservedVariables(gridIndex,nDimensions+2+O2) *                             &
                exp(- activationTemperature / temperature(gridIndex, 1)) *                   &
                activationTemperature / temperature(gridIndex, 1) ** 2

           F = - forcingFactor * targetMollifier(gridIndex, 1) * timeRampFactor

           source(gridIndex,:) = source(gridIndex,:) + F * temp * deltaTemperature

           temp = H * controlParameter * Damkohler *                                         &
                exp(- activationTemperature / temperature(gridIndex, 1))

           source(gridIndex,nDimensions+2+H2) = source(gridIndex,nDimensions+2+H2) +         &
                F * temp * conservedVariables(gridIndex,nDimensions+2+O2)

           source(gridIndex,nDimensions+2+O2) = source(gridIndex,nDimensions+2+O2) +         &
                F * temp * conservedVariables(gridIndex,nDimensions+2+H2)

        end do
     end do
  end do

  deallocate(deltaTemperature)

  return
end subroutine heat_release_adjoint_source
