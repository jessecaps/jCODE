module flame_temperature

  ! External modules
  use functional

  implicit none

  real(WP) :: burnRadius, rampPeak, rampWidth, rampSteepness, rampFraction
  real(WP), allocatable :: targetTemperature(:)
  logical :: weightBurnRegion, useTimeRamp, anchorFlame

contains

  subroutine verify_flame_temperature_patch(patch)

    ! External modules
    use combustion

    implicit none

    ! Arguments
    type(t_Patch), intent(in) :: patch

    ! Local variables
    integer :: i, n, extent(6)
    character(len = str_long) :: message

    if (chemistryModel .ne. ONE_STEP) call die('verify_flame_temperature_patch: &
         &flame temperature cost functional requires one-step chemistry!')

    extent(1) = patch%iMin; extent(2) = patch%iMax
    extent(3) = patch%jMin; extent(4) = patch%jMax
    extent(5) = patch%kMin; extent(6) = patch%kMax

    n = nDimensions
    do i = 1, nDimensions
       if (extent((i-1)*2+1) .eq. extent((i-1)*2+2)) n = n - 1
    end do

    if (n .ne. nDimensions) then
       write(message, '(2(A,I0.0),A)') 'verify_flame_temperature_patch: Expected a ',        &
            nDimensions, 'D patch, but extent represents a ', n, 'D patch!'
       call die(trim(message))
    end if

    return
  end subroutine verify_flame_temperature_patch

end module flame_temperature


! =========================================== !
! Setup the flame temperature cost functional !
! =========================================== !
subroutine flame_temperature_setup

  ! Internal modules
  use flame_temperature

  ! External modules
  use parser
  use solver_options
  use geometry
  use grid
  use onestep
  use state
  use equation_of_state

  implicit none

  ! Local variables
  integer :: i
  real(WP) :: flamePosition, referenceTemperature, flameTemperature

  ! Verify the patch type
  call verify_flame_temperature_patch(functionalPatch)

  call parser_read('flame temperature weight burn region', weightBurnRegion, .false.)
  if (weightBurnRegion) call parser_read('flame temperature burn radius', burnRadius)

  call parser_read('flame temperature use time ramp', useTimeRamp, .false.)
  if (useTimeRamp) then
     call parser_read('flame temperature ramp peak', rampPeak)
     call parser_read('flame temperature ramp width', rampWidth)
     call parser_read('flame temperature ramp steepness', rampSteepness, 10.0_WP)
     call parser_read('flame temperature ramp fraction', rampFraction, 0.2_WP)
  end if

  ! Allocate the target temperature
  allocate(targetTemperature(nGridPoints))

  call parser_read('flame temperature anchor flame', anchorFlame, .false.)
  if (anchorFlame) then
     ! Compute the target temperature to anchor a flame
     call parser_read('flame temperature anchor position', flamePosition)
     referenceTemperature = 1.0_WP / (ratioOfSpecificHeats - 1.0_WP)
     flameTemperature = referenceTemperature / (1.0_WP - heatRelease)
     do i = 1, nGridPoints
        targetTemperature(i) = referenceTemperature + 0.5_WP *                               &
             (flameTemperature - referenceTemperature) * (1.0_WP + tanh(coordinates(i,1) -   &
             flamePosition))
     end do
  else
     ! Get target temperature from target state
     if(.not.useTargetState) call die('flame_temperature_setup: target state needed!')
     call compute_dependent_variables(targetState, temperature = targetTemperature)
  end if

  return
end subroutine flame_temperature_setup


! ============================================= !
! Cleanup the flame temperature cost functional !
! ============================================= !
subroutine flame_temperature_cleanup

  ! Internal modules
  use flame_temperature

  implicit none

  if (allocated(targetTemperature)) deallocate(targetTemperature)

  return
end subroutine flame_temperature_cleanup


! ============================================= !
! Compute the flame temperature cost functional !
! ============================================= !
subroutine flame_temperature_compute

  ! Internal modules
  use flame_temperature

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
  integer :: i
  real(WP) :: YF0, YO0, s, Z, Zst, gaussianFactor, timeRampFactor, flameTemperature,         &
       currentTime
  real(WP), allocatable :: F(:), W(:)

  ! Flame temperature cost functional is based on one-step chemistry
  flameTemperature = 1.0_WP / (ratioOfSpecificHeats - 1.0_WP) / (1.0_WP - heatRelease)

  currentTime = time

  timeRampFactor = 1.0_WP
  if (useTimeRamp) timeRampFactor = exp(- 0.5_WP * (currentTime - rampPeak)**2 / rampWidth**2)
  !if (useTimeRamp) then
  !   tStart = rampPeak - 0.5_WP * rampWidth
  !   tEnd = rampPeak + 0.5_WP * rampWidth
  !   t = 2.0_WP * (currentTime - 0.5_WP * timeStepSize - tStart) / (tEnd - tStart) - 1.0_WP
  !   timeRampFactor = 0.5_WP * (tanh(rampSteepness * (t + 1.0_WP - 0.5_WP * rampFraction)) - &
  !        tanh(rampSteepness * (t - 1.0_WP + 0.5_WP * rampFraction)))
  !end if
  
  allocate(W(nGridPoints))

  if (weightBurnRegion) then

     YF0 = Y0(H2)
     YO0 = Y0(O2)
     s = 1.0_WP / stoichiometricRatio
     Zst = 1.0_WP / (1.0_WP + s * YF0 / YO0)
     gaussianFactor = -0.5_WP / burnRadius**2

     do i = 1, nGridPoints
        Z = (massFraction(i, H2) - massFraction(i, O2) * s  + YO0 * s) / (YF0 + YO0 * s)
        W(i) = targetMollifier(i, 1) * exp(gaussianFactor * (Z - Zst) **2)
     end do

  else

     W = targetMollifier(:, 1)

  end if

  allocate(F(nGridPoints))

  !F = (temperature(:,1) - targetTemperature)**2
  F = (temperature(:,1) - targetTemperature) / (flameTemperature - targetTemperature)
  instantaneousCostFunctional = inner_product(F, W)

  deallocate(F)
  deallocate(W)

  auxilaryCostFunctional = instantaneousCostFunctional
  instantaneousCostFunctional = instantaneousCostFunctional * timeRampFactor

  return
end subroutine flame_temperature_compute


! ============================================= !
! Compute the flame temperature adjoint forcing !
! ============================================= !
subroutine flame_temperature_adjoint_source(source)

  ! Internal modules
  use flame_temperature

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
  real(WP) :: forcingFactor, flameTemperature, F, W, Z, Zst, s, YF0, YO0, gaussianFactor,    &
       currentTime, timeRampFactor
  real(WP), dimension(:), allocatable :: deltaTemperature

  if (useContinuousAdjoint) then
     forcingFactor = 1.0_WP
  else
     forcingFactor = adjointForcingFactor
  end if

  currentTime = adjointCoefficientTime

  flameTemperature = 1.0_WP / (ratioOfSpecificHeats - 1.0_WP) / (1.0_WP - heatRelease)

  timeRampFactor = 1.0_WP
  if (useTimeRamp) timeRampFactor = exp(-0.5_WP * (currentTime - rampPeak)**2 / rampWidth **2)
  !if (useTimeRamp) then
  !   tStart = rampPeak - 0.5_WP * rampWidth
  !   tEnd = rampPeak + 0.5_WP * rampWidth
  !   t = 2.0_WP * (currentTime - 0.5_WP * timeStepSize - tStart) / (tEnd - tStart) - 1.0_WP
  !   timeRampFactor = 0.5_WP * (tanh(rampSteepness * (t + 1.0_WP - 0.5_WP * rampFraction)) - &
  !        tanh(rampSteepness * (t - 1.0_WP + 0.5_WP * rampFraction)))
  !end if

  if (weightBurnRegion) then

     gaussianFactor = -0.5_WP / burnRadius**2
     YF0 = Y0(H2)
     YO0 = Y0(O2)
     s = 1.0_WP / stoichiometricRatio
     Zst = 1.0_WP / (1.0_wp + YF0 / YO0 / s)

  end if

  allocate(deltaTemperature(nUnknowns))

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

           if (weightBurnRegion) then

              Z = (massFraction(gridIndex, H2) - massFraction(gridIndex, O2) * s  +          &
                   YO0 * s ) / (YF0 + YO0 * s)
              W = targetMollifier(gridIndex, 1) * exp(gaussianFactor * (Z - Zst) **2)

              ! First apply -W/(Tf - T0)*dT/dQ
              F = - forcingFactor * W * timeRampFactor /                                     &
                   (flameTemperature - targetTemperature(gridIndex))

              source(gridIndex,:) = source(gridIndex, :) + F * deltaTemperature

              ! Now apply -(T-T0)/(Tf-T0)*dW/dQ
              F =  forcingFactor * W * timeRampFactor *                                      &
                   (temperature(gridIndex, 1) - targetTemperature(gridIndex)) /              &
                   (flameTemperature - targetTemperature(gridIndex)) *                       &
                   ((Z - Zst) / burnRadius**2) * specificVolume(gridIndex,1) / (YF0 + YO0 * s)

              source(gridIndex, 1) = source(gridIndex, 1)  +                                 &
                   (massFraction(gridIndex, O2) * s - massFraction(gridIndex, H2)) * F
              source(gridIndex,nDimensions+2+H2) = source(gridIndex,nDimensions+2+H2) + F
              source(gridIndex,nDimensions+2+O2) = source(gridIndex,nDimensions+2+O2) - F * s

           else

              F = - forcingFactor * targetMollifier(gridIndex, 1) * timeRampFactor /         &
                   (flameTemperature - targetTemperature(gridIndex))
              !F = - 2.0_WP * forcingFactor * targetMollifier(gridIndex, 1) * timeRampFactor *&
              !     (temperature(gridIndex,1) - targetTemperature(gridIndex))

              source(gridIndex,:) = source(gridIndex,:) + F * deltaTemperature

           end if

        end do
     end do
  end do

  deallocate(deltaTemperature)

  return
end subroutine flame_temperature_adjoint_source
