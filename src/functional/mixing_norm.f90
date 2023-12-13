module mixing_norm

  ! External modules
  use functional

  implicit none

  ! Perturbation parameters
  real(WP) :: initialThicknessInverse, rampPeak, rampWidth, rampCenter, rampSteepness,       &
       rampStep
  character(len=str_medium) :: rampFunction
  logical  :: useTimeRamp

contains

  subroutine verify_mixing_norm_patch(patch)

    ! External modules
    use solver_options

    implicit none

    ! Arguments
    type(t_Patch), intent(in) :: patch

    ! Local variables
    integer :: i, n, extent(6)
    character(len = str_long) :: message

    if (nSpecies .ne. 1) call die('verify_mixing_norm_patch: nSpecies must = 1')

    extent(1) = patch%iMin; extent(2) = patch%iMax
    extent(3) = patch%jMin; extent(4) = patch%jMax
    extent(5) = patch%kMin; extent(6) = patch%kMax

    n = nDimensions

    do i = 1, nDimensions
       if (extent((i-1)*2+1) .eq. extent((i-1)*2+2)) n = n - 1
    end do

    n = nDimensions
    if (n .ne. nDimensions) then
       write(message, '(2(A,I0.0),A)') 'verify_mixing_norm_patch: Expected a ',              &
            nDimensions, 'D patch, but extent represents a ', n, 'D patch!'
       call die(trim(message))
    end if

    return
  end subroutine verify_mixing_norm_patch

  ! Compute the time ramp factor
  ! ----------------------------
  function get_time_ramp(currentTime)

    implicit none

    ! Arguments
    real(WP), intent(in) :: currentTime

    ! Result
    real(WP) :: get_time_ramp

    select case (rampFunction)

    case ('gaussian')
       get_time_ramp = exp(-0.5_WP * (currentTime - rampPeak)**2 / rampWidth**2)

    case ('tanh')
       get_time_ramp = 0.5_WP * (1.0_WP + tanh(rampSteepness * (currentTime - rampCenter)))

    case ('step')
       if (currentTime .lt. rampStep) then
          get_time_ramp = 0.0_WP
       else
          get_time_ramp = 1.0_WP
       end if

    end select

  end function get_time_ramp

end module mixing_norm


! ===================================== !
! Setup the mixing norm cost functional !
! ===================================== !
subroutine mixing_norm_setup

  ! Internal modules
  use mixing_norm

  ! External modules
  use parser
  use solver_options

  implicit none

  ! Verify the patch type
  call verify_mixing_norm_patch(functionalPatch)

  ! Read in initial diffusion thickness for Rayleigh-Taylor initialazing
  call parser_read('initial diffusion thickness', initialThicknessInverse)
  initialThicknessInverse = 1.0_WP / initialThicknessInverse

  call parser_read('use functional time ramp', useTimeRamp, .false.)
  if (useTimeRamp) then
     call parser_read('ramp function', rampFunction, '')

     select case (rampFunction)

     case ('gaussian')
        call parser_read('ramp peak', rampPeak)
        call parser_read('ramp width', rampWidth)

     case ('tanh')
        call parser_read('ramp center', rampCenter)
        call parser_read('ramp steepness', rampSteepness)

     case ('step')
        call parser_read('ramp step', rampStep)

     case default
        call die('mixing_norm_setup: Unknown ramp function: "' // trim(rampFunction) // '"')

     end select
  end if

  return
end subroutine mixing_norm_setup


! ======================================= !
! Cleanup the mixing norm cost functional !
! ======================================= !
subroutine mixing_norm_cleanup

  ! Internal modules
  use mixing_norm

  implicit none

  ! Nothing to do

  return
end subroutine mixing_norm_cleanup


! ======================================= !
! Compute the mixing norm cost functional !
! ======================================= !
subroutine mixing_norm_compute

  ! Internal modules
  use mixing_norm

  ! External modules
  use grid_functions, only : inner_product
  use state
  use onestep
  use time_info

  implicit none

  ! Local variables
  integer :: i
  real(WP) :: timeRampFactor, currentTime, molarFraction
  real(WP), allocatable :: F(:)

  allocate(F(nGridPoints))

  currentTime = time

  timeRampFactor = 1.0_WP
  if (useTimeRamp) timeRampFactor = get_time_ramp(currentTime)

  do i = 1, nGridPoints
     ! Compute the molar fraction
     molarFraction = massFraction(i,1) * mixtureMolecularWeight(i, 1) *                      &
          molecularWeightInverse(1)

     ! Define mixing norm as (X-X_0)^2 where X_0 is the reference molar fraction
     ! (an error fucntion distrinution) 
     F(i) = ( molarFraction - 0.5_WP * ( 1.0_WP +                                            &
          erf(coordinates(i, 2) * initialThicknessInverse) ) )**2
  end do

  instantaneousCostFunctional = inner_product(F, targetMollifier(:,1))

  deallocate(F)

  auxilaryCostFunctional = instantaneousCostFunctional
  instantaneousCostFunctional = instantaneousCostFunctional * timeRampFactor

  return
end subroutine mixing_norm_compute


! ======================================= !
! Compute the mixing norm adjoint forcing !
! ======================================= !
subroutine mixing_norm_adjoint_source(source)

  ! Internal modules
  use mixing_norm

  ! External modules
  use simulation_flags
  use solver_options
  use geometry
  use grid
  use state
  use time_info

  implicit none

  ! Arguments
  real(WP), dimension(nGridPoints, nUnknowns), intent(inout) :: source

  ! Local variables
  integer :: i, j, k, gridIndex, patchIndex
  real(WP) :: forcingFactor, F, currentTime, timeRampFactor, molarFraction

  if (useContinuousAdjoint) then
     forcingFactor = 1.0_WP
  else
     forcingFactor = adjointForcingFactor
  end if

  currentTime = adjointCoefficientTime

  timeRampFactor = 1.0_WP
  if (useTimeRamp) timeRampFactor = get_time_ramp(currentTime)

  do k = functionalPatch%iStart(3), functionalPatch%iEnd(3)
     do j = functionalPatch%iStart(2), functionalPatch%iEnd(2)
        do i = functionalPatch%iStart(1), functionalPatch%iEnd(1)
           gridIndex = i - gridOffset(1) + localGridSize(1) *                                &
                (j - 1 - gridOffset(2) + localGridSize(2) *                                  &
                (k - 1 - gridOffset(3)))
           patchIndex = i - functionalPatch%offset(1) +                                      &
                functionalPatch%localSize(1) * (j - 1 - functionalPatch%offset(2) +          &
                functionalPatch%localSize(2) * (k - 1 - functionalPatch%offset(3)))

           ! Compute the molar fraction
           molarFraction = massFraction(gridIndex,1) * mixtureMolecularWeight(gridIndex,1) * &
                molecularWeightInverse(1)

           F = - forcingFactor * targetMollifier(gridIndex, 1) * specificVolume(gridIndex, 1)&
                * mixtureMolecularWeight(gridIndex, 1) * 2.0_WP * (molarFraction -           &
                0.5_WP * (1.0_WP + erf(coordinates(gridIndex, 2) * initialThicknessInverse)))&
                * timeRampFactor

           source(gridIndex, 1) = source(gridIndex, 1) - F * molarFraction *                 &
                molecularWeightInverse(2)
           source(gridIndex, nDimensions+2+1) = source(gridIndex, nDimensions+2+1) + F *     &
                ( molecularWeightInverse(1) - molarFraction *                                &
                ( molecularWeightInverse(1) - molecularWeightInverse(2) ) )

        end do
     end do
  end do

  return
end subroutine mixing_norm_adjoint_source
