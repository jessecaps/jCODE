module functional

  ! External modules
  use precision
  use grid_patch

  implicit none

  ! Functional types
  integer, parameter ::                                                                      &
       NULL_COST              = 0,                                                           &
       SOUND_COST             = 1,                                                           &
       PRESSURE_DRAG_COST     = 2,                                                           &
       DRAG_COST              = 3,                                                           &
       REYNOLDS_STRESS_COST   = 4,                                                           &
       TEMPERATURE_COST       = 5,                                                           &
       HEAT_RELEASE_COST      = 6,                                                           &
       REACTANT_COST          = 7,                                                           &
       MIXING_COST            = 8,                                                           &
       VELOCITY_NORM_COST     = 9,                                                           &
       MIXING_NORM_COST       = 10,                                                          &
       DATA_ASSIMILATION_COST = 11,                                                          &
       SHADOWGRAPH_COST       = 12

  ! Functional data
  type(t_Patch), pointer :: functionalPatch
  integer :: costFunctionalType
  real(WP) :: instantaneousCostFunctional, currentCostFunctional,  auxilaryCostFunctional,   &
       adjointForcingFactor, adjointCorrection
  real(WP), allocatable :: targetMollifier(:,:)

contains

  ! Verify the functional patch is correct
  ! --------------------------------------
  subroutine verify_functional_patch(patch)

    ! External modules
    use geometry

    implicit none

    ! Arguments
    type(t_Patch), intent(in) :: patch

    ! Local variables
    integer :: i, extent(6)

    extent(1) = patch%iMin; extent(2) = patch%iMax
    extent(3) = patch%jMin; extent(4) = patch%jMax
    extent(5) = patch%kMin; extent(6) = patch%kMax

    do i = 1, nDimensions
       if (extent((i-1)*2+1) .lt. 0 .or. extent((i-1)*2+2) .gt. globalGridSize(i) .or.       &
            extent((i-1)*2+1) .gt. extent((i-1)*2+2))                                        &
            call die("verify_functional_patch: Invalid extent on '" //                       &
            trim(patch%name) // "'!")
    end do

    return
  end subroutine verify_functional_patch


  ! Normalize the target mollifier
  ! ------------------------------
  subroutine normalize_target_mollifier

    ! External modules
    use parallel

    implicit none

    ! Local variables
    real(WP) :: mollifierNorm
    logical :: hasNegativeMollifier

    hasNegativeMollifier = any(targetMollifier(:,1) .lt. 0.0_WP)
    call parallel_lor(hasNegativeMollifier)
    if (hasNegativeMollifier)                                                                &
         call die('Target mollifying support is not non-negative everywhere!')

    call patch_quadrature(functionalPatch, targetMollifier(:,1), mollifierNorm)

    if (mollifierNorm .gt. 0.0_WP) targetMollifier = targetMollifier / mollifierNorm

    return
  end subroutine normalize_target_mollifier

end module functional


! ========================= !
! Setup the cost functional !
! ========================= !
subroutine functional_setup

  ! Internal modules
  use functional

  ! External modules
  use string
  use parser
  use solver_options

  implicit none

  ! Local varianbles
  integer :: i, j, k, gridIndex
  real(WP) :: steepness, fraction
  character(len = str_medium) :: type
  character(len = str_short) :: xyz(3)

  ! Initialize values
  instantaneousCostFunctional = 0.0_WP
  currentCostFunctional = 0.0_WP
  auxilaryCostFunctional = 0.0_WP
  adjointForcingFactor = 1.0_WP

  ! Get the functional type
  call parser_read('cost functional type', type, '')

  select case (trim(type))

  case ('sound')
     costFunctionalType = SOUND_COST

  case ('pressure drag')
     costFunctionalType = PRESSURE_DRAG_COST

  case ('drag')
     costFunctionalType = DRAG_COST

  case ('reynolds stress')
     costFunctionalType = REYNOLDS_STRESS_COST

  case ('temperature')
     costFunctionalType = TEMPERATURE_COST

  case ('heat release')
     costFunctionalType = HEAT_RELEASE_COST

  case ('reactant')
     costFunctionalType = REACTANT_COST

  case ('binary mixing')
     costFunctionalType = MIXING_COST

  case ('velocity norm')
     costFunctionalType = VELOCITY_NORM_COST

  case ('mixing norm')
     costFunctionalType = MIXING_NORM_COST

  case ('shadowgraph')
     costFunctionalType = SHADOWGRAPH_COST

  case ('data assimilation')
     costFunctionalType = DATA_ASSIMILATION_COST

  case default

     call die("functional_setup: Unknown cost functional type '" // trim(type) // "'")

  end select

  ! Connect the functional patch
  do i = 1, nPatches
     if (patches(i)%patchType .eq. COST_TARGET) then
        allocate(functionalPatch)
        functionalPatch => patches(i)
        exit
     end if
  end do
  if (i .gt. nPatches) call die('functional_setup: no target patch defined')

  ! Ensure only one functional is specified
  do j = i + 1, nPatches
     if (patches(j)%patchType .eq. COST_TARGET)                                              &
          call die('functional_setup: only one cost functional patch can be specified!')
  end do

  ! Verify the cost functional patch
  call verify_functional_patch(functionalPatch)

  ! Setup the target mollifier
  allocate(targetMollifier(nGridPoints, 1))
  targetMollifier = 0.0_WP
  do k = functionalPatch%iStart(3), functionalPatch%iEnd(3)
     do j = functionalPatch%iStart(2), functionalPatch%iEnd(2)
        do i = functionalPatch%iStart(1), functionalPatch%iEnd(1)
           gridIndex = i - gridOffset(1) + localGridSize(1) *                                &
                (j - 1 - gridOffset(2) + localGridSize(2) *                                  &
                (k - 1 - gridOffset(3)))
           targetMollifier(gridIndex, 1) = 1.0_WP
        end do
     end do
  end do
  call parser_read('target mollifier steepness', steepness, 20.0_WP)
  call parser_read('target mollifier fraction', fraction, 0.1_WP)
  xyz(1) = 'x'; xyz(2) = 'y'; xyz(3) = 'z'
  do i = 1, nDimensions
     call parser_read('target support type in ' // trim(xyz(i)), type, 'none')
     select case (trim(type))
     case ('tanh')
        call patch_tanh_support(functionalPatch, i, targetMollifier(:,1), steepness, fraction)

     case ('cubic')
        call patch_cubic_bspline_support(functionalPatch, i, targetMollifier(:,1))

     case ('none')

     case default
        call die("Unknown target support type '" // trim(type) // "'")

     end select
  end do
  call normalize_target_mollifier

  ! Setup the specific cost functional type
  select case (costFunctionalType)

  case (SOUND_COST)
     call acoustic_noise_setup

  case (PRESSURE_DRAG_COST)
     call form_drag_setup

  case (DRAG_COST)
     call drag_force_setup

  case (REYNOLDS_STRESS_COST)
     call reynolds_stress_setup

  case (TEMPERATURE_COST)
     call flame_temperature_setup

  case (HEAT_RELEASE_COST)
     call heat_release_setup

  case (REACTANT_COST)
     call reactant_depletion_setup

  case (MIXING_COST)
     call binary_mixing_setup

  case (VELOCITY_NORM_COST)
     call velocity_norm_setup

  case (MIXING_NORM_COST)
     call mixing_norm_setup

  case (SHADOWGRAPH_COST)
     call shadowgraph_setup

  case (DATA_ASSIMILATION_COST)
     call data_assimilation_setup

  end select

  return
end subroutine functional_setup


! =========================== !
! Cleanup the cost functional !
! =========================== !
subroutine functional_cleanup

  ! Internal modules
  use functional

  implicit none

  select case (costFunctionalType)

  case (SOUND_COST)
     call acoustic_noise_cleanup

  case (PRESSURE_DRAG_COST)
     call form_drag_cleanup

  case (DRAG_COST)
     call drag_force_cleanup

  case (REYNOLDS_STRESS_COST)
     call reynolds_stress_cleanup

  case (TEMPERATURE_COST)
     call flame_temperature_cleanup

  case (HEAT_RELEASE_COST)
     call heat_release_cleanup

  case (REACTANT_COST)
     call reactant_depletion_cleanup

  case (MIXING_COST)
     call binary_mixing_cleanup

  case (VELOCITY_NORM_COST)
     call velocity_norm_cleanup

  case (MIXING_NORM_COST)
     call mixing_norm_cleanup

  case (SHADOWGRAPH_COST)
     call shadowgraph_cleanup

  case (DATA_ASSIMILATION_COST)
     call data_assimilation_cleanup

  end select

  if (associated(functionalPatch)) nullify(functionalPatch)
  if (allocated(targetMollifier)) deallocate(targetMollifier)

  return
end subroutine functional_cleanup


! =========================== !
! Compute the cost functional !
! =========================== !
subroutine functional_compute

  ! Internal modules
  use functional

  ! External modules
  use simulation_flags

  implicit none

  if (predictionOnly) return

  ! Start the functional timer
  call timing_start('functional')

  select case (costFunctionalType)

  case (SOUND_COST)
     call acoustic_noise_compute

  case (PRESSURE_DRAG_COST)
     call form_drag_compute

  case (DRAG_COST)
     call drag_force_compute

  case (REYNOLDS_STRESS_COST)
     call reynolds_stress_compute

  case (TEMPERATURE_COST)
     call flame_temperature_compute

  case (HEAT_RELEASE_COST)
     call heat_release_compute

  case (REACTANT_COST)
     call reactant_depletion_compute

  case (MIXING_COST)
     call binary_mixing_compute

  case (VELOCITY_NORM_COST)
     call velocity_norm_compute

  case (MIXING_NORM_COST)
     call mixing_norm_compute

  case (SHADOWGRAPH_COST)
     call shadowgraph_compute

  end select

  ! Stop the functional timer
  call timing_stop('functional')

  return
end subroutine functional_compute


! ================================================= !
! Compute the cost functional at terminal condition !
! ================================================= !
subroutine functional_compute_tc

  ! Internal modules
  use functional

  ! External modules
  use simulation_flags

  implicit none

  if (predictionOnly) return

  ! Start the functional timer
  call timing_start('functional')

  select case (costFunctionalType)

  case (DATA_ASSIMILATION_COST)
     call data_assimilation_compute

  end select

  ! Stop the functional timer
  call timing_stop('functional')

  return
end subroutine functional_compute_tc


! =========================== !
! Compute the adjoint forcing !
! =========================== !
subroutine functional_adjoint_source(source)

  ! Internal modules
  use functional

  ! External modules
  use simulation_flags
  use solver_options
  use geometry

  implicit none

  ! Arguments
  real(WP), dimension(nGridPoints, nUnknowns), intent(inout) :: source

  ! Start the functional timer
  call timing_start('functional')

  select case (costFunctionalType)

  case (SOUND_COST)
     call acoustic_noise_adjoint_source(source)

  case (PRESSURE_DRAG_COST)
     call form_drag_adjoint_source(source)

  case (DRAG_COST)
     call drag_force_adjoint_source(source)

  case (REYNOLDS_STRESS_COST)
     call reynolds_stress_adjoint_source(source)

  case (TEMPERATURE_COST)
     call flame_temperature_adjoint_source(source)

  case (HEAT_RELEASE_COST)
     call heat_release_adjoint_source(source)

  case (REACTANT_COST)
     call reactant_depletion_adjoint_source(source)

  case (MIXING_COST)
     call binary_mixing_adjoint_source(source)

  case (VELOCITY_NORM_COST)
     call velocity_norm_adjoint_source(source)

  case (MIXING_NORM_COST)
     call mixing_norm_adjoint_source(source)

  case (SHADOWGRAPH_COST)
     call shadowgraph_adjoint_source(source)

  end select

  ! Stop the functional timer
  call timing_stop('functional')

  return
end subroutine functional_adjoint_source


! ================================================= !
! Compute the adjoint forcing at terminal condition !
! ================================================= !
subroutine functional_adjoint_source_tc(source)

  ! Internal modules
  use functional

  ! External modules
  use simulation_flags
  use solver_options
  use geometry

  implicit none

  ! Arguments
  real(WP), dimension(nGridPoints, nUnknowns), intent(inout) :: source

  select case (costFunctionalType)

  case (DATA_ASSIMILATION_COST)
     call data_assimilation_adjoint_source(source)

  end select

  return
end subroutine functional_adjoint_source_tc
