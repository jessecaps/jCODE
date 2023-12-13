module reynolds_stress

  ! External modules
  use functional

  implicit none

  real(WP) :: firstDirection(3), secondDirection(3)
  real(WP), allocatable :: meanVelocity(:,:)

contains

  subroutine verify_reynolds_stress_patch(patch)

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
       write(message, '(2(A,I0.0),A)') 'verify_reynolds_stress_patch: Expected a ',          &
            nDimensions, 'D patch, but extent represents a ', n, 'D patch!'
       call die(trim(message))
    end if

    return
  end subroutine verify_reynolds_stress_patch

end module reynolds_stress


! ========================================= !
! Setup the Reynolds stress cost functional !
! ========================================= !
subroutine reynolds_stress_setup

  ! Internal modules
  use reynolds_stress

  ! External modules
  use parser
  use string
  use simulation_flags
  use state

  implicit none

  ! Local variables
  real(WP) :: constantVelocity
  character(len = str_short) :: dir

  ! Verify the patch type
  call verify_reynolds_stress_patch(functionalPatch)

  ! Get Reynolds stress components
  call parser_read('reynolds stress components', dir)
  if (len_trim(dir).ne.2)                                                                    &
       call die("reynolds_stress_setup: length of 'reynolds stress components' must be 2")

  ! Set the first direction
  firstDirection  = 0.0_WP
  select case (trim(dir(1:1)))
  case ('x')
     firstDirection(1) = 1.0_WP
  case ('y')
     firstDirection(2) = 1.0_WP
  case ('z')
     firstDirection(3) = 1.0_WP
  case default
     call die("reynolds_stress_setup: unknown components '" // trim(dir))
  end select

  ! Set the second direction
  secondDirection  = 0.0_WP
  select case (trim(dir(2:2)))
  case ('x')
     secondDirection(1) = 1.0_WP
  case ('y')
     secondDirection(2) = 1.0_WP
  case ('z')
     secondDirection(3) = 1.0_WP
  case default
     call die("reynolds_stress_setup: unknown components '" // trim(dir))
  end select

  ! Get the mean velocity
  allocate(meanVelocity(nGridPoints, nDimensions))
  meanVelocity = 0.0_WP
  if (index(dir,'x').ne.0) then
     call parser_read('mean x velocity', constantVelocity)
     meanVelocity(:,1) = constantVelocity
  end if
  if (index(dir,'y').ne.0) then
     call parser_read('mean y velocity', constantVelocity)
     meanVelocity(:,2) = constantVelocity
  end if
  if (index(dir,'z').ne.0) then
     call parser_read('mean z velocity', constantVelocity)
     meanVelocity(:,3) = constantVelocity
  end if

  return
end subroutine reynolds_stress_setup


! =========================================== !
! Cleanup the Reynolds stress cost functional !
! =========================================== !
subroutine reynolds_stress_cleanup

  ! Internal modules
  use reynolds_stress

  implicit none

  if (allocated(meanVelocity)) deallocate(meanVelocity)

  return
end subroutine reynolds_stress_cleanup


! =========================================== !
! Compute the Reynolds stress cost functional !
! =========================================== !
subroutine reynolds_stress_compute

  ! Internal modules
  use reynolds_stress

  ! External modules
  use geometry
  use simulation_flags
  use solver_options
  use grid
  use grid_functions, only : inner_product
  use state, only : velocity

  implicit none

  ! Local variables
  integer :: i
  real(WP), allocatable :: F(:,:)

  allocate(F(nGridPoints, 1))

  F(:,1) = 0.0_WP
  do i = 1, nGridPoints
     F(i,1) = 0.5_WP * dot_product(velocity(i,:) - meanVelocity(i,:),                        &
          firstDirection(1:nDimensions)) * dot_product(velocity(i,:) - meanVelocity(i,:),    &
          secondDirection(1:nDimensions))
  end do

  instantaneousCostFunctional = inner_product(F(:,1), targetMollifier(:,1))

  deallocate(F)

  return
end subroutine reynolds_stress_compute


! =========================================== !
! Compute the Reynolds stress adjoint forcing !
! =========================================== !
subroutine reynolds_stress_adjoint_source(source)

  ! Internal modules
  use reynolds_stress

  ! External modules
  use simulation_flags
  use solver_options
  use geometry
  use grid
  use grid_patch
  use state

  implicit none

  ! Arguments
  real(WP), dimension(nGridPoints, nUnknowns), intent(inout) :: source

  ! Local variables
  integer :: i, j, k, gridIndex, patchIndex
  real(WP) :: forcingFactor, F
  real(WP), allocatable :: meanVelocityOnPatch(:,:)

  if (useContinuousAdjoint) then
     forcingFactor = 1.0_WP
  else
     forcingFactor = adjointForcingFactor
  end if

  allocate(meanVelocityOnPatch(functionalPatch%nPatchPoints, nDimensions))
  call patch_collect(functionalPatch, meanVelocity, meanVelocityOnPatch)

  do k = functionalPatch%iStart(3), functionalPatch%iEnd(3)
     do j = functionalPatch%iStart(2), functionalPatch%iEnd(2)
        do i = functionalPatch%iStart(1), functionalPatch%iEnd(1)
           gridIndex = i - gridOffset(1) + localGridSize(1) *                                &
                (j - 1 - gridOffset(2) + localGridSize(2) *                                  &
                (k - 1 - gridOffset(3)))
           patchIndex = i - functionalPatch%offset(1) +                                      &
                functionalPatch%localSize(1) * (j - 1 - functionalPatch%offset(2) +          &
                functionalPatch%localSize(2) * (k - 1 - functionalPatch%offset(3)))

           F = - 0.5_WP * forcingFactor * targetMollifier(gridIndex, 1) *                    &
                specificVolume(gridIndex,1) * dot_product(firstDirection(1:nDimensions),     &
                velocity(gridIndex,:) - meanVelocityOnPatch(patchIndex,:))

           source(gridIndex,2:nDimensions+1) = source(gridIndex,2:nDimensions+1) +           &
                secondDirection(1:nDimensions) * F
           source(gridIndex,1) = source(gridIndex,1) -                                       &
                dot_product(secondDirection(1:nDimensions), velocity(gridIndex,:)) * F

           F = - 0.5_WP * forcingFactor * targetMollifier(gridIndex, 1) *                    &
                specificVolume(gridIndex,1) * dot_product(secondDirection(1:nDimensions),    &
                velocity(gridIndex,:) - meanVelocityOnPatch(patchIndex,:))

           source(gridIndex,2:nDimensions+1) = source(gridIndex,2:nDimensions+1) +           &
                firstDirection(1:nDimensions) * F
           source(gridIndex,1) = source(gridIndex,1) -                                       &
                dot_product(firstDirection(1:nDimensions), velocity(gridIndex,:)) * F

        end do
     end do
  end do

  deallocate(meanVelocityOnPatch)

  return
end subroutine reynolds_stress_adjoint_source
