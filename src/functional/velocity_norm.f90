module velocity_norm

  ! External modules
  use functional

  implicit none

  integer  :: power, velocityComponent
  real(WP) :: meanVelocity

contains

  subroutine verify_velocity_norm_patch(patch)

    ! External modules
    use solver_options

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

    n = nDimensions
    if (n .ne. nDimensions) then
       write(message, '(2(A,I0.0),A)') 'verify_velocity_norm_patch: Expected a ',         &
            nDimensions, 'D patch, but extent represents a ', n, 'D patch!'
       call die(trim(message))
    end if

    return
  end subroutine verify_velocity_norm_patch

end module velocity_norm


! ======================================= !
! Setup the velocity norm cost functional !
! ======================================= !
subroutine velocity_norm_setup

  ! Internal modules
  use velocity_norm

  ! External modules
  use parser
  use solver_options

  implicit none

  ! Verify the patch type
  call verify_velocity_norm_patch(functionalPatch)

  ! Get the norm power
  call parser_read('norm power', power)
  if (power .le. 1)                                                                          &
       call die('velocity_norm_setup: norm power must be > 1') 

  ! Get the velocity component
  call parser_read('velocity component', velocityComponent)
  if (velocityComponent .lt. 1)                                                              &
       call die('velocity_norm_setup: velocity component must be > 0')
  if (velocityComponent .gt. nDimensions)                                                    &
       call die('velocity_norm_setup: velocity component must be <= nDimensions')

  ! Get the mean velocity
  call parser_read('mean velocity', meanVelocity, 0.0_WP)

  return
end subroutine velocity_norm_setup


! ========================================= !
! Cleanup the velocity norm cost functional !
! ========================================= !
subroutine velocity_norm_cleanup

  ! Internal modules
  use velocity_norm

  implicit none

  ! Nothing to do

  return
end subroutine velocity_norm_cleanup


! ========================================= !
! Compute the velocity norm cost functional !
! ========================================= !
subroutine velocity_norm_compute

  ! Internal modules
  use velocity_norm

  ! External modules
  use grid_functions, only : inner_product
  use state
  use onestep

  implicit none

  ! Local variables
  real(WP), allocatable :: F(:)

  allocate(F(nGridPoints))

  ! Define the velocity norm as < (V - <V>)^power >
  F = (velocity(:, velocityComponent) - meanVelocity)**power

  instantaneousCostFunctional = inner_product(F, targetMollifier(:,1))

  deallocate(F)

  return
end subroutine velocity_norm_compute


! ========================================= !
! Compute the velocity norm adjoint forcing !
! ========================================= !
subroutine velocity_norm_adjoint_source(source)

  ! Internal modules
  use velocity_norm

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

           source(gridIndex, 1) = source(gridIndex, 1) + F * (-real(power,WP)) *             &
                velocity(gridIndex, velocityComponent) *                                     &   
                (velocity(gridIndex, velocityComponent) - meanVelocity)**(power-1)
           source(gridIndex, velocityComponent+1) = source(gridIndex, velocityComponent+1) + &
                F * real(power,WP) *                                                         &
                (velocity(gridIndex, velocityComponent) - meanVelocity)**(power-1)

        end do
     end do
  end do

  return
end subroutine velocity_norm_adjoint_source
