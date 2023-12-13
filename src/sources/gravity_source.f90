module gravity_source

  ! External modules
  use precision
  use solver_options
  use simulation_flags
  use geometry

  implicit none

  real(WP), allocatable, dimension(:) :: gravity

end module gravity_source


! ======================== !
! Setup the gravity source !
! ======================== !
subroutine gravity_setup

  ! Internal modules
  use gravity_source

  ! External modules
  use parser

  implicit none

  ! Local variables
  integer :: i

  if (.not. useGravity) return

  allocate(gravity(1:nDimensions))

  call parser_read('gravity norm', gravity)
  if (abs(sqrt(sum(gravity ** 2)) - 1.0_WP) .gt. epsilon(0.0_WP))                            &
       call die('gravity_source_setup: gravity norm must sum to 1!')
  do i = 1, nDimensions
     if (abs(gravity(i)) .le. epsilon(0.0_WP)) then
        gravity(i) = 0.0_WP
     else
        gravity(i) = froudeNumberInverse / gravity(i)
     end if
  end do

  return
end subroutine gravity_setup


! ========================== !
! Cleanup the gravity source !
! ========================== !
subroutine gravity_cleanup

  ! Internal modules
  use gravity_source

  implicit none

  if (allocated(gravity)) deallocate(gravity)

  return
end subroutine gravity_cleanup


! ================================== !
! Add gravity during the forward run !
! ================================== !
subroutine gravity_forward(source)

  ! Internal modules
  use gravity_source

  ! External modules
  use state

  implicit none

  ! Arguments
  real(WP), dimension(nGridPoints, nUnknowns), intent(inout) :: source

  ! Local variables
  integer :: i

  if (.not. useGravity) return

  do i = 1, nGridPoints
     source(i,2:nDimensions+1) = source(i,2:nDimensions+1) +                                 &
          conservedVariables(i,1) * gravity(1:nDimensions)

     source(i,nDimensions+2) = source(i,nDimensions+2) +                                     &
          dot_product(conservedVariables(i,2:nDimensions+1), gravity(1:nDimensions))
  end do

  return
end subroutine gravity_forward


! ================================== !
! Add gravity during the adjoint run !
! ================================== !
subroutine gravity_adjoint(source)

  ! Internal modules
  use gravity_source

  ! External modules
  use state

  implicit none

  ! Arguments
  real(WP), dimension(nGridPoints, nUnknowns), intent(inout) :: source

  ! Local variables
  integer :: i
  real(WP), allocatable :: localSourceJacobian(:,:), temp(:)

  if (.not. useGravity) return

  allocate(localSourceJacobian(nUnknowns, nUnknowns))
  allocate(temp(nUnknowns))

  localSourceJacobian = 0.0_WP 
  localSourceJacobian(2:nDimensions+1,1) = gravity(1:nDimensions)
  localSourceJacobian(nDimensions+2,2:nDimensions+1) = gravity(1:nDimensions)

  do i = 1, nGridPoints
     temp = matmul(transpose(localSourceJacobian), adjointVariables(i,:))
     source(i,:) = source(i,:) - temp
  end do

  deallocate(localSourceJacobian)
  deallocate(temp)

  return
end subroutine gravity_adjoint
