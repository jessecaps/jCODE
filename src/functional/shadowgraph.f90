module shadowgraph

  ! External modules
  use functional

  implicit none

  integer, parameter ::                                                                      &
       SHADOWGRAPH_GRAD = 1,                                                                 &
       SHADOWGRAPH_LAP  = 2

  integer :: shadowgraphMethod
  real(WP), allocatable :: buffer(:,:)

contains

  subroutine verify_shadowgraph_patch(patch)

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
       write(message, '(2(A,I0.0),A)') 'verify_shadowgraph_patch: Expected a ',              &
            nDimensions, 'D patch, but extent represents a ', n, 'D patch!'
       call die(trim(message))
    end if

    return
  end subroutine verify_shadowgraph_patch

end module shadowgraph


! ===================================== !
! Setup the shadowgraph cost functional !
! ===================================== !
subroutine shadowgraph_setup

  ! Internal modules
  use shadowgraph

  ! External modules
  use parser
  use state

  implicit none

  ! Local variables
  character(len = str_medium) :: key

  ! Verify the patch type
  call verify_shadowgraph_patch(functionalPatch)

  call parser_read('shadowgraph method', key, 'gradient')
  select case (trim(key))
  case ('gradient', 'grad', 'GRAD')
     shadowgraphMethod = SHADOWGRAPH_GRAD
     allocate(buffer(nGridPoints, nDimensions))
  case ('laplacian', 'Laplacian', 'lap')
     shadowgraphMethod = SHADOWGRAPH_LAP
     allocate(buffer(nGridPoints, 1))
  case default
     call die("shadowgraph_setup: Unknown shadowgraph method:" // trim(key))
  end select

  return
end subroutine shadowgraph_setup


! ======================================= !
! Cleanup the shadowgraph cost functional !
! ======================================= !
subroutine shadowgraph_cleanup

  ! Internal modules
  use shadowgraph

  implicit none

  deallocate(buffer)

  return
end subroutine shadowgraph_cleanup


! ======================================= !
! Compute the shadowgraph cost functional !
! ======================================= !
subroutine shadowgraph_compute

  ! Internal modules
  use shadowgraph

  ! External modules
  use geometry
  use grid_functions, only : inner_product, gradient, laplacian
  use state, only : conservedVariables!, specificVolume

  implicit none

  ! Local variables
  !integer :: i

  select case (shadowgraphMethod)

  case (SHADOWGRAPH_GRAD)

     call gradient(log(conservedVariables(:,1)), buffer)
     ! Divide the density gradient by density
     !do i = 1, nDimensions
     !   buffer(:,i) = specificVolume(:,1) * buffer(:,i)
     !end do

  case (SHADOWGRAPH_LAP)
     call laplacian(conservedVariables(:,1), buffer(:,1))

  end select

  instantaneousCostFunctional = inner_product(buffer, buffer, targetMollifier(:,1))

  return
end subroutine shadowgraph_compute


! ======================================= !
! Compute the shadowgraph adjoint forcing !
! ======================================= !
subroutine shadowgraph_adjoint_source(source)

  ! Internal modules
  use shadowgraph

  ! External modules
  use simulation_flags
  use solver_options
  use geometry
  use grid
  use grid_functions, only : gradient, laplacian
  use grid_patch
  use state, only : conservedVariables, specificVolume
  use first_derivative, only: adjoint_first_derivative_apply
  use second_derivative, only: adjoint_second_derivative_apply

  implicit none

  ! Arguments
  real(WP), dimension(nGridPoints, nUnknowns), intent(inout) :: source

  ! Local variables
  integer :: i
  real(WP) :: forcingFactor
  real(WP), allocatable :: temp(:,:), temp2(:,:)

  if (useContinuousAdjoint) then
     forcingFactor = 1.0_WP
  else
     forcingFactor = adjointForcingFactor
  end if

  allocate(temp(nGridPoints, nDimensions))

  select case (shadowgraphMethod)

  case (SHADOWGRAPH_GRAD)

     call gradient(log(conservedVariables(:,1)), buffer)
     do i = 1, nDimensions
        temp(:,i) = targetMollifier(:,1) * buffer(:,i)
        call adjoint_first_derivative_apply(i, temp(:,i:i))
     end do

     ! Multiply temp by metrics
     select case (nDimensions)

     case (1)
        temp(:,1) = metrics(:,1) * temp(:,1)

     case (2)
        if (isDomainCurvilinear) then
           allocate(temp2(nGridPoints, nDimensions - 1))
           temp2(:,1) = temp(:,1)
           temp(:,1) = metrics(:,1) * temp(:,1) + metrics(:,3) * temp(:,2)
           temp(:,2) = metrics(:,2) * temp2(:,1) + metrics(:,4) * temp(:,2)
           deallocate(temp2)
        else
           temp(:,1) = metrics(:,1) * temp(:,1)
           temp(:,2) = metrics(:,4) * temp(:,2)
        end if

     case (3)
        if (isDomainCurvilinear) then
           allocate(temp2(nGridPoints, nDimensions - 1))
           temp2(:,1:2) = temp(:,1:2)
           temp(:,1) = metrics(:,1) * temp(:,1) + metrics(:,4) * temp(:,2) +                    &
                metrics(:,7) * temp(:,3)
           temp(:,2) = metrics(:,2) * temp2(:,1) + metrics(:,5) * temp(:,2) +                   &
                metrics(:,8) * temp(:,3)
           temp(:,3) = metrics(:,3) * temp2(:,1) + metrics(:,6) * temp2(:,2) +                  &
                metrics(:,9) * temp(:,3)
           deallocate(temp2)
        else
           temp(:,1) = metrics(:,1) * temp(:,1)
           temp(:,2) = metrics(:,5) * temp(:,2)
           temp(:,3) = metrics(:,9) * temp(:,3)
        end if

     end select

     temp(:,1) = jacobian(:,1) * sum(temp, dim = 2) * specificVolume(:,1)

  case (SHADOWGRAPH_LAP)

     call laplacian(conservedVariables(:,1), buffer(:,1))
     do i = 1, nDimensions
        temp(:,i) = jacobian(:,1) * targetMollifier(:,1) * buffer(:,1)
        call adjoint_second_derivative_apply(i, temp(:,i:i))
     end do

     ! Multiply temp by metric^2
     select case (nDimensions)

     case (1)
        temp(:,1) = metrics(:,1)**2 * temp(:,1)

     case (2)
        if (isDomainCurvilinear) then
           temp(:,1) = metrics(:,1)**2 * temp(:,1) + metrics(:,3)**2 * temp(:,2) +           &
                metrics(:,2)**2 * temp(:,1) + metrics(:,4)**2 * temp(:,2)
        else
           temp(:,1) = metrics(:,1)**2 * temp(:,1) + metrics(:,4)**2 * temp(:,2)
        end if

     case (3)
        if (isDomainCurvilinear) then
           temp(:,1) =                                                                       &
                metrics(:,1)**2 * temp(:,1) + metrics(:,4)**2 * temp(:,2) +                  &
                metrics(:,7)**2 * temp(:,3) +                                                &
                metrics(:,2)**2 * temp(:,1) + metrics(:,5)**2 * temp(:,2) +                  &
                metrics(:,8)**2 * temp(:,3) +                                                &
                metrics(:,3)**2 * temp(:,1) + metrics(:,6)**2 * temp(:,2) +                  &
                metrics(:,9)**2 * temp(:,3)
        else
           temp(:,1) = metrics(:,1)**2 * temp(:,1) + metrics(:,5)**2 * temp(:,2) +           &
                metrics(:,9)**2 * temp(:,3)
        end if

     end select

     temp(:,1) = jacobian(:,1) * temp(:,1)

  end select

  source(:,1) = source(:,1) - 2.0_WP * forcingFactor * temp(:,1)

  deallocate(temp)

  return
end subroutine shadowgraph_adjoint_source
