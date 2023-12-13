module grid_functions

  ! External modules
  use grid

  implicit none

  interface inner_product
     module procedure scalar_inner_product
     module procedure vector_inner_product
  end interface inner_product

  interface gradient
     module procedure gradient_of_scalar
     module procedure gradient_of_vector
  end interface gradient

  interface divergence
     module procedure divergence_of_vector
     module procedure divergence_of_tensor
  end interface divergence

  interface laplacian
     module procedure laplacian_of_scalar
     module procedure laplacian_of_vector
  end interface laplacian

contains

  ! Compute inner product of a scalar
  ! ---------------------------------
  function scalar_inner_product(f, g, weight) result(innerProduct)

    ! External modules
    use parallel

    implicit none

    ! Arguments
    real(WP), intent(in) :: f(:), g(:)
    real(WP), intent(in), optional :: weight(:)

    ! Result
    real(WP) :: innerProduct

    ! Local variables
    real(WP) :: innerProduct_

    innerProduct_ = 0.0_WP

    if (present(weight)) then
       innerProduct_ = innerProduct_ + sum(f * gridNorm(:,1) * g * weight)
    else
       innerProduct_ = innerProduct_ + sum(f * gridNorm(:,1) * g)
    end if

    call parallel_sum(innerProduct_, innerProduct)

    return
  end function scalar_inner_product


  ! Compute inner product of a vector
  ! ---------------------------------
  function vector_inner_product(f, g, weight) result(innerProduct)

    ! External modules
    use parallel

    implicit none

    ! Arguments
    real(WP), intent(in) :: f(:,:), g(:,:)
    real(WP), intent(in), optional :: weight(:)

    ! Result
    real(WP) :: innerProduct

    ! Local variables
    integer :: i
    real(WP) :: innerProduct_

    innerProduct_ = 0.0_WP

    do i = 1, size(f, 2)
       if (present(weight)) then
          innerProduct_ = innerProduct_ + sum(f(:,i) * gridNorm(:,1) * g(:,i) * weight)
       else
          innerProduct_ = innerProduct_ + sum(f(:,i) * gridNorm(:,1) * g(:,i))
       end if
    end do

    call parallel_sum(innerProduct_, innerProduct)

    return
  end function vector_inner_product


  ! Compute gradient of a scalar
  ! ----------------------------
  subroutine gradient_of_scalar(f, gradF)

    ! External modules
    use parallel
    use simulation_flags
    use first_derivative

    ! Arguments
    real(WP), intent(in) :: f(:)
    real(WP), intent(out) :: gradF(:,:)

    ! Local variables
    real(WP), allocatable :: temp(:,:)

    if (size(f) .ne. nGridPoints)                                                            &
         call die ('gradient_of_scalar: size of `f` inconsistent with local grid size')

    if (isDomainCurvilinear .and. nDimensions .gt. 1)                                        &
         allocate(temp(nGridPoints, nDimensions - 1))

    select case (nDimensions)

    case (1)
       gradF(:,1) = f
       call first_derivative_apply(1, gradF(:,1:1))
       gradF(:,1) = jacobian(:,1) * metrics(:,1) * gradF(:,1)

    case (2)

       gradF(:,1) = f
       call first_derivative_apply(1, gradF(:,1:1))

       gradF(:,2) = f
       call first_derivative_apply(2, gradF(:,2:2))

       if (isDomainCurvilinear) then
          temp(:,1) = gradF(:,1)
          gradF(:,1) = jacobian(:,1) * (metrics(:,1) * gradF(:,1) + metrics(:,3) * gradF(:,2))
          gradF(:,2) = jacobian(:,1) * (metrics(:,2) * temp(:,1) + metrics(:,4) * gradF(:,2))
       else
          gradF(:,1) = jacobian(:,1) * metrics(:,1) * gradF(:,1)
          gradF(:,2) = jacobian(:,1) * metrics(:,4) * gradF(:,2)
       end if

    case (3)

       gradF(:,1) = f
       call first_derivative_apply(1, gradF(:,1:1))

       gradF(:,2) = f
       call first_derivative_apply(2, gradF(:,2:2))

       gradF(:,3) = f
       call first_derivative_apply(3, gradF(:,3:3))

       if (isDomainCurvilinear) then
          temp(:,1:2) = gradF(:,1:2)
          gradF(:,1) = jacobian(:,1) * (metrics(:,1) * gradF(:,1) +                          &
               metrics(:,4) * gradF(:,2) + metrics(:,7) * gradF(:,3))
          gradF(:,2) = jacobian(:,1) * (metrics(:,2) * temp(:,1) +                           &
               metrics(:,5) * gradF(:,2) + metrics(:,8) * gradF(:,3))
          gradF(:,3) = jacobian(:,1) * (metrics(:,3) * temp(:,1) +                           &
               metrics(:,6) * temp(:,2) + metrics(:,9) * gradF(:,3))
       else
          gradF(:,1) = jacobian(:,1) * metrics(:,1) * gradF(:,1)
          gradF(:,2) = jacobian(:,1) * metrics(:,5) * gradF(:,2)
          gradF(:,3) = jacobian(:,1) * metrics(:,9) * gradF(:,3)
       end if

    end select

    if (allocated(temp)) deallocate(temp)

  end subroutine gradient_of_scalar


  ! Compute gradient of a vector
  ! ----------------------------
  subroutine gradient_of_vector(f, gradF)

    ! External modules
    use parallel
    use simulation_flags
    use first_derivative

    ! Arguments
    real(WP), intent(in) :: f(:,:)
    real(WP), intent(out) :: gradF(:,:)

    ! Local variables
    real(WP), allocatable :: temp(:,:)

    if (size(f, 1) .ne. nGridPoints)                                                         &
         call die ('gradient_of_vector: size of `f` inconsistent with local grid size')

    select case (nDimensions)

    case (1)

       gradF = f
       call first_derivative_apply(1, gradF(:,1:1))
       gradF(:,1) = jacobian(:,1) * metrics(:,1) * gradF(:,1)

    case (2)

       gradF(:,1:2) = f(:,1:2)
       call first_derivative_apply(1, gradF(:,1:2))

       gradF(:,3:4) = f(:,1:2)
       call first_derivative_apply(2, gradF(:,3:4))

       if (isDomainCurvilinear) then

          allocate(temp(nGridPoints, 2))
          temp = gradF(:,1:2)

          gradF(:,1) = jacobian(:,1) * (metrics(:,1) * gradF(:,1) + metrics(:,3) * gradF(:,3))
          gradF(:,2) = jacobian(:,1) * (metrics(:,2) * temp(:,1) + metrics(:,4) * gradF(:,3))
          gradF(:,3) = jacobian(:,1) * (metrics(:,1) * temp(:,2) + metrics(:,3) * gradF(:,4))
          gradF(:,4) = jacobian(:,1) * (metrics(:,2) * temp(:,2) + metrics(:,4) * gradF(:,4))

       else

          allocate(temp(nGridPoints, 1))
          temp = gradF(:,2:2)

          gradF(:,1) = jacobian(:,1) * metrics(:,1) * gradF(:,1)
          gradF(:,2) = jacobian(:,1) * metrics(:,4) * gradF(:,3)
          gradF(:,3) = jacobian(:,1) * metrics(:,1) * temp(:,1)
          gradF(:,4) = jacobian(:,1) * metrics(:,4) * gradF(:,4)

       end if

    case (3)

       gradF(:,1:3) = f(:,1:3)
       call first_derivative_apply(1, gradF(:,1:3))

       gradF(:,4:6) = f(:,1:3)
       call first_derivative_apply(2, gradF(:,4:6))

       gradF(:,7:9) = f(:,1:3)
       call first_derivative_apply(3, gradF(:,7:9))

       if (isDomainCurvilinear) then

          allocate(temp(nGridPoints, 3))
          temp = gradF(:,1:3)

          gradF(:,1) = jacobian(:,1) * (metrics(:,1) * gradF(:,1) +                          &
               metrics(:,4) * gradF(:,4) + metrics(:,7) * gradF(:,7))
          gradF(:,2) = jacobian(:,1) * (metrics(:,2) * temp(:,1) +                           &
               metrics(:,5) * gradF(:,4) + metrics(:,8) * gradF(:,7))
          gradF(:,3) = jacobian(:,1) * (metrics(:,3) * temp(:,1) +                           &
               metrics(:,6) * gradF(:,4) + metrics(:,9) * gradF(:,7))
          gradF(:,4) = jacobian(:,1) * (metrics(:,1) * temp(:,2) +                           &
               metrics(:,4) * gradF(:,5) + metrics(:,7) * gradF(:,8))
          gradF(:,7) = jacobian(:,1) * (metrics(:,1) * temp(:,3) +                           &
               metrics(:,4) * gradF(:,6) + metrics(:,7) * gradF(:,9))

          temp(:,1) = gradF(:,8)
          gradF(:,8) = jacobian(:,1) * (metrics(:,2) * temp(:,3) +                           &
               metrics(:,5) * gradF(:,6) + metrics(:,8) * gradF(:,9))
          gradF(:,9) = jacobian(:,1) * (metrics(:,3) * temp(:,3) +                           &
               metrics(:,6) * gradF(:,6) + metrics(:,9) * gradF(:,9))

          temp(:,3) = gradF(:,5)
          gradF(:,5) = jacobian(:,1) * (metrics(:,2) * temp(:,2) +                           &
               metrics(:,5) * gradF(:,5) + metrics(:,8) * temp(:,1))
          gradF(:,6) = jacobian(:,1) * (metrics(:,3) * temp(:,2) +                           &
               metrics(:,6) * temp(:,3) + metrics(:,9) * temp(:,1))

       else

          allocate(temp(nGridPoints, 1))

          gradF(:,1) = jacobian(:,1) * metrics(:,1) * gradF(:,1)
          gradF(:,5) = jacobian(:,1) * metrics(:,5) * gradF(:,5)
          gradF(:,9) = jacobian(:,1) * metrics(:,9) * gradF(:,9)

          temp(:,1) = gradF(:,2)
          gradF(:,2) = jacobian(:,1) * metrics(:,5) * gradF(:,4)
          gradF(:,4) = jacobian(:,1) * metrics(:,1) * temp(:,1)

          temp(:,1) = gradF(:,3)
          gradF(:,3) = jacobian(:,1) * metrics(:,9) * gradF(:,7)
          gradF(:,7) = jacobian(:,1) * metrics(:,1) * temp(:,1)

          temp(:,1) = gradF(:,6)
          gradF(:,6) = jacobian(:,1) * metrics(:,9) * gradF(:,8)
          gradF(:,8) = jacobian(:,1) * metrics(:,5) * temp(:,1)

       end if

    end select

    if (allocated(temp)) deallocate(temp)

    return
  end subroutine gradient_of_vector


  ! Compute divergence of a vector
  ! ------------------------------
  subroutine divergence_of_vector(f, divF)

    ! External modules
    use first_derivative

    ! Arguments
    real(WP), intent(in) :: f(:,:)
    real(WP), intent(out) :: divF(:)

    ! Local variables
    real(WP), dimension(:,:), allocatable :: temp

    if (size(f, 1) .ne. nGridPoints)                                                         &
         call die ('divergence_of_vector: size of `f` inconsistent with local grid size')

    allocate(temp(size(f,1), size(f,2)))
    
    ! Transform the vector and take its derivative
    select case (nDimensions)

    case (1)
       temp(:,1) = metrics(:,1) * f(:,1)
       call first_derivative_apply(1, temp(:,1:1))
       divF = temp(:,1) * jacobian(:,1)

    case (2)
       if (isDomainCurvilinear) then
          temp(:,1) = metrics(:,1) * f(:,1) + metrics(:,2) * f(:,2)
          temp(:,2) = metrics(:,3) * f(:,1) + metrics(:,4) * f(:,2)
       else
          temp(:,1) = metrics(:,1) * f(:,1)
          temp(:,2) = metrics(:,4) * f(:,2)
       end if
       call first_derivative_apply(1, temp(:,1:1))
       call first_derivative_apply(2, temp(:,2:2))
       divF = (temp(:,1) + temp(:,2)) * jacobian(:,1)

    case (3)
       if (isDomainCurvilinear) then
          temp(:,1) = metrics(:,1) * f(:,1) + metrics(:,2) * f(:,2) + metrics(:,3) * f(:,3)
          temp(:,2) = metrics(:,4) * f(:,1) + metrics(:,5) * f(:,2) + metrics(:,6) * f(:,3)
          temp(:,3) = metrics(:,7) * f(:,1) + metrics(:,8) * f(:,2) + metrics(:,9) * f(:,3)
       else
          temp(:,1) = metrics(:,1) * f(:,1)
          temp(:,2) = metrics(:,5) * f(:,2)
          temp(:,3) = metrics(:,9) * f(:,3)
       end if
       call first_derivative_apply(1, temp(:,1:1))
       call first_derivative_apply(2, temp(:,2:2))
       call first_derivative_apply(3, temp(:,3:3))
       divF = (temp(:,1) + temp(:,2) + temp(:,3)) * jacobian(:,1)

    end select

    ! Cleanup
    deallocate(temp)

    return
  end subroutine divergence_of_vector


  ! Compute divergence of a tensor
  ! ------------------------------
  subroutine divergence_of_tensor(f, divF)

    ! External modules
    use first_derivative

    ! Arguments
    real(WP), intent(in) :: f(:,:,:)
    real(WP), intent(out) :: divF(:,:)

    ! Local variables
    integer :: i
    real(WP), dimension(:,:,:), allocatable :: temp

    if (size(f, 1) .ne. nGridPoints)                                                         &
         call die ('divergence_of_tensor: size of `f` inconsistent with local grid size')

    allocate(temp(size(f,1), size(f,2), size(f,3)))

    ! Transform the tensor and take its derivative
    select case (nDimensions)

    case (1)
       do i = 1, size(f, 2)
          temp(:,i,1) = metrics(:,1) * f(:,i,1)
       end do
       call first_derivative_apply(1, temp(:,:,1))
       do i = 1, size(f, 2)
          divF(:,i) = temp(:,i,1) * jacobian(:,1)
       end do

    case (2)
       if (isDomainCurvilinear) then
          do i = 1, size(f, 2)
             temp(:,i,1) = metrics(:,1) * f(:,i,1) + metrics(:,2) * f(:,i,2)
             temp(:,i,2) = metrics(:,3) * f(:,i,1) + metrics(:,4) * f(:,i,2)
          end do
       else
          do i = 1, size(f, 2)
             temp(:,i,1) = metrics(:,1) * f(:,i,1)
             temp(:,i,2) = metrics(:,4) * f(:,i,2)
          end do
       end if
       call first_derivative_apply(1, temp(:,:,1))
       call first_derivative_apply(2, temp(:,:,2))
       do i = 1, size(f, 2)
          divF(:,i) = (temp(:,i,1) + temp(:,i,2)) * jacobian(:,1)
       end do

    case (3)
       if (isDomainCurvilinear) then
          do i = 1, size(f, 2)
             temp(:,i,1) = metrics(:,1) * f(:,i,1) +                                         &
                  metrics(:,2) * f(:,i,2) + metrics(:,3) * f(:,i,3)
             temp(:,i,2) = metrics(:,4) * f(:,i,1) +                                         &
                  metrics(:,5) * f(:,i,2) + metrics(:,6) * f(:,i,3)
             temp(:,i,3) = metrics(:,7) * f(:,i,1) +                                         &
                  metrics(:,8) * f(:,i,2) + metrics(:,9) * f(:,i,3)
          end do
       else
          do i = 1, size(f, 2)
             temp(:,i,1) = metrics(:,1) * f(:,i,1)
             temp(:,i,2) = metrics(:,5) * f(:,i,2)
             temp(:,i,3) = metrics(:,9) * f(:,i,3)
          end do
       end if
       call first_derivative_apply(1, temp(:,:,1))
       call first_derivative_apply(2, temp(:,:,2))
       call first_derivative_apply(3, temp(:,:,3))
       do i = 1, size(f, 2)
          divF(:,i) = (temp(:,i,1) + temp(:,i,2) + temp(:,i,3)) * jacobian(:,1)
       end do

    end select

    ! Cleanup
    deallocate(temp)

    return
  end subroutine divergence_of_tensor


  ! Compute Laplacian of a scalar
  ! -----------------------------
  subroutine laplacian_of_scalar(f, lapF)

    ! External modules
    use parallel
    use simulation_flags
    use first_derivative
    use second_derivative

    ! Arguments
    real(WP), intent(in) :: f(:)
    real(WP), intent(out) :: lapF(:)

    ! Local variables
    real(WP), allocatable :: temp(:,:), temp2(:,:,:)

    if (size(f) .ne. nGridPoints)                                                            &
         call die ('laplacian_of_scalar: size of `f` inconsistent with local grid size')

    allocate(temp(size(f,1), nDimensions))
    allocate(temp2(size(f,1), 2, nDimensions))

    select case (nDimensions)

    case (1)
       
       temp(:,1) = f
       call second_derivative_apply(1, temp(:,1:1))

       temp2(:,1,1) = f
       temp2(:,2,1) = jacobian(:,1) * metrics(:,1)
       call first_derivative_apply(1, temp2(:,:,1))
       
       lapF(:) = jacobian(:,1)**2 * metrics(:,1)**2 * temp(:,1) +                            &
            jacobian(:,1) * metrics(:,1) * temp2(:,1,1) * temp2(:,2,1)

    case (2)

       temp(:,1) = f
       call second_derivative_apply(1, temp(:,1:1))
       
       temp(:,2) = f
       call second_derivative_apply(2, temp(:,2:2))

       if (isDomainCurvilinear) then
          
          ! Do we really want to do this?
          call die('laplacian_of_scalar: not yet implemented for curvilinear coordinates')
          
       else

          temp2(:,1,1) = f
          temp2(:,2,1) = jacobian(:,1) * metrics(:,1)
          call first_derivative_apply(1, temp2(:,:,1))

          temp2(:,1,2) = f
          temp2(:,2,2) = jacobian(:,1) * metrics(:,4)
          call first_derivative_apply(2, temp2(:,:,2))
          
          lapF(:) = jacobian(:,1)**2 * (                                                     &
               metrics(:,1)**2 * temp(:,1) + metrics(:,4)**2 * temp(:,2)) +                  &
               jacobian(:,1) * (                                                             &
               metrics(:,1) * temp2(:,1,1) * temp2(:,2,1) +                                  &
               metrics(:,4) * temp2(:,1,2) * temp2(:,2,2))
       end if

    case (3)

       temp(:,1) = f
       call second_derivative_apply(1, temp(:,1:1))

       temp(:,2) = f
       call second_derivative_apply(2, temp(:,2:2))

       temp(:,3) = f
       call second_derivative_apply(3, temp(:,3:3))

       if (isDomainCurvilinear) then
          
          ! Do we really want to do this?
          call die('laplacian_of_scalar: not yet implemented for curvilinear coordinates')
          
       else

          temp2(:,1,1) = f
          temp2(:,2,1) = jacobian(:,1) * metrics(:,1)
          call first_derivative_apply(1, temp2(:,:,1))

          temp2(:,1,2) = f
          temp2(:,2,2) = jacobian(:,1) * metrics(:,5)
          call first_derivative_apply(2, temp2(:,:,2))

          temp2(:,1,3) = f
          temp2(:,2,3) = jacobian(:,1) * metrics(:,9)
          call first_derivative_apply(3, temp2(:,:,3))
          
          lapF = jacobian(:,1)**2 * (                                                        &
               metrics(:,1)**2 * temp(:,1) + metrics(:,5)**2 * temp(:,2) +                   &
               metrics(:,9)**2 * temp(:,3)) + jacobian(:,1) * (                              &
               metrics(:,1) * temp2(:,1,1) * temp2(:,2,1) +                                  &
               metrics(:,5) * temp2(:,1,2) * temp2(:,2,2) +                                  &
               metrics(:,9) * temp2(:,1,3) * temp2(:,2,3))
       end if

    end select

    ! Cleanup
    deallocate(temp, temp2)

  end subroutine laplacian_of_scalar


  ! Compute Laplacian of a vector
  ! -----------------------------
  subroutine laplacian_of_vector(f, lapF)

    ! External modules
    use parallel
    use simulation_flags
    use first_derivative
    use second_derivative

    ! Arguments
    real(WP), intent(in) :: f(:,:)
    real(WP), intent(out) :: lapF(:,:)

    ! Local variables
    integer :: i
    real(WP), allocatable :: temp(:,:,:), temp2(:,:,:), temp3(:,:)

    if (size(f, 1) .ne. nGridPoints)                                                         &
         call die ('laplacian_of_vector: size of `f` inconsistent with local grid size')

    allocate(temp(size(f,1), size(f,2), nDimensions))
    allocate(temp2(size(f,1), size(f,2), nDimensions))
    allocate(temp3(size(f,1), nDimensions))
    
    select case (nDimensions)

    case (1)
       
       temp(:,:,1) = f
       call second_derivative_apply(1, temp(:,:,1))

       temp2(:,:,1) = f
       call first_derivative_apply(1, temp2(:,:,1))
       temp3(:,1) = jacobian(:,1) * metrics(:,1)
       call first_derivative_apply(1, temp3(:,1:1))
       
       do i = 1, size(f, 2)
          lapF(:,i) = jacobian(:,1)**2 * metrics(:,1)**2 * temp(:,i,1) +                     &
               jacobian(:,1) * metrics(:,1) * temp2(:,i,1) * temp3(:,1)
       end do

    case (2)

       temp(:,:,1) = f
       call second_derivative_apply(1, temp(:,:,1))
       
       temp(:,:,2) = f
       call second_derivative_apply(2, temp(:,:,2))

       if (isDomainCurvilinear) then

          ! Do we really want to do this?
          call die('laplacian_of_scalar: not yet implemented for curvilinear coordinates')
          
       else
          
          temp2(:,:,1) = f
          call first_derivative_apply(1, temp2(:,:,1))
          temp2(:,:,2) = f
          call first_derivative_apply(2, temp2(:,:,2))
       
          temp3(:,1) = jacobian(:,1) * metrics(:,1)
          call first_derivative_apply(1, temp3(:,1:1))
          temp3(:,2) = jacobian(:,1) * metrics(:,4)
          call first_derivative_apply(2, temp3(:,2:2))
       
          do i = 1, size(f, 2)
             lapF(:,i) = jacobian(:,1)**2 * (                                                &
                  metrics(:,1)**2 * temp(:,i,1) + metrics(:,4)**2 * temp(:,i,2)) +           &
                  jacobian(:,1) * (                                                          &
                  metrics(:,1) * temp2(:,i,1) * temp3(:,1) +                                 &
                  metrics(:,4) * temp2(:,i,2) * temp3(:,2))
          end do
       end if

    case (3)

       temp(:,:,1) = f
       call second_derivative_apply(1, temp(:,:,1))

       temp(:,:,2) = f
       call second_derivative_apply(2, temp(:,:,2))

       temp(:,:,3) = f
       call second_derivative_apply(3, temp(:,:,3))

       if (isDomainCurvilinear) then

          ! Do we really want to do this?
          call die('laplacian_of_scalar: not yet implemented for curvilinear coordinates')
         
       else

          temp2(:,:,1) = f
          call first_derivative_apply(1, temp2(:,:,1))
          temp2(:,:,2) = f
          call first_derivative_apply(2, temp2(:,:,2))
          temp2(:,:,3) = f
          call first_derivative_apply(3, temp2(:,:,3))
       
          temp3(:,1) = jacobian(:,1) * metrics(:,1)
          call first_derivative_apply(1, temp3(:,1:1))
          temp3(:,2) = jacobian(:,1) * metrics(:,5)
          call first_derivative_apply(2, temp3(:,2:2))
          temp3(:,3) = jacobian(:,1) * metrics(:,9)
          call first_derivative_apply(3, temp3(:,3:3))
          
          do i = 1, size(f, 2)
             lapF(:,i) = jacobian(:,1)**2 * (                                                &
                  metrics(:,1)**2 * temp(:,i,1) + metrics(:,5)**2 * temp(:,i,2) +            &
                  metrics(:,9)**2 * temp(:,i,3)) +                                           &
                  jacobian(:,1) * (                                                          &
                  metrics(:,1) * temp2(:,i,1) * temp3(:,1) +                                 &
                  metrics(:,5) * temp2(:,i,2) * temp3(:,2) +                                 &
                  metrics(:,9) * temp2(:,i,3) * temp3(:,3))
          end do
       end if

    end select

    ! Cleanup
    deallocate(temp, temp2, temp3)

  end subroutine laplacian_of_vector


  ! Solve diffusion implicitly using alternating direct implicit (ADI)
  ! (Second-order finite different with constant coefficient)
  ! ------------------------------------------------------------------
  subroutine implicit_diffusion(f, diffusionAmount)

    ! External modules
    use parallel
    use geometry
    use simulation_flags

    implicit none

    ! Arguments
    real(WP), intent(inout) :: f(:,:)
    real(WP), intent(in) :: diffusionAmount

    ! Local variables
    integer :: i, j, k, ii, jj, kk, nDiagonals, nVariables, no, gridIndex, dir,              &
         periodicOffset(2), nGhost(2)
    real(WP) :: dxi, dxiL, dxiR
    real(WP), dimension(:,:,:,:), allocatable :: A, dx, R1, R2, R3
    logical :: isPeriodicityOverlapping

    ! Determine number of diagonals
    nDiagonals = 3

    ! Number of overlapping points
    no = floor(0.5_WP * real(nDiagonals, WP))

    ! Number of variables
    nVariables = size(f,2)

    ! Allocate RHS arrays
    allocate(R1(iStart(2):iEnd(2),iStart(3):iEnd(3),iStart(1):iEnd(1),nVariables))
    allocate(R2(iStart(1):iEnd(1),iStart(3):iEnd(3),iStart(2):iEnd(2),nVariables))
    allocate(R3(iStart(1):iEnd(1),iStart(2):iEnd(2),iStart(3):iEnd(3),nVariables))

    ! >>>>> Solve in direction '1' <<<<<
    ! ----------------------------------
    dir = 1

    ! Set the number of ghost points
    nGhost = no
    if (.not. isPeriodic(dir) .and. procCoords(dir) .eq. 0) nGhost(1) = 0
    if (.not. isPeriodic(dir) .and. procCoords(dir) .eq. nProcsDir(dir) - 1) nGhost(2) = 0

    ! Set the periodic offset for overlap-type periodicity
    isPeriodicityOverlapping = (periodicityType(dir) .eq. OVERLAP)
    periodicOffset = 0
    if (isPeriodic(dir) .and. isPeriodicityOverlapping) then
       if (isPeriodicityOverlapping .and. procCoords(dir) .eq. 0)                            &
            periodicOffset(2) = 1
       if (isPeriodicityOverlapping .and. procCoords(dir) .eq. nProcsDir(dir) - 1)           &
            periodicOffset(1) = 1
    end if

    ! Allocate arrays
    allocate(A(iStart(2):iEnd(2),iStart(3):iEnd(3),iStart(1):iEnd(1),-no:+no))
    allocate(dx(localGridSize(1) + sum(nGhost), localGridSize(2), localGridSize(3),1))

    ! Store inverse grid spacing
    do k = 1, localGridSize(3)
       do j = 1, localGridSize(2)
          do i = 1, localGridSize(1)
             dx(i + nGhost(1), j, k, 1) = gridSpacing(i + localGridSize(1) *                &
                  (j - 1 + localGridSize(2) * (k - 1)), dir)
          end do
       end do
    end do
    call fill_ghost_points(dx, dir, nGhost, periodicOffset)

    ! Construct the operators
    do k = iStart(3), iEnd(3)
       do j = iStart(2), iEnd(2)
          do i = iStart(1), iEnd(1)
             gridIndex = i - gridOffset(1) + localGridSize(1) *                              &
                  (j - 1 - gridOffset(2) + localGridSize(2) *                                &
                  (k - 1 - gridOffset(3)))
             ii = i - gridOffset(1) + nGhost(1)
             jj = j - gridOffset(2)
             kk = k - gridOffset(3)
             ! Laplacian operator
             if (i.eq.1 .and. nGhost(1).eq.0) then
                dxi = 1.0_WP /  dx(ii,jj,kk,1)
                dxiR= 2.0_WP / (dx(ii,jj,kk,1) + dx(ii+1,jj,kk,1))
                A(j,k,i,-1) = 0.0_WP
                A(j,k,i, 0) = 1.0_WP + 1.0_WP * diffusionAmount * dxi * dxiR
                A(j,k,i,+1) = - diffusionAmount * dxi * dxiR
             else if (i.eq.globalGridSize(1) .and. nGhost(2).eq.0) then
                dxi = 1.0_WP /  dx(ii,jj,kk,1)
                dxiL= 2.0_WP / (dx(ii,jj,kk,1) + dx(ii-1,jj,kk,1))
                A(j,k,i,-1) = - diffusionAmount * dxi * dxiL
                A(j,k,i, 0) = 1.0_WP + 1.0_WP * diffusionAmount * dxi * dxiL
                A(j,k,i,+1) = 0.0_WP
             else
                dxi = 1.0_WP /  dx(ii,jj,kk,1)
                dxiL= 2.0_WP / (dx(ii,jj,kk,1) + dx(ii-1,jj,kk,1))
                dxiR= 2.0_WP / (dx(ii,jj,kk,1) + dx(ii+1,jj,kk,1))
                A(j,k,i,-1) = - diffusionAmount * dxi * dxiL
                A(j,k,i, 0) = 1.0_WP + 2.0_WP * diffusionAmount * dxiL * dxiR
                A(j,k,i,+1) = - diffusionAmount * dxi * dxiR
             end if
             ! RHS
             R1(j,k,i,:) = f(gridIndex,:)
          end do
       end do
    end do

    ! Invert the system
    select case (nDiagonals)
    case (3)
       do i = 1, nVariables
          call tridiagonal(                                                                  &
               A(iStart(2),iStart(3),iStart(1),-1),                                          &
               A(iStart(2),iStart(3),iStart(1), 0),                                          &
               A(iStart(2),iStart(3),iStart(1),+1),                                          &
               R1(:,:,:,i),localGridSize(1),localGridSize(2)*localGridSize(3),dir,isPeriodic)
       end do
    case (5)
       do i = 1, nVariables
          call pentadiagonal(                                                                &
               A(iStart(2),iStart(3),iStart(1),-2),                                          &
               A(iStart(2),iStart(3),iStart(1),-1),                                          &
               A(iStart(2),iStart(3),iStart(1), 0),                                          &
               A(iStart(2),iStart(3),iStart(1),+1),                                          &
               A(iStart(2),iStart(3),iStart(1),+2),                                          &
               R1(:,:,:,i),localGridSize(1),localGridSize(2)*localGridSize(3),dir,isPeriodic)
       end do
    end select

    ! Clean up
    deallocate(A, dx)

    ! >>>>> Solve in direction '2' <<<<<
    ! ----------------------------------
    dir = 2

    ! RHS
    do k = iStart(3), iEnd(3)
       do j = iStart(2), iEnd(2)
          do i = iStart(1), iEnd(1)
             R2(i,k,j,:) = R1(j,k,i,:)
          end do
       end do
    end do

    if (nDimensions.ge.2) then

       ! Set the number of ghost points
       nGhost = no
       if (.not. isPeriodic(dir) .and. procCoords(dir) .eq. 0) nGhost(1) = 0
       if (.not. isPeriodic(dir) .and. procCoords(dir) .eq. nProcsDir(dir) - 1) nGhost(2) = 0

       ! Set the periodic offset for overlap-type periodicity
       isPeriodicityOverlapping = (periodicityType(dir) .eq. OVERLAP)
       periodicOffset = 0
       if (isPeriodic(dir) .and. isPeriodicityOverlapping) then
          if (isPeriodicityOverlapping .and. procCoords(dir) .eq. 0)                         &
               periodicOffset(2) = 1
          if (isPeriodicityOverlapping .and. procCoords(dir) .eq. nProcsDir(dir) - 1)        &
               periodicOffset(1) = 1
       end if

       ! Allocate arrays
       allocate(A(iStart(1):iEnd(1),iStart(3):iEnd(3),iStart(2):iEnd(2),-no:+no))
       allocate(dx(localGridSize(1), localGridSize(2)+sum(nGhost), localGridSize(3), 1))

       ! Store inverse grid spacing
       do k = 1, localGridSize(3)
          do j = 1, localGridSize(2)
             do i = 1, localGridSize(1)
                dx(i, j + nGhost(1), k, 1) = gridSpacing(i + localGridSize(1) *              &
                     (j - 1 + localGridSize(2) * (k - 1)), dir)
             end do
          end do
       end do
       call fill_ghost_points(dx, dir, nGhost, periodicOffset)

       ! Construct the operators
       do k = iStart(3), iEnd(3)
          do j = iStart(2), iEnd(2)
             do i = iStart(1), iEnd(1)
                ii = i - gridOffset(1)
                jj = j - gridOffset(2) + nGhost(1)
                kk = k - gridOffset(3)
                ! Laplacian operator
                if (j.eq.1 .and. nGhost(1).eq.0) then
                   dxi = 1.0_WP /  dx(ii,jj,kk,1)
                   dxiR= 2.0_WP / (dx(ii,jj,kk,1) + dx(ii,jj+1,kk,1))
                   A(i,k,j,-1) = 0.0_WP
                   A(i,k,j, 0) = 1.0_WP + 1.0_WP * diffusionAmount * dxi * dxiR
                   A(i,k,j,+1) = - diffusionAmount * dxi * dxiR
                else if (j.eq.globalGridSize(2) .and. nGhost(2).eq.0) then
                   dxi = 1.0_WP /  dx(ii,jj,kk,1)
                   dxiL= 2.0_WP / (dx(ii,jj,kk,1) + dx(ii,jj-1,kk,1))
                   A(i,k,j,-1) = - diffusionAmount * dxi * dxiL
                   A(i,k,j, 0) = 1.0_WP + 1.0_WP * diffusionAmount * dxi * dxiL
                   A(i,k,j,+1) = 0.0_WP
                else
                   dxi = 1.0_WP /  dx(ii,jj,kk,1)
                   dxiL= 2.0_WP / (dx(ii,jj,kk,1) + dx(ii,jj-1,kk,1))
                   dxiR= 2.0_WP / (dx(ii,jj,kk,1) + dx(ii,jj+1,kk,1))
                   A(i,k,j,-1) = - diffusionAmount * dxi * dxiL
                   A(i,k,j, 0) = 1.0_WP + 2.0_WP * diffusionAmount * dxiL * dxiR
                   A(i,k,j,+1) = - diffusionAmount * dxi * dxiR
                end if
             end do
          end do
       end do

       ! Invert the system
       select case (nDiagonals)
       case (3)
          do i = 1, nVariables
             call tridiagonal(                                                               &
                  A(iStart(1),iStart(3),iStart(2),-1),                                       &
                  A(iStart(1),iStart(3),iStart(2), 0),                                       &
                  A(iStart(1),iStart(3),iStart(2),+1),                                       &
                  R2(:,:,:,i),localGridSize(2),localGridSize(1)*localGridSize(3),dir,        &
                  isPeriodic)
          end do
       case (5)
          do i = 1, nVariables
             call pentadiagonal(                                                             &
                  A(iStart(1),iStart(3),iStart(2),-2),                                       &
                  A(iStart(1),iStart(3),iStart(2),-1),                                       &
                  A(iStart(1),iStart(3),iStart(2), 0),                                       &
                  A(iStart(1),iStart(3),iStart(2),+1),                                       &
                  A(iStart(1),iStart(3),iStart(2),+2),                                       &
                  R2(:,:,:,i),localGridSize(2),localGridSize(1)*localGridSize(3),dir,        &
                  isPeriodic)
          end do
       end select

       ! Clean up
       deallocate(A, dx)

    end if

    ! >>>>> Solve in direction '3' <<<<<
    ! ----------------------------------
    dir = 3

    ! RHS
    do k = iStart(3), iEnd(3)
       do j = iStart(2), iEnd(2)
          do i = iStart(1), iEnd(1)
             R3(i,j,k,:) = R2(i,k,j,:)
          end do
       end do
    end do

    if (nDimensions.ge.3) then

       ! Set the number of ghost points
       nGhost = no
       if (.not. isPeriodic(dir) .and. procCoords(dir) .eq. 0) nGhost(1) = 0
       if (.not. isPeriodic(dir) .and. procCoords(dir) .eq. nProcsDir(dir) - 1) nGhost(2) = 0

       ! Set the periodic offset for overlap-type periodicity
       isPeriodicityOverlapping = (periodicityType(dir) .eq. OVERLAP)
       periodicOffset = 0
       if (isPeriodic(dir) .and. isPeriodicityOverlapping) then
          if (isPeriodicityOverlapping .and. procCoords(dir) .eq. 0)                         &
               periodicOffset(2) = 1
          if (isPeriodicityOverlapping .and. procCoords(dir) .eq. nProcsDir(dir) - 1)        &
               periodicOffset(1) = 1
       end if

       ! Allocate arrays
       allocate(A(iStart(1):iEnd(1),iStart(2):iEnd(2),iStart(3):iEnd(3),-no:+no))
       allocate(dx(localGridSize(1), localGridSize(2), localGridSize(3) + sum(nGhost), 1))

       ! Store inverse grid spacing
       do k = 1, localGridSize(3)
          do j = 1, localGridSize(2)
             do i = 1, localGridSize(1)
                dx(i, j, k + nGhost(1), 1) = gridSpacing(i + localGridSize(1) *              &
                     (j - 1 + localGridSize(2) * (k - 1)), dir)
             end do
          end do
       end do
       call fill_ghost_points(dx, dir, nGhost, periodicOffset)

       ! Construct the operators
       do k = iStart(3), iEnd(3)
          do j = iStart(2), iEnd(2)
             do i = iStart(1), iEnd(1)
                ii = i - gridOffset(1)
                jj = j - gridOffset(2)
                kk = k - gridOffset(3) + nGhost(1)
                ! Laplacian operator
                if (k.eq.1 .and. nGhost(1).eq.0) then
                   dxi = 1.0_WP /  dx(ii,jj,kk,1)
                   dxiR= 2.0_WP / (dx(ii,jj,kk,1) + dx(ii,jj,kk+1,1))
                   A(i,j,k,-1) = 0.0_WP
                   A(i,j,k, 0) = 1.0_WP + 1.0_WP * diffusionAmount * dxi * dxiR
                   A(i,j,k,+1) = - diffusionAmount * dxi * dxiR
                else if (k.eq.globalGridsize(3) .and. nGhost(2).eq.0) then
                   dxi = 1.0_WP /  dx(ii,jj,kk,1)
                   dxiL= 2.0_WP / (dx(ii,jj,kk,1) + dx(ii,jj,kk-1,1))
                   A(i,j,k,-1) = - diffusionAmount * dxi * dxiL
                   A(i,j,k, 0) = 1.0_WP + 1.0_WP * diffusionAmount * dxi * dxiL
                   A(i,j,k,+1) = - 0.0_WP
                else
                   dxi = 1.0_WP /  dx(ii,jj,kk,1)
                   dxiL= 2.0_WP / (dx(ii,jj,kk,1) + dx(ii,jj,kk-1,1))
                   dxiR= 2.0_WP / (dx(ii,jj,kk,1) + dx(ii,jj,kk+1,1))
                   A(i,j,k,-1) = - diffusionAmount * dxi * dxiL
                   A(i,j,k, 0) = 1.0_WP + 2.0_WP * diffusionAmount * dxiL * dxiR
                   A(i,j,k,+1) = - diffusionAmount * dxi * dxiR
                end if
             end do
          end do
       end do

       ! Invert the system
       select case (nDiagonals)
       case (3)
          do i = 1, nVariables
             call tridiagonal(                                                               &
                  A(iStart(1),iStart(2),iStart(3),-1),                                       &
                  A(iStart(1),iStart(2),iStart(3), 0),                                       &
                  A(iStart(1),iStart(2),iStart(3),+1),                                       &
                  R3(:,:,:,i),localGridSize(3),localGridSize(1)*localGridSize(2),dir,        &
                  isPeriodic)
          end do
       case (5)
          do i = 1, nVariables
             call pentadiagonal(                                                             &
                  A(iStart(1),iStart(2),iStart(3),-2),                                       &
                  A(iStart(1),iStart(2),iStart(3),-1),                                       &
                  A(iStart(1),iStart(2),iStart(3), 0),                                       &
                  A(iStart(1),iStart(2),iStart(3),+1),                                       &
                  A(iStart(1),iStart(2),iStart(3),+2),                                       &
                  R3(:,:,:,i),localGridSize(3),localGridSize(1)*localGridSize(2),dir,        &
                  isPeriodic)
          end do
       end select

       ! Clean up
       deallocate(A, dx)

    end if

    ! Update f
    do k = iStart(3), iEnd(3)
       do j = iStart(2), iEnd(2)
          do i = iStart(1), iEnd(1)
             gridIndex = i - gridOffset(1) + localGridSize(1) *                              &
                  (j - 1 - gridOffset(2) + localGridSize(2) *                                &
                  (k - 1 - gridOffset(3)))
             f(gridIndex,:)=R3(i,j,k,:)
          end do
       end do
    end do

    ! Clean up
    if (allocated(R1)) deallocate(R1)
    if (allocated(R2)) deallocate(R2)
    if (allocated(R3)) deallocate(R3)

    return
  end subroutine implicit_diffusion


  ! Compute a coordinate derivative
  ! -------------------------------
  subroutine compute_coordinate_derivative(direction, coordinateDerivatives)

    ! External modules
    use parallel
    use first_derivative

    implicit none

    ! Arguments
    integer, intent(in) :: direction
    real(WP), intent(inout) :: coordinateDerivatives(:,:)

    ! Local variables
    integer :: i, j, k, l, gridSizeWithGhostPoints(3), numGhostPointsBegin(3)
    real(WP), allocatable :: xWithGhostPoints(:,:,:,:)

    ! Straightforward if periodicity type is not `PLANE`
    if (periodicityType(direction) .ne. PLANE) then
       coordinateDerivatives = coordinates
       call first_derivative_apply(direction, coordinateDerivatives)
       return
    end if

    ! Special hack require for periodicity type `PLANE`:

    ! Find the grid size including ghost points, and the number of ghost points at the
    ! beginning (which is also the offset of the physical grid points)
    gridSizeWithGhostPoints = localGridSize
    gridSizeWithGhostPoints(direction) = gridSizeWithGhostPoints(direction) +                &
         sum(firstDerivative(direction)%nGhost)
    numGhostPointsBegin = 0
    numGhostPointsBegin(direction) = firstDerivative(direction)%nGhost(1)

    ! Allocate an array to hold both physical and ghost point coordinates
    allocate(xWithGhostPoints(gridSizeWithGhostPoints(1), gridSizeWithGhostPoints(2),        &
         gridSizeWithGhostPoints(3), nDimensions))

    ! Copy coordinates at physical grid points to the ghost array
    do l = 1, nDimensions
       do k = 1, localGridSize(3)
          do j = 1, localGridSize(2)
             do i = 1, localGridSize(1)
                xWithGhostPoints(i + numGhostPointsBegin(1),                                 &
                     j + numGhostPointsBegin(2), k + numGhostPointsBegin(3), l) =            &
                     coordinates(i + localGridSize(1) * (j - 1 +                             &
                     localGridSize(2) * (k - 1)), l)
             end do
          end do
       end do
    end do

    ! Exchange ghost points between MPI processes
    call fill_ghost_points(xWithGhostPoints, direction,                                      &
         firstDerivative(direction)%nGhost,                                                  &
         firstDerivative(direction)%periodicOffset)

    ! At the first process, subtract the periodic length from coordinates received from
    ! the last process
    if (procCoords(direction) .eq. 0) then
       select case (direction)
       case (1)
          xWithGhostPoints(1:numGhostPointsBegin(1),:,:,1) =                                 &
               xWithGhostPoints(1:numGhostPointsBegin(1),:,:,1) -                            &
               periodicLength(1)
       case (2)
          xWithGhostPoints(:,1:numGhostPointsBegin(2),:,2) =                                 &
               xWithGhostPoints(:,1:numGhostPointsBegin(2),:,2) -                            &
               periodicLength(2)
       case (3)
          xWithGhostPoints(:,:,1:numGhostPointsBegin(3),3) =                                 &
               xWithGhostPoints(:,:,1:numGhostPointsBegin(3),3) -                            &
               periodicLength(3)
       end select
    end if

    ! At the last process, add the periodic length from coordinates received from the
    ! first process
    if (procCoords(direction) .eq. nProcsDir(direction) - 1) then
       select case (direction)
       case (1)
          xWithGhostPoints(gridSizeWithGhostPoints(1) + 1 -                                  &
               numGhostPointsBegin(1) : gridSizeWithGhostPoints(1),:,:,1)                    &
               = xWithGhostPoints(gridSizeWithGhostPoints(1) + 1 -                           &
               numGhostPointsBegin(1) : gridSizeWithGhostPoints(1),:,:,1) +                  &
               periodicLength(1)
       case (2)
          xWithGhostPoints(:,gridSizeWithGhostPoints(2) + 1 -                                &
               numGhostPointsBegin(2) : gridSizeWithGhostPoints(2),:,2)                      &
               = xWithGhostPoints(:,gridSizeWithGhostPoints(2) + 1 -                         &
               numGhostPointsBegin(2) : gridSizeWithGhostPoints(2),:,2) +                    &
               periodicLength(2)
       case (3)
          xWithGhostPoints(:,:,gridSizeWithGhostPoints(3) + 1 -                              &
               numGhostPointsBegin(3) : gridSizeWithGhostPoints(3),3)                        &
               = xWithGhostPoints(:,:,gridSizeWithGhostPoints(3) + 1 -                       &
               numGhostPointsBegin(3) : gridSizeWithGhostPoints(3),3) +                      &
               periodicLength(3)
       end select
    end if

    ! Stencil needs to be applied only at interior points
    call first_derivative_apply_interior(direction, xWithGhostPoints, coordinateDerivatives)

    deallocate(xWithGhostPoints)

    return
  end subroutine compute_coordinate_derivative


  ! Find the minimum value and return its location
  ! ----------------------------------------------
  subroutine find_minimum(f, fMin, iMin, jMin, kMin)

    ! External modules
    use parallel

    implicit none

    ! Arguments
    real(WP), intent(in) :: f(:)
    real(WP), intent(out) :: fMin
    integer, intent(out), optional :: iMin, jMin, kMin

    ! Local variables
    integer :: i, j, k, minIndex(3), ierror
    real(WP), allocatable :: minValues(:)
    integer, allocatable :: minIndices(:,:)
    real(WP) :: a
    real(WP) :: minValue

    allocate(minValues(nProcs), minIndices(3, nProcs))

    minValue = huge(0.0_WP)
    do k = iStart(3), iEnd(3)
       do j = iStart(2), iEnd(2)
          do i = iStart(1), iEnd(1)
             a = f(i - gridOffset(1) + localGridSize(1) * (j - 1 - gridOffset(2) +           &
                  localGridSize(2) * (k - 1 - gridOffset(3))))
             if (a < minValue) then
                minIndex = (/ i, j, k /)
                minValue = a
             end if
          end do
       end do
    end do

    if (present(iMin) .or. present(jMin) .or. present(kMin)) then
       call MPI_Allgather(minIndex, 3, MPI_INTEGER, minIndices, 3, MPI_INTEGER, comm, ierror)
    end if
    call MPI_Allgather(minValue, 1, MPI_REAL_WP, minValues, 1, MPI_REAL_WP, comm, ierror)

    i = minloc(minValues, 1)
    fMin = minValues(i)
    if (present(iMin)) iMin = minIndices(1,i)
    if (present(jMin)) jMin = minIndices(2,i)
    if (present(kMin)) kMin = minIndices(3,i)

    deallocate(minValues)
    deallocate(minIndices)

    return
  end subroutine find_minimum


  ! Find the maximum value and return its location
  ! ----------------------------------------------
  subroutine find_maximum(f, fMax, iMax, jMax, kMax)

    ! External modules
    use parallel

    implicit none

    ! Arguments
    real(WP), intent(in) :: f(:)
    real(WP), intent(out) :: fMax
    integer, intent(out), optional :: iMax, jMax, kMax

    ! Local variables
    integer :: i, j, k, maxIndex(3), ierror
    real(WP), allocatable :: maxValues(:)
    integer, allocatable :: maxIndices(:,:)
    real(WP) :: a
    real(WP) :: maxValue

    allocate(maxValues(nProcs), maxIndices(3, nProcs))

    maxValue = - huge(0.0_WP)
    do k = iStart(3), iEnd(3)
       do j = iStart(2), iEnd(2)
          do i = iStart(1), iEnd(1)
             a = f(i - gridOffset(1) + localGridSize(1) * (j - 1 - gridOffset(2) +           &
                  localGridSize(2) * (k - 1 - gridOffset(3))))
             if (a > maxValue) then
                maxIndex = (/ i, j, k /)
                maxValue = a
             end if
          end do
       end do
    end do

    if (present(iMax) .or. present(jMax) .or. present(kMax)) then
       call MPI_Allgather(maxIndex, 3, MPI_INTEGER, maxIndices, 3, MPI_INTEGER, comm, ierror)
    end if
    call MPI_Allgather(maxValue, 1, MPI_REAL_WP, maxValues, 1, MPI_REAL_WP, comm, ierror)

    i = maxloc(maxValues, 1)
    fMax = maxValues(i)
    if (present(iMax)) iMax = maxIndices(1,i)
    if (present(jMax)) jMax = maxIndices(2,i)
    if (present(kMax)) kMax = maxIndices(3,i)

    deallocate(maxValues)
    deallocate(maxIndices)

    return
  end subroutine find_maximum


  ! Determine if variable is within some range
  ! ------------------------------------------
  function isVariableWithinRange(f, fOutsideRange, iOutsideRange, jOutsideRange,          &
       kOutsideRange, minValue, maxValue)

    implicit none

    ! Arguments
    real(WP), intent(in) :: f(:)
    real(WP), intent(out) :: fOutsideRange
    integer, intent(out) :: iOutsideRange, jOutsideRange, kOutsideRange
    real(WP), intent(in), optional :: minValue, maxValue
    logical :: isVariableWithinRange

    ! Local variables
    real(WP) :: fMin, fMax

    isVariableWithinRange = .true.
    if (.not. present(minValue) .and. .not. present(maxValue)) return

    if (present(minValue)) then
       call find_minimum(f, fMin, iOutsideRange, jOutsideRange, kOutsideRange)
       if (real(fMin, WP) .le. minValue) then
          fOutsideRange = fMin
          isVariableWithinRange = .false.
          return
       end if
    end if

    if (present(maxValue)) then
       call find_maximum(f, fMax, iOutsideRange, jOutsideRange, kOutsideRange)
       if (real(fMax, WP) .ge. maxValue) then
          fOutsideRange = fMax
          isVariableWithinRange = .false.
       end if
    end if

    return
  end function isVariableWithinRange


  ! Correct distance in periodic directions
  ! ---------------------------------------
  subroutine correct_periodic_distance(distance)

    implicit none

    ! Arguments
    real(WP), intent(inout) :: distance(:)

    ! Local variables
    integer :: i
    real(WP) :: buf

    do i = 1, nDimensions
       if (isPeriodic(i)) then
          buf = distance(i)
          if (abs(distance(i)) .gt. abs(distance(i) - periodicLength(i)))                    &
               distance(i) = distance(i) - periodicLength(i)
          if (abs(distance(i)) .gt. abs(buf + periodicLength(i)))                            &
               distance(i) = buf + periodicLength(i)
       end if
    end do

    return
  end subroutine correct_periodic_distance


  ! Transform fluxes from Cartesian to computational coordinates
  ! ------------------------------------------------------------
  subroutine transform_fluxes(fluxes, transformedFluxes)

    ! External modules
    use simulation_flags

    implicit none

    ! Arguments
    real(WP), intent(in) :: fluxes(:,:,:)
    real(WP), intent(out) :: transformedFluxes(:,:,:)

    ! Local variables
    integer :: i

    select case (nDimensions)

    case (1)
       do i = 1, size(fluxes, 2)
          transformedFluxes(:,i,1) = metrics(:,1) * fluxes(:,i,1)
       end do

    case (2)
       if (isDomainCurvilinear) then
          do i = 1, size(fluxes, 2)
             transformedFluxes(:,i,1) = metrics(:,1) * fluxes(:,i,1) +                         &
                  metrics(:,2) * fluxes(:,i,2)
             transformedFluxes(:,i,2) = metrics(:,3) * fluxes(:,i,1) +                         &
                  metrics(:,4) * fluxes(:,i,2)
          end do
       else
          do i = 1, size(fluxes, 2)
             transformedFluxes(:,i,1) = metrics(:,1) * fluxes(:,i,1)
             transformedFluxes(:,i,2) = metrics(:,4) * fluxes(:,i,2)
          end do
       end if

    case (3)
       if (isDomainCurvilinear) then
          do i = 1, size(fluxes, 2)
             transformedFluxes(:,i,1) = metrics(:,1) * fluxes(:,i,1) +                         &
                  metrics(:,2) * fluxes(:,i,2) + metrics(:,3) * fluxes(:,i,3)
             transformedFluxes(:,i,2) = metrics(:,4) * fluxes(:,i,1) +                         &
                  metrics(:,5) * fluxes(:,i,2) + metrics(:,6) * fluxes(:,i,3)
             transformedFluxes(:,i,3) = metrics(:,7) * fluxes(:,i,1) +                         &
                  metrics(:,8) * fluxes(:,i,2) + metrics(:,9) * fluxes(:,i,3)
          end do
       else
          do i = 1, size(fluxes, 2)
             transformedFluxes(:,i,1) = metrics(:,1) * fluxes(:,i,1)
             transformedFluxes(:,i,2) = metrics(:,5) * fluxes(:,i,2)
             transformedFluxes(:,i,3) = metrics(:,9) * fluxes(:,i,3)
          end do
       end if

    end select

    return
  end subroutine transform_fluxes


  ! Get grid index in lexicographical ordering
  ! ------------------------------------------
  function grid_index(i, j, k) result(l)

    ! Internal modules
    use geometry

    implicit none

    ! Arguments
    integer, intent(in) :: i, j, k

    ! Result
    integer :: l

    l = i - gridOffset(1) + localGridSize(1) *                                                 &
         (j - 1 - gridOffset(2) + localGridSize(2) *                                           &
         (k - 1 - gridOffset(3)))

    return
  end function grid_index

end module grid_functions


! ===================================== !
! Compute directional distances to wall !
! ===================================== !
subroutine get_distance_to_wall(distanceToWall)

  ! Internal modules
  use grid

  ! External modules
  use parallel

  implicit none

  ! Arguments
  integer, dimension(nGridPoints, nDimensions, 2), intent(out) :: distanceToWall

  ! Local parameters
  integer :: i, j, k, n, dir, gridIndex, ierror
  integer, dimension(:), allocatable :: iblank1D

  ! Direction '1'
  dir = 1
  allocate(iblank1D(GlobalGridSize(dir)))

  do k = iStart(3), iEnd(3)
     do j = iStart(2), iEnd(2)
        ! Initialize distance and populate 1D iblank
        iblank1D = 0
        do i = iStart(1), iEnd(1)
           gridIndex = i - gridOffset(1) + localGridSize(1) *                                &
                (j - 1 - gridOffset(2) + localGridSize(2) *                                  &
                (k - 1 - gridOffset(3)))
           distanceToWall(gridIndex, dir, :) = globalGridSize(dir)
           iblank1D(i) = iblank(gridIndex)
        end do
        ! Communicate 1D iblank
        call MPI_ALLREDUCE(MPI_IN_PLACE, iblank1D, globalGridSize(dir), MPI_INTEGER,         &
             MPI_SUM, commDir(dir), ierror)
        ! <<< Left <<<
        do i = iStart(1), iEnd(1)
           loopxL:do n=0,i-1
              if (iblank1D(i-n).eq.0) exit loopxL
           end do loopxL
           gridIndex = i - gridOffset(1) + localGridSize(1) *                                &
                (j - 1 - gridOffset(2) + localGridSize(2) *                                  &
                (k - 1 - gridOffset(3)))
           if (n.ne.i) then
              ! Distance to iblank
              distanceToWall(gridIndex,dir,1) = n
           else if (.not. isPeriodic(dir)) then
              ! Distance to domain boundary
              distanceToWall(gridIndex,dir,1) = i
           end if
        end do
        ! Correct for periodicity
        if (isPeriodic(dir)) then
           do i = iStart(1), iEnd(1)
              loopxLp:do n=globalGridSize(dir),i+1,-1
                 if (iblank1D(n).eq.0) exit loopxLp
              end do loopxLp
              gridIndex = i - gridOffset(1) + localGridSize(1) *                             &
                   (j - 1 - gridOffset(2) + localGridSize(2) *                               &
                   (k - 1 - gridOffset(3)))
              if (n.ne.i) distanceToWall(gridIndex,dir,1) =                                  &
                   min(distanceToWall(gridIndex,dir,1), globalGridSize(dir) + i - n)
           end do
        end if

        ! >>> Right >>>
        do i = iStart(1), iEnd(1)
           loopxR:do n=0,globalGridSize(dir)-i
              if (iblank1D(i+n).eq.0) exit loopxR
           end do loopxR
           gridIndex = i - gridOffset(1) + localGridSize(1) *                                &
                (j - 1 - gridOffset(2) + localGridSize(2) *                                  &
                (k - 1 - gridOffset(3)))
           if (n.ne.globalGridSize(dir)-i+1) then
              ! Distance to iblank
              distanceToWall(gridIndex,dir,2) = n
           else if (.not. isPeriodic(dir)) then
              ! Distance to domain boundary
              distanceToWall(gridIndex,dir,2) = globalGridSize(dir)-i+1
           end if
        end do
        ! Correct for periodicity
        if (isPeriodic(dir)) then
           do i = iStart(1), iEnd(1)
              loopxRp:do n=1,i-1
                 if (iblank1D(n).eq.0) exit loopxRp
              end do loopxRp
              gridIndex = i - gridOffset(1) + localGridSize(1) *                             &
                   (j - 1 - gridOffset(2) + localGridSize(2) *                               &
                   (k - 1 - gridOffset(3)))
              if (n.ne.i) distanceToWall(gridIndex,dir,2) =                                  &
                   min(distanceToWall(gridIndex,dir,2), globalGridSize(dir) + n - i)
           end do
        end if
     end do
  end do

  ! Cleanup
  deallocate(iblank1D)

  ! Direction '2'
  if (nDimensions.ge.2) then
     dir = 2
     allocate(iblank1D(GlobalGridSize(dir)))

     do k = iStart(3), iEnd(3)
        do i = iStart(1), iEnd(1)
           ! Initialize distance and populate 1D iblank
           iblank1D = 0
           do j = iStart(2), iEnd(2)
              gridIndex = i - gridOffset(1) + localGridSize(1) *                             &
                   (j - 1 - gridOffset(2) + localGridSize(2) *                               &
                   (k - 1 - gridOffset(3)))
              distanceToWall(gridIndex, dir, :) = globalGridSize(dir)
              iblank1D(j) = iblank(gridIndex)
           end do
           ! Communicate 1D iblank
           call MPI_ALLREDUCE(MPI_IN_PLACE, iblank1D, globalGridSize(dir), MPI_INTEGER,      &
                MPI_SUM, commDir(dir), ierror)
           ! <<< Left <<<
           do j = iStart(2), iEnd(2)
              loopyL:do n=0,j-1
                 if (iblank1D(j-n).eq.0) exit loopyL
              end do loopyL
              gridIndex = i - gridOffset(1) + localGridSize(1) *                             &
                   (j - 1 - gridOffset(2) + localGridSize(2) *                               &
                   (k - 1 - gridOffset(3)))
              if (n.ne.j) then
                 ! Distance to iblank
                 distanceToWall(gridIndex,dir,1) = n
              else if (.not. isPeriodic(dir)) then
                 ! Distance to domain boundary
                 distanceToWall(gridIndex,dir,1) = j
              end if
           end do
           ! Correct for periodicity
           if (isPeriodic(dir)) then
              do j = iStart(2), iEnd(2)
                 loopyLp:do n=globalGridSize(dir),j+1,-1
                    if (iblank1D(n).eq.0) exit loopyLp
                 end do loopyLp
                 gridIndex = i - gridOffset(1) + localGridSize(1) *                          &
                      (j - 1 - gridOffset(2) + localGridSize(2) *                            &
                      (k - 1 - gridOffset(3)))
                 if (n.ne.j) distanceToWall(gridIndex,dir,1) =                               &
                      min(distanceToWall(gridIndex,dir,1), globalGridSize(dir) + j - n)
              end do
           end if

           ! >>> Right >>>
           do j = iStart(2), iEnd(2)
              loopyR:do n=0,globalGridSize(dir)-j
                 if (iblank1D(j+n).eq.0) exit loopyR
              end do loopyR
              gridIndex = i - gridOffset(1) + localGridSize(1) *                             &
                   (j - 1 - gridOffset(2) + localGridSize(2) *                               &
                   (k - 1 - gridOffset(3)))
              if (n.ne.globalGridSize(dir)-j+1) then
                 ! Distance to iblank
                 distanceToWall(gridIndex,dir,2) = n
              else if (.not. isPeriodic(dir)) then
                 ! Distance to domain boundary
                 distanceToWall(gridIndex,dir,2) = globalGridSize(dir)-j+1
              end if
           end do
           ! Correct for periodicity
           if (isPeriodic(dir)) then
              do j = iStart(2), iEnd(2)
                 loopyRp:do n=1,j-1
                    if (iblank1D(n).eq.0) exit loopyRp
                 end do loopyRp
                 gridIndex = i - gridOffset(1) + localGridSize(1) *                          &
                      (j - 1 - gridOffset(2) + localGridSize(2) *                            &
                      (k - 1 - gridOffset(3)))
                 if (n.ne.j) distanceToWall(gridIndex,dir,2) =                               &
                      min(distanceToWall(gridIndex,dir,2), globalGridSize(dir) + n - j)
              end do
           end if
        end do
     end do

     ! Cleanup
     deallocate(iblank1D)

  end if

  ! Direction '3'
  if (nDimensions.ge.3) then
     dir = 3
     allocate(iblank1D(GlobalGridSize(dir)))

     do j = iStart(2), iEnd(2)
        do i = iStart(1), iEnd(1)
           ! Initialize distance and populate 1D iblank
           iblank1D = 0
           do k = iStart(3), iEnd(3)
              gridIndex = i - gridOffset(1) + localGridSize(1) *                             &
                   (j - 1 - gridOffset(2) + localGridSize(2) *                               &
                   (k - 1 - gridOffset(3)))
              distanceToWall(gridIndex, dir, :) = globalGridSize(dir)
              iblank1D(k) = iblank(gridIndex)
           end do
           ! Communicate 1D iblank
           call MPI_ALLREDUCE(MPI_IN_PLACE, iblank1D, globalGridSize(dir), MPI_INTEGER,      &
                MPI_SUM, commDir(dir), ierror)
           ! <<< Left <<<
           do k = iStart(3), iEnd(3)
              loopzL:do n=0,k-1
                 if (iblank1D(k-n).eq.0) exit loopzL
              end do loopzL
              gridIndex = i - gridOffset(1) + localGridSize(1) *                             &
                   (j - 1 - gridOffset(2) + localGridSize(2) *                               &
                   (k - 1 - gridOffset(3)))
              if (n.ne.k) then
                 ! Distance to iblank
                 distanceToWall(gridIndex,dir,1) = n
              else if (.not. isPeriodic(dir)) then
                 ! Distance to domain boundary
                 distanceToWall(gridIndex,dir,1) = k
              end if
           end do
           ! Correct for periodicity
           if (isPeriodic(dir)) then
              do k = iStart(3), iEnd(3)
                 loopzLp:do n=globalGridSize(dir),k+1,-1
                    if (iblank1D(n).eq.0) exit loopzLp
                 end do loopzLp
                 gridIndex = i - gridOffset(1) + localGridSize(1) *                          &
                      (j - 1 - gridOffset(2) + localGridSize(2) *                            &
                      (k - 1 - gridOffset(3)))
                 if (n.ne.k) distanceToWall(gridIndex,dir,1) =                               &
                      min(distanceToWall(gridIndex,dir,1), globalGridSize(dir) + k - n)
              end do
           end if

           ! >>> Right >>>
           do k = iStart(3), iEnd(3)
              loopzR:do n=0,globalGridSize(dir)-k
                 if (iblank1D(k+n).eq.0) exit loopzR
              end do loopzR
              gridIndex = i - gridOffset(1) + localGridSize(1) *                             &
                   (j - 1 - gridOffset(2) + localGridSize(2) *                               &
                   (k - 1 - gridOffset(3)))
              if (n.ne.globalGridSize(dir)-k+1) then
                 ! Distance to iblank
                 distanceToWall(gridIndex,dir,2) = n
              else if (.not. isPeriodic(dir)) then
                 ! Distance to domain boundary
                 distanceToWall(gridIndex,dir,2) = globalGridSize(dir)-k+1
              end if
           end do
           ! Correct for periodicity
           if (isPeriodic(dir)) then
              do k = iStart(3), iEnd(3)
                 loopzRp:do n=1,k-1
                    if (iblank1D(n).eq.0) exit loopzRp
                 end do loopzRp
                 gridIndex = i - gridOffset(1) + localGridSize(1) *                          &
                      (j - 1 - gridOffset(2) + localGridSize(2) *                            &
                      (k - 1 - gridOffset(3)))
                 if (n.ne.k) distanceToWall(gridIndex,dir,2) =                               &
                      min(distanceToWall(gridIndex,dir,2), globalGridSize(dir) + n - k)
              end do
           end if
        end do
     end do

     ! Cleanup
     deallocate(iblank1D)
  end if

  return
end subroutine get_distance_to_wall


! ========================== !
! Compute volume of the grid !
! ========================== !
subroutine compute_grid_volume(volume)

  ! Internal modules
  use grid_functions

  ! External modules
  use parallel
  
  implicit none

  ! Arguments
  real(WP), intent(out) :: volume

  ! Local variables
  real(WP), dimension(nGridPoints) :: F

  F = 1.0_WP
  volume = inner_product(F, F)
  
  return
end subroutine compute_grid_volume
