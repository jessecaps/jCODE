module operator

  ! External modules
  use precision
  use string
  use geometry

  implicit none

  integer, parameter ::                                                                      &
       SYMMETRIC      = 0,                                                                   &
       SKEW_SYMMETRIC = 1,                                                                   &
       ASYMMETRIC     = 2

  character(len = str_medium) :: discretizationType
  integer, allocatable :: distanceToWall(:,:,:)

  type :: t_StencilOperator

     integer :: symmetryType, interiorWidth, boundaryWidth, boundaryDepth,                   &
          nGhost(2), periodicOffset(2), nDiagonals
     logical :: hasDomainBoundary(2), implicit
     real(WP), allocatable :: normBoundary(:), rhsInterior(:), rhsBoundary1(:,:),            &
          rhsBoundary2(:,:), lhsInterior(:), lhsBoundary1(:,:), lhsBoundary2(:,:)

  end type t_StencilOperator

contains

  subroutine allocate_operator(stencilOperator)

    implicit none

    ! Arguments
    type(t_StencilOperator), intent(inout) :: stencilOperator

    if (stencilOperator%symmetryType .ne. ASYMMETRIC) then
       allocate(stencilOperator%rhsInterior(                                                 &
            -stencilOperator%interiorWidth/2:stencilOperator%interiorWidth/2))
       stencilOperator%rhsInterior = 0.0_WP
    end if
    allocate(stencilOperator%rhsBoundary1(stencilOperator%boundaryWidth,                     &
         stencilOperator%boundaryDepth))
    stencilOperator%rhsBoundary1 = 0.0_WP
    allocate(stencilOperator%rhsBoundary2(stencilOperator%boundaryWidth,                     &
         stencilOperator%boundaryDepth))
    stencilOperator%rhsBoundary2 = 0.0_WP
    allocate(stencilOperator%normBoundary(stencilOperator%boundaryDepth))
    stencilOperator%normBoundary = 1.0_WP
    if (stencilOperator%implicit) then
       allocate(stencilOperator%lhsinterior(                                                 &
            -stencilOperator%nDiagonals/2:stencilOperator%nDiagonals/2))
       stencilOperator%lhsInterior = 0.0_WP
       allocate(stencilOperator%lhsBoundary1(                                                &
            -stencilOperator%nDiagonals/2:stencilOperator%nDiagonals/2,                      &
            stencilOperator%boundaryDepth))
       stencilOperator%lhsBoundary1 = 0.0_WP
       allocate(stencilOperator%lhsBoundary2(                                                &
            -stencilOperator%nDiagonals/2:stencilOperator%nDiagonals/2,                      &
            stencilOperator%boundaryDepth))
       stencilOperator%lhsBoundary2 = 0.0_WP
    else
       stencilOperator%nDiagonals = 0
    end if

    return
  end subroutine allocate_operator


  subroutine deallocate_operator(stencilOperator)

    implicit none

    ! Arguments
    type(t_StencilOperator), intent(inout) :: stencilOperator

    if (allocated(stencilOperator%rhsInterior))  deallocate(stencilOperator%rhsInterior)
    if (allocated(stencilOperator%rhsBoundary1)) deallocate(stencilOperator%rhsBoundary1)
    if (allocated(stencilOperator%rhsBoundary2)) deallocate(stencilOperator%rhsBoundary2)
    if (allocated(stencilOperator%normBoundary)) deallocate(stencilOperator%normBoundary)
    if (allocated(stencilOperator%lhsInterior))  deallocate(stencilOperator%lhsInterior)
    if (allocated(stencilOperator%lhsBoundary1)) deallocate(stencilOperator%lhsBoundary1)
    if (allocated(stencilOperator%lhsBoundary2)) deallocate(stencilOperator%lhsBoundary2)

    return
  end subroutine deallocate_operator


  subroutine operator_update(stencilOperator, dir)

    ! External modules
    use parallel

    implicit none

    ! Arguments
    type(t_StencilOperator), intent(inout) :: stencilOperator
    integer, intent(in) :: dir

    ! Local variables
    logical :: isPeriodicityOverlapping

    ! Does this process contain the boundary of the computational domain?
    stencilOperator%hasDomainBoundary(1) = (procCoords(dir) .eq. 0)
    stencilOperator%hasDomainBoundary(2) = (procCoords(dir) .eq. nProcsDir(dir) - 1)
    stencilOperator%hasDomainBoundary(1) = stencilOperator%hasDomainBoundary(1) .and.        &
         .not. isPeriodic(dir)
    stencilOperator%hasDomainBoundary(2) = stencilOperator%hasDomainBoundary(2) .and.        &
         .not. isPeriodic(dir)

    ! Set the number of ghost points
    stencilOperator%nGhost = max(abs(lbound(stencilOperator%rhsInterior,1)),                 &
         abs(ubound(stencilOperator%rhsInterior,1)))
   
    if (.not. isPeriodic(dir) .and. procCoords(dir) .eq. 0) stencilOperator%nGhost(1) = 0
    if (.not. isPeriodic(dir) .and. procCoords(dir) .eq. nProcsDir(dir) - 1)                 &
         stencilOperator%nGhost(2) = 0

    ! Set the periodic offset for overlap-type periodicity
    isPeriodicityOverlapping = (periodicityType(dir) .eq. OVERLAP)
    stencilOperator%periodicOffset = 0
    if (isPeriodic(dir) .and. isPeriodicityOverlapping) then
       if (isPeriodicityOverlapping .and. procCoords(dir) .eq. 0)                            &
            stencilOperator%periodicOffset(2) = 1
       if (isPeriodicityOverlapping .and. procCoords(dir) .eq. nProcsDir(dir) - 1)           &
            stencilOperator%periodicOffset(1) = 1
    end if
    ! Hack to handle parallel communication in fill_ghost_points routine
    if (periodicityType(dir) .eq. POLAR) then
       if (procCoords(dir) .eq. 0) stencilOperator%periodicOffset(1) = 2
       if (procCoords(dir) .eq. nProcsDir(dir) - 1) stencilOperator%periodicOffset(2) = 2
    end if

    return
  end subroutine operator_update


  subroutine operator_adjoint(stencilOperator, adjointOperator)

    implicit none

    ! Arguments
    type(t_StencilOperator), intent(in) :: stencilOperator
    type(t_StencilOperator), intent(inout) :: adjointOperator

    ! Local variables
    integer :: i, j

    ! Copy basic information
    adjointOperator%implicit = stencilOperator%implicit
    adjointOperator%symmetryType = stencilOperator%symmetryType
    adjointOperator%interiorWidth = stencilOperator%interiorWidth
    adjointOperator%boundaryDepth = stencilOperator%boundaryWidth
    adjointOperator%boundaryWidth = stencilOperator%boundaryWidth +                          &
         stencilOperator%interiorWidth / 2
    call allocate_operator(adjointOperator)

    ! Reverse the interior stencil
    do i = - stencilOperator%interiorWidth / 2, stencilOperator%interiorWidth / 2
       adjointOperator%rhsInterior(i) = stencilOperator%rhsInterior(-i)
    end do

    ! Copy `normBoundary`
    if (allocated(adjointOperator%normBoundary)) deallocate(adjointOperator%normBoundary)
    allocate(adjointOperator%normBoundary(stencilOperator%boundaryDepth))
    adjointOperator%normBoundary = stencilOperator%normBoundary

    ! Transpose the left-boundary coefficients
    adjointOperator%rhsBoundary1 = 0.0_WP
    adjointOperator%rhsBoundary1(1:stencilOperator%boundaryDepth,:) =                        &
         transpose(stencilOperator%rhsBoundary1)
    do i = stencilOperator%boundaryDepth + 1, stencilOperator%boundaryWidth +                &
         stencilOperator%interiorWidth / 2
       do j = - stencilOperator%interiorWidth / 2, stencilOperator%interiorWidth / 2
          if (i + j > stencilOperator%boundaryWidth) exit
          adjointOperator%rhsBoundary1(i,i+j) = stencilOperator%rhsInterior(j)
       end do
    end do

    ! Pre-multiply by the inverse of the norm matrix
    do i = 1, adjointOperator%boundaryWidth
       adjointOperator%rhsBoundary1(i,1:size(adjointOperator%normBoundary)) =                &
            adjointOperator%rhsBoundary1(i,1:size(adjointOperator%normBoundary)) /           &
            adjointOperator%normBoundary
    end do

    ! Post-multiply by the norm matrix
    do i = 1, adjointOperator%boundaryDepth
       adjointOperator%rhsBoundary1(1:size(adjointOperator%normBoundary),i) =                &
            adjointOperator%rhsBoundary1(1:size(adjointOperator%normBoundary),i) *           &
            adjointOperator%normBoundary
    end do

    ! Fill the right-boundary coefficients
    select case (adjointOperator%symmetryType)
    case (SYMMETRIC)
       adjointOperator%rhsBoundary2(1:adjointOperator%boundaryWidth,:) =                     &
            +adjointOperator%rhsBoundary1(adjointOperator%boundaryWidth:1:-1,:)
    case (SKEW_SYMMETRIC)
       adjointOperator%rhsBoundary2(1:adjointOperator%boundaryWidth,:) =                     &
            -adjointOperator%rhsBoundary1(adjointOperator%boundaryWidth:1:-1,:)
    end select
    if (allocated(adjointOperator%lhsBoundary1) .and.                                        &
         allocated(adjointOperator%lhsBoundary2)) then
       adjointOperator%lhsBoundary2(                                                         &
            -adjointOperator%nDiagonals/2:adjointOperator%ndiagonals/2,:) =                  &
            adjointOperator%lhsBoundary1(                                                    &
            adjointOperator%nDiagonals/2:-adjointOperator%nDiagonals/2:-1,:)
    end if

    return
  end subroutine operator_adjoint


  subroutine operator_transpose(stencilOperator, transposeOperator, isQ)

    implicit none

    ! Arguments
    type(t_StencilOperator), intent(in) :: stencilOperator
    type(t_StencilOperator), intent(inout) :: transposeOperator
    logical, intent(in), optional :: isQ

    ! Local variables
    integer :: i
    logical :: isQ_ = .false.

    if (present(isQ)) isQ_ = isQ

    ! Copy basic information
    transposeOperator%implicit = .false.
    transposeOperator%symmetryType = stencilOperator%symmetryType
    transposeOperator%interiorWidth = stencilOperator%interiorWidth
    transposeOperator%boundaryDepth = stencilOperator%boundaryDepth
    transposeOperator%boundaryWidth = stencilOperator%boundaryWidth

    call allocate_operator(transposeOperator)

    ! Reverse the interior stencil
    if (allocated(transposeOperator%rhsInterior)) deallocate(transposeOperator%rhsInterior)
    allocate(transposeOperator%rhsInterior(-ubound(stencilOperator%rhsInterior,1):           &
         -lbound(stencilOperator%rhsInterior,1)))

    do i = lbound(stencilOperator%rhsInterior,1), ubound(stencilOperator%rhsInterior,1)
       transposeOperator%rhsInterior(-i) =  stencilOperator%rhsInterior(i)
    end do

    ! Copy normBoundary
    transposeOperator%normBoundary = stencilOperator%normBoundary

    ! Transpose the left-boundary coefficients (only square block)
    transposeOperator%rhsBoundary1 = 0.0_WP
    transposeOperator%rhsBoundary1(1:stencilOperator%boundaryDepth,                          &
         1:stencilOperator%boundaryDepth) =                                                  &
         transpose(stencilOperator%rhsBoundary1(1:stencilOperator%boundaryDepth,             &
         1:stencilOperator%boundaryDepth))

    ! Fill the boundary coefficients
    call operator_fill_rhsBoundary(transposeOperator, isQ_)

    return
  end subroutine operator_transpose


  subroutine operator_fill_rhsBoundary(stencilOperator, isQ)

    implicit none

    ! Arguments
    type(t_StencilOperator), intent(inout) :: stencilOperator
    logical, intent(in), optional :: isQ

    ! Local variables
    integer :: i, s
    logical :: isQ_ = .false.

    if (present(isQ)) isQ_ = isQ

    select case (stencilOperator%symmetryType)

    case (SYMMETRIC)
       stencilOperator%rhsBoundary2(1:stencilOperator%boundaryWidth,:) =                     &
            +stencilOperator%rhsBoundary1(stencilOperator%boundaryWidth:1:-1,:)

    case (SKEW_SYMMETRIC)
       stencilOperator%rhsBoundary2(1:stencilOperator%boundaryWidth,:) =                     &
            -stencilOperator%rhsBoundary1(stencilOperator%boundaryWidth:1:-1,:)

    case (ASYMMETRIC)
       s = stencilOperator%boundaryWidth - stencilOperator%boundaryDepth

       ! fill the square block of rhsBoundary2
       stencilOperator%rhsBoundary2(s+1:stencilOperator%boundaryWidth,                       &
            1:stencilOperator%boundaryDepth) =                                               &
            stencilOperator%rhsBoundary1(1:stencilOperator%boundaryDepth,                    &
            stencilOperator%boundaryDepth:1:-1)

       stencilOperator%rhsBoundary2(s+1:stencilOperator%boundaryWidth,                       &
            1:stencilOperator%boundaryDepth) =                                               &
            transpose( stencilOperator%rhsBoundary2(s+1:stencilOperator%boundaryWidth,       &
            1:stencilOperator%boundaryDepth) )

       if (isQ_) then

          ! complete rhsBoundary1 with interior diagonal terms
          do i = 1, ubound(stencilOperator%rhsInterior,1)
             stencilOperator%rhsBoundary1(stencilOperator%boundaryDepth+1:                   &
                  stencilOperator%boundaryDepth+i,                                           &
                  stencilOperator%boundaryDepth-ubound(stencilOperator%rhsInterior,1)+i) =   &
                  stencilOperator%rhsInterior(ubound(stencilOperator%rhsInterior,1)-i+1:     &
                  ubound(stencilOperator%rhsInterior,1))
          end do

          ! complete rhsBoundary2 with interior diagonal terms
          do i = 1, abs(lbound(stencilOperator%rhsInterior,1))
             stencilOperator%rhsBoundary2(s+lbound(stencilOperator%rhsInterior,1)+i:s,       &
                  stencilOperator%boundaryDepth+1-i) =                                       &
                  stencilOperator%rhsInterior(lbound(stencilOperator%rhsInterior,1):-i)
          end do

       else
          ! if operator D=H^(-1) Q

          ! complete rhsBoundary1 with interior diagonal terms
          do i = 1, ubound(stencilOperator%rhsInterior,1)
             stencilOperator%rhsBoundary1(stencilOperator%boundaryDepth+1:                   &
                  stencilOperator%boundaryDepth+i,                                           &
                  stencilOperator%boundaryDepth-ubound(stencilOperator%rhsInterior,1)+i) =   &
                  stencilOperator%rhsInterior(ubound(stencilOperator%rhsInterior,1)-i+1:     &
                  ubound(stencilOperator%rhsInterior,1)) /                                   &
                  stencilOperator%normBoundary(i)
          end do

          ! complete rhsBoundary2 with interior diagonal terms
          do i = 1, abs(lbound(stencilOperator%rhsInterior,1))
             stencilOperator%rhsBoundary2(s+lbound(stencilOperator%rhsInterior,1)+i:s,       &
                  stencilOperator%boundaryDepth+1-i) =                                       &
                  stencilOperator%rhsInterior(lbound(stencilOperator%rhsInterior,1):-i) /    &
                  stencilOperator%normBoundary(i)
          end do

       end if

    end select

    return
  end subroutine operator_fill_rhsBoundary


  subroutine operator_multiply_inverse_diagonal_norm(stencilOperator, factor)

    implicit none

    ! Arguments
    type(t_StencilOperator), intent(inout) :: stencilOperator
    real(WP), intent(in), optional :: factor

    ! Local variables
    integer :: i
    real(WP) :: factor_

    if (present(factor)) then
       factor_ = factor
    else
       factor_ = 1.0_WP
    end if

    ! Multiply rhsBoundary1
    do i = 1, stencilOperator%boundaryDepth
       stencilOperator%rhsBoundary1(:,i) = stencilOperator%rhsBoundary1(:,i) /               &
            stencilOperator%normBoundary(i) * factor_
    end do

    ! Multiply rhsBoundary2
    do i = 1, stencilOperator%boundaryDepth
       stencilOperator%rhsBoundary2(:,i) = stencilOperator%rhsBoundary2(:,i) /               &
            stencilOperator%normBoundary(i) * factor_
    end do

    ! Multiply interior (if factor is provided)
    if (present(factor)) then
       stencilOperator%rhsInterior(:) = stencilOperator%rhsInterior(:) * factor_
    end if

  end subroutine operator_multiply_inverse_diagonal_norm


  subroutine operator_apply_1(stencilOperator, x)

    ! External modules
    use parallel

    implicit none

    ! Arguments
    type(t_StencilOperator), intent(in) :: stencilOperator
    real(WP), intent(inout) :: x(:,:)

    ! Local variables
    integer :: i, j, k, l, m, n
    real(WP), allocatable :: xWithGhostPoints(:,:,:,:)

    allocate(xWithGhostPoints(localGridSize(1) + sum(stencilOperator%nGhost),                &
         localGridSize(2), localGridSize(3), size(x, 2)))

    do l = 1, size(x, 2)
       do k = 1, localGridSize(3)
          do j = 1, localGridSize(2)
             do i = 1, localGridSize(1)
                xWithGhostPoints(i + stencilOperator%nGhost(1), j, k, l) =                   &
                     x(i + localGridSize(1) * (j - 1 + localGridSize(2) * (k - 1)), l)
             end do
          end do
       end do
    end do

    call fill_ghost_points(xWithGhostPoints, 1, stencilOperator%nGhost,                      &
         stencilOperator%periodicOffset)

    call operator_apply_interior_1(stencilOperator, xWithGhostPoints, x)

    n = stencilOperator%boundaryWidth

    ! Left boundary
    if (stencilOperator%hasDomainBoundary(1)) then
       i = stencilOperator%nGhost(1)
       do l = 1, size(x, 2)
          do k = 1, localGridSize(3)
             do j = 1, localGridSize(2)
                do m = 1, stencilOperator%boundaryDepth
                   x(i + m - stencilOperator%nGhost(1) +                                     &
                        localGridSize(1) * (j - 1 + localGridSize(2) * (k - 1)), l) =        &
                        sum(stencilOperator%rhsBoundary1(:,m) *                              &
                        xWithGhostPoints(i+1:i+n,j,k,l))
                end do
             end do
          end do
       end do
    end if

    ! Right boundary
    if (stencilOperator%hasDomainBoundary(2)) then
       i = stencilOperator%nGhost(1) + localGridSize(1) + 1
       do l = 1, size(x, 2)
          do k = 1, localGridSize(3)
             do j = 1, localGridSize(2)
                do m = stencilOperator%boundaryDepth, 1, -1
                   x(i - m - stencilOperator%nGhost(1) +                                     &
                        localGridSize(1) * (j - 1 + localGridSize(2) * (k - 1)), l) =        &
                        sum(stencilOperator%rhsBoundary2(:,m) *                              &
                        xWithGhostPoints(i-n:i-1,j,k,l))
                end do
             end do
          end do
       end do
    end if

    deallocate(xWithGhostPoints)

    return
  end subroutine operator_apply_1


  subroutine operator_apply_2(stencilOperator, x)

    ! External modules
    use parallel

    implicit none

    ! Arguments
    type(t_StencilOperator), intent(in) :: stencilOperator
    real(WP), intent(inout) :: x(:,:)

    ! Local variables
    integer :: i, j, k, l, m, n
    real(WP), allocatable :: xWithGhostPoints(:,:,:,:)

    allocate(xWithGhostPoints(localGridSize(1), localGridSize(2) +                           &
         sum(stencilOperator%nGhost), localGridSize(3), size(x, 2)))

    do l = 1, size(x, 2)
       do k = 1, localGridSize(3)
          do j = 1, localGridSize(2)
             do i = 1, localGridSize(1)
                xWithGhostPoints(i, j + stencilOperator%nGhost(1), k, l) =                   &
                     x(i + localGridSize(1) * (j - 1 + localGridSize(2) * (k - 1)), l)
             end do
          end do
       end do
    end do

    call fill_ghost_points(xWithGhostPoints, 2, stencilOperator%nGhost,                      &
         stencilOperator%periodicOffset)

    call operator_apply_interior_2(stencilOperator, xWithGhostPoints, x)

    n = stencilOperator%boundaryWidth

    ! Left boundary
    if (stencilOperator%hasDomainBoundary(1)) then
       j = stencilOperator%nGhost(1)
       do l = 1, size(x, 2)
          do k = 1, localGridSize(3)
             do i = 1, localGridSize(1)
                do m = 1, stencilOperator%boundaryDepth
                   x(i + localGridSize(1) * (j - 1 + m - stencilOperator%nGhost(1) +         &
                        localGridSize(2) * (k - 1)), l) =                                    &
                        sum(stencilOperator%rhsBoundary1(:,m) *                              &
                        xWithGhostPoints(i,j+1:j+n,k,l))
                end do
             end do
          end do
       end do
    end if

    ! Right boundary
    if (stencilOperator%hasDomainBoundary(2)) then
       j = stencilOperator%nGhost(1) + localGridSize(2) + 1
       do l = 1, size(x, 2)
          do k = 1, localGridSize(3)
             do i = 1, localGridSize(1)
                do m = stencilOperator%boundaryDepth, 1, -1
                   x(i + localGridSize(1) * (j - 1 - m - stencilOperator%nGhost(1) +         &
                        localGridSize(2) * (k - 1)), l) =                                    &
                        sum(stencilOperator%rhsBoundary2(:,m) *                              &
                        xWithGhostPoints(i,j-n:j-1,k,l))
                end do
             end do
          end do
       end do
    end if

    deallocate(xWithGhostPoints)

    return
  end subroutine operator_apply_2


  subroutine operator_apply_3(stencilOperator, x)

    ! External modules
    use parallel

    implicit none

    ! Arguments
    type(t_StencilOperator), intent(in) :: stencilOperator
    real(WP), intent(inout) :: x(:,:)

    ! Local variables
    integer :: i, j, k, l, m, n
    real(WP), allocatable :: xWithGhostPoints(:,:,:,:)

    allocate(xWithGhostPoints(localGridSize(1), localGridSize(2),                            &
         localGridSize(3) + sum(stencilOperator%nGhost), size(x, 2)))

    do l = 1, size(x, 2)
       do k = 1, localGridSize(3)
          do j = 1, localGridSize(2)
             do i = 1, localGridSize(1)
                xWithGhostPoints(i, j, k + stencilOperator%nGhost(1), l) =                   &
                     x(i + localGridSize(1) * (j - 1 + localGridSize(2) * (k - 1)), l)
             end do
          end do
       end do
    end do

    call fill_ghost_points(xWithGhostPoints, 3, stencilOperator%nGhost,                        &
         stencilOperator%periodicOffset)

    call operator_apply_interior_3(stencilOperator, xWithGhostPoints, x)

    n = stencilOperator%boundaryWidth

    ! Left boundary
    if (stencilOperator%hasDomainBoundary(1)) then
       k = stencilOperator%nGhost(1)
       do l = 1, size(x, 2)
          do j = 1, localGridSize(2)
             do i = 1, localGridSize(1)
                do m = 1, stencilOperator%boundaryDepth
                   x(i + localGridSize(1) * (j - 1 +                                         &
                        localGridSize(2) * (k - 1 + m - stencilOperator%nGhost(1))), l) =    &
                        sum(stencilOperator%rhsBoundary1(:,m) *                              &
                        xWithGhostPoints(i,j,k+1:k+n,l))
                end do
             end do
          end do
       end do
    end if

    ! Right boundary
    if (stencilOperator%hasDomainBoundary(2)) then
       k = stencilOperator%nGhost(1) + localGridSize(3) + 1
       do l = 1, size(x, 2)
          do j = 1, localGridSize(2)
             do i = 1, localGridSize(1)
                do m = stencilOperator%boundaryDepth, 1, -1
                   x(i + localGridSize(1) * (j - 1 +                                         &
                        localGridSize(2) * (k - 1 - m - stencilOperator%nGhost(1))), l) =    &
                        sum(stencilOperator%rhsBoundary2(:,m) *                              &
                        xWithGhostPoints(i,j,k-n:k-1,l))
                end do
             end do
          end do
       end do
    end if

    deallocate(xWithGhostPoints)

    return
  end subroutine operator_apply_3


  subroutine operator_apply_interior_1(stencilOperator, xWithGhostPoints, x)

    implicit none

    ! Arguments
    type(t_StencilOperator), intent(in) :: stencilOperator
    real(WP), intent(in) :: xWithGhostPoints(:,:,:,:)
    real(WP), intent(out) :: x(:,:)

    ! Local variables
    integer :: i, j, k, l, n, is, ie

    n = stencilOperator%interiorWidth / 2
    is = 1 + stencilOperator%nGhost(1)
    ie = localGridSize(1) + stencilOperator%nGhost(1)
    if (stencilOperator%nGhost(1) .eq. 0) is = is + stencilOperator%boundaryDepth
    if (stencilOperator%nGhost(2) .eq. 0) ie = ie - stencilOperator%boundaryDepth

    ! Save FLOPS based on symmetry of interior stencil
    select case (stencilOperator%symmetryType)

    case (SKEW_SYMMETRIC)
       do l = 1, size(x, 2)
          do k = 1, localGridSize(3)
             do j = 1, localGridSize(2)
                do i = is, ie
                   x(i - stencilOperator%nGhost(1) + localGridSize(1) *                      &
                        (j - 1 + localGridSize(2) * (k - 1)), l) =                           &
                        sum(stencilOperator%rhsInterior(1:n) *                               &
                        (xWithGhostPoints(i+1:i+n,j,k,l) -                                   &
                        xWithGhostPoints(i-1:i-n:-1,j,k,l)))
                end do
             end do
          end do
       end do

    case (SYMMETRIC)
       do l = 1, size(x, 2)
          do k = 1, localGridSize(3)
             do j = 1, localGridSize(2)
                do i = is, ie
                   x(i - stencilOperator%nGhost(1) + localGridSize(1) *                      &
                        (j - 1 + localGridSize(2) * (k - 1)), l) =                           &
                        sum(stencilOperator%rhsInterior(1:n) *                               &
                        (xWithGhostPoints(i+1:i+n,j,k,l) +                                   &
                        xWithGhostPoints(i-1:i-n:-1,j,k,l))) +                               &
                        stencilOperator%rhsInterior(0) * xWithGhostPoints(i,j,k,l)
                end do
             end do
          end do
       end do

    case default
       do l = 1, size(x, 2)
          do k = 1, localGridSize(3)
             do j = 1, localGridSize(2)
                do i = is, ie
                   x(i - stencilOperator%nGhost(1) + localGridSize(1) *                      &
                        (j - 1 + localGridSize(2) * (k - 1)), l) =                           &
                        sum(stencilOperator%rhsInterior * xWithGhostPoints(i +               &
                        lbound(stencilOperator%rhsInterior,1) : i +                          &
                        ubound(stencilOperator%rhsInterior,1), j, k, l))
                end do
             end do
          end do
       end do

    end select

    return
  end subroutine operator_apply_interior_1


  subroutine operator_apply_interior_2(stencilOperator, xWithGhostPoints, x)

    implicit none

    ! Arguments
    type(t_StencilOperator), intent(in) :: stencilOperator
    real(WP), intent(in) :: xWithGhostPoints(:,:,:,:)
    real(WP), intent(out) :: x(:,:)

    ! Local variables
    integer :: i, j, k, l, n, js, je

    n = stencilOperator%interiorWidth / 2
    js = 1 + stencilOperator%nGhost(1)
    je = localGridSize(2) + stencilOperator%nGhost(1)
    if (stencilOperator%nGhost(1) .eq. 0) js = js + stencilOperator%boundaryDepth
    if (stencilOperator%nGhost(2) .eq. 0) je = je - stencilOperator%boundaryDepth

    ! Save FLOPS based on symmetry of interior stencil
    select case (stencilOperator%symmetryType)

    case (SKEW_SYMMETRIC)
       do l = 1, size(x, 2)
          do k = 1, localGridSize(3)
             do j = js, je
                do i = 1, localGridSize(1)
                   x(i + localGridSize(1) * (j - 1 - stencilOperator%nGhost(1) +             &
                        localGridSize(2) * (k - 1)), l) =                                    &
                        sum(stencilOperator%rhsInterior(1:n) *                               &
                        (xWithGhostPoints(i,j+1:j+n,k,l) -                                   &
                        xWithGhostPoints(i,j-1:j-n:-1,k,l)))
                end do
             end do
          end do
       end do

    case (SYMMETRIC)
       do l = 1, size(x, 2)
          do k = 1, localGridSize(3)
             do j = js, je
                do i = 1, localGridSize(1)
                   x(i + localGridSize(1) * (j - 1 - stencilOperator%nGhost(1) +             &
                        localGridSize(2) * (k - 1)), l) =                                    &
                        sum(stencilOperator%rhsInterior(1:n) *                               &
                        (xWithGhostPoints(i,j+1:j+n,k,l) +                                   &
                        xWithGhostPoints(i,j-1:j-n:-1,k,l))) +                               &
                        stencilOperator%rhsInterior(0) * xWithGhostPoints(i,j,k,l)
                end do
             end do
          end do
       end do

    case default
       do l = 1, size(x, 2)
          do k = 1, localGridSize(3)
             do j = js, je
                do i = 1, localGridSize(1)
                   x(i + localGridSize(1) * (j - 1 - stencilOperator%nGhost(1) +             &
                        localGridSize(2) * (k - 1)), l) =                                    &
                        sum(stencilOperator%rhsInterior * xWithGhostPoints(i, j +            &
                        lbound(stencilOperator%rhsInterior,1) : j +                          &
                        ubound(stencilOperator%rhsInterior,1), k, l))
                end do
             end do
          end do
       end do

    end select

    return
  end subroutine operator_apply_interior_2


  subroutine operator_apply_interior_3(stencilOperator, xWithGhostPoints, x)

    implicit none

    ! Arguments
    type(t_StencilOperator), intent(in) :: stencilOperator
    real(WP), intent(in) :: xWithGhostPoints(:,:,:,:)
    real(WP), intent(out) :: x(:,:)

    ! Local variables
    integer :: i, j, k, l, n, ks, ke

    n = stencilOperator%interiorWidth / 2
    ks = 1 + stencilOperator%nGhost(1)
    ke = localGridSize(3) + stencilOperator%nGhost(1)
    if (stencilOperator%nGhost(1) .eq. 0) ks = ks + stencilOperator%boundaryDepth
    if (stencilOperator%nGhost(2) .eq. 0) ke = ke - stencilOperator%boundaryDepth

    ! Save FLOPS based on symmetry of interior stencil
    select case (stencilOperator%symmetryType)

    case (SKEW_SYMMETRIC)
       do l = 1, size(x, 2)
          do k = ks, ke
             do j = 1, localGridSize(2)
                do i = 1, localGridSize(1)
                   x(i + localGridSize(1) * (j - 1 + localGridSize(2) * (k - 1 -             &
                        stencilOperator%nGhost(1))), l) =                                    &
                        sum(stencilOperator%rhsInterior(1:n) *                               &
                        (xWithGhostPoints(i,j,k+1:k+n,l) -                                   &
                        xWithGhostPoints(i,j,k-1:k-n:-1,l)))
                end do
             end do
          end do
       end do

    case (SYMMETRIC)
       do l = 1, size(x, 2)
          do k = ks, ke
             do j = 1, localGridSize(2)
                do i = 1, localGridSize(1)
                   x(i + localGridSize(1) * (j - 1 + localGridSize(2) * (k - 1 -             &
                        stencilOperator%nGhost(1))), l) =                                    &
                        sum(stencilOperator%rhsInterior(1:n) *                               &
                        (xWithGhostPoints(i,j,k+1:k+n,l) +                                   &
                        xWithGhostPoints(i,j,k-1:k-n:-1,l))) +                               &
                        stencilOperator%rhsInterior(0) * xWithGhostPoints(i,j,k,l)
                end do
             end do
          end do
       end do

    case default
       do l = 1, size(x, 2)
          do k = ks, ke
             do j = 1, localGridSize(2)
                do i = 1, localGridSize(1)
                   x(i + localGridSize(1) * (j - 1 + localGridSize(2) * (k - 1 -             &
                        stencilOperator%nGhost(1))), l) = sum(stencilOperator%rhsInterior *  &
                        xWithGhostPoints(i, j, k + lbound(stencilOperator%rhsInterior, 1) :  &
                        k + ubound(stencilOperator%rhsInterior, 1), l))
                end do
             end do
          end do
       end do

    end select

    return
  end subroutine operator_apply_interior_3


  subroutine operator_hybrid_apply_1(stencilOperator, stencilOperator2, phi, x)

    ! External modules
    use parallel

    implicit none

    ! Arguments
    type(t_StencilOperator), intent(in) :: stencilOperator, stencilOperator2
    real(WP), intent(in) :: phi(:)
    real(WP), intent(inout) :: x(:,:)

    ! Local variables
    integer :: i, j, k, l, m, n, n2, is, ie, gi, lb, ub
    real(WP), allocatable :: xWithGhostPoints(:,:,:,:)

    allocate(xWithGhostPoints(localGridSize(1) + sum(stencilOperator%nGhost),                &
            localGridSize(2), localGridSize(3), size(x, 2)))

    do l = 1, size(x, 2)
       do k = 1, localGridSize(3)
          do j = 1, localGridSize(2)
             do i = 1, localGridSize(1)
                xWithGhostPoints(i + stencilOperator%nGhost(1), j, k, l) =                   &
                     x(i + localGridSize(1) * (j - 1 + localGridSize(2) * (k - 1)), l)
             end do
          end do
       end do
    end do

    call fill_ghost_points(xWithGhostPoints, 1, stencilOperator%nGhost,                      &
         stencilOperator%periodicOffset)

    ! Interior
    n = stencilOperator%interiorWidth / 2
    n2 = stencilOperator2%interiorWidth / 2
    is = 1 + stencilOperator%nGhost(1)
    ie = localGridSize(1) + stencilOperator%nGhost(1)
    if (stencilOperator%nGhost(1) .eq. 0) is = is + stencilOperator%boundaryDepth
    if (stencilOperator%nGhost(2) .eq. 0) ie = ie - stencilOperator%boundaryDepth

    ! Save FLOPS based on symmetry of interior stencil
    ! Assumption : Either stencilOperator%symmetryType = stencilOperator2%symmetryType OR    &
    ! stencilOperator%symmetryType = ASYMMETRIC
    select case (stencilOperator%symmetryType)

    case (SKEW_SYMMETRIC)
       do l = 1, size(x, 2)
          do k = 1, localGridSize(3)
             do j = 1, localGridSize(2)
                do i = is, ie
                   gi = i - stencilOperator%nGhost(1) + localGridSize(1) *                   &
                        (j - 1 + localGridSize(2) * (k - 1))
                   x(gi, l) =                                                                &
                        (1.0_WP - phi(gi)) *                                                 &
                        sum(stencilOperator%rhsInterior(1:n) *                               &
                        (xWithGhostPoints(i+1:i+n,j,k,l) -                                   &
                        xWithGhostPoints(i-1:i-n:-1,j,k,l)))                                 &
                        +                                                                    &
                        phi(gi) *                                                            &
                        sum(stencilOperator2%rhsInterior(1:n2) *                             &
                        (xWithGhostPoints(i+1:i+n2,j,k,l) -                                  &
                        xWithGhostPoints(i-1:i-n2:-1,j,k,l)))
                end do
             end do
          end do
       end do

    case (SYMMETRIC)
       do l = 1, size(x, 2)
          do k = 1, localGridSize(3)
             do j = 1, localGridSize(2)
                do i = is, ie
                   gi = i - stencilOperator%nGhost(1) + localGridSize(1) *                   &
                        (j - 1 + localGridSize(2) * (k - 1))
                   x(gi, l) =                                                                &
                        (1.0_WP - phi(gi)) *                                                 &
                        (sum(stencilOperator%rhsInterior(1:n) *                              &
                        (xWithGhostPoints(i+1:i+n,j,k,l) +                                   &
                        xWithGhostPoints(i-1:i-n:-1,j,k,l))) +                               &
                        stencilOperator%rhsInterior(0) * xWithGhostPoints(i,j,k,l))          &
                        +                                                                    &
                        phi(gi) *                                                            &
                        (sum(stencilOperator2%rhsInterior(1:n2) *                            &
                        (xWithGhostPoints(i+1:i+n2,j,k,l) +                                  &
                        xWithGhostPoints(i-1:i-n2:-1,j,k,l))) +                              &
                        stencilOperator2%rhsInterior(0) * xWithGhostPoints(i,j,k,l))
                end do
             end do
          end do
       end do

    case default
       do l = 1, size(x, 2)
          do k = 1, localGridSize(3)
             do j = 1, localGridSize(2)
                do i = is, ie
                   gi = i - stencilOperator%nGhost(1) + localGridSize(1) *                   &
                        (j - 1 + localGridSize(2) * (k - 1))
                   x(gi, l) =                                                                &
                        (1.0_WP - phi(gi)) *                                                 &
                        sum(stencilOperator%rhsInterior * xWithGhostPoints(                  &
                        i + lbound(stencilOperator%rhsInterior,1) :                          &
                        i + ubound(stencilOperator%rhsInterior,1), j, k, l))                 &
                        +                                                                    &
                        phi(gi) *                                                            &
                        sum(stencilOperator2%rhsInterior * xWithGhostPoints(                 &
                        i + lbound(stencilOperator2%rhsInterior,1) :                         &
                        i + ubound(stencilOperator2%rhsInterior,1), j, k, l))
                end do
             end do
          end do
       end do

    end select

    n = stencilOperator%boundaryWidth
    n2 = stencilOperator2%boundaryWidth
    lb = lbound(stencilOperator2%rhsInterior,1)
    ub = ubound(stencilOperator2%rhsInterior,1)

    ! Left boundary 
    ! Assumption: stencilOperator2%boundaryDepth <= stencilOperator%boundaryDepth
    if (stencilOperator%hasDomainBoundary(1)) then
       i = stencilOperator%nGhost(1)
       do l = 1, size(x, 2)
          do k = 1, localGridSize(3)
             do j = 1, localGridSize(2)
                do m = 1, stencilOperator2%boundaryDepth
                   gi = i + m - stencilOperator%nGhost(1) + localGridSize(1) *               &
                        (j - 1 + localGridSize(2) * (k - 1))
                   x(gi, l) =                                                                &
                        (1.0_WP - phi(gi)) *                                                 &
                        sum(stencilOperator%rhsBoundary1(:,m) *                              &
                        xWithGhostPoints(i+1:i+n,j,k,l))                                     &
                        +                                                                    &
                        phi(gi) *                                                            &
                        sum(stencilOperator2%rhsBoundary1(:,m) *                             &
                        xWithGhostPoints(i+1:i+n2,j,k,l))
                end do
                do m = stencilOperator2%boundaryDepth+1, stencilOperator%boundaryDepth
                   gi = i + m - stencilOperator%nGhost(1) + localGridSize(1) *               &
                        (j - 1 + localGridSize(2) * (k - 1))
                   x(gi, l) =                                                                &
                        (1.0_WP - phi(gi)) *                                                 &
                        sum(stencilOperator%rhsBoundary1(:,m) *                              &
                        xWithGhostPoints(i+1:i+n,j,k,l))                                     &
                        +                                                                    &
                        phi(gi) *                                                            &
                        sum(stencilOperator2%rhsInterior *                                   &
                        xWithGhostPoints(m+lb:m+ub,j,k,l))
                end do
             end do
          end do
       end do
    end if

    ! Right boundary
    ! Assumption: stencilOperator2%boundaryDepth <= stencilOperator%boundaryDepth
    if (stencilOperator%hasDomainBoundary(2)) then
       i = stencilOperator%nGhost(1) + localGridSize(1) + 1
       do l = 1, size(x, 2)
          do k = 1, localGridSize(3)
             do j = 1, localGridSize(2)
                do m = stencilOperator%boundaryDepth, stencilOperator2%boundaryDepth+1, -1
                   gi = i - m - stencilOperator%nGhost(1) + localGridSize(1) *               &
                        (j - 1 + localGridSize(2) * (k - 1))
                   x(gi, l) =                                                                &
                        (1.0_WP - phi(gi)) *                                                 &
                        sum(stencilOperator%rhsBoundary2(:,m) *                              &
                        xWithGhostPoints(i-n:i-1,j,k,l))                                     &
                        +                                                                    &
                        phi(gi) *                                                            &
                        sum(stencilOperator2%rhsInterior *                                   &
                        xWithGhostPoints(i-m+lb:i-m+ub,j,k,l))
                end do
                do m = stencilOperator2%boundaryDepth, 1, -1
                   gi = i - m - stencilOperator%nGhost(1) + localGridSize(1) *               &
                        (j - 1 + localGridSize(2) * (k - 1))
                   x(gi, l) =                                                                &
                        (1.0_WP - phi(gi)) *                                                 &
                        sum(stencilOperator%rhsBoundary2(:,m) *                              &
                        xWithGhostPoints(i-n:i-1,j,k,l))                                     &
                        +                                                                    &
                        phi(gi) *                                                            &
                        sum(stencilOperator2%rhsBoundary2(:,m) *                             &
                        xWithGhostPoints(i-n2:i-1,j,k,l))
                end do
             end do
          end do
       end do
    end if

    deallocate(xWithGhostPoints)

    return
  end subroutine operator_hybrid_apply_1


  subroutine operator_hybrid_apply_2(stencilOperator, stencilOperator2, phi, x)

    ! External modules
    use parallel

    implicit none

    ! Arguments
    type(t_StencilOperator), intent(in) :: stencilOperator, stencilOperator2
    real(WP), intent(in) :: phi(:)
    real(WP), intent(inout) :: x(:,:)

    ! Local variables
    integer :: i, j, k, l, m, n, n2, js, je, gi, lb, ub
    real(WP), allocatable :: xWithGhostPoints(:,:,:,:)

    allocate(xWithGhostPoints(localGridSize(1), localGridSize(2) +                           &
         sum(stencilOperator%nGhost), localGridSize(3), size(x, 2)))

    do l = 1, size(x, 2)
       do k = 1, localGridSize(3)
          do j = 1, localGridSize(2)
             do i = 1, localGridSize(1)
                xWithGhostPoints(i, j + stencilOperator%nGhost(1), k, l) =                   &
                     x(i + localGridSize(1) * (j - 1 + localGridSize(2) * (k - 1)), l)
             end do
          end do
       end do
    end do

    call fill_ghost_points(xWithGhostPoints, 2, stencilOperator%nGhost,                      &
         stencilOperator%periodicOffset)

    ! Interior
    n = stencilOperator%interiorWidth / 2
    n2 = stencilOperator2%interiorWidth / 2
    js = 1 + stencilOperator%nGhost(1)
    je = localGridSize(2) + stencilOperator%nGhost(1)
    if (stencilOperator%nGhost(1) .eq. 0) js = js + stencilOperator%boundaryDepth
    if (stencilOperator%nGhost(2) .eq. 0) je = je - stencilOperator%boundaryDepth

    ! Save FLOPS based on symmetry of interior stencil
    ! Assumption : Either stencilOperator%symmetryType = stencilOperator2%symmetryType OR    &
    ! stencilOperator%symmetryType = ASYMMETRIC
    select case (stencilOperator%symmetryType)

    case (SKEW_SYMMETRIC)
       do l = 1, size(x, 2)
          do k = 1, localGridSize(3)
             do j = js, je
                do i = 1, localGridSize(1)
                   gi = i + localGridSize(1) * (j - 1 - stencilOperator%nGhost(1) +          &
                        localGridSize(2) * (k - 1))
                   x(gi, l) =                                                                &
                        (1.0_WP - phi(gi)) *                                                 &
                        sum(stencilOperator%rhsInterior(1:n) *                               &
                        (xWithGhostPoints(i,j+1:j+n,k,l) -                                   &
                        xWithGhostPoints(i,j-1:j-n:-1,k,l)))                                 &
                        +                                                                    &
                        phi(gi) *                                                            &
                        sum(stencilOperator2%rhsInterior(1:n2) *                             &
                        (xWithGhostPoints(i,j+1:j+n2,k,l) -                                  &
                        xWithGhostPoints(i,j-1:j-n2:-1,k,l)))
                end do
             end do
          end do
       end do

    case (SYMMETRIC)
       do l = 1, size(x, 2)
          do k = 1, localGridSize(3)
             do j = js, je
                do i = 1, localGridSize(1)
                   gi = i + localGridSize(1) * (j - 1 - stencilOperator%nGhost(1) +          &
                        localGridSize(2) * (k - 1))
                   x(gi, l) =                                                                &
                        (1.0_WP - phi(gi)) *                                                 &
                        (sum(stencilOperator%rhsInterior(1:n) *                              &
                        (xWithGhostPoints(i,j+1:j+n,k,l) +                                   &
                        xWithGhostPoints(i,j-1:j-n:-1,k,l))) +                               &
                        stencilOperator%rhsInterior(0) * xWithGhostPoints(i,j,k,l))          &
                        +                                                                    &
                        phi(gi) *                                                            &
                        (sum(stencilOperator2%rhsInterior(1:n2) *                            &
                        (xWithGhostPoints(i,j+1:j+n2,k,l) +                                  &
                        xWithGhostPoints(i,j-1:j-n2:-1,k,l))) +                              &
                        stencilOperator2%rhsInterior(0) * xWithGhostPoints(i,j,k,l))
                end do
             end do
          end do
       end do

    case default
       do l = 1, size(x, 2)
          do k = 1, localGridSize(3)
             do j = js, je
                do i = 1, localGridSize(1)
                   gi = i + localGridSize(1) * (j - 1 - stencilOperator%nGhost(1) +          &
                        localGridSize(2) * (k - 1))
                   x(gi, l) =                                                                &
                        (1.0_WP - phi(gi)) *                                                 &
                        sum(stencilOperator%rhsInterior * xWithGhostPoints(i,                &
                        j + lbound(stencilOperator%rhsInterior,1) :                          &
                        j + ubound(stencilOperator%rhsInterior,1), k, l))                    &
                        +                                                                    &
                        phi(gi) *                                                            &
                        sum(stencilOperator2%rhsInterior * xWithGhostPoints(i,               &
                        j + lbound(stencilOperator2%rhsInterior,1) :                         &
                        j + ubound(stencilOperator2%rhsInterior,1), k, l))
                end do
             end do
          end do
       end do

    end select

    n = stencilOperator%boundaryWidth
    n2 = stencilOperator2%boundaryWidth
    lb = lbound(stencilOperator2%rhsInterior,1)
    ub = ubound(stencilOperator2%rhsInterior,1)

    ! Left boundary 
    ! Assumption: stencilOperator2%boundaryDepth <= stencilOperator%boundaryDepth
    if (stencilOperator%hasDomainBoundary(1)) then
       j = stencilOperator%nGhost(1)
       do l = 1, size(x, 2)
          do k = 1, localGridSize(3)
             do i = 1, localGridSize(1)
                do m = 1, stencilOperator2%boundaryDepth
                   gi = i + localGridSize(1) * (j - 1 + m - stencilOperator%nGhost(1) +      &
                        localGridSize(2) * (k - 1))
                   x(gi, l) =                                                                &
                        (1.0_WP - phi(gi)) *                                                 &
                        sum(stencilOperator%rhsBoundary1(:,m) *                              &
                        xWithGhostPoints(i,j+1:j+n,k,l))                                     &
                        +                                                                    &
                        phi(gi) *                                                            &
                        sum(stencilOperator2%rhsBoundary1(:,m) *                             &
                        xWithGhostPoints(i,j+1:j+n2,k,l))
                end do
                do m = stencilOperator2%boundaryDepth+1, stencilOperator%boundaryDepth
                   gi = i + localGridSize(1) * (j - 1 + m - stencilOperator%nGhost(1) +      &
                        localGridSize(2) * (k - 1))
                   x(gi, l) =                                                                &
                        (1.0_WP - phi(gi)) *                                                 &
                        sum(stencilOperator%rhsBoundary1(:,m) *                              &
                        xWithGhostPoints(i,j+1:j+n,k,l))                                     &
                        +                                                                    &
                        phi(gi) *                                                            &
                        sum(stencilOperator2%rhsInterior *                                   &
                        xWithGhostPoints(i,m+lb:m+ub,k,l))
                end do
             end do
          end do
       end do
    end if

    ! Right boundary
    ! Assumption: stencilOperator2%boundaryDepth <= stencilOperator%boundaryDepth
    if (stencilOperator%hasDomainBoundary(2)) then
       j = stencilOperator%nGhost(1) + localGridSize(2) + 1
       do l = 1, size(x, 2)
          do k = 1, localGridSize(3)
             do i = 1, localGridSize(1)
                do m = stencilOperator%boundaryDepth, stencilOperator2%boundaryDepth+1, -1
                   gi = i + localGridSize(1) * (j - 1 - m - stencilOperator%nGhost(1) +      &
                        localGridSize(2) * (k - 1))
                   x(gi, l) =                                                                &
                        (1.0_WP - phi(gi)) *                                                 &
                        sum(stencilOperator%rhsBoundary2(:,m) *                              &
                        xWithGhostPoints(i,j-n:j-1,k,l))                                     &
                        +                                                                    &
                        phi(gi) *                                                            &
                        sum(stencilOperator2%rhsInterior *                                   &
                        xWithGhostPoints(i,j-m+lb:j-m+ub,k,l))
                end do
                do m = stencilOperator2%boundaryDepth, 1, -1
                   gi = i + localGridSize(1) * (j - 1 - m - stencilOperator%nGhost(1) +      &
                        localGridSize(2) * (k - 1))
                   x(gi, l) =                                                                &
                        (1.0_WP - phi(gi)) *                                                 &
                        sum(stencilOperator%rhsBoundary2(:,m) *                              &
                        xWithGhostPoints(i,j-n:j-1,k,l))                                     &
                        +                                                                    &
                        phi(gi) *                                                            &
                        sum(stencilOperator2%rhsBoundary2(:,m) *                             &
                        xWithGhostPoints(i,j-n2:j-1,k,l))
                end do
             end do
          end do
       end do
    end if

    deallocate(xWithGhostPoints)

    return
  end subroutine operator_hybrid_apply_2


  subroutine operator_hybrid_apply_3(stencilOperator, stencilOperator2, phi, x)

    ! External modules
    use parallel

    implicit none

    ! Arguments
    type(t_StencilOperator), intent(in) :: stencilOperator, stencilOperator2
    real(WP), intent(in) :: phi(:)
    real(WP), intent(inout) :: x(:,:)

    ! Local variables
    integer :: i, j, k, l, m, n, n2, ks, ke, gi, lb, ub
    real(WP), allocatable :: xWithGhostPoints(:,:,:,:)

    allocate(xWithGhostPoints(localGridSize(1), localGridSize(2),                            &
         localGridSize(3) + sum(stencilOperator%nGhost), size(x, 2)))

    do l = 1, size(x, 2)
       do k = 1, localGridSize(3)
          do j = 1, localGridSize(2)
             do i = 1, localGridSize(1)
                xWithGhostPoints(i, j, k + stencilOperator%nGhost(1), l) =                   &
                     x(i + localGridSize(1) * (j - 1 + localGridSize(2) * (k - 1)), l)
             end do
          end do
       end do
    end do

    call fill_ghost_points(xWithGhostPoints, 3, stencilOperator%nGhost,                      &
         stencilOperator%periodicOffset)

    ! Interior
    n = stencilOperator%interiorWidth / 2
    n2 = stencilOperator2%interiorWidth / 2
    ks = 1 + stencilOperator%nGhost(1)
    ke = localGridSize(3) + stencilOperator%nGhost(1)
    if (stencilOperator%nGhost(1) .eq. 0) ks = ks + stencilOperator%boundaryDepth
    if (stencilOperator%nGhost(2) .eq. 0) ke = ke - stencilOperator%boundaryDepth

    ! Save FLOPS based on symmetry of interior stencil
    ! Assumption : Either stencilOperator%symmetryType = stencilOperator2%symmetryType OR    &
    ! stencilOperator%symmetryType = ASYMMETRIC
    select case (stencilOperator%symmetryType)

    case (SKEW_SYMMETRIC)
       do l = 1, size(x, 2)
          do k = ks, ke
             do j = 1, localGridSize(2)
                do i = 1, localGridSize(1)
                   gi = i + localGridSize(1) * (j - 1 + localGridSize(2) *                   &
                        (k - 1 - stencilOperator%nGhost(1)))
                   x(gi, l) =                                                                &
                        (1.0_WP - phi(gi)) *                                                 &
                        sum(stencilOperator%rhsInterior(1:n) *                               &
                        (xWithGhostPoints(i,j,k+1:k+n,l) -                                   &
                        xWithGhostPoints(i,j,k-1:k-n:-1,l)))                                 &
                        +                                                                    &
                        phi(gi) *                                                            &
                        sum(stencilOperator2%rhsInterior(1:n2) *                             &
                        (xWithGhostPoints(i,j,k+1:k+n2,l) -                                  &
                        xWithGhostPoints(i,j,k-1:k-n2:-1,l)))
                end do
             end do
          end do
       end do

    case (SYMMETRIC)
       do l = 1, size(x, 2)
          do k = ks, ke
             do j = 1, localGridSize(2)
                do i = 1, localGridSize(1)
                   gi = i + localGridSize(1) * (j - 1 + localGridSize(2) *                   &
                        (k - 1 - stencilOperator%nGhost(1)))
                   x(gi, l) =                                                                &
                        (1.0_WP - phi(gi)) *                                                 &
                        (sum(stencilOperator%rhsInterior(1:n) *                              &
                        (xWithGhostPoints(i,j,k+1:k+n,l) +                                   &
                        xWithGhostPoints(i,j,k-1:k-n:-1,l))) +                               &
                        stencilOperator%rhsInterior(0) * xWithGhostPoints(i,j,k,l))          &
                        +                                                                    &
                        phi(gi) *                                                            &
                        (sum(stencilOperator2%rhsInterior(1:n2) *                            &
                        (xWithGhostPoints(i,j,k+1:k+n2,l) +                                  &
                        xWithGhostPoints(i,j,k-1:k-n2:-1,l))) +                              &
                        stencilOperator2%rhsInterior(0) * xWithGhostPoints(i,j,k,l))
                end do
             end do
          end do
       end do

    case default
       do l = 1, size(x, 2)
          do k = ks, ke
             do j = 1, localGridSize(2)
                do i = 1, localGridSize(1)
                   gi = i + localGridSize(1) * (j - 1 + localGridSize(2) *                   &
                        (k - 1 - stencilOperator%nGhost(1)))
                   x(gi, l) =                                                                &
                        (1.0_WP - phi(gi)) *                                                 &
                        sum(stencilOperator%rhsInterior * xWithGhostPoints(i, j,             &
                        k + lbound(stencilOperator%rhsInterior, 1) :                         &
                        k + ubound(stencilOperator%rhsInterior, 1), l))                      &
                        +                                                                    &
                        phi(gi) *                                                            &
                        sum(stencilOperator2%rhsInterior * xWithGhostPoints(i, j,            &
                        k + lbound(stencilOperator2%rhsInterior, 1) :                        &
                        k + ubound(stencilOperator2%rhsInterior, 1), l))
                end do
             end do
          end do
       end do

    end select

    n = stencilOperator%boundaryWidth
    n2 = stencilOperator2%boundaryWidth
    lb = lbound(stencilOperator2%rhsInterior,1)
    ub = ubound(stencilOperator2%rhsInterior,1)

    ! Left boundary 
    ! Assumption: stencilOperator2%boundaryDepth <= stencilOperator%boundaryDepth
    if (stencilOperator%hasDomainBoundary(1)) then
       k = stencilOperator%nGhost(1)
       do l = 1, size(x, 2)
          do j = 1, localGridSize(2)
             do i = 1, localGridSize(1)
                do m = 1, stencilOperator2%boundaryDepth
                   gi = i + localGridSize(1) * (j - 1 + localGridSize(2) *                   &
                        (k - 1 + m - stencilOperator%nGhost(1)))
                   x(gi, l) =                                                                &
                        (1.0_WP - phi(gi)) *                                                 &
                        sum(stencilOperator%rhsBoundary1(:,m) *                              &
                        xWithGhostPoints(i,j,k+1:k+n,l))                                     &
                        +                                                                    &
                        phi(gi) *                                                            &
                        sum(stencilOperator2%rhsBoundary1(:,m) *                             &
                        xWithGhostPoints(i,j,k+1:k+n2,l))
                end do
                do m = stencilOperator2%boundaryDepth+1, stencilOperator%boundaryDepth
                   gi = i + localGridSize(1) * (j - 1 + localGridSize(2) *                   &
                        (k - 1 + m - stencilOperator%nGhost(1)))
                   x(gi, l) =                                                                &
                        (1.0_WP - phi(gi)) *                                                 &
                        sum(stencilOperator%rhsBoundary1(:,m) *                              &
                        xWithGhostPoints(i,j,k+1:k+n,l))                                     &
                        +                                                                    &
                        phi(gi) *                                                            &
                        sum(stencilOperator2%rhsInterior *                                   &
                        xWithGhostPoints(i,j,m+lb:m+ub,l))
                end do
             end do
          end do
       end do
    end if

    ! Right boundary
    ! Assumption: stencilOperator2%boundaryDepth <= stencilOperator%boundaryDepth
    if (stencilOperator%hasDomainBoundary(2)) then
       k = stencilOperator%nGhost(1) + localGridSize(3) + 1
       do l = 1, size(x, 2)
          do j = 1, localGridSize(2)
             do i = 1, localGridSize(1)
                do m = stencilOperator%boundaryDepth, stencilOperator2%boundaryDepth+1, -1
                   gi = i + localGridSize(1) * (j - 1 + localGridSize(2) *                   &
                        (k - 1 - m - stencilOperator%nGhost(1)))
                   x(gi, l) =                                                                &
                        (1.0_WP - phi(gi)) *                                                 &
                        sum(stencilOperator%rhsBoundary2(:,m) *                              &
                        xWithGhostPoints(i,j,k-n:k-1,l))                                     &
                        +                                                                    &
                        phi(gi) *                                                            &
                        sum(stencilOperator2%rhsInterior *                                   &
                        xWithGhostPoints(i,j,k-m+lb:k-m+ub,l))
                end do
                do m = stencilOperator2%boundaryDepth, 1, -1
                   gi = i + localGridSize(1) * (j - 1 + localGridSize(2) *                   &
                        (k - 1 - m - stencilOperator%nGhost(1)))
                   x(gi, l) =                                                                &
                        (1.0_WP - phi(gi)) *                                                 &
                        sum(stencilOperator%rhsBoundary2(:,m) *                              &
                        xWithGhostPoints(i,j,k-n:k-1,l))                                     &
                        +                                                                    &
                        phi(gi) *                                                            &
                        sum(stencilOperator2%rhsBoundary2(:,m) *                             &
                        xWithGhostPoints(i,j,k-n2:k-1,l))
                end do
             end do
          end do
       end do
    end if

    deallocate(xWithGhostPoints)

    return
  end subroutine operator_hybrid_apply_3


  subroutine operator_adjoint_hybrid_apply_1(stencilOperator, stencilOperator2, phi, x)

    ! External modules
    use parallel

    implicit none

    ! Arguments
    type(t_StencilOperator), intent(in) :: stencilOperator, stencilOperator2
    real(WP), intent(in) :: phi(:)
    real(WP), intent(inout) :: x(:,:)

    ! Local variables
    integer :: i, j, k, l, m, n, n2, p, is, ie
    real(WP), allocatable :: xWithGhostPoints(:,:,:,:)

    allocate(xWithGhostPoints(localGridSize(1) + sum(stencilOperator%nGhost),                &
            localGridSize(2), localGridSize(3), size(x, 2) + 1)) !... include phi

    p = size(xWithGhostPoints,4) !... index of phi

    do k = 1, localGridSize(3)
       do j = 1, localGridSize(2)
          do i = 1, localGridSize(1)
             do l = 1, size(x, 2)
                xWithGhostPoints(i + stencilOperator%nGhost(1), j, k, l) =                   &
                     x(i + localGridSize(1) * (j - 1 + localGridSize(2) * (k - 1)), l)
             end do
             ! x = stenciOpertor * (1-phi) * x + stenciOpertor2 * phi * x
             ! Communicare (1 - phi) since stenciOpertor has a wider stencil size than 
             ! sencilOperator2 --> x = stenciOpertor * phi * x + stenciOpertor2 * (1-phi) * x
             xWithGhostPoints(i + stencilOperator%nGhost(1), j, k, p) =                      &
                  1.0_WP - phi(i + localGridSize(1) * (j - 1 + localGridSize(2) * (k - 1)))
          end do
       end do
    end do

    call fill_ghost_points(xWithGhostPoints, 1, stencilOperator%nGhost,                      &
         stencilOperator%periodicOffset)

    ! Interior
    n = stencilOperator%interiorWidth / 2
    n2 = stencilOperator2%interiorWidth / 2
    is = 1 + stencilOperator%nGhost(1)
    ie = localGridSize(1) + stencilOperator%nGhost(1)
    if (stencilOperator%nGhost(1) .eq. 0) is = is + stencilOperator%boundaryDepth
    if (stencilOperator%nGhost(2) .eq. 0) ie = ie - stencilOperator%boundaryDepth

    ! Save FLOPS based on symmetry of interior stencil
    ! Assumption : Either stencilOperator%symmetryType = stencilOperator2%symmetryType OR    &
    ! stencilOperator%symmetryType = ASYMMETRIC
    select case (stencilOperator%symmetryType)

    case (SKEW_SYMMETRIC)
       do l = 1, size(x, 2)
          do k = 1, localGridSize(3)
             do j = 1, localGridSize(2)
                do i = is, ie
                   x(i - stencilOperator%nGhost(1) + localGridSize(1) *                      &
                        (j - 1 + localGridSize(2) * (k - 1)), l) =                           &
                        sum(stencilOperator%rhsInterior(1:n) *                               &
                        (xWithGhostPoints(i+1:i+n,j,k,l) *                                   &
                        xWithGhostPoints(i+1:i+n,j,k,p) -                                    &
                        xWithGhostPoints(i-1:i-n:-1,j,k,l) *                                 &
                        xWithGhostPoints(i-1:i-n:-1,j,k,p)))                                 &
                        +                                                                    &
                        sum(stencilOperator2%rhsInterior(1:n2) *                             &
                        (xWithGhostPoints(i+1:i+n2,j,k,l) *                                  &
                        (1.0_WP - xWithGhostPoints(i+1:i+n2,j,k,p)) -                        &
                        xWithGhostPoints(i-1:i-n2:-1,j,k,l) *                                &
                        (1.0_WP - xWithGhostPoints(i-1:i-n2:-1,j,k,p))))
                end do
             end do
          end do
       end do

    case (SYMMETRIC)
       do l = 1, size(x, 2)
          do k = 1, localGridSize(3)
             do j = 1, localGridSize(2)
                do i = is, ie
                   x(i - stencilOperator%nGhost(1) + localGridSize(1) *                      &
                        (j - 1 + localGridSize(2) * (k - 1)), l) =                           &
                        sum(stencilOperator%rhsInterior(1:n) *                               &
                        (xWithGhostPoints(i+1:i+n,j,k,l) *                                   &
                        xWithGhostPoints(i+1:i+n,j,k,p) +                                    &
                        xWithGhostPoints(i-1:i-n:-1,j,k,l) *                                 &
                        xWithGhostPoints(i-1:i-n:-1,j,k,p))) +                               &
                        stencilOperator%rhsInterior(0) * xWithGhostPoints(i,j,k,l) *         &
                        xWithGhostPoints(i,j,k,p)                                            &
                        +                                                                    &
                        sum(stencilOperator2%rhsInterior(1:n2) *                             &
                        (xWithGhostPoints(i+1:i+n2,j,k,l) *                                  &
                        (1.0_WP - xWithGhostPoints(i+1:i+n2,j,k,p)) +                        &
                        xWithGhostPoints(i-1:i-n2:-1,j,k,l) *                                &
                        (1.0_WP - xWithGhostPoints(i-1:i-n2:-1,j,k,p)))) +                   &
                        stencilOperator2%rhsInterior(0) * xWithGhostPoints(i,j,k,l) *        &
                        (1.0_WP - xWithGhostPoints(i,j,k,p))
                end do
             end do
          end do
       end do

    case default
       do l = 1, size(x, 2)
          do k = 1, localGridSize(3)
             do j = 1, localGridSize(2)
                do i = is, ie
                   x(i - stencilOperator%nGhost(1) + localGridSize(1) *                      &
                        (j - 1 + localGridSize(2) * (k - 1)), l) =                           &
                        sum(stencilOperator%rhsInterior * xWithGhostPoints(                  &
                        i + lbound(stencilOperator%rhsInterior,1) :                          &
                        i + ubound(stencilOperator%rhsInterior,1), j, k, l) *                &
                        xWithGhostPoints(                                                    &
                        i + lbound(stencilOperator%rhsInterior,1) :                          &
                        i + ubound(stencilOperator%rhsInterior,1), j, k, p))                 &
                        +                                                                    &
                        sum(stencilOperator2%rhsInterior * xWithGhostPoints(                 &
                        i + lbound(stencilOperator2%rhsInterior,1) :                         &
                        i + ubound(stencilOperator2%rhsInterior,1), j, k, l) *               &
                        (1.0_WP - xWithGhostPoints(                                          &
                        i + lbound(stencilOperator2%rhsInterior,1) :                         &
                        i + ubound(stencilOperator2%rhsInterior,1), j, k, p)))
                end do
             end do
          end do
       end do

    end select

    n = stencilOperator%boundaryWidth
    n2 = stencilOperator2%boundaryWidth

    ! Left boundary
    if (stencilOperator%hasDomainBoundary(1)) then
       i = stencilOperator%nGhost(1)
       do l = 1, size(x, 2)
          do k = 1, localGridSize(3)
             do j = 1, localGridSize(2)
                do m = 1, stencilOperator%boundaryDepth
                   x(i + m - stencilOperator%nGhost(1) +                                     &
                        localGridSize(1) * (j - 1 + localGridSize(2) * (k - 1)), l) =        &
                        sum(stencilOperator%rhsBoundary1(:,m) *                              &
                        xWithGhostPoints(i+1:i+n,j,k,l) *                                    &
                        xWithGhostPoints(i+1:i+n,j,k,p))                                     &
                        +                                                                    &
                        sum(stencilOperator2%rhsBoundary1(:,m) *                             &
                        xWithGhostPoints(i+1:i+n2,j,k,l) *                                   &
                        (1.0_WP - xWithGhostPoints(i+1:i+n2,j,k,p)))
                end do
             end do
          end do
       end do
    end if

    ! Right boundary
    if (stencilOperator%hasDomainBoundary(2)) then
       i = stencilOperator%nGhost(1) + localGridSize(1) + 1
       do l = 1, size(x, 2)
          do k = 1, localGridSize(3)
             do j = 1, localGridSize(2)
                do m = stencilOperator%boundaryDepth, 1, -1
                   x(i - m - stencilOperator%nGhost(1) +                                     &
                        localGridSize(1) * (j - 1 + localGridSize(2) * (k - 1)), l) =        &
                        sum(stencilOperator%rhsBoundary2(:,m) *                              &
                        xWithGhostPoints(i-n:i-1,j,k,l) *                                    &
                        xWithGhostPoints(i-n:i-1,j,k,p))                                     &
                        +                                                                    &
                        sum(stencilOperator2%rhsBoundary2(:,m) *                             &
                        xWithGhostPoints(i-n2:i-1,j,k,l) *                                   &
                        (1.0_WP - xWithGhostPoints(i-n2:i-1,j,k,p)))
                end do
             end do
          end do
       end do
    end if

    deallocate(xWithGhostPoints)

    return
  end subroutine operator_adjoint_hybrid_apply_1


  subroutine operator_adjoint_hybrid_apply_2(stencilOperator, stencilOperator2, phi, x)

    ! External modules
    use parallel

    implicit none

    ! Arguments
    type(t_StencilOperator), intent(in) :: stencilOperator, stencilOperator2
    real(WP), intent(in) :: phi(:)
    real(WP), intent(inout) :: x(:,:)

    ! Local variables
    integer :: i, j, k, l, m, n, n2, p, js, je
    real(WP), allocatable :: xWithGhostPoints(:,:,:,:)

    allocate(xWithGhostPoints(localGridSize(1), localGridSize(2) +                           &
         sum(stencilOperator%nGhost), localGridSize(3), size(x, 2) + 1)) !... include phi

    p = size(xWithGhostPoints,4) !... index of phi

    do k = 1, localGridSize(3)
       do j = 1, localGridSize(2)
          do i = 1, localGridSize(1)
             do l = 1, size(x, 2)
                xWithGhostPoints(i, j + stencilOperator%nGhost(1), k, l) =                   &
                     x(i + localGridSize(1) * (j - 1 + localGridSize(2) * (k - 1)), l)
             end do
             ! x = stenciOpertor * (1-phi) * x + stenciOpertor2 * phi * x
             ! Communicare (1 - phi) since stenciOpertor has a wider stencil size than 
             ! sencilOperator2 --> x = stenciOpertor * phi * x + stenciOpertor2 * (1-phi) * x
             xWithGhostPoints(i, j + stencilOperator%nGhost(1), k, p) =                      &
                  1.0_WP - phi(i + localGridSize(1) * (j - 1 + localGridSize(2) * (k - 1)))
          end do
       end do
    end do

    call fill_ghost_points(xWithGhostPoints, 2, stencilOperator%nGhost,                      &
         stencilOperator%periodicOffset)

    ! Interior
    n = stencilOperator%interiorWidth / 2
    n2 = stencilOperator2%interiorWidth / 2
    js = 1 + stencilOperator%nGhost(1)
    je = localGridSize(2) + stencilOperator%nGhost(1)
    if (stencilOperator%nGhost(1) .eq. 0) js = js + stencilOperator%boundaryDepth
    if (stencilOperator%nGhost(2) .eq. 0) je = je - stencilOperator%boundaryDepth

    ! Save FLOPS based on symmetry of interior stencil
    ! Assumption : Either stencilOperator%symmetryType = stencilOperator2%symmetryType OR    &
    ! stencilOperator%symmetryType = ASYMMETRIC
    select case (stencilOperator%symmetryType)

    case (SKEW_SYMMETRIC)
       do l = 1, size(x, 2)
          do k = 1, localGridSize(3)
             do j = js, je
                do i = 1, localGridSize(1)
                   x(i + localGridSize(1) * (j - 1 - stencilOperator%nGhost(1) +             &
                        localGridSize(2) * (k - 1)), l) =                                    &
                        sum(stencilOperator%rhsInterior(1:n) *                               &
                        (xWithGhostPoints(i,j+1:j+n,k,l) *                                   &
                        xWithGhostPoints(i,j+1:j+n,k,p) -                                    &
                        xWithGhostPoints(i,j-1:j-n:-1,k,l) *                                 &
                        xWithGhostPoints(i,j-1:j-n:-1,k,p)))                                 &
                        +                                                                    &
                        sum(stencilOperator2%rhsInterior(1:n2) *                             &
                        (xWithGhostPoints(i,j+1:j+n2,k,l) *                                  &
                        (1.0_WP - xWithGhostPoints(i,j+1:j+n2,k,p)) -                        &
                        xWithGhostPoints(i,j-1:j-n2:-1,k,l)  *                               &
                        (1.0_WP - xWithGhostPoints(i,j-1:j-n2:-1,k,p))))
                end do
             end do
          end do
       end do

    case (SYMMETRIC)
       do l = 1, size(x, 2)
          do k = 1, localGridSize(3)
             do j = js, je
                do i = 1, localGridSize(1)
                   x(i + localGridSize(1) * (j - 1 - stencilOperator%nGhost(1) +             &
                        localGridSize(2) * (k - 1)), l) =                                    &
                        sum(stencilOperator%rhsInterior(1:n) *                               &
                        (xWithGhostPoints(i,j+1:j+n,k,l) *                                   &
                        xWithGhostPoints(i,j+1:j+n,k,p) +                                    &
                        xWithGhostPoints(i,j-1:j-n:-1,k,l) *                                 &
                        xWithGhostPoints(i,j-1:j-n:-1,k,p))) +                               &
                        stencilOperator%rhsInterior(0) * xWithGhostPoints(i,j,k,l) *         &
                        xWithGhostPoints(i,j,k,p)                                            &
                        +                                                                    &
                        sum(stencilOperator2%rhsInterior(1:n2) *                             &
                        (xWithGhostPoints(i,j+1:j+n2,k,l) *                                  &
                        (1.0_WP - xWithGhostPoints(i,j+1:j+n2,k,p)) +                        &
                        xWithGhostPoints(i,j-1:j-n2:-1,k,l) *                                &
                        (1.0_WP - xWithGhostPoints(i,j-1:j-n2:-1,k,p)))) +                   &
                        stencilOperator2%rhsInterior(0) * xWithGhostPoints(i,j,k,l) *        &
                        (1.0_WP - xWithGhostPoints(i,j,k,p))
                end do
             end do
          end do
       end do

    case default
       do l = 1, size(x, 2)
          do k = 1, localGridSize(3)
             do j = js, je
                do i = 1, localGridSize(1)
                   x(i + localGridSize(1) * (j - 1 - stencilOperator%nGhost(1) +             &
                        localGridSize(2) * (k - 1)), l) =                                    &
                        sum(stencilOperator%rhsInterior * xWithGhostPoints(i,                &
                        j + lbound(stencilOperator%rhsInterior,1) :                          &
                        j + ubound(stencilOperator%rhsInterior,1), k, l) *                   &
                        xWithGhostPoints(i,                                                  &
                        j + lbound(stencilOperator%rhsInterior,1) :                          &
                        j + ubound(stencilOperator%rhsInterior,1), k, p))                    &
                        +                                                                    &
                        sum(stencilOperator2%rhsInterior * xWithGhostPoints(i,               &
                        j + lbound(stencilOperator2%rhsInterior,1) :                         &
                        j + ubound(stencilOperator2%rhsInterior,1), k, l) *                  &
                        (1.0_WP - xWithGhostPoints(i,                                        &
                        j + lbound(stencilOperator2%rhsInterior,1) :                         &
                        j + ubound(stencilOperator2%rhsInterior,1), k, p)))
                end do
             end do
          end do
       end do

    end select

    n = stencilOperator%boundaryWidth
    n2 = stencilOperator2%boundaryWidth

    ! Left boundary
    if (stencilOperator%hasDomainBoundary(1)) then
       j = stencilOperator%nGhost(1)
       do l = 1, size(x, 2)
          do k = 1, localGridSize(3)
             do i = 1, localGridSize(1)
                do m = 1, stencilOperator%boundaryDepth
                   x(i + localGridSize(1) * (j - 1 + m - stencilOperator%nGhost(1) +         &
                        localGridSize(2) * (k - 1)), l) =                                    &
                        sum(stencilOperator%rhsBoundary1(:,m) *                              &
                        xWithGhostPoints(i,j+1:j+n,k,l) *                                    &
                        xWithGhostPoints(i,j+1:j+n,k,p))                                     &
                        +                                                                    &
                        sum(stencilOperator2%rhsBoundary1(:,m) *                             &
                        xWithGhostPoints(i,j+1:j+n2,k,l) *                                   &
                        (1.0_WP - xWithGhostPoints(i,j+1:j+n2,k,p)))
                end do
             end do
          end do
       end do
    end if

    ! Right boundary
    if (stencilOperator%hasDomainBoundary(2)) then
       j = stencilOperator%nGhost(1) + localGridSize(2) + 1
       do l = 1, size(x, 2)
          do k = 1, localGridSize(3)
             do i = 1, localGridSize(1)
                do m = stencilOperator%boundaryDepth, 1, -1
                   x(i + localGridSize(1) * (j - 1 - m - stencilOperator%nGhost(1) +         &
                        localGridSize(2) * (k - 1)), l) =                                    &
                        sum(stencilOperator%rhsBoundary2(:,m) *                              &
                        xWithGhostPoints(i,j-n:j-1,k,l) *                                    &
                        xWithGhostPoints(i,j-n:j-1,k,p))                                     &
                        +                                                                    &
                        sum(stencilOperator2%rhsBoundary2(:,m) *                             &
                        xWithGhostPoints(i,j-n2:j-1,k,l) *                                   &
                        (1.0_WP - xWithGhostPoints(i,j-n2:j-1,k,p)))
                end do
             end do
          end do
       end do
    end if

    deallocate(xWithGhostPoints)

    return
  end subroutine operator_adjoint_hybrid_apply_2


  subroutine operator_adjoint_hybrid_apply_3(stencilOperator, stencilOperator2, phi, x)

    ! External modules
    use parallel

    implicit none

    ! Arguments
    type(t_StencilOperator), intent(in) :: stencilOperator, stencilOperator2
    real(WP), intent(in) :: phi(:)
    real(WP), intent(inout) :: x(:,:)

    ! Local variables
    integer :: i, j, k, l, m, n, n2, p, ks, ke
    real(WP), allocatable :: xWithGhostPoints(:,:,:,:)

    allocate(xWithGhostPoints(localGridSize(1), localGridSize(2),                            &
         localGridSize(3) + sum(stencilOperator%nGhost), size(x, 2) + 1)) !... include phi

    p = size(xWithGhostPoints,4) !... index of phi

    do k = 1, localGridSize(3)
       do j = 1, localGridSize(2)
          do i = 1, localGridSize(1)
             do l = 1, size(x, 2)
                xWithGhostPoints(i, j, k + stencilOperator%nGhost(1), l) =                   &
                     x(i + localGridSize(1) * (j - 1 + localGridSize(2) * (k - 1)), l)
             end do
             ! x = stenciOpertor * (1-phi) * x + stenciOpertor2 * phi * x
             ! Communicare (1 - phi) since stenciOpertor has a wider stencil size than 
             ! sencilOperator2 --> x = stenciOpertor * phi * x + stenciOpertor2 * (1-phi) * x
             xWithGhostPoints(i, j, k + stencilOperator%nGhost(1), p) =                      &
                  1.0_WP - phi(i + localGridSize(1) * (j - 1 + localGridSize(2) * (k - 1)))
          end do
       end do
    end do

    call fill_ghost_points(xWithGhostPoints, 3, stencilOperator%nGhost,                      &
         stencilOperator%periodicOffset)

    ! Interior
    n = stencilOperator%interiorWidth / 2
    n2 = stencilOperator2%interiorWidth / 2
    ks = 1 + stencilOperator%nGhost(1)
    ke = localGridSize(3) + stencilOperator%nGhost(1)
    if (stencilOperator%nGhost(1) .eq. 0) ks = ks + stencilOperator%boundaryDepth
    if (stencilOperator%nGhost(2) .eq. 0) ke = ke - stencilOperator%boundaryDepth

    ! Save FLOPS based on symmetry of interior stencil
    ! Assumption : Either stencilOperator%symmetryType = stencilOperator2%symmetryType OR    &
    ! stencilOperator%symmetryType = ASYMMETRIC
    select case (stencilOperator%symmetryType)

    case (SKEW_SYMMETRIC)
       do l = 1, size(x, 2)
          do k = ks, ke
             do j = 1, localGridSize(2)
                do i = 1, localGridSize(1)
                   x(i + localGridSize(1) * (j - 1 + localGridSize(2) * (k - 1 -             &
                        stencilOperator%nGhost(1))), l) =                                    &
                        sum(stencilOperator%rhsInterior(1:n) *                               &
                        (xWithGhostPoints(i,j,k+1:k+n,l) *                                   &
                        xWithGhostPoints(i,j,k+1:k+n,p) -                                    &
                        xWithGhostPoints(i,j,k-1:k-n:-1,l) *                                 &
                        xWithGhostPoints(i,j,k-1:k-n:-1,p)))                                 &
                        +                                                                    &
                        sum(stencilOperator2%rhsInterior(1:n2) *                             &
                        (xWithGhostPoints(i,j,k+1:k+n2,l) *                                  &
                        (1.0_WP - xWithGhostPoints(i,j,k+1:k+n2,p)) -                        &
                        xWithGhostPoints(i,j,k-1:k-n2:-1,l) *                                &
                        (1.0_WP - xWithGhostPoints(i,j,k-1:k-n2:-1,p))))
                end do
             end do
          end do
       end do

    case (SYMMETRIC)
       do l = 1, size(x, 2)
          do k = ks, ke
             do j = 1, localGridSize(2)
                do i = 1, localGridSize(1)
                   x(i + localGridSize(1) * (j - 1 + localGridSize(2) * (k - 1 -             &
                        stencilOperator%nGhost(1))), l) =                                    &
                        sum(stencilOperator%rhsInterior(1:n) *                               &
                        (xWithGhostPoints(i,j,k+1:k+n,l) *                                   &
                        xWithGhostPoints(i,j,k+1:k+n,p) +                                    &
                        xWithGhostPoints(i,j,k-1:k-n:-1,l) *                                 &
                        xWithGhostPoints(i,j,k-1:k-n:-1,p))) +                               &
                        stencilOperator%rhsInterior(0) * xWithGhostPoints(i,j,k,l) *         &
                        xWithGhostPoints(i,j,k,p)                                            &
                        +                                                                    &
                        sum(stencilOperator2%rhsInterior(1:n2) *                             &
                        (xWithGhostPoints(i,j,k+1:k+n2,l) *                                  &
                        (1.0_WP - xWithGhostPoints(i,j,k+1:k+n2,p)) +                        &
                        xWithGhostPoints(i,j,k-1:k-n2:-1,l) *                                &
                        (1.0_WP - xWithGhostPoints(i,j,k-1:k-n2:-1,p)))) +                   &
                        stencilOperator2%rhsInterior(0) * xWithGhostPoints(i,j,k,l) *        &
                        (1.0_WP - xWithGhostPoints(i,j,k,p))
                end do
             end do
          end do
       end do

    case default
       do l = 1, size(x, 2)
          do k = ks, ke
             do j = 1, localGridSize(2)
                do i = 1, localGridSize(1)
                   x(i + localGridSize(1) * (j - 1 + localGridSize(2) * (k - 1 -             &
                        stencilOperator%nGhost(1))), l) =                                    &
                        sum(stencilOperator%rhsInterior * xWithGhostPoints(i, j,             &
                        k + lbound(stencilOperator%rhsInterior, 1) :                         &
                        k + ubound(stencilOperator%rhsInterior, 1), l) *                     &
                        xWithGhostPoints(i, j,                                               &
                        k + lbound(stencilOperator%rhsInterior, 1) :                         &
                        k + ubound(stencilOperator%rhsInterior, 1), p))                      &
                        +                                                                    &
                        sum(stencilOperator2%rhsInterior * xWithGhostPoints(i, j,            &
                        k + lbound(stencilOperator2%rhsInterior, 1) :                        &
                        k + ubound(stencilOperator2%rhsInterior, 1), l) *                    &
                        (1.0_WP - xWithGhostPoints(i, j,                                     &
                        k + lbound(stencilOperator2%rhsInterior, 1) :                        &
                        k + ubound(stencilOperator2%rhsInterior, 1), p)))
                end do
             end do
          end do
       end do

    end select

    n = stencilOperator%boundaryWidth
    n2 = stencilOperator2%boundaryWidth

    ! Left boundary
    if (stencilOperator%hasDomainBoundary(1)) then
       k = stencilOperator%nGhost(1)
       do l = 1, size(x, 2)
          do j = 1, localGridSize(2)
             do i = 1, localGridSize(1)
                do m = 1, stencilOperator%boundaryDepth
                   x(i + localGridSize(1) * (j - 1 +                                         &
                        localGridSize(2) * (k - 1 + m - stencilOperator%nGhost(1))), l) =    &
                        sum(stencilOperator%rhsBoundary1(:,m) *                              &
                        xWithGhostPoints(i,j,k+1:k+n,l) *                                    &
                        xWithGhostPoints(i,j,k+1:k+n,p))                                     &
                        +                                                                    &
                        sum(stencilOperator2%rhsBoundary1(:,m) *                             &
                        xWithGhostPoints(i,j,k+1:k+n2,l) *                                   &
                        (1.0_WP - xWithGhostPoints(i,j,k+1:k+n2,p)))
                end do
             end do
          end do
       end do
    end if

    ! Right boundary
    if (stencilOperator%hasDomainBoundary(2)) then
       k = stencilOperator%nGhost(1) + localGridSize(3) + 1
       do l = 1, size(x, 2)
          do j = 1, localGridSize(2)
             do i = 1, localGridSize(1)
                do m = stencilOperator%boundaryDepth, 1, -1
                   x(i + localGridSize(1) * (j - 1 +                                         &
                        localGridSize(2) * (k - 1 - m - stencilOperator%nGhost(1))), l) =    &
                        sum(stencilOperator%rhsBoundary2(:,m) *                              &
                        xWithGhostPoints(i,j,k-n:k-1,l) *                                    &
                        xWithGhostPoints(i,j,k-n:k-1,p))                                     &
                        +                                                                    &
                        sum(stencilOperator2%rhsBoundary2(:,m) *                             &
                        xWithGhostPoints(i,j,k-n2:k-1,l) *                                   &
                        (1.0_WP - xWithGhostPoints(i,j,k-n2:k-1,p)))
                end do
             end do
          end do
       end do
    end if

    deallocate(xWithGhostPoints)

    return
  end subroutine operator_adjoint_hybrid_apply_3


  subroutine operator_apply_and_project_boundary_1(stencilOperator, x, faceOrientation)

    implicit none

    ! Arguments
    type(t_StencilOperator), intent(in) :: stencilOperator
    real(WP), intent(inout) :: x(:,:)
    integer, intent(in) :: faceOrientation

    ! Local variables
    integer :: i, j, k, l, n

    n = stencilOperator%boundaryWidth

    if (faceOrientation .gt. 0) then

       ! Left boundary
       if (stencilOperator%hasDomainBoundary(1)) then
          do l = 1, size(x, 2)
             do k = 1, localGridSize(3)
                do j = 1, localGridSize(2)
                   i = 1 + localGridSize(1) * (j - 1 + localGridSize(2) * (k - 1))
                   x(i,l) = sum(stencilOperator%rhsBoundary1(:,1) * x(i:i+n-1,l))
                   x(i+1:i+localGridSize(1)-1,l) = 0.0_WP
                end do
             end do
          end do
       else
          x = 0.0_WP
       end if

    else if (faceOrientation .lt. 0) then

       ! Right boundary
       if (stencilOperator%hasDomainBoundary(2)) then
          do l = 1, size(x, 2)
             do k = 1, localGridSize(3)
                do j = 1, localGridSize(2)
                   i = localGridSize(1) * (j + localGridSize(2) * (k - 1))
                   x(i,l) = sum(stencilOperator%rhsBoundary2(:,1) * x(i-n+1:i,l))
                   x(i-localGridSize(1)+1:i-1,l) = 0.0_WP
                end do
             end do
          end do
       else
          x = 0.0_WP
       end if

    end if

    return
  end subroutine operator_apply_and_project_boundary_1


  subroutine operator_apply_and_project_boundary_2(stencilOperator, x, faceOrientation)

    implicit none

    ! Arguments
    type(t_StencilOperator), intent(in) :: stencilOperator
    real(WP), intent(inout) :: x(:,:)
    integer, intent(in) :: faceOrientation

    ! Local variables
    integer :: i, j, k, l, m, n

    n = stencilOperator%boundaryWidth
    m = localGridSize(1)

    if (faceOrientation .gt. 0) then

       ! Left boundary
       if (stencilOperator%hasDomainBoundary(1)) then
          do l = 1, size(x, 2)
             do k = 1, localGridSize(3)
                do i = 1, localGridSize(1)
                   j = i + localGridSize(1) * localGridSize(2) * (k - 1)
                   x(j,l) = sum(stencilOperator%rhsBoundary1(:,1) *                          &
                        x(j:j+(n-1)*m:m,l))
                   x(j+m:j+(localGridSize(2)-1)*m:m,l) = 0.0_WP
                end do
             end do
          end do
       else
          x = 0.0_WP
       end if

    else if (faceOrientation .lt. 0) then

       ! Right boundary
       if (stencilOperator%hasDomainBoundary(2)) then
          do l = 1, size(x, 2)
             do k = 1, localGridSize(3)
                do i = 1, localGridSize(1)
                   j = i + localGridSize(1) * (localGridSize(2) * k - 1)
                   x(j,l) = sum(stencilOperator%rhsBoundary2(:,1) * x(j-(n-1)*m:j:m,l))
                   x(j-(localGridSize(2)-1)*m:j-m:m,l) = 0.0_WP
                end do
             end do
          end do
       else
          x = 0.0_WP
       end if

    end if

    return
  end subroutine operator_apply_and_project_boundary_2


  subroutine operator_apply_and_project_boundary_3(stencilOperator, x, faceOrientation)

    implicit none

    ! Arguments
    type(t_StencilOperator), intent(in) :: stencilOperator
    real(WP), intent(inout) :: x(:,:)
    integer, intent(in) :: faceOrientation

    ! Local variables
    integer :: i, j, k, l, m, n

    n = stencilOperator%boundaryWidth
    m = localGridSize(1) * localGridSize(2)

    if (faceOrientation .gt. 0) then

       ! Left boundary
       if (stencilOperator%hasDomainBoundary(1)) then
          do l = 1, size(x, 2)
             do j = 1, localGridSize(2)
                do i = 1, localGridSize(1)
                   k = i + localGridSize(1) * (j - 1)
                   x(k,l) = sum(stencilOperator%rhsBoundary1(:,1) *                          &
                        x(k:k+(n-1)*m:m,l))
                   x(k+m:k+(localGridSize(3)-1)*m:m,l) = 0.0_WP
                end do
             end do
          end do
       else
          x = 0.0_WP
       end if

    else if (faceOrientation .lt. 0) then

       ! Right boundary
       if (stencilOperator%hasDomainBoundary(2)) then
          do l = 1, size(x, 2)
             do j = 1, localGridSize(2)
                do i = 1, localGridSize(1)
                   k = i + localGridSize(1) * (j - 1 + localGridSize(2) *                    &
                        (localGridSize(3) - 1))
                   x(k,l) = sum(stencilOperator%rhsBoundary2(:,1) * x(k-(n-1)*m:k:m,l))
                   x(k-(localGridSize(3)-1)*m:k-m:m,l) = 0.0_WP
                end do
             end do
          end do
       else
          x = 0.0_WP
       end if

    end if

    return
  end subroutine operator_apply_and_project_boundary_3


  subroutine operator_project_boundary_and_apply_1(stencilOperator, x, faceOrientation)

    implicit none

    ! Arguments
    type(t_StencilOperator), intent(in) :: stencilOperator
    real(WP), intent(inout) :: x(:,:)
    integer, intent(in) :: faceOrientation

    ! Local variables
    integer :: i, j, k, l, m, n
    real(WP), allocatable :: x_(:,:,:,:)

    allocate(x_(localGridSize(1), localGridSize(2), localGridSize(3), size(x, 2)))

    do l = 1, size(x, 2)
       do k = 1, localGridSize(3)
          do j = 1, localGridSize(2)
             do i = 1, localGridSize(1)
                x_(i,j,k,l) = x(i + localGridSize(1) * (j - 1 + localGridSize(2) *           &
                     (k - 1)), l)
             end do
          end do
       end do
    end do

    n = stencilOperator%boundaryWidth

    x = 0.0_WP

    if (faceOrientation .gt. 0) then

       ! Left boundary
       if (stencilOperator%hasDomainBoundary(1)) then
          i = 1
          do l = 1, size(x, 2)
             do k = 1, localGridSize(3)
                do j = 1, localGridSize(2)
                   do m = 1, stencilOperator%boundaryDepth
                      x(i + m - 1 + localGridSize(1) * (j - 1 + localGridSize(2) *           &
                           (k - 1)), l) = stencilOperator%rhsBoundary1(1,m) * x_(i,j,k,l)
                   end do
                end do
             end do
          end do
       end if

    else

       ! Right boundary
       if (stencilOperator%hasDomainBoundary(2)) then
          i = localGridSize(1)
          do l = 1, size(x, 2)
             do k = 1, localGridSize(3)
                do j = 1, localGridSize(2)
                   do m = stencilOperator%boundaryDepth, 1, -1
                      x(i - m + 1 + localGridSize(1) * (j - 1 + localGridSize(2) *           &
                           (k - 1)), l) = stencilOperator%rhsBoundary2(n,m) * x_(i,j,k,l)
                   end do
                end do
             end do
          end do
       end if

    end if

    deallocate(x_)

    return
  end subroutine operator_project_boundary_and_apply_1


  subroutine operator_project_boundary_and_apply_2(stencilOperator, x, faceOrientation)

    implicit none

    ! Arguments
    type(t_StencilOperator), intent(in) :: stencilOperator
    real(WP), intent(inout) :: x(:,:)
    integer, intent(in) :: faceOrientation

    ! Local variables
    integer :: i, j, k, l, m, n
    real(WP), allocatable :: x_(:,:,:,:)

    allocate(x_(localGridSize(1), localGridSize(2), localGridSize(3), size(x, 2)))

    do l = 1, size(x, 2)
       do k = 1, localGridSize(3)
          do j = 1, localGridSize(2)
             do i = 1, localGridSize(1)
                x_(i,j,k,l) = x(i + localGridSize(1) * (j - 1 + localGridSize(2) *           &
                     (k - 1)), l)
             end do
          end do
       end do
    end do

    n = stencilOperator%boundaryWidth

    x = 0.0_WP

    if (faceOrientation .gt. 0) then

       ! Left boundary
       if (stencilOperator%hasDomainBoundary(1)) then
          j = 1
          do l = 1, size(x, 2)
             do k = 1, localGridSize(3)
                do i = 1, localGridSize(1)
                   do m = 1, stencilOperator%boundaryDepth
                      x(i + localGridSize(1) * (j + m - 2 + localGridSize(2) *               &
                           (k - 1)), l) = stencilOperator%rhsBoundary1(1,m) * x_(i,j,k,l)
                   end do
                end do
             end do
          end do
       end if

    else

       ! Right boundary
       if (stencilOperator%hasDomainBoundary(2)) then
          j = localGridSize(2)
          do l = 1, size(x, 2)
             do k = 1, localGridSize(3)
                do i = 1, localGridSize(1)
                   do m = stencilOperator%boundaryDepth, 1, -1
                      x(i + localGridSize(1) * (j - m + localGridSize(2) *                   &
                           (k - 1)), l) = stencilOperator%rhsBoundary2(n,m) * x_(i,j,k,l)
                   end do
                end do
             end do
          end do
       end if

    end if

    deallocate(x_)

    return
  end subroutine operator_project_boundary_and_apply_2


  subroutine operator_project_boundary_and_apply_3(stencilOperator, x, faceOrientation)

    implicit none

    ! Arguments
    type(t_StencilOperator), intent(in) :: stencilOperator
    real(WP), intent(inout) :: x(:,:)
    integer, intent(in) :: faceOrientation

    ! Local variables
    integer :: i, j, k, l, m, n
    real(WP), allocatable :: x_(:,:,:,:)

    allocate(x_(localGridSize(1), localGridSize(2), localGridSize(3), size(x, 2)))

    do l = 1, size(x, 2)
       do k = 1, localGridSize(3)
          do j = 1, localGridSize(2)
             do i = 1, localGridSize(1)
                x_(i,j,k,l) = x(i + localGridSize(1) * (j - 1 + localGridSize(2) *           &
                     (k - 1)), l)
             end do
          end do
       end do
    end do

    n = stencilOperator%boundaryWidth

    x = 0.0_WP

    if (faceOrientation .gt. 0) then

       ! Left boundary
       if (stencilOperator%hasDomainBoundary(1)) then
          k = 1
          do l = 1, size(x, 2)
             do j = 1, localGridSize(2)
                do i = 1, localGridSize(1)
                   do m = 1, stencilOperator%boundaryDepth
                      x(i + localGridSize(1) * (j - 1 + localGridSize(2) *                   &
                           (k + m - 2)), l) = stencilOperator%rhsBoundary1(1,m) * x_(i,j,k,l)
                   end do
                end do
             end do
          end do
       end if

    else

       ! Right boundary
       if (stencilOperator%hasDomainBoundary(2)) then
          k = localGridSize(3)
          do l = 1, size(x, 2)
             do j = 1, localGridSize(2)
                do i = 1, localGridSize(1)
                   do m = stencilOperator%boundaryDepth, 1, -1
                      x(i + localGridSize(1) * (j - 1 + localGridSize(2) * (k - m)), l) =    &
                           stencilOperator%rhsBoundary2(n,m) * x_(i,j,k,l)
                   end do
                end do
             end do
          end do
       end if

    end if

    deallocate(x_)

    return
  end subroutine operator_project_boundary_and_apply_3


  subroutine operator_apply_norm_1(stencilOperator, x)

    implicit none

    ! Arguments
    type(t_StencilOperator), intent(in) :: stencilOperator
    real(WP), intent(inout) :: x(:,:)

    ! Local variables
    integer :: i, j, k, l

    ! Left boundary
    if (stencilOperator%hasDomainBoundary(1)) then
       do l = 1, size(x, 2)
          do k = 1, localGridSize(3)
             do j = 1, localGridSize(2)
                do i = 1, stencilOperator%boundaryDepth
                   x(i + localGridSize(1) * (j - 1 + localGridSize(2) * (k - 1)), l) =       &
                        stencilOperator%normBoundary(i) *                                    &
                        x(i + localGridSize(1) * (j - 1 + localGridSize(2) * (k - 1)), l)
                end do
             end do
          end do
       end do
    end if

    ! Right boundary
    if (stencilOperator%hasDomainBoundary(2)) then
       do l = 1, size(x, 2)
          do k = 1, localGridSize(3)
             do j = 1, localGridSize(2)
                do i = localGridSize(1) + 1 - stencilOperator%boundaryDepth, localGridSize(1)
                   x(i + localGridSize(1) * (j - 1 + localGridSize(2) * (k - 1)), l) =       &
                        stencilOperator%normBoundary(localGridSize(1) + 1 - i) *             &
                        x(i + localGridSize(1) * (j - 1 + localGridSize(2) * (k - 1)), l)
                end do
             end do
          end do
       end do
    end if

    return
  end subroutine operator_apply_norm_1


  subroutine operator_apply_norm_2(stencilOperator, x)

    implicit none

    ! Arguments
    type(t_StencilOperator), intent(in) :: stencilOperator
    real(WP), intent(inout) :: x(:,:)

    ! Local variables
    integer :: i, j, k, l

    ! Left boundary
    if (stencilOperator%hasDomainBoundary(1)) then
       do l = 1, size(x, 2)
          do k = 1, localGridSize(3)
             do j = 1, stencilOperator%boundaryDepth
                do i = 1, localGridSize(1)
                   x(i + localGridSize(1) * (j - 1 + localGridSize(2) * (k - 1)), l) =       &
                        stencilOperator%normBoundary(j) *                                    &
                        x(i + localGridSize(1) * (j - 1 + localGridSize(2) * (k - 1)), l)
                end do
             end do
          end do
       end do
    end if

    ! Right boundary
    if (stencilOperator%hasDomainBoundary(2)) then
       do l = 1, size(x, 2)
          do k = 1, localGridSize(3)
             do j = localGridSize(2) + 1 - stencilOperator%boundaryDepth, localGridSize(2)
                do i = 1, localGridSize(1)
                   x(i + localGridSize(1) * (j - 1 + localGridSize(2) * (k - 1)), l) =       &
                        stencilOperator%normBoundary(localGridSize(2) + 1 - j) *             &
                        x(i + localGridSize(1) * (j - 1 + localGridSize(2) * (k - 1)), l)
                end do
             end do
          end do
       end do
    end if

    return
  end subroutine operator_apply_norm_2


  subroutine operator_apply_norm_3(stencilOperator, x)

    implicit none

    ! Arguments
    type(t_StencilOperator), intent(in) :: stencilOperator
    real(WP), intent(inout) :: x(:,:)

    ! Local variables
    integer :: i, j, k, l

    ! Left boundary
    if (stencilOperator%hasDomainBoundary(1)) then
       do l = 1, size(x, 2)
          do k = 1, stencilOperator%boundaryDepth
             do j = 1, localGridSize(2)
                do i = 1, localGridSize(1)
                   x(i + localGridSize(1) * (j - 1 + localGridSize(2) * (k - 1)), l) =       &
                        stencilOperator%normBoundary(k) *                                    &
                        x(i + localGridSize(1) * (j - 1 + localGridSize(2) * (k - 1)), l)
                end do
             end do
          end do
       end do
    end if

    ! Right boundary
    if (stencilOperator%hasDomainBoundary(2)) then
       do l = 1, size(x, 2)
          do k = localGridSize(3) + 1 - stencilOperator%boundaryDepth, localGridSize(3)
             do j = 1, localGridSize(2)
                do i = 1, localGridSize(1)
                   x(i + localGridSize(1) * (j - 1 + localGridSize(2) * (k - 1)), l) =       &
                        stencilOperator%normBoundary(localGridSize(3) + 1 - k) *             &
                        x(i + localGridSize(1) * (j - 1 + localGridSize(2) * (k - 1)), l)
                end do
             end do
          end do
       end do
    end if

    return
  end subroutine operator_apply_norm_3


  subroutine operator_apply_norm_inverse_1(stencilOperator, x)

    implicit none

    ! Arguments
    type(t_StencilOperator), intent(in) :: stencilOperator
    real(WP), intent(inout) :: x(:,:)

    ! Local variables
    integer :: i, j, k, l

    ! Left boundary
    if (stencilOperator%hasDomainBoundary(1)) then
       do l = 1, size(x, 2)
          do k = 1, localGridSize(3)
             do j = 1, localGridSize(2)
                do i = 1, stencilOperator%boundaryDepth
                   x(i + localGridSize(1) * (j - 1 + localGridSize(2) * (k - 1)), l) =       &
                        x(i + localGridSize(1) * (j - 1 + localGridSize(2) * (k - 1)), l) /  &
                        stencilOperator%normBoundary(i)
                end do
             end do
          end do
       end do
    end if

    ! Right boundary
    if (stencilOperator%hasDomainBoundary(2)) then
       do l = 1, size(x, 2)
          do k = 1, localGridSize(3)
             do j = 1, localGridSize(2)
                do i = localGridSize(1) + 1 - stencilOperator%boundaryDepth, localGridSize(1)
                   x(i + localGridSize(1) * (j - 1 + localGridSize(2) * (k - 1)), l) =       &
                        x(i + localGridSize(1) * (j - 1 + localGridSize(2) * (k - 1)), l) /  &
                        stencilOperator%normBoundary(localGridSize(1) + 1 - i)

                end do
             end do
          end do
       end do
    end if

    return
  end subroutine operator_apply_norm_inverse_1


  subroutine operator_apply_norm_inverse_2(stencilOperator, x)

    implicit none

    ! Arguments
    type(t_StencilOperator), intent(in) :: stencilOperator
    real(WP), intent(inout) :: x(:,:)

    ! Local variables
    integer :: i, j, k, l

    ! Left boundary
    if (stencilOperator%hasDomainBoundary(1)) then
       do l = 1, size(x, 2)
          do k = 1, localGridSize(3)
             do j = 1, stencilOperator%boundaryDepth
                do i = 1, localGridSize(1)
                   x(i + localGridSize(1) * (j - 1 + localGridSize(2) * (k - 1)), l) =       &
                        x(i + localGridSize(1) * (j - 1 + localGridSize(2) * (k - 1)), l) /  &
                        stencilOperator%normBoundary(j)
                end do
             end do
          end do
       end do
    end if

    ! Right boundary
    if (stencilOperator%hasDomainBoundary(2)) then
       do l = 1, size(x, 2)
          do k = 1, localGridSize(3)
             do j = localGridSize(2) + 1 - stencilOperator%boundaryDepth, localGridSize(2)
                do i = 1, localGridSize(1)
                   x(i + localGridSize(1) * (j - 1 + localGridSize(2) * (k - 1)), l) =       &
                        x(i + localGridSize(1) * (j - 1 + localGridSize(2) * (k - 1)), l) /  &
                        stencilOperator%normBoundary(localGridSize(2) + 1 - j)
                end do
             end do
          end do
       end do
    end if

    return
  end subroutine operator_apply_norm_inverse_2


  subroutine operator_apply_norm_inverse_3(stencilOperator, x)

    implicit none

    ! Arguments
    type(t_StencilOperator), intent(in) :: stencilOperator
    real(WP), intent(inout) :: x(:,:)

    ! Local variables
    integer :: i, j, k, l

    ! Left boundary
    if (stencilOperator%hasDomainBoundary(1)) then
       do l = 1, size(x, 2)
          do k = 1, stencilOperator%boundaryDepth
             do j = 1, localGridSize(2)
                do i = 1, localGridSize(1)
                   x(i + localGridSize(1) * (j - 1 + localGridSize(2) * (k - 1)), l) =       &
                        x(i + localGridSize(1) * (j - 1 + localGridSize(2) * (k - 1)), l) /  &
                        stencilOperator%normBoundary(k)
                end do
             end do
          end do
       end do
    end if

    ! Right boundary
    if (stencilOperator%hasDomainBoundary(2)) then
       do l = 1, size(x, 2)
          do k = localGridSize(3) + 1 - stencilOperator%boundaryDepth, localGridSize(3)
             do j = 1, localGridSize(2)
                do i = 1, localGridSize(1)
                   x(i + localGridSize(1) * (j - 1 + localGridSize(2) * (k - 1)), l) =       &
                        x(i + localGridSize(1) * (j - 1 + localGridSize(2) * (k - 1)), l) /  &
                        stencilOperator%normBoundary(localGridSize(3) + 1 - k)
                end do
             end do
          end do
       end do
    end if

    return
  end subroutine operator_apply_norm_inverse_3


  subroutine operator_apply_implicit_1(stencilOperator, x)

    ! External modules
    use parallel

    implicit none

    ! Arguments
    type(t_StencilOperator), intent(in) :: stencilOperator
    real(WP), intent(inout) :: x(:,:)

    ! Local variables
    integer :: i, j, k, l, n, gridIndex
    real(WP), allocatable :: RHS(:,:,:), A(:,:,:,:)

    if (.not. stencilOperator%implicit) return

    if (stencilOperator%nDiagonals .ne. 3)                                                   &
         call die('implicit operator not yet implemented for non-tridiagonal systems')

    n = stencilOperator%nDiagonals / 2    
    allocate(RHS(iStart(2):iEnd(2),iStart(3):iEnd(3),iStart(1):iEnd(1)))
    allocate(A(iStart(2):iEnd(2),iStart(3):iEnd(3),iStart(1):iEnd(1),-n:n))
    
    do l = 1, size(x,2)
       do i = iStart(1), iEnd(1)
          do k = iStart(3), iEnd(3)
             do j = iStart(2), iEnd(2)
                gridIndex = i - gridOffset(1) + localGridSize(1) *                           &
                     (j - 1 - gridOffset(2) + localGridSize(2) *                             &
                     (k - 1 - gridOffset(3)))
                RHS(j,k,i) = x(gridIndex, l)
                if (stencilOperator%nGhost(1).eq.0 .and.                                     &
                     i.le.stencilOperator%boundaryDepth) then
                   A(j,k,i,-n:+n) = stencilOperator%lhsBoundary1(-n:+n,i)
                else if (stencilOperator%nGhost(2).eq.0 .and.                                &
                     i.ge.globalGridSize(1)-stencilOperator%boundaryDepth+1) then
                   A(j,k,i,-n:+n) = stencilOperator%lhsBoundary2(-n:+n,globalGridSize(1)-i+1)
                else
                   A(j,k,i,-n:+n) = stencilOperator%lhsInterior
                end if
             end do
          end do
       end do
       call tridiagonal(                                                                     &
            A(iStart(2),iStart(3),iStart(1),-1),                                             &
            A(iStart(2),iStart(3),iStart(1), 0),                                             &
            A(iStart(2),iStart(3),iStart(1),+1),                                             &    
            RHS,localGridSize(1),localGridSize(2)*localGridSize(3),1,isPeriodic)
       do i = iStart(1), iEnd(1)
          do k = iStart(3), iEnd(3)
             do j = iStart(2), iEnd(2)
                gridIndex = i - gridOffset(1) + localGridSize(1) *                           &
                     (j - 1 - gridOffset(2) + localGridSize(2) *                             &
                     (k - 1 - gridOffset(3)))
                x(gridIndex, l) = RHS(j, k, i)
             end do
          end do
       end do
    end do

    ! Cleanup
    deallocate(RHS, A)

    return
  end subroutine operator_apply_implicit_1


  subroutine operator_apply_implicit_2(stencilOperator, x)

    ! External modules
    use parallel

    implicit none

    ! Arguments
    type(t_StencilOperator), intent(in) :: stencilOperator
    real(WP), intent(inout) :: x(:,:)

    ! Local variables
    integer :: i, j, k, l, n, gridIndex
    real(WP), allocatable :: RHS(:,:,:), A(:,:,:,:)
    
    if (.not. stencilOperator%implicit) return

    if (stencilOperator%nDiagonals .ne. 3)                                                   &
         call die('implicit operator not yet implemented for non-tridiagonal systems')

    n = stencilOperator%nDiagonals / 2    
    allocate(RHS(iStart(1):iEnd(1),iStart(3):iEnd(3),iStart(2):iEnd(2)))
    allocate(A(iStart(1):iEnd(1),iStart(3):iEnd(3),iStart(2):iEnd(2),-n:n))

    do l = 1, size(x,2)
       do j = iStart(2), iEnd(2)
          do k = iStart(3), iEnd(3)
             do i = iStart(1), iEnd(1)
                gridIndex = i - gridOffset(1) + localGridSize(1) *                           &
                     (j - 1 - gridOffset(2) + localGridSize(2) *                             &
                     (k - 1 - gridOffset(3)))
                RHS(i,k,j) = x(gridIndex, l)
                if (stencilOperator%nGhost(1).eq.0 .and.                                     &
                     j.le.stencilOperator%boundaryDepth) then
                   A(i,k,j,-n:+n) = stencilOperator%lhsBoundary1(-n:+n,j)
                else if (stencilOperator%nGhost(2).eq.0 .and.                                &
                     j.ge.globalGridSize(2)-stencilOperator%boundaryDepth+1) then
                   A(i,k,j,-n:+n) = stencilOperator%lhsBoundary2(-n:+n,globalGridSize(2)-j+1)
                else
                   A(i,k,j,-n:+n) = stencilOperator%lhsInterior
                end if

             end do
          end do
       end do
       call tridiagonal(                                                                     &
            A(iStart(1),iStart(3),iStart(2),-1),                                             &
            A(iStart(1),iStart(3),iStart(2), 0),                                             &
            A(iStart(1),iStart(3),iStart(2),+1),                                             &
            RHS,localGridSize(2),localGridSize(1)*localGridSize(3),2,isPeriodic)
       do j = iStart(2), iEnd(2)
          do k = iStart(3), iEnd(3)
             do i = iStart(1), iEnd(1)
                gridIndex = i - gridOffset(1) + localGridSize(1) *                           &
                     (j - 1 - gridOffset(2) + localGridSize(2) *                             &
                     (k - 1 - gridOffset(3)))
                x(gridIndex, l) = RHS(i, k, j)
             end do
          end do
       end do
    end do

    ! Cleanup
    deallocate(RHS, A)

    return
  end subroutine operator_apply_implicit_2


  subroutine operator_apply_implicit_3(stencilOperator, x)

    ! External modules
    use parallel

    implicit none

    ! Arguments
    type(t_StencilOperator), intent(in) :: stencilOperator
    real(WP), intent(inout) :: x(:,:)

    ! Local variables
    integer :: i, j, k, l, n, gridIndex
    real(WP), allocatable :: RHS(:,:,:), A(:,:,:,:)

    if (.not. stencilOperator%implicit) return

    if (stencilOperator%nDiagonals .ne. 3)                                                   &
         call die('implicit operator not yet implemented for non-tridiagonal systems')

    n = stencilOperator%nDiagonals / 2    
    allocate(RHS(iStart(1):iEnd(1),iStart(2):iEnd(2),iStart(3):iEnd(3)))
    allocate(A(iStart(1):iEnd(1),iStart(2):iEnd(2),iStart(3):iEnd(3),-n:n))

    do l = 1, size(x,2)
       do k = iStart(3), iEnd(3)
          do j = iStart(2), iEnd(2)
             do i = iStart(1), iEnd(1)
                gridIndex = i - gridOffset(1) + localGridSize(1) *                           &
                     (j - 1 - gridOffset(2) + localGridSize(2) *                             &
                     (k - 1 - gridOffset(3)))
                RHS(i,j,k) = x(gridIndex, l)
                if (stencilOperator%nGhost(1).eq.0 .and.                                     &
                     k.le.stencilOperator%boundaryDepth) then
                   A(i,j,k,-n:+n) = stencilOperator%lhsBoundary1(-n:+n,k)
                else if (stencilOperator%nGhost(2).eq.0 .and.                                &
                     k.ge.globalGridSize(3)-stencilOperator%boundaryDepth+1) then
                   A(i,j,k,-n:+n) = stencilOperator%lhsBoundary2(-n:+n,globalGridSize(3)-k+1)
                else
                   A(i,j,k,-n:+n) = stencilOperator%lhsInterior
                end if


             end do
          end do
       end do
       call tridiagonal(                                                                     &
            A(iStart(1),iStart(2),iStart(3),-1),                                             &
            A(iStart(1),iStart(2),iStart(3), 0),                                             &
            A(iStart(1),iStart(2),iStart(3),+1),                                             &
            RHS,localGridSize(3),localGridSize(1)*localGridSize(2),3,isPeriodic)
       do k = iStart(3), iEnd(3)
          do j = iStart(2), iEnd(2)
             do i = iStart(1), iEnd(1)
                gridIndex = i - gridOffset(1) + localGridSize(1) *                           &
                     (j - 1 - gridOffset(2) + localGridSize(2) *                             &
                     (k - 1 - gridOffset(3)))
                x(gridIndex, l) = RHS(i, j, k)
             end do
          end do
       end do
    end do

    ! Cleanup
    deallocate(RHS, A)

    return
  end subroutine operator_apply_implicit_3


  subroutine print_operator(stencilOperator)

    implicit none

    ! Arguments
    type(t_StencilOperator), intent(in) :: stencilOperator

    ! Local variables
    integer :: j

    print *
    print *,'---------------------------------------------'
    print *,'    PRINT STENCIL OPERATOR IN MATRIX FORM    '

    print *

    print *, 'STENCIL OPERATOR TYPE: ', stencilOperator%symmetryType

    print *,'LEFT BOUNDARY'
    do j = 1, stencilOperator%boundaryDepth
       print *, 'ROW', j, 'STENCIL',                                                        &
            stencilOperator%rhsBoundary1(1:stencilOperator%boundaryWidth,j)
    end do

    print *
    print *,'INTERIOR'
    print *,stencilOperator%rhsInterior

    print *
    print *,'RIGHT BOUNDARY'
    do j = stencilOperator%boundaryDepth, 1, -1
       print *, 'ROW', j, 'STENCIL',                                                        &
            stencilOperator%rhsBoundary2(1:stencilOperator%boundaryWidth,j)
    end do

    print *
    print *,'---------------------------------------------'
    print *

  end subroutine print_operator

end module operator


! =========================== !
! Setup the stencil operators !
! =========================== !
subroutine operator_setup

  ! Internal modules
  use operator

  ! External modules
  use parser
  use simulation_flags

  implicit none

  ! Read the default discretization scheme from input
  call parser_read('default discretization scheme', discretizationType, 'SBP 4-8')

  ! Get distance to walls
  if (useIblank) then
     allocate(distanceToWall(nGridPoints, nDimensions, 2))
     call get_distance_to_wall(distanceToWall)
  end if

  ! First derivative and adjoint derivative operators
  call first_derivative_setup

  ! Second derivative operators
  call second_derivative_setup

  ! Fourth derivative operators
  call fourth_derivative_setup

  ! Artificial dissipation operators
  call dissipation_setup

  ! Filtering
  call filter_setup

  return
end subroutine operator_setup


! ============================= !
! Cleanup the stencil operators !
! ============================= !
subroutine operator_cleanup

  ! Internal modules
  use operator

  implicit none

  call first_derivative_cleanup
  call second_derivative_cleanup
  call fourth_derivative_cleanup
  call dissipation_cleanup
  call filter_cleanup

  if (allocated(distanceToWall)) deallocate(distanceToWall)

  return
end subroutine operator_cleanup
