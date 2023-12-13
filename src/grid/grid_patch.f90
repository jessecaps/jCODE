module grid_patch

  ! External modules
  use precision
  use string
  use parallel
  use grid

  implicit none

  ! Patch types
  integer, parameter ::                                                                      &
       SPONGE_PATCH          = 0,                                                            &
       SAT_FAR_FIELD         = 1,                                                            &
       SAT_SLIP_WALL         = 2,                                                            &
       SAT_ISOTHERMAL_WALL   = 3,                                                            &
       SAT_ADIABATIC_WALL    = 4,                                                            &
       ACTUATOR              = 5,                                                            &
       COST_TARGET           = 6,                                                            &
       VISUALIZATION         = 7,                                                            &
       IBLANK_PATCH          = 8,                                                            &
       SAT_OUTFLOW           = 9,                                                            &
       SAT_INFLOW            = 10,                                                           &
       DIRICHLET_BC          = 11,                                                           &
       NEUMANN_BC            = 12,                                                           &
       SLIP_BC               = 13,                                                           &
       ISOTHERMAL_BC         = 14,                                                           &
       EXCITATION_PATCH      = 15,                                                           &
       STATISTICS_PATCH      = 16

  type :: t_Patch
     character(len = str_medium) :: name
     integer :: patchType, normalDirection, iMin, iMax, jMin, jMax, kMin, kMax,              &
          globalSize(3), localSize(3), offset(3), iStart(3), iEnd(3), nPatchPoints = 0,      &
          comm = MPI_COMM_NULL, mpiRealSubarrayType = MPI_DATATYPE_NULL, masterRank
     integer, allocatable :: hole(:)
     real(WP), allocatable :: norm(:,:)
  end type t_Patch

  integer :: nPatches
  type(t_Patch), pointer :: patches(:)

  interface patch_collect
     module procedure patch_collect_scalar
     module procedure patch_collect_vector
     module procedure patch_collect_tensor
     module procedure patch_collect_integer
  end interface patch_collect

  interface patch_disperse
     module procedure patch_disperse_scalar
     module procedure patch_disperse_vector
     module procedure patch_disperse_tensor
  end interface patch_disperse

contains

  ! Quick-sort routine to sort patches by type
  ! ------------------------------------------
  recursive subroutine patch_sort(patch, first, last)

    implicit none

    ! Arguments
    type(t_Patch), intent(inout) :: patch(:)
    integer, intent(in) :: first, last

    ! Local variables
    type(t_Patch) :: patch1, patch2
    integer :: i, j

    patch1 = patch( (first + last) / 2 )
    i = first
    j = last
    do
       do while (patch(i)%patchType .lt. patch1%patchType)
          i = i + 1
       end do
       do while (patch1%patchType .lt. patch(j)%patchType)
          j = j - 1
       end do
       if (i .ge. j) exit
       patch2 = patch(i)
       patch(i) = patch(j)
       patch(j) = patch2
       i = i + 1
       j = j - 1
    end do

    if (first .lt. i - 1) call patch_sort(patch, first, i - 1)
    if (j + 1 .lt. last)  call patch_sort(patch, j + 1, last)

    return
  end subroutine patch_sort


  ! Collect scalar from grid to patch
  ! ---------------------------------
  subroutine patch_collect_scalar(patch, gridArray, patchArray)

    implicit none

    ! Arguments
    type(t_Patch), intent(in) :: patch
    real(WP), intent(in) :: gridArray(:)
    real(WP), intent(out) :: patchArray(:)

    ! Local variables
    integer :: i, j, k, patchIndex, localIndex

    if (patch%nPatchPoints .eq. 0) return

    do k = patch%iStart(3), patch%iEnd(3)
       do j = patch%iStart(2), patch%iEnd(2)
          do i = patch%iStart(1), patch%iEnd(1)
             patchIndex = i - patch%offset(1) +                                              &
                  patch%localSize(1) * (j - 1 - patch%offset(2) +                            &
                  patch%localSize(2) * (k - 1 - patch%offset(3)))
             localIndex = i - gridOffset(1) +                                                &
                  localGridSize(1) * (j - 1 - gridOffset(2) +                                &
                  localGridSize(2) * (k - 1 - gridOffset(3)))
             patchArray(patchIndex) = gridArray(localIndex)
          end do
       end do
    end do

    return
  end subroutine patch_collect_scalar


  ! Collect vector from grid to patch
  ! ---------------------------------
  subroutine patch_collect_vector(patch, gridArray, patchArray)

    implicit none

    ! Arguments
    type(t_Patch), intent(in) :: patch
    real(WP), intent(in) :: gridArray(:,:)
    real(WP), intent(out) :: patchArray(:,:)

    ! Local variables
    integer :: i, j, k, l, patchIndex, localIndex

    if (patch%nPatchPoints .eq. 0) return

    do l = 1, size(patchArray, 2)
       do k = patch%iStart(3), patch%iEnd(3)
          do j = patch%iStart(2), patch%iEnd(2)
             do i = patch%iStart(1), patch%iEnd(1)
                patchIndex = i - patch%offset(1) +                                           &
                     patch%localSize(1) * (j - 1 - patch%offset(2) +                         &
                     patch%localSize(2) * (k - 1 - patch%offset(3)))
                localIndex = i - gridOffset(1) +                                             &
                     localGridSize(1) * (j - 1 - gridOffset(2) +                             &
                     localGridSize(2) * (k - 1 - gridOffset(3)))
                patchArray(patchIndex,l) = gridArray(localIndex,l)
             end do
          end do
       end do
    end do

    return
  end subroutine patch_collect_vector


  ! Collect tensor from grid to patch
  ! ---------------------------------
  subroutine patch_collect_tensor(patch, gridArray, patchArray)

    implicit none

    ! Arguments
    type(t_Patch), intent(in) :: patch
    real(WP), intent(in) :: gridArray(:,:,:)
    real(WP), intent(out) :: patchArray(:,:,:)

    ! Local variables
    integer :: i, j, k, l, m, patchIndex, localIndex

    if (patch%nPatchPoints .eq. 0) return

    do m = 1, size(patchArray, 3)
       do l = 1, size(patchArray, 2)
          do k = patch%iStart(3), patch%iEnd(3)
             do j = patch%iStart(2), patch%iEnd(2)
                do i = patch%iStart(1), patch%iEnd(1)
                   patchIndex = i - patch%offset(1) +                                        &
                        patch%localSize(1) * (j - 1 - patch%offset(2) +                      &
                        patch%localSize(2) * (k - 1 - patch%offset(3)))
                   localIndex = i - gridOffset(1) +                                          &
                        localGridSize(1) * (j - 1 - gridOffset(2) +                          &
                        localGridSize(2) * (k - 1 - gridOffset(3)))
                   patchArray(patchIndex,l,m) = gridArray(localIndex,l,m)
                end do
             end do
          end do
       end do
    end do

    return
  end subroutine patch_collect_tensor


  ! Collect integer from grid to patch
  ! ----------------------------------
  subroutine patch_collect_integer(patch, gridArray, patchArray)

    implicit none

    ! Arguments
    type(t_Patch), intent(in) :: patch
    integer, intent(in) :: gridArray(:)
    integer, intent(out) :: patchArray(:)

    ! Local variables
    integer :: i, j, k, patchIndex, localIndex

    if (patch%nPatchPoints .eq. 0) return

    do k = patch%iStart(3), patch%iEnd(3)
       do j = patch%iStart(2), patch%iEnd(2)
          do i = patch%iStart(1), patch%iEnd(1)
             patchIndex = i - patch%offset(1) +                                              &
                  patch%localSize(1) * (j - 1 - patch%offset(2) +                            &
                  patch%localSize(2) * (k - 1 - patch%offset(3)))
             localIndex = i - gridOffset(1) +                                                &
                  localGridSize(1) * (j - 1 - gridOffset(2) +                                &
                  localGridSize(2) * (k - 1 - gridOffset(3)))
             patchArray(patchIndex) = gridArray(localIndex)
          end do
       end do
    end do

    return
  end subroutine patch_collect_integer


  ! Disperse a scalar from patch to grid
  ! ------------------------------------
  subroutine patch_disperse_scalar(patch, patchArray, gridArray)

    implicit none

    ! Arguments
    type(t_Patch), intent(in) :: patch
    real(WP), intent(in) :: patchArray(:)
    real(WP), intent(out) :: gridArray(:)

    ! Local variables
    integer :: i, j, k, patchIndex, localIndex

    gridArray = 0.0_WP

    if (patch%nPatchPoints .eq. 0) return

    do k = patch%iStart(3), patch%iEnd(3)
       do j = patch%iStart(2), patch%iEnd(2)
          do i = patch%iStart(1), patch%iEnd(1)
             patchIndex = i - patch%offset(1) +                                              &
                  patch%localSize(1) * (j - 1 - patch%offset(2) +                            &
                  patch%localSize(2) * (k - 1 - patch%offset(3)))
             localIndex = i - gridOffset(1) +                                                &
                  localGridSize(1) * (j - 1 - gridOffset(2) +                                &
                  localGridSize(2) * (k - 1 - gridOffset(3)))
             gridArray(localIndex) = patchArray(patchIndex)
          end do
       end do
    end do

    return
  end subroutine patch_disperse_scalar


  ! Disperse a vector from patch to grid
  ! ------------------------------------
  subroutine patch_disperse_vector(patch, patchArray, gridArray)

    implicit none

    ! Arguments
    type(t_Patch), intent(in) :: patch
    real(WP), intent(in) :: patchArray(:,:)
    real(WP), intent(out) :: gridArray(:,:)

    ! Local variables
    integer :: i, j, k, l, patchIndex, localIndex

    gridArray = 0.0_WP

    if (patch%nPatchPoints .eq. 0) return

    do l = 1, size(patchArray, 2)
       do k = patch%iStart(3), patch%iEnd(3)
          do j = patch%iStart(2), patch%iEnd(2)
             do i = patch%iStart(1), patch%iEnd(1)
                patchIndex = i - patch%offset(1) +                                           &
                     patch%localSize(1) * (j - 1 - patch%offset(2) +                         &
                     patch%localSize(2) * (k - 1 - patch%offset(3)))
                localIndex = i - gridOffset(1) +                                             &
                     localGridSize(1) * (j - 1 - gridOffset(2) +                             &
                     localGridSize(2) * (k - 1 - gridOffset(3)))
                gridArray(localIndex,l) = patchArray(patchIndex,l)
             end do
          end do
       end do
    end do

    return
  end subroutine patch_disperse_vector


  ! Disperse a tensor from patch to grid
  ! ------------------------------------
  subroutine patch_disperse_tensor(patch, patchArray, gridArray)

    implicit none

    ! Arguments
    type(t_Patch), intent(in) :: patch
    real(WP), intent(in) :: patchArray(:,:,:)
    real(WP), intent(out) :: gridArray(:,:,:)

    ! Local variables
    integer :: i, j, k, l, m, patchIndex, localIndex

    gridArray = 0.0_WP

    if (patch%nPatchPoints .eq. 0) return

    do m = 1, size(patchArray, 3)
       do l = 1, size(patchArray, 2)
          do k = patch%iStart(3), patch%iEnd(3)
             do j = patch%iStart(2), patch%iEnd(2)
                do i = patch%iStart(1), patch%iEnd(1)
                   patchIndex = i - patch%offset(1) +                                        &
                        patch%localSize(1) * (j - 1 - patch%offset(2) +                      &
                        patch%localSize(2) * (k - 1 - patch%offset(3)))
                   localIndex = i - gridOffset(1) +                                          &
                        localGridSize(1) * (j - 1 - gridOffset(2) +                          &
                        localGridSize(2) * (k - 1 - gridOffset(3)))
                   gridArray(localIndex,l,m) = patchArray(patchIndex,l,m)
                end do
             end do
          end do
       end do
    end do

    return
  end subroutine patch_disperse_tensor


  ! Compute quadrature on a patch
  ! -----------------------------
  subroutine patch_quadrature(patch, integrand, integral)

    ! External modules
    use grid_functions

    implicit none

    ! Arguments
    type(t_Patch), intent(in) :: patch
    real(WP), intent(in) :: integrand(:)
    real(WP), intent(out) :: integral

    ! Local variables
    integer :: i, j, k, gridIndex
    real(WP), allocatable :: mask(:)

    allocate(mask(nGridPoints))
    mask = 0.0_WP

    if (patch%nPatchPoints .gt. 0) then
       do k = patch%iStart(3), patch%iEnd(3)
          do j = patch%iStart(2), patch%iEnd(2)
             do i = patch%iStart(1), patch%iEnd(1)
                gridIndex = i - gridOffset(1) +                                              &
                     localGridSize(1) * (j - 1 - gridOffset(2) +                             &
                     localGridSize(2) * (k - 1 - gridOffset(3)))
                mask(gridIndex) = 1.0_WP
             end do
          end do
       end do
    end if

    where (iblank .eq. 0)
       mask = 0.0_WP
    end where

    integral = inner_product(mask, integrand)

    deallocate(mask)

    return
  end subroutine patch_quadrature


  ! Create hyperbolic tangent support on a patch
  ! --------------------------------------------
  subroutine patch_tanh_support(patch, direction, gridArray, steepness, fraction)

    ! External modules
    use parallel

    implicit none

    ! Arguments
    type(t_Patch), intent(in) :: patch
    integer, intent(in) :: direction
    real(WP), intent(inout) :: gridArray(:)
    real(WP), intent(in), optional :: steepness, fraction

    ! Local variables
    integer :: i, j, k, gridIndex
    real(WP) :: steepness_, fraction_, x, y, z

    if (patch%globalSize(direction) .ge. globalGridSize(direction)) return

    if (present(steepness)) then
       steepness_ = steepness
    else
       steepness_ = 20.0_WP
    end if

    if (present(fraction)) then
       fraction_ = fraction
    else
       fraction_ = 0.1_WP
    end if

    select case (direction)

    case (1)

       do k = iStart(3), iEnd(3)
          do j = iStart(2), iEnd(2)
             do i = iStart(1), iEnd(1)
                gridIndex = i - gridOffset(1) +                                              &
                     localGridSize(1) * (j - 1 - gridOffset(2) +                             &
                     localGridSize(2) * (k - 1 - gridOffset(3)))
                if (i .ge. patch%imin .and. i .le. patch%imax) then
                   x = 2.0_WP * real(i - patch%imin, WP) /                                   &
                        real(patch%imax - patch%imin, WP) - 1.0_WP
                   gridArray(gridIndex) = gridArray(gridIndex) * 0.5_WP * (                  &
                        tanh(steepness_ * (x + 1.0_WP - 0.5_WP * fraction_)) -               &
                        tanh(steepness_ * (x - 1.0_WP + 0.5_WP * fraction_)))
                else
                   gridArray(gridIndex) = 0.0_WP
                end if
             end do
          end do
       end do

    case (2)

       do k = iStart(3), iEnd(3)
          do j = iStart(2), iEnd(2)
             do i = iStart(1), iEnd(1)
                gridIndex = i - gridOffset(1) +                                              &
                     localGridSize(1) * (j - 1 - gridOffset(2) +                             &
                     localGridSize(2) * (k - 1 - gridOffset(3)))
                if (j .ge. patch%jmin .and. j .le. patch%jmax) then
                   y = 2.0_WP * real(j - patch%jmin, WP) /                                   &
                        real(patch%jmax - patch%jmin, WP) - 1.0_WP
                   gridArray(gridIndex) = gridArray(gridIndex) * 0.5_WP * (                  &
                        tanh(steepness_ * (y + 1.0_WP - 0.5_WP * fraction_)) -               &
                        tanh(steepness_ * (y - 1.0_WP + 0.5_WP * fraction_)))
                else
                   gridArray(gridIndex) = 0.0_WP
                end if
             end do
          end do
       end do

    case (3)

       do k = iStart(3), iEnd(3)
          do j = iStart(2), iEnd(2)
             do i = iStart(1), iEnd(1)
                gridIndex = i - gridOffset(1) +                                              &
                     localGridSize(1) * (j - 1 - gridOffset(2) +                             &
                     localGridSize(2) * (k - 1 - gridOffset(3)))
                if (k .ge. patch%kmin .and. k .le. patch%kmax) then
                   z = 2.0_WP * real(k - patch%kmin, WP) /                                   &
                        real(patch%kmax - patch%kmin, WP) - 1.0_WP
                   gridArray(gridIndex) = gridArray(gridIndex) * 0.5_WP * (                  &
                        tanh(steepness_ * (z + 1.0_WP - 0.5_WP * fraction_)) -               &
                        tanh(steepness_ * (z - 1.0_WP + 0.5_WP * fraction_)))
                else
                   gridArray(gridIndex) = 0.0_WP
                end if
             end do
          end do
       end do

    end select

    return
  end subroutine patch_tanh_support


  ! Create compact support on a patch using cubic b-splines
  ! -------------------------------------------------------
  subroutine patch_cubic_bspline_support(patch, direction, gridArray)

    ! External modules
    use parallel
    use math, only : bspline2

    implicit none

    ! Arguments
    type(t_Patch), intent(in) :: patch
    integer, intent(in) :: direction
    real(WP), intent(inout) :: gridArray(:)

    ! Local variables
    integer :: i, j, k, gridIndex
    real(WP) :: x, y, z

    if (patch%globalSize(direction) .ge. globalGridSize(direction)) return

    select case (direction)

    case (1)

       do k = iStart(3), iEnd(3)
          do j = iStart(2), iEnd(2)
             do i = iStart(1), iEnd(1)
                gridIndex = i - gridOffset(1) +                                              &
                     localGridSize(1) * (j - 1 - gridOffset(2) +                             &
                     localGridSize(2) * (k - 1 - gridOffset(3)))
                if (i .ge. patch%imin .and. i .le. patch%imax) then
                   x = real(i - patch%imin, WP) / real(patch%imax - patch%imin, WP)
                   gridArray(gridIndex) = gridArray(gridIndex) * bspline2(x)
                else
                   gridArray(gridIndex) = 0.0_WP
                end if
             end do
          end do
       end do

    case (2)

       do k = iStart(3), iEnd(3)
          do j = iStart(2), iEnd(2)
             do i = iStart(1), iEnd(1)
                gridIndex = i - gridOffset(1) +                                              &
                     localGridSize(1) * (j - 1 - gridOffset(2) +                             &
                     localGridSize(2) * (k - 1 - gridOffset(3)))
                if (j .ge. patch%jmin .and. j .le. patch%jmax) then
                   y = real(j - patch%jmin, WP) / real(patch%jmax - patch%jmin, WP)
                   gridArray(gridIndex) = gridArray(gridIndex) * bspline2(y)
                else
                   gridArray(gridIndex) = 0.0_WP
                end if
             end do
          end do
       end do

    case (3)

       do k = iStart(3), iEnd(3)
          do j = iStart(2), iEnd(2)
             do i = iStart(1), iEnd(1)
                gridIndex = i - gridOffset(1) +                                              &
                     localGridSize(1) * (j - 1 - gridOffset(2) +                             &
                     localGridSize(2) * (k - 1 - gridOffset(3)))
                if (k .ge. patch%kmin .and. k .le. patch%kmax) then
                   z = real(k - patch%kmin, WP) / real(patch%kmax - patch%kmin, WP)
                   gridArray(gridIndex) = gridArray(gridIndex) * bspline2(z)
                else
                   gridArray(gridIndex) = 0.0_WP
                end if
             end do
          end do
       end do
    end select

    return
  end subroutine patch_cubic_bspline_support


  ! Setup an individual patch
  ! -------------------------
  subroutine patch_setup(patch)

    implicit none

    ! Arguments
    type(t_Patch), intent(inout) :: patch

    ! Local variables
    integer :: i, extent(6), offset(3), color, ierror
    logical :: isExtentValid
    character(len = str_long) :: message

    ! Temporarily store the patch extents
    extent = (/ patch%iMin, patch%iMax, patch%jMin, patch%jMax, patch%kMin, patch%kMax /)

    ! Extents must be non-zero
    isExtentValid = all(extent .ne. 0)

    ! Negative values indicate counting backwards from the end
    do i = 1, nDimensions
       if (extent(1+2*(i-1)) .lt. 0) extent(1+2*(i-1)) = extent(1+2*(i-1))                   &
            + globalGridSize(i) + 1
       if (extent(2+2*(i-1)) .lt. 0) extent(2+2*(i-1)) = extent(2+2*(i-1))                   &
            + globalGridSize(i) + 1
    end do

    ! Check that extent describes a part of the grid
    do i = 1, nDimensions
       isExtentValid = isExtentValid .and. (extent(2+2*(i-1)) .ge. extent(1+2*(i-1)))
       isExtentValid = isExtentValid .and. (extent(2+2*(i-1))                                &
            - extent(1+2*(i-1)) + 1 .le. globalGridSize(i))
    end do
    do while (i .le. 3) ! ... reset for direction > number of dimensions
       extent(1+2*(i-1)) = 1
       extent(2+2*(i-1)) = 1
       i = i + 1
    end do

    ! Fail if the extent is invalid
    if (.not. isExtentValid) then
       write(message, '(3A,6(1X,I0.0),A)') "Patch '", trim(patch%name),                      &
            "' has an invalid extent:", patch%iMin, patch%iMax, patch%jMin, patch%jMax,      &
            patch%kMin, patch%kMax, "!"
       call die(trim(message))
    end if

    ! Copy over extent to iMin, iMax, etc
    patch%iMin = extent(1)
    patch%iMax = extent(2)
    patch%jMin = extent(3)
    patch%jMax = extent(4)
    patch%kMin = extent(5)
    patch%kMax = extent(6)

    ! Global patch size
    patch%globalSize(1) = patch%iMax - patch%iMin + 1
    patch%globalSize(2) = patch%jMax - patch%jMin + 1
    patch%globalSize(3) = patch%kMax - patch%kMin + 1

    ! Zero-based index of first point on the patch belonging to the ``current'' process (this
    ! value has no meaning if the patch lies outside the part of the grid belonging to the
    ! ``current'' process)
    patch%offset(1) = max(patch%iMin, gridOffset(1) + 1)
    patch%offset(2) = max(patch%jMin, gridOffset(2) + 1)
    patch%offset(3) = max(patch%kMin, gridOffset(3) + 1)
    patch%offset = min(patch%offset, gridOffset + localGridSize) - 1

    ! Extent of the patch belonging to the ``current'' process (this value has no meaning if
    ! the patch lies outside the part of the grid belonging to the ``current'' process)
    patch%localSize(1) = max(patch%iMax, gridOffset(1) + 1)
    patch%localSize(2) = max(patch%jMax, gridOffset(2) + 1)
    patch%localSize(3) = max(patch%kMax, gridOffset(3) + 1)
    patch%localSize = min(patch%localSize, gridOffset + localGridSize)
    patch%localSize = patch%localSize - patch%offset

    ! Reset size and offset if the patch lies outside the part of the grid belonging to the
    ! ``current'' process
    if (any(patch%localSize .lt. 0 .or.                                                      &
         patch%iMax .lt. gridOffset(1) + 1 .or.                                              &
         patch%iMin .gt. gridOffset(1) + localGridSize(1) .or.                               &
         patch%jMax .lt. gridOffset(2) + 1 .or.                                              &
         patch%jMin .gt. gridOffset(2) + localGridSize(2) .or.                               &
         patch%kMax .lt. gridOffset(3) + 1 .or.                                              &
         patch%kMin .gt. gridOffset(3) + localGridSize(3))) then
       patch%offset = 0
       patch%localSize = 0
    end if

    ! Store the patch partition
    patch%iStart = patch%offset + 1
    patch%iEnd = patch%offset + patch%localSize

    ! Get the total number of local patch points
    patch%nPatchPoints = product(patch%localSize)

    ! Patch data for processes belonging to this patch
    if (patch%nPatchPoints .gt. 0) then
       offset = patch%offset - extent(1::2) + 1
       call MPI_Type_create_subarray(3, patch%globalSize, patch%localSize, offset,              &
            MPI_ORDER_FORTRAN, MPI_REAL_WP, patch%mpiRealSubarrayType, ierror)
       call MPI_Type_commit(patch%mpiRealSubarrayType, ierror)

       if (allocated(gridNorm)) then
          allocate(patch%norm(patch%nPatchPoints, 1))
          call patch_collect(patch, gridNorm, patch%norm)
       end if
    end if

    ! Create a patch communicator
    color = MPI_UNDEFINED
    if (patch%nPatchPoints .gt. 0) color = 1
    call MPI_Comm_split(comm, color, iRank, patch%comm, ierror)

    ! Find the master rank on this patch
    patch%masterRank = huge(1)
    if (patch%nPatchPoints .gt. 0) patch%masterRank = iRank
    call parallel_min(patch%masterRank)

    ! Setup the iblank region
    if (patch%patchType .eq. IBLANK_PATCH) call setup_iblank_patch(patch)

    return
  end subroutine patch_setup


  ! Cleanup an individual patch
  ! ---------------------------
  subroutine patch_cleanup(patch)

    implicit none

    ! Arguments
    type(t_Patch), intent(inout) :: patch

    if (allocated(patch%norm)) deallocate(patch%norm)

    return
  end subroutine patch_cleanup


  ! Setup the iblank patch
  ! ----------------------
  subroutine setup_iblank_patch(patch)

    ! External modules
    use geometry

    implicit none

    ! Arguments
    type(t_Patch), intent(inout) :: patch

    ! Local variables
    integer :: i, j, k, gridIndex, extent(6)

    ! Verify the patch type
    if (patch%patchType .ne. IBLANK_PATCH)                                                   &
         call die('verify_iblank_patch: Patch type mismatch!')

    if (patch%normalDirection .ne. 0)                                                        &
         call die("verify_iblank_patch: Normal direction normal direction /= 0!")

    extent(1) = patch%iMin; extent(2) = patch%iMax
    extent(3) = patch%jMin; extent(4) = patch%jMax
    extent(5) = patch%kMin; extent(6) = patch%kMax

    do i = 1, nDimensions
       if (extent((i-1)*2+1) .lt. 0 .or. extent((i-1)*2+2) .gt. globalGridSize(i) .or.       &
            extent((i-1)*2+1) .gt. extent((i-1)*2+2))                                        &
            call die("verify_iblank_patch: Invalid extent on '" // trim(patch%name) // "'!")
    end do

    ! Turn on iblank
    useIblank = .true.

    ! Set iblank value within the patch
    if (patch%nPatchPoints .gt. 0) then
       do k = patch%iStart(3), patch%iEnd(3)
          do j = patch%iStart(2), patch%iEnd(2)
             do i = patch%iStart(1), patch%iEnd(1)
                gridIndex = i - gridOffset(1) + localGridSize(1) *                              &
                     (j - 1 - gridOffset(2) + localGridSize(2) *                                &
                     (k - 1 - gridOffset(3)))
                iblank(gridIndex) = 0
             end do
          end do
       end do
    end if

    return
  end subroutine setup_iblank_patch

end module grid_patch


! ====================== !
! Setup the grid patches !
! ====================== !
subroutine grid_patch_setup

  ! Internal modules
  use grid_patch

  ! External modules
  use parser
  use fileio

  implicit none

  ! Local variables
  integer :: i, istat
  integer, parameter :: nPatchEntries = 9
  logical :: isDefined
  character(len = str_medium), allocatable :: patchData(:)
  character(len = str_medium) :: message, patchType

  ! Determine if patches are used
  call parser_fieldfortag('patches', i, isDefined)
  if (.not. isDefined) then
     nPatches = 0
     return
  end if

  ! Determine the number of patches
  call parser_getsize('patches', i)
  if (mod(i,nPatchEntries) .ne. 0)                                                           &
       call die('grid_patch_setup: Invalid number of patch entries!')
  nPatches = i / nPatchEntries

  if (nPatches .eq. 1) then
     write(message, "(A,I0.0,A)") 'Found ', nPatches, ' patch'
  else
     write(message, "(A,I0.0,A)") 'Found ', nPatches, ' patches'
  end if
  call monitor_log(trim(message))

  ! Prepare the patch array
  allocate(patches(nPatches))
  allocate(patchData(nPatches * nPatchEntries))
  call parser_read('patches', patchData)
  do i = 1, nPatches
     ! Read the patch data
     read(patchData((i-1) * nPatchEntries + 1 : i * nPatchEntries), *, iostat = istat)       &
          patches(i)%name, patchType, patches(i)%normalDirection,                            &
          patches(i)%iMin, patches(i)%iMax, patches(i)%jMin, patches(i)%jMax,                &
          patches(i)%kMin, patches(i)%kMax

     ! Get the patch type
     select case (trim(patchType))

     case ('SPONGE')
        patches(i)%patchType = SPONGE_PATCH

     case ('SAT_FAR_FIELD')
        patches(i)%patchType = SAT_FAR_FIELD

     case ('SAT_OUTFLOW')
        patches(i)%patchType = SAT_OUTFLOW

     case ('SAT_INFLOW')
        patches(i)%patchType = SAT_INFLOW

     case ('SAT_SLIP_WALL')
        patches(i)%patchType = SAT_SLIP_WALL

     case ('SAT_ISOTHERMAL_WALL')
        patches(i)%patchType = SAT_ISOTHERMAL_WALL

     case ('SAT_ADIABATIC_WALL')
        patches(i)%patchType = SAT_ADIABATIC_WALL

     case ('DIRICHLET')
        patches(i)%patchType = DIRICHLET_BC

     case ('NEUMANN')
        patches(i)%patchType = NEUMANN_BC

     case ('SLIP_WALL')
        patches(i)%patchType = SLIP_BC

     case ('ISOTHERMAL_WALL')
        patches(i)%patchType = ISOTHERMAL_BC

     case ('ACTUATOR')
        patches(i)%patchType = ACTUATOR

     case ('COST_TARGET')
        patches(i)%patchType = COST_TARGET

     case ('VISUALIZATION')
        patches(i)%patchType = VISUALIZATION

     case ('IBLANK')
        patches(i)%patchType = IBLANK_PATCH

     case ('EXCITATION')
        patches(i)%patchType = EXCITATION_PATCH

     case ('STATISTICS')
        patches(i)%patchType = STATISTICS_PATCH

     case default
        call die("grid_patch_setup: Unknown patch type '" //trim(patchType) // "'!")

     end select

     ! Check for duplicate names
     if (any(patches(:i-1)%name .eq. patches(i)%name))                                       &
          call die("grid_patch_setup: A patch with name '" // trim(patches(i)%name) //       &
          "' already exists!")
  end do
  deallocate(patchData)

  ! Reorder the patches
  call patch_sort(patches, 1, nPatches)

  ! Configure the individual patches
  do i = 1, nPatches
     call patch_setup(patches(i))
  end do

  return
end subroutine grid_patch_setup


! ======================== !
! Cleanup the grid patches !
! ======================== !
subroutine grid_patch_cleanup

  ! Internal modules
  use grid_patch

  implicit none

  ! Local variables
  integer :: i

  if (associated(patches)) then
     do i = 1, size(patches)
        call patch_cleanup(patches(i))
     end do
     deallocate(patches)
  end if

  nPatches = 0

  return
end subroutine grid_patch_cleanup
