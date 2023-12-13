module sponge

  ! External modules
  use boundary
  use precision

  implicit none

  integer :: nSponges, defaultSpongeExponent
  real(WP) :: defaultSpongeAmount

  type, private :: t_Sponge
     integer :: spongeExponent
     real(WP) :: spongeAmount
     real(WP), allocatable :: spongeStrength(:)
  end type t_Sponge

  type(t_Sponge), allocatable :: spongeData(:)
  type(t_Patch), pointer :: spongePatch(:)

contains

  subroutine setup_sponge_patch(patch, bc)

    ! External modules
    use parser
    use geometry
    use grid_functions

    implicit none

    ! Arguments
    type(t_Patch), intent(inout) :: patch
    type(t_Sponge), intent(inout) :: bc

    ! Local variables
    integer :: direction, i, j, k
    real(WP), dimension(:,:), allocatable :: coordinateDerivatives, arcLength,               &
         globalArcLengthsAlongDirection
    real(WP), allocatable :: curveLengthIntegrand(:)

    ! Verify the patch type
    call verify_sponge_patch(patch)

    ! Read in the sponge parameters
    call parser_read(trim(patch%name) //  ' sponge amount', bc%spongeAmount,                 &
         defaultspongeAmount)
    call parser_read(trim(patch%name) //  ' sponge exponent', bc%spongeExponent,             &
         defaultspongeExponent)

    ! Get the direction
    direction = abs(patch%normalDirection)

    ! Allocate sponge strength
    if (patch%nPatchPoints .gt. 0) then
       allocate(bc%spongeStrength(patch%nPatchPoints))
       bc%spongeStrength = 0.0_WP
    end if

    ! Compute local arc length
    allocate(arcLength(nGridPoints, 1))
    allocate(coordinateDerivatives(nGridPoints, nDimensions))
    call compute_coordinate_derivative(direction, coordinateDerivatives)
    arcLength(:,1) = sqrt(sum(coordinateDerivatives ** 2, dim = 2))
    deallocate(coordinateDerivatives) !... no longer needed

    ! Gather arc length along normal direction
    allocate(globalArcLengthsAlongDirection(nGridPoints / localGridSize(direction) *         &
         globalGridSize(direction), 1))
    call gather_along_direction(arcLength, localGridSize, direction, gridOffset(direction),  &
         globalArcLengthsAlongDirection)
    deallocate(arcLength) !... no longer needed

    allocate(curveLengthIntegrand(globalGridSize(direction)))

    if (patch%nPatchPoints .gt. 0) then

       select case (direction)

       case (1)

          do k = patch%iStart(3), patch%iEnd(3)
             do j = patch%iStart(2), patch%iEnd(2)

                do i = 1, globalGridSize(1)
                   curveLengthIntegrand(i) = globalArcLengthsAlongDirection(i +              &
                        globalGridSize(1) * (j - 1 - gridOffset(2) +                         &
                        localGridSize(2) * (k - 1 - gridOffset(3))), 1)
                end do

                if (patch%normalDirection .gt. 0) then
                   do i = patch%iStart(1), patch%iEnd(1)
                      bc%spongeStrength(i - patch%offset(1) +                                &
                           patch%localSize(1) * (j - 1 - patch%offset(2) +                   &
                           patch%localSize(2) * (k - 1 - patch%offset(3)))) =                &
                           sum(curveLengthIntegrand(patch%iMin : i - 1)) /                   &
                           sum(curveLengthIntegrand(patch%iMin : patch%iMax - 1))
                   end do
                else
                   do i = patch%iStart(1), patch%iEnd(1)
                      bc%spongeStrength(i - patch%offset(1) +                                &
                           patch%localSize(1) * (j - 1 - patch%offset(2) +                   &
                           patch%localSize(2) * (k - 1 - patch%offset(3)))) =                &
                           sum(curveLengthIntegrand(i + 1 : patch%iMax)) /                   &
                           sum(curveLengthIntegrand(patch%iMin + 1 : patch%iMax))
                   end do
                end if

             end do
          end do

       case (2)

          do k = patch%iStart(3), patch%iEnd(3)
             do i =  patch%iStart(1), patch%iEnd(1)

                do j = 1, globalGridSize(2)
                   curveLengthIntegrand(j) =                                                 &
                        globalArcLengthsAlongDirection(i - gridOffset(1) +                   &
                        localGridSize(1) * (j - 1 +                                          &
                        globalGridSize(2) * (k - 1 - gridOffset(3))), 1)
                end do

                if (patch%normalDirection .gt. 0) then
                   do j = patch%iStart(2), patch%iEnd(2)
                      bc%spongeStrength(i - patch%offset(1) +                                &
                           patch%localSize(1) * (j - 1 - patch%offset(2) +                   &
                           patch%localSize(2) * (k - 1 - patch%offset(3)))) =                &
                           sum(curveLengthIntegrand(patch%jMin : j - 1)) /                   &
                           sum(curveLengthIntegrand(patch%jMin : patch%jMax - 1))
                   end do
                else
                   do j = patch%iStart(2), patch%iEnd(2)
                      bc%spongeStrength(i - patch%offset(1) +                                &
                           patch%localSize(1) * (j - 1 - patch%offset(2) +                   &
                           patch%localSize(2) * (k - 1 - patch%offset(3)))) =                &
                           sum(curveLengthIntegrand(j + 1 : patch%jMax)) /                   &
                           sum(curveLengthIntegrand(patch%jMin + 1 : patch%jMax))
                   end do
                end if

             end do
          end do

       case (3)

          do j = patch%iStart(2), patch%iEnd(2)
             do i = patch%iStart(1), patch%iEnd(1)

                do k = 1, globalGridSize(3)
                   curveLengthIntegrand(k) =                                                 &
                        globalArcLengthsAlongDirection(i - gridOffset(1) +                   &
                        localGridSize(1) * (j - 1 - gridOffset(2) +                          &
                        localGridSize(2) * (k - 1)), 1)
                end do

                if (patch%normalDirection .gt. 0) then
                   do k = patch%iStart(3), patch%iEnd(3)
                      bc%spongeStrength(i - patch%offset(1) +                                &
                           patch%localSize(1) * (j - 1 - patch%offset(2) +                   &
                           patch%localSize(2) * (k - 1 - patch%offset(3)))) =                &
                           sum(curveLengthIntegrand(patch%kMin : k - 1)) /                   &
                           sum(curveLengthIntegrand(patch%kMin : patch%kMax - 1))
                   end do
                else
                   do k = patch%iStart(3), patch%iEnd(3)
                      bc%spongeStrength(i - patch%offset(1) +                                &
                           patch%localSize(1) * (j - 1 - patch%offset(2) +                   &
                           patch%localSize(2) * (k - 1 - patch%offset(3)))) =                &
                           sum(curveLengthIntegrand(k + 1 : patch%kMax)) /                   &
                           sum(curveLengthIntegrand(patch%kMin + 1 : patch%kMax))
                   end do
                end if

             end do
          end do

       end select

       bc%spongeStrength = bc%spongeAmount *                                                 &
            (1.0_WP - bc%spongeStrength) ** real(bc%spongeExponent, WP)

    end if

    deallocate(curveLengthIntegrand)
    deallocate(globalArcLengthsAlongDirection)

    return
  end subroutine setup_sponge_patch


  subroutine cleanup_sponge_patch(patch, bc)

    implicit none

    ! Arguments
    type(t_Patch), intent(inout) :: patch
    type(t_Sponge), intent(inout) :: bc

    if (allocated(bc%spongeStrength)) deallocate(bc%spongeStrength)

    return
  end subroutine cleanup_sponge_patch


  subroutine verify_sponge_patch(patch)

    ! External modules
    use simulation_flags
    use geometry

    implicit none

    ! Arguments
    type(t_Patch), intent(in) :: patch

    ! Local variables
    integer :: i, n, extent(6)
    character(len = str_long) :: message

    if (patch%patchType .ne. SPONGE_PATCH)                                                   &
         call die('verify_sponge_patch: Patch type mismatch!')

    if (patch%normalDirection .gt. nDimensions .or. patch%normalDirection .eq. 0)            &
         call die("verify_sponge_patch: Normal direction is invalid for '" //                &
         trim(patch%name) // "'!")

    if (.not. useTargetState)                                                                &
         call die('verify_sponge_patch: No target state available for sponge damping!')

    extent(1) = patch%iMin; extent(2) = patch%iMax
    extent(3) = patch%jMin; extent(4) = patch%jMax
    extent(5) = patch%kMin; extent(6) = patch%kMax

    n = nDimensions
    do i = 1, nDimensions
       if (extent((i-1)*2+1) .lt. 0 .or. extent((i-1)*2+2) .gt. globalGridSize(i) .or.       &
            extent((i-1)*2+1) .gt. extent((i-1)*2+2))                                        &
            call die("verify_sponge_patch: Invalid extent on '" // trim(patch%name) // "'!")
       if (extent((i-1)*2+1) .eq. extent((i-1)*2+2)) n = n - 1
    end do


    if (n .ne. nDimensions) then
       write(message, '(2(A,I0.0),3A)') "verify_sponge_patch: Expected a ", nDimensions,     &
            "D patch, but extent represents a ", n, "D patch on '", trim(patch%name), "'!"
       call die(trim(message))
    end if

    i = abs(patch%normalDirection)
    if ((patch%normalDirection .gt. 0 .and. extent((i-1)*2+1) .ne. 1) .or.                   &
         (patch%normalDirection .lt. 0 .and. extent((i-1)*2+2) .ne. globalGridSize(i)))      &
         call die("verify_sponge_patch: '" // trim(patch%name) //                            &
         "' not aligned wiht a computational boundary!")

    return
  end subroutine verify_sponge_patch


  subroutine sponge_patch_forward(patch, bc, source)

    ! External modules
    use solver_options
    use geometry
    use grid
    use state

    implicit none

    ! Arguments
    type(t_Patch), intent(inout) :: patch
    type(t_Sponge), intent(inout) :: bc
    real(WP), dimension(:,:), intent(inout) :: source

    ! Local variables
    integer :: i, j, k, l, gridIndex, patchIndex
    real(WP) :: localConservedVariable, vfCorrection

    do l = 1, nUnknowns
       do k = patch%iStart(3), patch%iEnd(3)
          do j = patch%iStart(2), patch%iEnd(2)
             do i = patch%iStart(1), patch%iEnd(1)
                gridIndex = i - gridOffset(1) + localGridSize(1) *                           &
                     (j - 1 - gridOffset(2) + localGridSize(2) *                             &
                     (k - 1 - gridOffset(3)))
                patchIndex = i - patch%offset(1) + patch%localSize(1) *                      &
                     (j - 1 - patch%offset(2) + patch%localSize(2) *                         &
                     (k - 1 - patch%offset(3)))

                if (twoWayCoupling) then
                   localConservedVariable = conservedVariables(gridIndex,l) /                &
                        volumeFraction(gridIndex,1)
                   vfCorrection = volumeFraction(gridIndex,1)
                else
                   localConservedVariable = conservedVariables(gridIndex,l)
                   vfCorrection = 1.0_WP
                end if

                source(gridIndex, l) = source(gridIndex, l) - bc%spongeStrength(patchIndex) *&
                     (localConservedVariable - targetState(gridIndex, l)) * vfCorrection

             end do
          end do
       end do
    end do

    return
  end subroutine sponge_patch_forward


  subroutine sponge_patch_adjoint(patch, bc, source)

    ! External modules
    use simulation_flags
    use solver_options
    use geometry
    use grid
    use state

    implicit none

    ! Arguments
    type(t_Patch), intent(inout) :: patch
    type(t_Sponge), intent(inout) :: bc
    real(WP), dimension(:,:), intent(inout) :: source

    ! Local variables
    integer :: i, j, k, l, gridIndex, patchIndex

    do l = 1, nUnknowns
       do k = patch%iStart(3), patch%iEnd(3)
          do j = patch%iStart(2), patch%iEnd(2)
             do i = patch%iStart(1), patch%iEnd(1)
                gridIndex = i - gridOffset(1) + localGridSize(1) *                                &
                     (j - 1 - gridOffset(2) + localGridSize(2) *                                  &
                     (k - 1 - gridOffset(3)))
                patchIndex = i - patch%offset(1) + patch%localSize(1) *                           &
                     (j - 1 - patch%offset(2) + patch%localSize(2) *                              &
                     (k - 1 - patch%offset(3)))

                source(gridIndex, l) = source(gridIndex, l) + bc%spongeStrength(patchIndex) *  &
                     adjointVariables(gridIndex, l)

             end do
          end do
       end do
    end do

    return
  end subroutine sponge_patch_adjoint

end module sponge


! ======================== !
! Setup the sponge patches !
! ======================== !
subroutine sponge_setup

  ! Internal modules
  use sponge

  ! External modules
  use parser

  implicit none

  ! Local variables
  integer :: i, j

  ! Find the number of boundaries of this type
  nSponges = 0; j = 0
  do i = 1, nPatches
     if (patches(i)%patchType .eq. SPONGE_PATCH) then
        nSponges = nSponges + 1
        if (j .eq. 0) j = i
     end if
  end do
  if (nSponges .eq. 0) return

  ! Allocate the sponge type
  allocate(spongeData(nSponges))

  ! Read in default parameters
  call parser_read('default sponge amount', defaultSpongeAmount, 1.0_WP)
  call parser_read('default sponge exponent', defaultSpongeExponent, 2)

  ! Connect the boundary patch
  spongePatch => patches(j:j+nSponges-1)

  ! Setup the boundary conditions
  do i = 1, nSponges
     call setup_sponge_patch(spongePatch(i), spongeData(i))
  end do

  return
end subroutine sponge_setup


! ========================== !
! Cleanup the sponge patches !
! ========================== !
subroutine sponge_cleanup

  ! Internal modules
  use sponge

  implicit none

  ! Local variables
  integer :: i

  if (nSponges .gt. 0) then
     do i = 1, nSponges
        call cleanup_sponge_patch(spongePatch(i), spongeData(i))
     end do
     deallocate(spongeData)
     nullify(spongePatch)
  end if

  nSponges = 0

  return
end subroutine sponge_cleanup


! ========================================= !
! Add sponge sources during the forward run !
! ========================================= !
subroutine sponge_forward(source)

  ! Internal modules
  use sponge

  ! External modules
  use solver_options
  use geometry

  implicit none

  ! Arguments
  real(WP), dimension(nGridPoints, nUnknowns), intent(inout) :: source

  ! Local variables
  integer :: i

  do i = 1, nSponges
     call sponge_patch_forward(spongePatch(i), spongeData(i), source)
  end do

  return
end subroutine sponge_forward


! ========================================= !
! Add sponge sources during the adjoint run !
! ========================================= !
subroutine sponge_adjoint(source)

  ! Internal modules
  use sponge

  ! External modules
  use solver_options
  use geometry

  implicit none

  ! Arguments
  real(WP), dimension(nGridPoints, nUnknowns), intent(inout) :: source

  ! Local variables
  integer :: i

  do i = 1, nSponges
     call sponge_patch_adjoint(spongePatch(i), spongeData(i), source)
  end do

  return
end subroutine sponge_adjoint
