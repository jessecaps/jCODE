module boundary_strong

  ! External modules
  use boundary
  use precision

  implicit none

  integer :: nBoundaries

  type, private :: t_Boundary
     real(WP), allocatable :: temperature(:)
  end type t_Boundary

  type(t_Boundary), allocatable :: boundaryData(:)
  type(t_Patch), pointer :: boundaryPatch(:)

contains

  ! Dirichlet boundary treatment based on targetState
  ! -------------------------------------------------
  subroutine dirichlet_patch_inject(patch, stateVector)

    ! External modules
    use state_jacobian
    use solver_options
    use geometry
    use grid
    use grid_levelset
    use state

    implicit none

    ! Arguments
    type(t_Patch), intent(inout) :: patch
    real(WP), dimension(:,:), intent(inout) :: stateVector

    ! Local variables
    integer :: i, j, k, gridIndex

    do k = patch%iStart(3), patch%iEnd(3)
       do j = patch%iStart(2), patch%iEnd(2)
          do i = patch%iStart(1), patch%iEnd(1)
             gridIndex = i - gridOffset(1) + localGridSize(1) *                              &
                  (j - 1 - gridOffset(2) + localGridSize(2) *                                &
                  (k - 1 - gridOffset(3)))

             stateVector(gridIndex,:) = targetState(gridIndex,:)

          end do
       end do
    end do

    return
  end subroutine dirichlet_patch_inject


  ! Neumann boundary treatment using neighboring grid point
  ! -------------------------------------------------------
   subroutine neumann_patch_inject(patch, stateVector)

    ! External modules
    use state_jacobian
    use solver_options
    use geometry
    use grid
    use grid_levelset
    use state

    implicit none

    ! Arguments
    type(t_Patch), intent(inout) :: patch
    real(WP), dimension(:,:), intent(inout) :: stateVector

    ! Local variables
    integer :: i, j, k, ijk(3), direction, gridIndex, gridIndex2

    direction = abs(patch%normalDirection)

    do k = patch%iStart(3), patch%iEnd(3)
       do j = patch%iStart(2), patch%iEnd(2)
          do i = patch%iStart(1), patch%iEnd(1)
             gridIndex = i - gridOffset(1) + localGridSize(1) *                              &
                  (j - 1 - gridOffset(2) + localGridSize(2) *                                &
                  (k - 1 - gridOffset(3)))

             ! Get the index of the neighboring grid point
             ijk(1) = i; ijk(2) = j; ijk(3) = k
             ijk(direction) = ijk(direction) + sign(1, patch%normalDirection)
             gridIndex2 = ijk(1) - gridOffset(1) + localGridSize(1) *                        &
                  (ijk(2) - 1 - gridOffset(2) + localGridSize(2) *                           &
                  (ijk(3) - 1 - gridOffset(3)))

             stateVector(gridIndex,:) = stateVector(gridIndex2,:)

          end do
       end do
    end do

    return
  end subroutine neumann_patch_inject


  ! No-penetration slip wall treatment
  ! ----------------------------------
  subroutine slip_patch_inject(patch, stateVector)

    ! External modules
    use state_jacobian
    use solver_options
    use geometry
    use grid
    use grid_levelset
    use state

    implicit none

    ! Arguments
    type(t_Patch), intent(inout) :: patch
    real(WP), dimension(:,:), intent(inout) :: stateVector

    ! Local variables
    integer :: i, j, k, ijk(3), direction, gridIndex, gridIndex2

    direction = abs(patch%normalDirection)

    do k = patch%iStart(3), patch%iEnd(3)
       do j = patch%iStart(2), patch%iEnd(2)
          do i = patch%iStart(1), patch%iEnd(1)
             gridIndex = i - gridOffset(1) + localGridSize(1) *                              &
                  (j - 1 - gridOffset(2) + localGridSize(2) *                                &
                  (k - 1 - gridOffset(3)))

             ! Get the index of the neighboring grid point
             ijk(1) = i; ijk(2) = j; ijk(3) = k
             ijk(direction) = ijk(direction) + sign(1, patch%normalDirection)
             gridIndex2 = ijk(1) - gridOffset(1) + localGridSize(1) *                        &
                  (ijk(2) - 1 - gridOffset(2) + localGridSize(2) *                           &
                  (ijk(3) - 1 - gridOffset(3)))

             ! Neumann condition
             stateVector(gridIndex,:) = stateVector(gridIndex2,:)

             ! No-penetration
             stateVector(gridIndex,direction+1) = 0.0_WP

          end do
       end do
    end do

    return
  end subroutine slip_patch_inject


  ! No-slip isothermal wall treatment
  ! ---------------------------------
  subroutine isothermal_patch_inject(patch, bc, stateVector)

    ! External modules
    use state_jacobian
    use solver_options
    use geometry
    use grid
    use grid_levelset
    use state

    implicit none

    ! Arguments
    type(t_Patch), intent(inout) :: patch
    type(t_Boundary), intent(inout) :: bc
    real(WP), dimension(:,:), intent(inout) :: stateVector
    
    ! Local variables
    integer :: i, j, k, direction, gridIndex, patchIndex

    direction = abs(patch%normalDirection)

    do k = patch%iStart(3), patch%iEnd(3)
       do j = patch%iStart(2), patch%iEnd(2)
          do i = patch%iStart(1), patch%iEnd(1)
             gridIndex = i - gridOffset(1) + localGridSize(1) *                              &
                  (j - 1 - gridOffset(2) + localGridSize(2) *                                &
                  (k - 1 - gridOffset(3)))
             patchIndex = i - patch%offset(1) + patch%localSize(1) *                         &
                  (j - 1 - patch%offset(2) + patch%localSize(2) *                            &
                  (k - 1 - patch%offset(3)))

             ! No-slip treatment
             stateVector(gridIndex,2:nDimensions+1) = 0.0_WP

             ! Isothermal treatment
             if (equationOfState .eq. IDEAL_GAS_MIXTURE) then
                stateVector(gridIndex,nDimensions+2) = stateVector(gridIndex,1) *            &
                     bc%temperature(patchIndex) / ratioOfSpecificHeats /                     &
                     mixtureMolecularWeight(gridIndex, 1)
             else
                stateVector(gridIndex,nDimensions+2) = stateVector(gridIndex,1) *            &
                     bc%temperature(patchIndex) / ratioOfSpecificHeats
             end if

          end do
       end do
    end do

    return
  end subroutine isothermal_patch_inject


  ! Setup the individual boundary patch
  ! -----------------------------------
  subroutine setup_boundary_patch(patch, bc)

    ! External modules
    use parser
    use simulation_flags
    use solver_options
    use first_derivative
    use state
    use equation_of_state
    use impenetrable, only : setup_impenetrable_patch

    implicit none

    ! Arguments
    type(t_Patch), intent(inout) :: patch
    type(t_Boundary), intent(inout) :: bc

    ! Local variables
    real(WP) :: wallTemperature
    real(WP), allocatable :: targetTemperature(:)

    ! Verify the patch type
    call verify_boundary_patch(patch)

    ! Store temperature if this is a no-slip patch
    if (patch%patchType.eq.ISOTHERMAL_BC .and. patch%nPatchPoints.gt.0 .and. useViscosity) then

       allocate(bc%temperature(patch%nPatchPoints))

       ! Assign the wall temperature & mass fraction
       if (useTargetState) then ! ... get values from target state
          allocate(targetTemperature(nGridPoints))
          call compute_dependent_variables(targetState, temperature = targetTemperature)
          call patch_collect(patch, targetTemperature, bc%temperature)
          deallocate(targetTemperature)

       else ! ... read in values from the input file
          call parser_read(trim(patch%name) // ' temperature', wallTemperature,              &
               1.0_WP / (ratioOfSpecificHeats - 1.0_WP))
          bc%temperature(:) = wallTemperature
       end if
    end if

    ! Dirichlet requires a target state
    if (patch%patchType.eq.DIRICHLET_BC .and. .not.useTargetState)                           &
         call die('verify_boundary_patch: No target state available for Dirichlet patch!')

    return
  end subroutine setup_boundary_patch


  ! Cleanup the individual boundary partch
  ! --------------------------------------
  subroutine cleanup_boundary_patch(patch, bc)

    implicit none

    ! Arguments
    type(t_Patch), intent(inout) :: patch
    type(t_Boundary), intent(inout) :: bc

    if (allocated(bc%temperature)) deallocate(bc%temperature)

    return
  end subroutine cleanup_boundary_patch
  

  ! Verify the individual boundary patch is defined correctly
  ! ---------------------------------------------------------
  subroutine verify_boundary_patch(patch)

    ! External modules
    use simulation_flags
    use geometry

    implicit none

    ! Arguments
    type(t_Patch), intent(in) :: patch

    ! Local variables
    integer :: i, extent(6)

    if (patch%normalDirection .gt. nDimensions .or. patch%normalDirection .eq. 0)            &
         call die("verify_dirichlet_patch: Normal direction is invalid for '" //             &
         trim(patch%name) // "'!")

    extent(1) = patch%iMin; extent(2) = patch%iMax
    extent(3) = patch%jMin; extent(4) = patch%jMax
    extent(5) = patch%kMin; extent(6) = patch%kMax

    do i = 1, nDimensions
       if (extent((i-1)*2+1) .lt. 0 .or. extent((i-1)*2+2) .gt. globalGridSize(i) .or.       &
            extent((i-1)*2+1) .gt. extent((i-1)*2+2))                                        &
            call die("verify_dirichlet_patch: Invalid extent on '" // trim(patch%name) //    &
            "'!")
    end do

    i = abs(patch%normalDirection)
    if (extent((i-1)*2+1) .ne. extent((i-1)*2+2)) call die("verify_dirichlet_patch: '" //    &
         trim(patch%name) //  "' extends more than 1 grid point along normal direction!")

    return
  end subroutine verify_boundary_patch

end module boundary_strong


! ================================= !
! Setup the strong boundary patches !
! ================================= !
subroutine boundary_strong_setup

  ! Internal modules
  use boundary_strong

  ! External modules
  use simulation_flags
  use solver_options
  use geometry
  use state_functions
  use state

  implicit none

  ! Local variables
  integer :: i, j

  ! Find the number of boundaries of this type
  nBoundaries = 0; j = 0
  do i = 1, nPatches
     if (patches(i)%patchType .eq. DIRICHLET_BC .or.                                         &
          patches(i)%patchType .eq. NEUMANN_BC .or.                                          &
          patches(i)%patchType .eq. SLIP_BC .or.                                             &
          patches(i)%patchType .eq. ISOTHERMAL_BC) then
        nBoundaries = nBoundaries + 1
        if (j .eq. 0) j = i
     end if
  end do
  if (nBoundaries .eq. 0) return

  ! Allocate the boundary type
  allocate(boundaryData(nBoundaries))

  ! Connect the boundary patch
  boundaryPatch => patches(j:j+nBoundaries-1)

  ! Setup the boundary conditions
  do i = 1, nBoundaries
     call setup_boundary_patch(boundaryPatch(i), boundaryData(i))
  end do

  return
end subroutine boundary_strong_setup


! =================================== !
! Cleanup the strong boundary patches !
! =================================== !
subroutine boundary_strong_cleanup

  ! Internal modules
  use boundary_strong

  implicit none

  ! Local variables
  integer :: i

  if (nBoundaries .gt. 0) then
     do i = 1, nBoundaries
        call cleanup_boundary_patch(boundaryPatch(i), boundaryData(i))
        nullify(boundaryPatch)
     end do
     deallocate(boundaryData)
  end if
 

  nBoundaries = 0

  return
end subroutine boundary_strong_cleanup


! ================================================= !
! Overwrite the state vector at the boundary points !
! Strongly-imposed (non-SAT) boundary conditions    !
! ================================================= !
subroutine boundary_correct_state(stateVector)

  ! Internal modules
  use boundary_strong

  ! External modules
  use solver_options
  use geometry

  implicit none

  ! Arguments
  real(WP), dimension(nGridPoints,nUnknowns), intent(inout) :: stateVector

  ! Local variables
  integer :: i

  do i = 1, nBoundaries
     
     select case (boundaryPatch(i)%patchType)
        
     case (DIRICHLET_BC)
        call dirichlet_patch_inject(boundaryPatch(i), stateVector)
        
     case (NEUMANN_BC)
        call neumann_patch_inject(boundaryPatch(i), stateVector)
        
     case (SLIP_BC)
        call slip_patch_inject(boundaryPatch(i), stateVector)

     case (ISOTHERMAL_BC)
        call isothermal_patch_inject(boundaryPatch(i), boundaryData(i), stateVector)
        
     end select
     
  end do

  return
end subroutine boundary_correct_state
