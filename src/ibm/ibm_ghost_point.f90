module ibm_ghost_point

  ! External modules
  use precision

  implicit none

  ! Global arrays
  integer :: nOverlapGhost(3,2), nGhostPoints, nInteriorPoints, nLayers
  integer, allocatable :: interiorIndex(:)
  real(WP) :: ibmPenalty
  real(WP), dimension(:,:,:,:), allocatable :: ghostTemperature, ghostPressure,              &
       ghostVelocity, ghostMassFraction, ghostLevelset, ghostBuffer
  real(WP), dimension(:,:), allocatable :: cartesianCoordinates
  logical :: filterIBM, useWeakForm

  ! Ghost point type
  type :: t_GhostPoint
     integer :: gridIndex
     integer :: ghostIndex(3)
     integer :: imageIndex(3)
     real(WP) :: imagePoint(3)
  end type t_GhostPoint
  type(t_GhostPoint), dimension(:), allocatable :: ghostPoint

contains

  ! Interpolation routine used for ghost point method
  ! -------------------------------------------------
  subroutine interpolate_image_point(ghostIndex, imageIndex, imagePoint, T_IP, P_IP, U_IP,   &
       Y_IP, flag)

    ! External modules
    use geometry
    use solver_options
    use grid_functions

    implicit none

    ! Arguments
    integer, intent(in) :: ghostIndex(3), imageIndex(3)
    real(WP), intent(in) :: imagePoint(3)
    real(WP), intent(out) :: T_IP, P_IP, U_IP(3), Y_IP(nSpecies)
    logical, intent(out) :: flag

    ! Local variables
    integer :: i, j, k, i1, j1, k1, i2, j2, k2, n
    real(WP) :: interpCoeff(3,3,3), alpha(3,3,3), dist(2,2,2), eta(2,2,2), buf
    real(WP), parameter :: eps = 1.0e-10_WP

    select case (nDimensions)

    case (2)

       ! Get the interpolation points
       i1 = imageIndex(1); i2 = i1 + 1
       j1 = imageIndex(2); j2 = j1 + 1

       ! Correct for walls
       if (nOverlapGhost(1,1) .eq. 0) i1 = max(i1, iStart(1))
       if (nOverlapGhost(1,2) .eq. 0) i2 = min(i2, iEnd(1))
       if (nOverlapGhost(2,1) .eq. 0) j1 = max(j1, iStart(2))
       if (nOverlapGhost(2,2) .eq. 0) j2 = min(j2, iEnd(2))
       

!!$       ! Compute the interpolation coefficients
!!$       call lagrange_basis(2, cartesianCoordinates(i1:i2,1), imagePoint(1), weight(1:2,1))
!!$       call lagrange_basis(2, cartesianCoordinates(j1:j2,2), imagePoint(2), weight(1:2,2))
!!$
!!$       ! Combine the interpolation coefficients
!!$       interpCoeff = 0.0_WP
!!$       interpCoeff(1,1,1) = weight(1,1) * weight(1,2)
!!$       interpCoeff(2,1,1) = weight(2,1) * weight(1,2)
!!$       interpCoeff(1,2,1) = weight(1,1) * weight(2,2)
!!$       interpCoeff(2,2,1) = weight(2,1) * weight(2,2)
!!$
!!$       alpha = 1.0_WP
!!$       if (ghostLevelset(i1,j1,1,1).lt.0.0_WP) alpha(1,1,1) = 0.0_WP
!!$       if (ghostLevelset(i2,j1,1,1).lt.0.0_WP) alpha(2,1,1) = 0.0_WP
!!$       if (ghostLevelset(i1,j2,1,1).lt.0.0_WP) alpha(1,2,1) = 0.0_WP
!!$       if (ghostLevelset(i2,j2,1,1).lt.0.0_WP) alpha(2,2,1) = 0.0_WP
!!$       buf = sum(alpha(:,:,1) * interpCoeff(:,:,1))
!!$       if (buf .gt. 0.0_WP) then
!!$          interpCoeff(:,:,1) = alpha(:,:,1) * interpCoeff(:,:,1) / buf
!!$          flag = .false.

       ! Check if the image point is inside the levelset
       alpha = 1.0_WP
       if (ghostLevelset(i1,j1,1,1).lt.0.0_WP) alpha(1,1,1) = 0.0_WP
       if (ghostLevelset(i2,j1,1,1).lt.0.0_WP) alpha(2,1,1) = 0.0_WP
       if (ghostLevelset(i1,j2,1,1).lt.0.0_WP) alpha(1,2,1) = 0.0_WP
       if (ghostLevelset(i2,j2,1,1).lt.0.0_WP) alpha(2,2,1) = 0.0_WP
       buf = sum(alpha(1:2,1:2,1))

       if (buf .gt. 0.0_WP) then
          ! Get interpolation weights at image point (Chaudhuri et al. 2011, JCP)
          flag = .false.
          dist = 0.0_WP
          dist(1,1,1) = sqrt(                                                                   &
               (cartesianCoordinates(i1, 1) - imagePoint(1))**2 +                               &
               (cartesianCoordinates(j1, 2) - imagePoint(2))**2 )
          dist(2,1,1) = sqrt(                                                                   &
               (cartesianCoordinates(i2, 1) - imagePoint(1))**2 +                               &
               (cartesianCoordinates(j1, 2) - imagePoint(2))**2 )
          dist(1,2,1) = sqrt(                                                                   &
               (cartesianCoordinates(i1, 1) - imagePoint(1))**2 +                               &
               (cartesianCoordinates(j2, 2) - imagePoint(2))**2 )
          dist(2,2,1) = sqrt(                                                                   &
               (cartesianCoordinates(i2, 1) - imagePoint(1))**2 +                               &
               (cartesianCoordinates(j2, 2) - imagePoint(2))**2 )
          interpCoeff = 0.0_WP
          if (dist(1,1,1) .le. eps * minGridSpacing) then
             interpCoeff(1,1,1) = 1.0_WP
          else if (dist(2,1,1) .le. eps * minGridSpacing) then
             interpCoeff(2,1,1) = 1.0_WP
          else if (dist(1,2,1) .le. eps * minGridSpacing) then
             interpCoeff(1,2,1) = 1.0_WP
          else if (dist(2,2,1) .le. eps * minGridSpacing) then
             interpCoeff(2,2,1) = 1.0_WP
          else
             eta(:,:,1) = 1.0_WP / dist(:,:,1)**2
             buf = sum(eta(:,:,1) * alpha(1:2,1:2,1))
             interpCoeff(1:2,1:2,1) = alpha(1:2,1:2,1) * eta(:,:,1) / buf
          end if
          
       else

          ! Image point resides within the levelset, average the neighbors instead
          flag = .true.
          i1 = ghostIndex(1) - 1; i2 = i1 + 2
          j1 = ghostIndex(2) - 1; j2 = j1 + 2

          if (nOverlapGhost(1,1) .eq. 0) i1 = max(i1, iStart(1))
          if (nOverlapGhost(1,2) .eq. 0) i2 = min(i2, iEnd(1))
          if (nOverlapGhost(2,1) .eq. 0) j1 = max(j1, iStart(2))
          if (nOverlapGhost(2,2) .eq. 0) j2 = min(j2, iEnd(2))

          interpCoeff = 1.0_WP
          interpCoeff(2,2,1) = 0.0_WP
          alpha = 1.0_WP
          do j = j1, j2
             do i = i1, i2
                ! Remove points inside fluid
                if (ghostLevelset(i,j,1,1) .gt. 0.0_WP) alpha(i-i1+1,j-j1+1,1) = 0.0_WP
             end do
          end do
          buf = sum(alpha(:,:,1) * interpCoeff(:,:,1))
          interpCoeff(:,:,1) = alpha(:,:,1) * interpCoeff(:,:,1) / buf
       end if

       ! Interpolate fluid quantities
       T_IP = sum(interpCoeff(1:i2-i1+1,1:j2-j1+1,1) * ghostTemperature(i1:i2, j1:j2, 1, 1))
       P_IP = sum(interpCoeff(1:i2-i1+1,1:j2-j1+1,1) * ghostPressure(i1:i2, j1:j2, 1, 1))
       U_IP = 0.0_WP
       do n = 1, nDimensions
          U_IP(n) = sum(interpCoeff(1:i2-i1+1,1:j2-j1+1,1) *                                 &
               ghostVelocity(i1:i2, j1:j2, 1, n))
       end do
       do n = 1, nSpecies
          Y_IP(n) = sum(interpCoeff(1:i2-i1+1,1:j2-j1+1,1) *                                 &
               ghostMassFraction(i1:i2, j1:j2, 1, n))
       end do

    case (3)

       ! Get the interpolation points
       i1 = imageIndex(1); i2 = i1 + 1
       j1 = imageIndex(2); j2 = j1 + 1
       k1 = imageIndex(3); k2 = k1 + 1

       ! Correct for walls
       if (nOverlapGhost(1,1) .eq. 0) i1 = max(i1, iStart(1))
       if (nOverlapGhost(1,2) .eq. 0) i2 = min(i2, iEnd(1))
       if (nOverlapGhost(2,1) .eq. 0) j1 = max(j1, iStart(2))
       if (nOverlapGhost(2,2) .eq. 0) j2 = min(j2, iEnd(2))
       if (nOverlapGhost(3,1) .eq. 0) k1 = max(k1, iStart(3))
       if (nOverlapGhost(3,2) .eq. 0) k2 = min(k2, iEnd(3))

       ! Check if the image point is inside the levelset
       alpha = 1.0_WP
       if (ghostLevelset(i1,j1,k1,1).lt.0.0_WP) alpha(1,1,1) = 0.0_WP
       if (ghostLevelset(i2,j1,k1,1).lt.0.0_WP) alpha(2,1,1) = 0.0_WP
       if (ghostLevelset(i1,j2,k1,1).lt.0.0_WP) alpha(1,2,1) = 0.0_WP
       if (ghostLevelset(i2,j2,k1,1).lt.0.0_WP) alpha(2,2,1) = 0.0_WP
       if (ghostLevelset(i1,j1,k2,1).lt.0.0_WP) alpha(1,1,2) = 0.0_WP
       if (ghostLevelset(i2,j1,k2,1).lt.0.0_WP) alpha(2,1,2) = 0.0_WP
       if (ghostLevelset(i1,j2,k2,1).lt.0.0_WP) alpha(1,2,2) = 0.0_WP
       if (ghostLevelset(i2,j2,k2,1).lt.0.0_WP) alpha(2,2,2) = 0.0_WP
       buf = sum(alpha(1:2,1:2,1:2))

       if (buf .gt. 0.0_WP) then
          ! Get interpolation weights at image point (Chaudhuri et al. 2011, JCP)
          flag = .false.
          dist = 0.0_WP
          dist(1,1,1) = sqrt(                                                                &
               (cartesianCoordinates(i1, 1) - imagePoint(1))**2 +                            &
               (cartesianCoordinates(j1, 2) - imagePoint(2))**2 +                            &
               (cartesianCoordinates(k1, 3) - imagePoint(3))**2 )
          dist(2,1,1) = sqrt(                                                                &
               (cartesianCoordinates(i2, 1) - imagePoint(1))**2 +                            &
               (cartesianCoordinates(j1, 2) - imagePoint(2))**2 +                            &
               (cartesianCoordinates(k1, 3) - imagePoint(3))**2 )
          dist(1,2,1) = sqrt(                                                                &
               (cartesianCoordinates(i1, 1) - imagePoint(1))**2 +                            &
               (cartesianCoordinates(j2, 2) - imagePoint(2))**2 +                            &
               (cartesianCoordinates(k1, 3) - imagePoint(3))**2 )
          dist(2,2,1) = sqrt(                                                                &
               (cartesianCoordinates(i2, 1) - imagePoint(1))**2 +                            &
               (cartesianCoordinates(j2, 2) - imagePoint(2))**2 +                            &
               (cartesianCoordinates(k1, 3) - imagePoint(3))**2 )
          dist(1,1,2) = sqrt(                                                                &
               (cartesianCoordinates(i1, 1) - imagePoint(1))**2 +                            &
               (cartesianCoordinates(j1, 2) - imagePoint(2))**2 +                            &
               (cartesianCoordinates(k2, 3) - imagePoint(3))**2 )
          dist(2,1,2) = sqrt(                                                                &
               (cartesianCoordinates(i2, 1) - imagePoint(1))**2 +                            &
               (cartesianCoordinates(j1, 2) - imagePoint(2))**2 +                            &
               (cartesianCoordinates(k2, 3) - imagePoint(3))**2 )
          dist(1,2,2) = sqrt(                                                                &
               (cartesianCoordinates(i1, 1) - imagePoint(1))**2 +                            &
               (cartesianCoordinates(j2, 2) - imagePoint(2))**2 +                            &
               (cartesianCoordinates(k2, 3) - imagePoint(3))**2 )
          dist(2,2,2) = sqrt(                                                                &
               (cartesianCoordinates(i2, 1) - imagePoint(1))**2 +                            &
               (cartesianCoordinates(j2, 2) - imagePoint(2))**2 +                            &
               (cartesianCoordinates(k2, 3) - imagePoint(3))**2 )
          interpCoeff = 0.0_WP
          if (dist(1,1,1) .le. eps * minGridSpacing) then
             interpCoeff(1,1,1) = 1.0_WP
          else if (dist(2,1,1) .le. eps * minGridSpacing) then
             interpCoeff(2,1,1) = 1.0_WP
          else if (dist(1,2,1) .le. eps * minGridSpacing) then
             interpCoeff(1,2,1) = 1.0_WP
          else if (dist(2,2,1) .le. eps * minGridSpacing) then
             interpCoeff(2,2,1) = 1.0_WP
          else if (dist(1,1,2) .le. eps * minGridSpacing) then
             interpCoeff(1,1,2) = 1.0_WP
          else if (dist(2,1,2) .le. eps * minGridSpacing) then
             interpCoeff(2,1,2) = 1.0_WP
          else if (dist(1,2,2) .le. eps * minGridSpacing) then
             interpCoeff(1,2,2) = 1.0_WP
          else if (dist(2,2,2) .le. eps * minGridSpacing) then
             interpCoeff(2,2,2) = 1.0_WP
          else
             eta = 1.0_WP / dist**2
             buf = sum(eta * alpha(1:2,1:2,1:2))
             interpCoeff(1:2,1:2,1:2) = alpha(1:2,1:2,1:2) * eta / buf
          end if

       else

          ! Image point resides within the levelset, average the neighbors instead
          flag = .true.
          i1 = ghostIndex(1) - 1; i2 = i1 + 2
          j1 = ghostIndex(2) - 1; j2 = j1 + 2
          k1 = ghostIndex(3) - 1; k2 = k1 + 2

          ! Correct for walls
          if (nOverlapGhost(1,1) .eq. 0) i1 = max(i1, iStart(1))
          if (nOverlapGhost(1,2) .eq. 0) i2 = min(i2, iEnd(1))
          if (nOverlapGhost(2,1) .eq. 0) j1 = max(j1, iStart(2))
          if (nOverlapGhost(2,2) .eq. 0) j2 = min(j2, iEnd(2))
          if (nOverlapGhost(3,1) .eq. 0) k1 = max(k1, iStart(3))
          if (nOverlapGhost(3,2) .eq. 0) k2 = min(k2, iEnd(3))

          interpCoeff = 1.0_WP
          interpCoeff(2,2,2) = 0.0_WP
          alpha = 1.0_WP
          do k = k1, k2
             do j = j1, j2
                do i = i1, i2
                   ! Remove points inside fluid
                   if (ghostLevelset(i,j,k,1) .gt. 0.0_WP) alpha(i-i1+1,j-j1+1,k-k1+1) =     &
                        0.0_WP
                end do
             end do
          end do
          buf = sum(alpha * interpCoeff)
          interpCoeff = alpha * interpCoeff / buf
       end if

       ! Interpolate fluid quantities
       T_IP = sum(interpCoeff(1:i2-i1+1,1:j2-j1+1,1:k2-k1+1) *                               &
            ghostTemperature(i1:i2, j1:j2, k1:k2, 1))
       P_IP = sum(interpCoeff(1:i2-i1+1,1:j2-j1+1,1:k2-k1+1) *                               &
            ghostPressure(i1:i2, j1:j2, k1:k2, 1))
       do n = 1, nDimensions
          U_IP(n) = sum(interpCoeff(1:i2-i1+1,1:j2-j1+1,1:k2-k1+1) *                         &
               ghostVelocity(i1:i2, j1:j2, k1:k2, n))
       end do
       do n = 1, nSpecies
          Y_IP(n) = sum(interpCoeff(1:i2-i1+1,1:j2-j1+1,1:k2-k1+1) *                         &
               ghostMassFraction(i1:i2, j1:j2, k1:k2, n))
       end do

    end select

    return
  end subroutine interpolate_image_point


  ! Prepare the arrays for interpolation of image points
  ! ----------------------------------------------------
  subroutine prepare_ghost_arrays

    ! External modules
    use parallel
    use state
    use grid_levelset

    implicit none

    ! Local variables
    integer :: gridIndex, i, j, k

    ! Get data belonging to local processor and pack in temporary array
    do k = iStart(3), iEnd(3)
       do j = iStart(2), iEnd(2)
          do i = iStart(1), iEnd(1)
             gridIndex = i - gridOffset(1) + localGridSize(1) *                              &
                  (j - 1 - gridOffset(2) + localGridSize(2) *                                &
                  (k - 1 - gridOffset(3)))
             ghostBuffer(i,j,k,1) = temperature(gridIndex, 1) !... temperature
             ghostBuffer(i,j,k,2) = pressure(gridIndex, 1) !... pressure
             ghostBuffer(i,j,k,3) = levelset(gridIndex,1) !... levelset
             ghostBuffer(i,j,k,4:nDimensions+3) = velocity(gridIndex,:) !... velocity
          end do
       end do
    end do
    if (nSpecies .gt. 0) then
       do k = iStart(3), iEnd(3)
          do j = iStart(2), iEnd(2)
             do i = iStart(1), iEnd(1)
                gridIndex = i - gridOffset(1) + localGridSize(1) *                           &
                     (j - 1 - gridOffset(2) + localGridSize(2) *                             &
                     (k - 1 - gridOffset(3)))
                ghostBuffer(i,j,k,nDimensions+4:nDimensions+3+nSpecies) =                    &
                     massFraction(gridIndex, 1:nSpecies) !... mass fraction
             end do
          end do
       end do
    end if

    ! Exhange data in direction `1`
    if (nDimensions.ge.1) call fill_ghost_points(ghostBuffer(:, iStart(2) : iEnd(2),         &
         iStart(3) : iEnd(3), 1:nDimensions+3+nSpecies), 1,                                  &
         (/nOverlapGhost(1,1), nOverlapGhost(1,2) /))

    ! Exhange data in direction `2`
    if (nDimensions.ge.2) call fill_ghost_points(ghostBuffer(iStart(1) - nOverlapGhost(1,1) :&
         iEnd(1) + nOverlapGhost(1,2), :, iStart(3):iEnd(3), 1:nDimensions+3+nSpecies), 2,   &
         (/nOverlapGhost(2,1), nOverlapGhost(2,2) /))

    ! Exhange data in direction `3`
    if (nDimensions.ge.3) call fill_ghost_points(ghostBuffer(iStart(1) - nOverlapGhost(1,1) :&
         iEnd(1) + nOverlapGhost(1,2), iStart(2) - nOverlapGhost(2,1) : iEnd(2) +            &
         nOverlapGhost(2,2), :, 1:nDimensions+3+nSpecies), 3, (/nOverlapGhost(3,1),          &
         nOverlapGhost(3,2) /))

    ! Unpack the temporary array
    ghostTemperature(:,:,:,1) = ghostBuffer(:,:,:,1)
    ghostPressure(:,:,:,1) = ghostBuffer(:,:,:,2)
    ghostLevelset(:,:,:,1) = ghostBuffer(:,:,:,3)
    ghostVelocity(:,:,:,1:nDimensions) = ghostBuffer(:,:,:,4:nDimensions+3)
    if (nSpecies .gt. 0) ghostMassFraction(:,:,:,1:nSpecies) =                               &
         ghostBuffer(:,:,:,nDimensions+4:nDimensions+3+nSpecies)

    return
  end subroutine prepare_ghost_arrays

end module ibm_ghost_point


! ======================================== !
! Setup ghost point information for more   !
! efficient implementation during run time !
! ======================================== !
subroutine ibm_ghost_point_setup

  ! Internal modules
  use ibm_ghost_point

  ! External modules
  use math
  use parser
  use parallel
  use geometry
  use simulation_flags
  use solver_options
  use first_derivative
  use grid
  use ibm

  implicit none
  
  ! Local variables
  integer :: i, j, n, imin, imax, jmin, jmax, kmin, kmax, gridIndex, ijk(3), ierror

  ! Return if not used
  if (.not.useIBM .or. .not.use_ghost_points) return

  if (.not. useTargetState) call die('IBM ghost points requires target state!')

  ! Determine whether interior IBM points should be filtered
  call parser_read('filter ibm', filterIBM, .false.)

  ! Determine if ghost points should be enforced weakly
  call parser_read('ibm weak form', useWeakForm, .false.)
  call parser_read('ibm penalty amount', ibmPenalty, 1.0_WP)

  ! Determine depth ghost points act over based on the discretization scheme
  nLayers = floor(0.5_WP * real(maxval(firstDerivative(1:nDimensions)%interiorWidth), WP))

  ! Adjust for upwinding
  if (useUpwinding) nLayers = nLayers + 1

  ! Adjust for repeated first derivatives
  if (useViscosity .and. .not.useSplitViscosity) nLayers = min(2*nLayers, nLayers + 2)

  ! Initialize to zero for now
  nInteriorPoints = 0

  ! Need sufficient overlapping grid points to accomodate interpolation to image point
  nOverlapGhost = 3 * nLayers
  
  ! Adjust the ghost points at the boundaries
  do i = 1, 3
     if (.not. isPeriodic(i)) then
        if (procCoords(i) .eq. 0) nOverlapGhost(i,1) = 0
        if (procCoords(i) .eq. nProcsDir(i) - 1) nOverlapGhost(i,2) = 0
     end if
  end do

  ! Get local indices of grid arrays with overlap
  imin = gridOffset(1) + 1 - nOverlapGhost(1,1)
  imax = gridOffset(1) + localGridSize(1) + nOverlapGhost(1,2)
  if (nDimensions.gt.1) then
     jmin = gridOffset(2) + 1 - nOverlapGhost(2,1)
     jmax = gridOffset(2) + localGridSize(2) + nOverlapGhost(2,2)
  else
     jmin = 1
     jmax = 1
  end if
  if (nDimensions.gt.2) then
     kmin = gridOffset(3) + 1 - nOverlapGhost(3,1)
     kmax = gridOffset(3) + localGridSize(3) + nOverlapGhost(3,2)
  else
     kmin = 1
     kmax = 1
  end if

  allocate(ghostTemperature(imin:imax, jmin:jmax, kmin:kmax, 1))
  allocate(ghostPressure(imin:imax, jmin:jmax, kmin:kmax, 1))
  allocate(ghostLevelset(imin:imax, jmin:jmax, kmin:kmax, 1))
  allocate(ghostVelocity(imin:imax, jmin:jmax, kmin:kmax, nDimensions))
  if (nSpecies.gt.0) allocate(ghostMassFraction(imin:imax, jmin:jmax, kmin:kmax, nSpecies))
  allocate(ghostBuffer(imin:imax, jmin:jmax, kmin:kmax, nDimensions+3+nSpecies))

  ! Store Cartesian coordinates for interpolation
  allocate(cartesianCoordinates(1 - maxval(nOverlapGhost) : maxval(globalGridSize) +         &
       maxval(nOverlapGhost), nDimensions))
  cartesianCoordinates = -huge(1.0_WP)
  do j = 1, nDimensions
     do i = iStart(j), iEnd(j)
        ijk = iStart; ijk(j) = i
        gridIndex = ijk(1) - gridOffset(1) + localGridSize(1) *                              &
             (ijk(2) - 1 - gridOffset(2) + localGridSize(2) *                                &
             (ijk(3) - 1 - gridOffset(3)))
        cartesianCoordinates(i, j) = coordinates(gridIndex, j)
     end do
     call MPI_Allreduce(MPI_IN_PLACE, cartesianCoordinates(:, j),                            &
          size(cartesianCoordinates, 1), MPI_REAL_WP, MPI_MAX,                               &
          commDir(j), ierror)
     if (isPeriodic(j) .and. procCoords(j) .eq. 0) then
        do n = 1, nOverlapGhost(j, 1)
           cartesianCoordinates(1 - n, j) = - periodicLength(j) +                            &
                cartesianCoordinates(globalGridSize(j) - n + 1, j)
        end do
     end if
     if (isPeriodic(j) .and. procCoords(j) .eq. nProcsDir(j) - 1) then
        do n = 1, nOverlapGhost(j, 2)
           cartesianCoordinates(globalGridSize(j) + n, j) = periodicLength(j) +              &
                cartesianCoordinates(n, j)
        end do
     end if
     call MPI_Allreduce(MPI_IN_PLACE, cartesianCoordinates(:, j),                            &
          size(cartesianCoordinates, 1), MPI_REAL_WP, MPI_MAX, commDir(j), ierror)
  end do

  return
end subroutine ibm_ghost_point_setup


! ================================== !
! Cleanup the IBM ghost point module !
! ================================== !
subroutine ibm_ghost_point_cleanup
  
  ! Internal modules
  use ibm_ghost_point

  implicit none

  if (allocated(ghostTemperature)) deallocate(ghostTemperature)
  if (allocated(ghostPressure)) deallocate(ghostPressure)
  if (allocated(ghostMassFraction)) deallocate(ghostMassFraction)
  if (allocated(ghostLevelset)) deallocate(ghostLevelset)
  if (allocated(ghostBuffer)) deallocate(ghostBuffer)
  if (allocated(ghostPoint)) deallocate(ghostPoint)
  if (allocated(cartesianCoordinates)) deallocate(cartesianCoordinates)

  return
end subroutine ibm_ghost_point_cleanup


! ======================================== !
! Update ghost point information for more  !
! efficient implementation during run time !
! ======================================== !
subroutine ibm_ghost_point_update

  ! Internal modules
  use ibm_ghost_point

  ! External modules
  use parallel
  use math, only : bisection
  use simulation_flags
  use filter
  use grid
  use grid_levelset
  use state
  use ibm

  implicit none

  ! Local variables
  integer :: i, j, k, l, n, i1, i2, j1, j2, k1, k2, ii, jj, kk, gridIndex, nn(3,2)
  integer, dimension(nGridPoints) :: tempInteriorIndex
  real(WP) :: dist
  real(WP), dimension(:,:), allocatable :: filtered
  logical :: success
  type(t_GhostPoint), dimension(:), allocatable :: tempGhostPoint

  if (.not.useIBM .or. .not.use_ghost_points) return

  ! Communicate the ghost point arrays
  call prepare_ghost_arrays

  ! Identify `freshly emerged fluid points' and set velocity to filtered velocity field
!!$  if (ibm_move) then
!!$     allocate(filtered(nGridPoints, nUnknowns))
!!$     filtered = conservedVariables
!!$     do i = 1, nDimensions
!!$        call gaussian_filter_apply(i, filtered)
!!$     end do
!!$     do j = 1, nInteriorPoints
!!$        i = interiorIndex(j)
!!$        if (levelset(i,1) .gt. 0.0_WP) then
!!$           ! These points were previously inside the IB, now they are outside
!!$           conservedVariables(i,:) = filtered(i,:)
!!$        end if
!!$     end do
!!$     deallocate(filtered)
!!$  end if

  ! Get extents for localizing image index
  nn = 0
  do i = 1, nDimensions
     if (isPeriodic(i)) then
        if (procCoords(i) .eq. 0) nn(i,1) = nOverlapGhost(i,1)
        if (procCoords(i) .eq. nProcsDir(i) - 1) nn(i,2) = nOverlapGhost(i,2)
     end if
  end do

  ! Store ghost point information
  allocate(tempGhostPoint(nGridPoints))
  nInteriorPoints = 0
  n = 0
  do k = iStart(3), iEnd(3)
     do j = iStart(2), iEnd(2)
        do i = iStart(1), iEnd(1)

           gridIndex = i - gridOffset(1) + localGridSize(1) *                                &
                (j - 1 - gridOffset(2) + localGridSize(2) *                                  &
                (k - 1 - gridOffset(3)))

           ! Cycle if inside the fluid
           if (levelset(gridIndex, 1) .gt. 0.0_WP) cycle

           nInteriorPoints = nInteriorPoints + 1
           tempInteriorIndex(nInteriorPoints) = gridIndex

           ! Get the neighboring points
           i1=i-nLayers; i2=i+nLayers
           j1=j-nLayers; j2=j+nLayers
           k1=k-nLayers; k2=k+nLayers

           ! Correct for walls
           select case (nDimensions)
           case (2)
              k1 = k; k2 = k
              if (nOverlapGhost(1,1) .eq. 0) i1 = max(i1, iStart(1))
              if (nOverlapGhost(1,2) .eq. 0) i2 = min(i2, iEnd(1))
              if (nOverlapGhost(2,1) .eq. 0) j1 = max(j1, iStart(2))
              if (nOverlapGhost(2,2) .eq. 0) j2 = min(j2, iEnd(2))
           case (3)
              if (nOverlapGhost(1,1) .eq. 0) i1 = max(i1, iStart(1))
              if (nOverlapGhost(1,2) .eq. 0) i2 = min(i2, iEnd(1))
              if (nOverlapGhost(2,1) .eq. 0) j1 = max(j1, iStart(2))
              if (nOverlapGhost(2,2) .eq. 0) j2 = min(j2, iEnd(2))
              if (nOverlapGhost(3,1) .eq. 0) k1 = max(k1, iStart(3))
              if (nOverlapGhost(3,2) .eq. 0) k2 = min(k2, iEnd(3))
           end select

           ! Determine if neighbors are inside the fluid
           success = .false.
           do kk = k1, k2
              do jj = j1, j2
                 do ii = i1, i2
                    if (ghostLevelset(ii,jj,kk,1) .gt. 0.0_WP) success = .true.
                 end do
              end do
           end do

           if (success) then
              n = n + 1
              
              ! Store the local grid index
              tempGhostPoint(n)%gridIndex = gridIndex
              tempGhostPoint(n)%ghostIndex(1) = i
              tempGhostPoint(n)%ghostIndex(2) = j
              tempGhostPoint(n)%ghostIndex(3) = k

              ! Get distance to surface
              dist = abs(ghostLevelset(i,j,k,1))

              ! Get image point location and index
              tempGhostpoint(n)%imagePoint(1:nDimensions) =                                  &
                   coordinates(gridIndex, 1:nDimensions) +                                   &
                   2.0_WP * dist * levelsetNormal(gridIndex, 1:nDimensions)
              tempGhostPoint(n)%imageIndex = 1
              do l = 1, nDimensions
                 call bisection(tempGhostPoint(n)%imagePoint(l),                             &
                      tempGhostPoint(n)%imageIndex(l), cartesianCoordinates(                 &
                      1 - nn(l,1):globalGridSize(l) + nn(l,2), l),                           &
                      1 - nn(l,1), globalGridSize(l) + nn(l,2))
                 ! Ensure image points are inside the domain
                 if (nn(l,1) .eq. 0) then
                    tempGhostPoint(n)%imageIndex(l) = max(tempGhostPoint(n)%imageIndex(l), 1)
                    tempGhostpoint(n)%imagePoint(l) = max(tempGhostpoint(n)%imagePoint(l),   &
                         cartesianCoordinates(1, l))
                 end if
                 if (nn(l,2) .eq. 0) then
                    tempGhostPoint(n)%imageIndex(l) = min(tempGhostPoint(n)%imageIndex(l),   &
                         globalGridSize(l))
                    tempGhostpoint(n)%imagePoint(l) = min(tempGhostpoint(n)%imagePoint(l),   &
                         cartesianCoordinates(globalGridSize(l), l))
                 end if
              end do

           end if

        end do
     end do
  end do

  ! Store grid index associated with interior points
  if (allocated(interiorIndex)) deallocate(interiorIndex)
  allocate(interiorIndex(nInteriorPoints))
  do i = 1, nInteriorPoints
     interiorIndex(i) = tempInteriorIndex(i)
  end do

  ! Assign ghost points and cleanup
  if (allocated(ghostPoint)) deallocate(ghostPoint)
  nGhostPoints = n
  allocate(ghostPoint(nGhostPoints))
  do i = 1, nGhostPoints
     ghostPoint(i) = tempGhostPoint(i)
  end do
  deallocate(tempGhostPoint)

  return
end subroutine ibm_ghost_point_update


! ============================= !
! Add ghost point data weakly   !
! This subroutine is not called !
! ============================= !
subroutine ibm_ghost_point_source(source)

  ! Internal modules
  use ibm_ghost_point

  ! External modules
  use math
  use parallel
  use simulation_flags
  use geometry
  use grid
  use grid_functions
  use grid_levelset
  use state
  use time_info
  use ibm

  implicit none

  ! Arguments
  real(WP), intent(inout) :: source(nGridPoints, nUnknowns)

  ! Local variables
  integer :: i, j, n
  real(WP) :: T_IP, P_IP, U_IP(3), Un_IP(3), Y_IP(nSpecies), Ug(3), Tg, Pg, rhog, Eg,        &
        Yg(nSpecies), objectVelocity(nDimensions), r(3), omegaR(3), dti
  logical :: flag, isGhost(nGridPoints)

  ! Return of not used
  if (.not.useIBM .or. .not.use_ghost_points .or. .not.useWeakForm) return

  call timing_start('ibm')

  dti = ibmPenalty / timeStepSize

  ! Prepare arrays for interpolation
  call prepare_ghost_arrays

  ! Store object velocity
  if (ibm_move) then
     ibmVelocity = 0.0_WP
     r = 0.0_WP
     do j = 1, nInteriorPoints
        i = interiorIndex(j)
        n = objectIndex(i)
        r(1:nDimensions) = coordinates(i,:) - object(n)%position(1:nDimensions)
        call correct_periodic_distance(r)
        omegaR = cross_product(object(n)%angularVelocity, r)
        ibmVelocity(i,:) = object(n)%velocity(1:nDimensions) + omegaR(1:nDimensions)
     end do
  end if
  
  ! Ghost point treatment for scalars
  isGhost = .false.
  do j = 1, nGhostPoints
     ! Get the grid index
     i = ghostPoint(j)%gridIndex
     isGhost(i) = .true.

     ! Interpolate to image point
     call interpolate_image_point(ghostPoint(j)%ghostIndex, ghostPoint(j)%imageIndex,        &
          ghostPoint(j)%imagePoint, T_IP, P_IP, U_IP, Y_IP, flag)

     ! Velocity of this object
     if (ibm_move) then
        objectVelocity = ibmVelocity(i,:)
     else
        objectVelocity = 0.0_WP
     end if

     ! Enforce boundary conditions at the ghost points
     if (flag) then
        ! If image point lies within the solid, use neighboring points instead
        Ug(1:nDimensions) = objectVelocity! U_IP(1:nDimensions)
        Tg = T_IP
        Pg = P_IP
     else
        ! Velocity treatment
        Ug(1:nDimensions) = 2.0_WP * objectVelocity - U_IP(1:nDimensions)

        ! Temperature treatment
        if (ibm_isothermal) then
           Tg = 2.0_WP * ibmTemperature - T_IP
        else
           Tg = T_IP
        end if

        ! Pressure treatment
        if (ibm_move) then
           n = objectIndex(i)
           Pg = P_IP! / (1.0_WP - 2.0_WP * abs(levelset(i, 1)) * ratioOfSpecificHeats /       &
               ! (ratioOfSpecificHeats - 1.0_WP) * sum(object(n)%dudt(1:nDimensions) *        &
               ! levelsetNormal(i,1:nDimensions)) / Tg)
           !rhog = conservedVariables(i,1)
           !Pg = P_IP + 2.0_WP * abs(levelset(i,1)) * rhog *                                  &
           !     sum(object(n)%dudt(1:nDimensions) * levelsetNormal(i,1:nDimensions))
        else
           Pg = P_IP
        end if
     end if

     ! Use ideal gas law to get density and energy at the ghost point
     rhog = ratioOfSpecificHeats / (ratioOfSpecificHeats - 1.0_WP) * Pg / Tg
     Eg = Pg / (ratioOfSpecificHeats - 1.0_WP) + 0.5_WP * rhog * sum(Ug(1:nDimensions)**2)

     ! Overwrite the state vector
     !source(i,:) = 0.0_WP
     source(i,1) = source(i,1) + (rhog - conservedVariables(i,1)) * dti
     source(i,2:nDimensions+1) = source(i,2:nDimensions+1) + rhog * (Ug(1:nDimensions) - velocity(i,:)) * dti
     source(i,nDimensions+2) = source(i,nDimensions+2) + (Eg - conservedVariables(i,nDimensions+2)) * dti
  end do

  ! Penalize non-ghost point interior points to target state
  do j = 1, nInteriorPoints
     i = interiorIndex(j)
     if (.not.isGhost(i)) then
        rhog = 1.0_WP
        Ug = 0.0_WP
        Pg = 1.0_WP / ratioOfSpecificHeats
        if (ibm_move) Ug(1:nDimensions) = ibmVelocity(i,1:nDimensions)
        Eg = Pg / (ratioOfSpecificHeats - 1.0_WP) + 0.5_WP * rhog * sum(Ug(1:nDimensions)**2)
        source(i,1) = source(i,1) + (rhog - conservedVariables(i,1)) * dti
        source(i,2:nDimensions+1) = source(i,2:nDimensions+1) + rhog * &
             (Ug(1:nDimensions) - velocity(i,1:nDimensions)) * dti
        source(i,nDimensions+2) = source(i,nDimensions+2) + (Eg - conservedVariables(i,nDimensions+2)) * dti
     end if
  end do

  call timing_stop('ibm')

  return
end subroutine ibm_ghost_point_source


! ======================================================= !
! Correct interior points using the ghost point treatment !
! ======================================================= !
subroutine ibm_ghost_point_correct_state(stateVector)

  ! Internal modules
  use ibm_ghost_point

  ! External modules
  use parallel
  use math
  use simulation_flags
  use solver_options
  use geometry
  use grid
  use grid_functions
  use grid_levelset
  use filter
  use state
  use ibm
  
  implicit none

  ! Arguments
  real(WP), intent(inout) :: stateVector(nGridPoints, nUnknowns)

  ! Local variables
  integer :: i, j, n
  real(WP) :: T_IP, P_IP, U_IP(3), Un_IP(3), Ug(3), Tg, Pg, rhog, Yg(nSpecies), Y_IP(nSpecies),&
       objectVelocity(nDimensions), r(3), omegaR(3), kappa
  logical :: flag

  ! Return if not used
  if (.not.useIBM .or. .not.use_ghost_points .or.useWeakForm) return

  call timing_start('ibm')

  ! Prepare arrays for interpolation
  call update_state
  call prepare_ghost_arrays

  ! Reset interior points
  do j = 1, nInteriorPoints
     i = interiorIndex(j)
     stateVector(i,:) = targetState(i,:)
  end do

  ! Store object velocity
  if (ibm_move) then
     ibmVelocity = 0.0_WP
     r = 0.0_WP
     do j = 1, nInteriorPoints
        i = interiorIndex(j)
        n = objectIndex(i)
        r(1:nDimensions) = coordinates(i,:) - object(n)%position(1:nDimensions)
        call correct_periodic_distance(r)
        omegaR = cross_product(object(n)%angularVelocity, r)
        ibmVelocity(i,:) = object(n)%velocity(1:nDimensions) + omegaR(1:nDimensions)
     end do
  end if
  
  ! Ghost point treatment for scalars
  do j = 1, nGhostPoints
     ! Get the grid index
     i = ghostPoint(j)%gridIndex

     ! Interpolate to image point
     call interpolate_image_point(ghostPoint(j)%ghostIndex, ghostPoint(j)%imageIndex,        &
          ghostPoint(j)%imagePoint, T_IP, P_IP, U_IP, Y_IP, flag)

     ! Velocity of this object
     if (ibm_move) then
        objectVelocity = ibmVelocity(i,:)
     else
        objectVelocity = 0.0_WP
     end if

     ! Enforce boundary conditions at the ghost points
     if (flag) then
        ! If image point lies within the solid, use neighboring points instead
        Ug(1:nDimensions) = objectVelocity! U_IP(1:nDimensions)
        Tg = T_IP
        Pg = P_IP
        if (nSpecies .gt. 0) Yg = Y_IP
     else
        ! Velocity treatment
        Ug(1:nDimensions) = objectVelocity!2.0_WP * objectVelocity - U_IP(1:nDimensions)

        ! Temperature treatment
        if (ibm_isothermal) then
           Tg = ibmTemperature! 2.0_WP * ibmTemperature - T_IP
        else
           Tg = T_IP
        end if

        ! Pressure treatment
        if (ibm_move) then
           n = objectIndex(i)
           Pg = P_IP / (1.0_WP - 2.0_WP * abs(levelset(i, 1)) * ratioOfSpecificHeats /       &
                (ratioOfSpecificHeats - 1.0_WP) * sum(object(n)%dudt(1:nDimensions) *        &
                levelsetNormal(i,1:nDimensions)) / Tg)
        else
           Pg = P_IP
        end if

        ! Scalar treatment
        if (nSpecies .gt. 0) Yg = Y_IP
     end if

     ! Use ideal gas law to get density at the ghost point
     rhog = ratioOfSpecificHeats / (ratioOfSpecificHeats - 1.0_WP) * Pg / Tg

     ! Overwrite the state vector
     stateVector(i, 1) = rhog
     stateVector(i, 2:nDimensions+1) = stateVector(i, 1) * Ug(1:nDimensions)
     stateVector(i, nDimensions+2) = Pg / (ratioOfSpecificHeats - 1.0_WP) +                  &
          0.5_WP * sum(stateVector(i, 2:nDimensions+1)**2) / stateVector(i, 1)
     if (nSpecies .gt. 0) stateVector(i, nDimensions+3:nDimensions+2+nSpecies) =             *
     stateVector(i, 1) * Yg

     ! Correct for two-way coupling
     if (twoWayCoupling) then
        do n = 1, nUnknowns
           stateVector(i, n) = stateVector(i, n) * volumeFraction(i,1)
        end do
     end if
  end do

  call timing_stop('ibm')

  return
end subroutine ibm_ghost_point_correct_state


! ======================================================= !
! Correct interior points using the ghost point treatment !
! ======================================================= !
subroutine ibm_filter_state(stateVector)

  ! Internal modules
  use ibm_ghost_point

  ! External modules
  use simulation_flags
  use solver_options
  use geometry
  use filter
  
  implicit none

  ! Arguments
  real(WP), intent(inout) :: stateVector(nGridPoints, nUnknowns)

  ! Local variables
  integer :: i, j
  real(WP), dimension(:,:), allocatable :: filtered
  
  ! Return of not used
  if (.not.useIBM .or. .not.filterIBM) return

  call timing_start('ibm')

  allocate(filtered(nGridPoints, nUnknowns))
  filtered = stateVector
  do i = 1, nDimensions
     call gaussian_filter_apply(i, filtered)
  end do
  do j = 1, nInteriorPoints
     i = interiorIndex(j)
     stateVector(i,:) = filtered(i,:)
  end do
  deallocate(filtered)

  call timing_stop('ibm')

  return
end subroutine ibm_filter_state
