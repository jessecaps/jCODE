module ibm

  ! External modules
  use precision

  implicit none

  ! IBM type
  integer :: ibmType
  integer, parameter ::                                                                      &
       IBM_GENERIC   = 1,                                                                    &
       IBM_PARTICLE  = 2

  ! IBM objects
  integer ::  nObjects
  type t_Object
     real(WP)               :: volume          !... Object volume
     real(WP), dimension(3) :: position        !... Center of mass
     real(WP), dimension(3) :: velocity        !... Translational velocity
     real(WP), dimension(3) :: angularVelocity !... Angular velocity
     real(WP), dimension(3) :: pForce          !... Force due to pressure
     real(WP), dimension(3) :: vForce          !... Force due to viscous stress
     real(WP), dimension(3) :: cForce          !... Force due to collisions
     real(WP), dimension(3) :: hTorque         !... Torque due to hydrodynamic stress
     real(WP), dimension(3) :: cTorque         !... Torque due to collisions
     real(WP), dimension(3) :: dudt            !... RHS for advancing velocity
     real(WP), dimension(3) :: dwdt            !... RHS for advancing ang. velocity
     logical                :: remove          !... Flag to remove object
  end type t_Object
  type(t_Object), dimension(:), allocatable :: object
  integer :: MPI_OBJECT, SIZE_MPI_OBJECT

  ! IBM parameters
  integer :: rigidBodyID
  logical :: ibm_move, use_ghost_points, ibm_slip, ibm_isothermal, ibm_rigid_motion,         &
       correctCurvature
  integer, dimension(:), allocatable :: objectIndex
  real(WP) :: ibmTemperature, ibmEpsilon
  real(WP), dimension(:,:), allocatable :: ibmVelocity, baselineLevelset

contains

  ! Send `removed` objects to the end and resize
  ! ---------------------------------------------
  subroutine recycle_objects(obj)

    implicit none

    ! Arguments
    type(t_Object), allocatable, intent(inout) :: obj(:)

    ! Local variables
    integer :: i, newSize

    if (.not.ibm_move) return

    ! Compact real objects at the beginning of the array
    newSize = 0
    if (allocated(obj)) then
       do i = 1, size(obj)
          if (.not.obj(i)%remove) then
             newSize = newSize + 1
             if (i .ne. newSize) then
                obj(newSize) = obj(i)
                obj(i)%remove = .true.
             end if
          end if
       end do
    end if

    ! Resize array to newSize
    call resize_objects(obj, newSize)

    ! Update object size
    nObjects = size(obj)

    return
  end subroutine recycle_objects


  ! Resize the object based on new size `n`
  ! ---------------------------------------
  subroutine resize_objects(obj, n)

    implicit none

    ! Arguments
    type(t_Object), allocatable, intent(inout) :: obj(:)
    integer, intent(in) :: n

    ! Local variables
    integer :: i, nObjectsOld
    type(t_Object), allocatable :: tempObject(:)

    ! Resize part array to size n
    if (.not. allocated(obj)) then
       ! Object is of size 0
       if (n .eq. 0) then
          ! Nothing to do, that's what we want
       else
          ! Allocate directly of size n
          allocate(obj(n))
          obj(1:n)%remove = .true.
       end if
    else if (n .eq. 0) then
       ! Empty the object array
       deallocate(obj); allocate(obj(0))
    else
       ! Update non zero size to another non zero size
       nObjectsOld = size(obj)
       if (n .gt. nObjectsOld) then
          ! Increase from `nObjectsOld` to `n`
          allocate(tempObject(n))
          do i = 1, nObjectsOld
             tempObject(i) = obj(i)
          end do
          deallocate(obj)
          allocate(obj(size(tempObject)))
          obj = tempObject
          obj(nObjectsOld + 1:n)%remove = .true.
       else if (n .lt. nObjectsOld) then
          ! Decrease from `nObjectsOld` to `n`
          allocate(tempObject(n))
          do i = 1, n
             tempObject(i) = obj(i)
          end do
          deallocate(obj)
          allocate(obj(size(tempObject)))
          obj = tempObject
       end if
    end if

    return
  end subroutine resize_objects

end module ibm

! ======================== !
! Setup MPI object for i/o !
! ======================== !
subroutine prepare_mpi_object

  ! Internal modules
  use ibm

  ! External modules
  use parallel

  implicit none

  ! Local variables
  integer, dimension(32) :: types, lengths, displacement ! ... size based on t_Object
  integer :: i, ierror

  ! Create the MPI structure to read IBM objects
  types(1:31) = MPI_REAL_WP
  types(32)   = MPI_LOGICAL
  lengths(:)  = 1

  ! Hard-code IBM object type displacement here
  displacement(1) = 0
  do i = 1, 31
     displacement(i+1) = displacement(i) + WP
  end do

  ! Finalize by creating and commiting the new type
  call MPI_Type_struct(32, lengths, displacement, types, MPI_OBJECT, ierror)
  call MPI_Type_commit(MPI_OBJECT, ierror)

  ! If problem, say it
  if (ierror .ne. 0) call die('Problem with MPI_OBJECT!')

  ! Get the size of this type
  call MPI_type_size(MPI_OBJECT, SIZE_MPI_OBJECT, ierror)

  return
end subroutine prepare_mpi_object


! ==================== !
! Setup the IBM module !
! ==================== !
subroutine ibm_setup

  ! Internal modules
  use ibm

  ! External modules
  use string
  use parallel
  use parser
  use simulation_flags
  use geometry
  use solver_options

  implicit none

  ! Local variables
  character(len = str_medium) :: val

  if (.not. useIBM) return

  ! IBM not yet implemented for curvilinear grids
  if (isDomainCurvilinear) call die('IBM currently implemented for Cartesian grids only!')

  ! Determine IBM type
  call parser_read('ibm type', val, '')
  select case (trim(val))

  case('PARTICLE', 'particle', 'Particle')
     ! Spheres (3D) or cylinders (2D)
     ibmType = IBM_PARTICLE

  case default
     ! Generic shape
     ibmType = IBM_GENERIC
  end select

  ! Initialize IBM routines
  nObjects = 0 !... To be determined after IBM file is read
  allocate(objectIndex(nGridPoints)); objectIndex = 1

  ! Regularization parameter
  call parser_read('regularization parameter', ibmEpsilon, 0.0_WP)

  ! Enforce slip condition?
  call parser_read('ibm slip', ibm_slip, .false.)
  if (ibm_slip) call parser_read('ibm curvature correction', correctCurvature, .false.)

  ! Temperature treatment
  call parser_read('ibm isothermal', ibm_isothermal, .false.)
  if (ibm_isothermal) call parser_read('ibm temperature', ibmTemperature,                    &
       1.0_WP / (ratioOfSpecificHeats - 1.0_WP))

  ! Moving IBM?
  call parser_read('ibm move', ibm_move, .false.)
  if (ibm_move) then
     allocate(ibmVelocity(nGridPoints, nDimensions))
     call parser_read('ibm rigid body motion', ibm_rigid_motion, .false.)
     if (ibm_rigid_motion) then
        call parser_read('ibm rigid body id', rigidBodyID, 1)
     else
        call ibm_solver_setup
     end if
  end if

  ! Store baseline levelset in case multiple levelsets are used (from file and IBM objects)
  allocate(baselineLevelset(nGridPoints, 1)); baselineLevelset = huge(1.0_WP)

  ! Setup the ghost point routine
  call parser_read('ibm ghost points', use_ghost_points, .true.)
  if (use_ghost_points) call ibm_ghost_point_setup

  return
end subroutine ibm_setup


! ====================== !
! Cleanup the IBM module !
! ====================== !
subroutine ibm_cleanup

  ! Internal modules
  use ibm

  implicit none

  call ibm_solver_cleanup
  call ibm_ghost_point_cleanup

  if (allocated(ibmVelocity)) deallocate(ibmVelocity)
  if (allocated(baselineLevelset)) deallocate(baselineLevelset)

  return
end subroutine ibm_cleanup


! ===================================== !
! Setup the objects associated with IBM !
! ===================================== !
subroutine ibm_setup_objects

  ! Internal modules
  use ibm

  ! External modules
  use string
  use parser
  use simulation_flags
  use grid_levelset

  implicit none

  ! Local parameters
  integer :: i
  logical :: readLevelset
  character(len=str_medium) :: levelsetFile

  ! Return if not used
  if (.not. useIBM) return

  ! Start the IBM timer
  call timing_start('ibm')

  ! Zero-out object properties and remove objects that are out-of-bounds
  do i = 1, nObjects
     object(i)%pForce = 0.0_WP
     object(i)%vForce = 0.0_WP
     object(i)%cForce = 0.0_WP
     object(i)%cTorque = 0.0_WP
     object(i)%dudt = 0.0_WP
     object(i)%dwdt = 0.0_WP
     call check_object_bounds(object(i))
  end do
  call recycle_objects(object)
  
  ! Get the levelset
  call parser_is_defined('levelset file to read', readLevelset)
  if (readLevelset) then
     call parser_read('levelset file to read', levelsetFile)
     ! Assume levelset was already read and setup
     call ibm_ghost_point_update
     ! Set baseline levelset to levelset read from file
     if (ibmType .eq. IBM_PARTICLE) then
        baselineLevelset = levelset
        call ibm_compute_levelset
     end if
  else
     call allocate_levelset
     call ibm_compute_levelset
  end if

  ! Stop the IBM timer
  call timing_stop('ibm')

  return
end subroutine ibm_setup_objects


! ==================== !
! Compute the levelset !
! ==================== !
subroutine ibm_compute_levelset

  ! Internal modules
  use ibm

  ! External modules
  use math
  use parallel
  use parser
  use simulation_flags
  use solver_options
  use geometry
  use grid
  use grid_levelset
  use grid_functions

  implicit none

  ! Local variables
  integer :: i, j, n, nObjectsLocal
  integer, dimension(nObjects) :: localObjectIndex
  real(WP) :: distance, radius, px, pxp, pxm, x1, x2
  real(WP), dimension(nDimensions) :: dist, minCoords, maxCoords
  logical, dimension(nDimensions) :: success

  if (.not. useIBM) return

  select case (ibmType)

  case(IBM_PARTICLE)

     if (ibm_move) then
        ! Localize IBM objects to processor to speed up levelset computation
        minCoords = minval(coordinates, 1)
        maxCoords = maxval(coordinates, 1)
        nObjectsLocal = 0
        do n = 1, nObjects

           ! Ignore objects that will be removed
           if (object(n)%remove) cycle

           ! Get radius from volume
           radius = 0.5_WP * (2.0_WP * real(nDimensions, WP) * object(n)%volume / pi)        &
                ** (1.0_WP / real(nDimensions, WP))

           ! Check if object is within local grid extent
           do i = 1, nDimensions
              success(i) = .false.
              px = object(n)%position(i)
              pxp = px; pxm = px
              if (isPeriodic(i)) then
                 pxp = pxp + periodicLength(i)
                 pxm = pxm - periodicLength(i)
              end if
              x1 = minCoords(i) - radius - 4.0_WP * minGridSpacing
              x2 = maxCoords(i) + radius + 4.0_WP * minGridSpacing
              if ((px.ge.x1.and.px.le.x2) .or. (pxp.ge.x1.and.pxp.le.x2) .or.                &
                   (pxm.ge.x1.and.pxm.le.x2)) success(i) = .true.
           end do

           ! Update list
           if (all(success)) then
              nObjectsLocal = nObjectsLocal + 1
              localObjectIndex(nObjectsLocal) = n
           end if

        end do

     else

        ! Revert to full list
        nObjectsLocal = nObjects
        do n = 1, nObjects
           localObjectIndex(n) = n
        end do

     end if

     do i = 1, nGridpoints
        ! Initialize distance       
        distance = baselineLevelset(i,1)

        ! Loop over all particles
        do j = 1, nObjectsLocal

           ! Get index
           n = localObjectIndex(j)
           
           ! Get radius from volume
           radius = 0.5_WP * (2.0_WP * real(nDimensions, WP) * object(n)%volume / pi)        &
                ** (1.0_WP / real(nDimensions, WP))

           ! Get distance to particle center
           dist = coordinates(i,:) - object(n)%position(1:nDimensions)
           call correct_periodic_distance(dist)

           ! Get signed distance to particle interface and compare to current
           distance = min(distance, sqrt(sum(dist**2)) - radius)

           ! Store object index
           if (sqrt(sum(dist**2)) .le. radius) objectIndex(i) = n
        end do

        ! Store min distance
        levelset(i,1) = distance
     end do

  case default

     call die('ibm_compute_levelset: levelset file must be specified in the input!')

  end select

  ! Compute levelset dependent variables
  call grid_levelset_dependent_variables

  ! Update ghost point
  call ibm_ghost_point_update

  return
end subroutine ibm_compute_levelset


! ==================================== !
! IBM volume penalization source terms !
! ==================================== !
subroutine ibm_cbvp_source(source)

  ! Internal modules
  use ibm

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

  implicit none

  ! Arguments
  real(WP), intent(inout) :: source(nGridPoints, nUnknowns)

  ! Local variables
  integer :: i, n
  real(WP) :: rho_, nDotGradRho, uDotGradRHo, dti, diffusionAmount, maxSpeed, kappa,         &
       weight, buf, velocitySquared, source_(nUnknowns)
  real(WP), dimension(nDimensions) :: velocityPenalty, normalVelocity, tangentVelocity,      &
       objectVelocity
  real(WP), dimension(nGridPoints, nDimensions) :: densityGradient, energyGradient,          &
       normalVelocityGradient
  real(WP), dimension(nGridPoints, nUnknowns) :: laplacianConservedVariables

  ! Return of not used
  if (.not.useIBM .or. use_ghost_points) return

  ! Start the IBM timer
  call timing_start('ibm')
  
  ! Set the max speed based on the reference sound speed
  maxSpeed = 1.0_WP

  ! Store inverse timestep size
  dti = 1.0_WP / timeStepSize

  ! Set the diffusion amount based on the stability limit
  diffusionAmount = 0.05_WP * minGridSpacing**2 * dti / real(nDimensions, WP)

  ! Compute derivatives
  call laplacian(conservedVariables, laplacianConservedVariables)
  call gradient(conservedVariables(:,1), densityGradient)
  call gradient(conservedVariables(:,1) * temperature(:, 1), energyGradient)
  energyGradient = energyGradient / ratioOfSpecificHeats

  ! Compute normal velocity gradient
  if (ibm_slip) then
     do n = 1, nDimensions
        normalVelocityGradient(:,n) = -sum(levelsetNormal *                                  &
             velocityGradient(:, 1+nDimensions*(n-1):nDimensions*n), dim = 2)
     end do
  else
     normalVelocityGradient = 0.0_WP
  end if

  if (ibm_move) ibmVelocity = 0.0_WP

  ! Update the source terms and IBM forcing
  do i = 1, nGridPoints

     ! Get local density
     rho_ = conservedVariables(i, 1)

     ! Get velocity of associated object
     if (ibm_move) then
        n = objectIndex(i)
        objectVelocity = object(n)%velocity(1:nDimensions)
        ibmVelocity(i,:) = objectVelocity
     else
        objectVelocity = 0.0_WP
     end if

     ! Compute the velocity penalty
     if (ibm_slip) then
        normalVelocity = sum((velocity(i,:) - objectVelocity) * levelsetNormal(i,:)) *       &
             levelsetNormal(i,:)
        velocityPenalty = - normalVelocity
     else
        velocityPenalty = objectVelocity - velocity(i,:)
     end if

     ! Store dot products
     nDotGradRho = sum(levelsetNormal(i,:) * densityGradient(i,:))
     uDotGradRho = sum(objectVelocity * densityGradient(i,:))
     
     ! Zero-out the local source terms
     source_ = 0.0_WP

     ! Density treatment
     source_(1) = maxSpeed * nDotGradRho - uDotGradRho +                                     &
          diffusionAmount * laplacianConservedVariables(i, 1)

     ! Momentum treatment
     do n = 1, nDimensions
        source_(n+1) = rho_ * velocityPenalty(n) * dti +                                     &
             velocity(i, n) * (maxSpeed * nDotGradRho - uDotGradRho) +                       &
             diffusionAmount * laplacianConservedVariables(i, n+1)
     end do

     ! Energy treatment
     velocitySquared = sum(velocity(i,:)**2)
     source_(nDimensions+2) = 0.5_WP * velocitySquared * (maxSpeed * nDotGradRho -           &
          uDotGradRho) + diffusionAmount * laplacianConservedVariables(i, nDimensions+2) +   &
          sum(conservedVariables(i, 2:nDimensions+1) * velocityPenalty(1:nDimensions)) * dti

     if (ibm_isothermal) then
        ! Isothermal
        source_(nDimensions+2) = source_(nDimensions+2) +                                    &
             rho_ * (ibmTemperature - temperature(i,1)) / ratioOfSpecificHeats * dti
     else
        ! Adiabatic
        source_(nDimensions+2) = source_(nDimensions+2) + maxSpeed *                         &
             sum(levelsetNormal(i,:) * energyGradient(i,:))
     end if

     ! Slip treatment
     if (ibm_slip) then
        do n = 1, nDimensions
           buf = - rho_ * maxSpeed * normalVelocityGradient(i,n)
           source_(n+1) = source_(n+1) + buf
           source_(nDimensions+2) = source_(nDimensions+2) + velocity(i, n) * buf
        end do
        ! Curvature correction (need to update how we compute curvature)
        if (correctCurvature .and. abs(levelset(i, 1)).lt.8.0_WP * minGridSpacing) then
           tangentVelocity = velocity(i,:) - normalVelocity
           !kappa = sign(1.0_WP, levelsetCurvature(i, 1)) * min(abs(levelsetCurvature(i, 1)), &
           !     0.5_WP / primitiveGridNorm(i,1)**(1.0_WP / real(nDimensions,WP)))
           buf = -maxSpeed * kappa * rho_**2 / ratioOfSpecificHeats *                        &
                sum(tangentVelocity**2) / pressure(i, 1)
           source_(1) = source_(1) + buf
           source_(2:nDimensions+1) = source_(2:nDimensions+1) + buf * velocity(i,:) +       &
                maxSpeed * rho_ * kappa * tangentVelocity
           source_(nDimensions+2) = source_(nDimensions+2) +                                 &
                0.5_WP * buf * velocitySquared + maxSpeed * rho_ *                           &
                sum(velocity(i,:) * tangentVelocity) * kappa *                               &
                (1.0_WP - 1.0_WP / (ratioOfSpecificHeats - 1.0_WP))
        end if
     end if

     ! Add the IBM contribution
     buf = ibmEpsilon * sqrt(sum((levelsetNormal(i,:) * gridSpacing(i, :))**2))
     weight = 1.0_WP - regularize_heaviside(levelset(i,1), buf)
     if (levelset(i, 1) .le. 0.0_WP) then
        source(i,:) = source_ * weight
     else
        source(i,:) = source(i,:) + source_ * weight
     end if
  end do

  ! Stop the IBM timer
  call timing_stop('ibm')

  return
end subroutine ibm_cbvp_source


! ===================================================================== !
! Compute drag by IBM by integrading divergence of stresses on the grid !
! ===================================================================== !
subroutine ibm_integrate_forces

  ! Internal modules
  use ibm

  ! External modules
  use parallel
  use math
  use first_derivative
  use grid_functions
  use grid_levelset
  use grid
  use state

  implicit none

  ! Local variables
  integer :: i, n
  real(WP) :: r(3), stress(3), torque(3)
  real(WP), dimension(nGridPoints, nDimensions, nDimensions) :: pressureStress,              &
       viscousStress, temp1, temp2

  if (nObjects.eq.0) return

  ! Zero out forces
  do i = 1, nObjects
     object(i)%pForce  = 0.0_WP
     object(i)%vForce  = 0.0_WP
     object(i)%hTorque = 0.0_WP
  end do

  ! Zero out the array
  pressureStress = 0.0_WP
  viscousStress = 0.0_WP

  ! Compute inviscid & viscous fluxes
  select case (nDimensions)

  case (1)
     pressureStress(:,1,1) = -pressure(:,1)
     if (useViscosity) then
        viscousStress(:,1,1) = stressTensor(:,1)
     end if

  case (2)
     pressureStress(:,1,1) = -pressure(:,1)
     pressureStress(:,2,2) = -pressure(:,1)
     if (useViscosity) then
        viscousStress(:,1,1) = stressTensor(:,1)
        viscousStress(:,1,2) = stressTensor(:,2)
        viscousStress(:,2,1) = stressTensor(:,3)
        viscousStress(:,2,2) = stressTensor(:,4)
     end if

  case (3)
     pressureStress(:,1,1) = -pressure(:,1)
     pressureStress(:,2,2) = -pressure(:,1)
     pressureStress(:,3,3) = -pressure(:,1)
     if (useViscosity) then
        viscousStress(:,1,1) = stressTensor(:,1)
        viscousStress(:,1,2) = stressTensor(:,2)
        viscousStress(:,1,3) = stressTensor(:,3)
        viscousStress(:,2,1) = stressTensor(:,4)
        viscousStress(:,2,2) = stressTensor(:,5)
        viscousStress(:,2,3) = stressTensor(:,6)
        viscousStress(:,3,1) = stressTensor(:,7)
        viscousStress(:,3,2) = stressTensor(:,8)
        viscousStress(:,3,3) = stressTensor(:,9)
     end if

  end select

  ! Transform stress tensor from Cartesian to contravariant form
  call transform_fluxes(pressureStress, temp1)
  call transform_fluxes(viscousStress, temp2)

  ! Take derivatives of stress tensor
  do i = 1, nDimensions
     call first_derivative_apply(i, temp1(:,:,i))
     call first_derivative_apply(i, temp2(:,:,i))
  end do

  r = 0.0_WP
  stress = 0.0_WP
  do i = 1, nGridPoints

     ! Only consider interior points
     if (levelset(i, 1) .gt. 0.0_WP) cycle

     ! Get the object index
     n = objectIndex(i)

     ! Take divergence of stress and sum
     object(n)%pForce(1:nDimensions) = object(n)%pForce(1:nDimensions) +                     &
          sum(temp1(i,:,:), dim = 2) * jacobian(i,1) * primitiveGridNorm(i, 1)
     object(n)%vForce(1:nDimensions) = object(n)%vForce(1:nDimensions) +                     &
          sum(temp2(i,:,:), dim = 2) * jacobian(i,1) * primitiveGridNorm(i, 1)

     ! Compute hydrodynamic torque and sum
     r(1:nDimensions) = coordinates(i,:) - object(n)%position(1:nDimensions)
     call correct_periodic_distance(r)
     stress(1:nDimensions) = sum(temp1(i,:,:) + temp2(i,:,:), dim = 2) * jacobian(i,1)
     torque = cross_product(r, stress)
     object(n)%hTorque(1:3) = object(n)%hTorque(1:3) + torque(1:3) * primitiveGridNorm(i, 1)
  end do

  ! Sum up forces and torques
  do n = 1, nObjects
     do i = 1, nDimensions
        call parallel_sum(object(n)%pForce(i))
        call parallel_sum(object(n)%vForce(i))
     end do
     do i = 1, 3
        call parallel_sum(object(n)%hTorque(i))
     end do
  end do

  return
end subroutine ibm_integrate_forces
