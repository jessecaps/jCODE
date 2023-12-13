module controller

  ! External modules
  use precision
  use string
  use grid_patch

  implicit none

  ! Controller types
  integer, parameter ::                                                                      &
       NULL_ACTUATOR         = 0,                                                            &
       THERMAL_ACTUATOR      = 1,                                                            &
       MOMENTUM_ACTUATOR     = 2,                                                            &
       FUEL_ACTUATOR         = 3,                                                            &
       IGNITION_ACTUATOR     = 4,                                                            &
       CHEMICAL_ACTUATOR     = 5,                                                            &
       PERTURBATION_ACTUATOR = 6,                                                            &
       GENERIC_IC_ACTUATOR   = 7

  ! Actuator data
  type(t_Patch), pointer :: controllerPatch
  integer :: controllerType, nControlParameters, iGradientBuffer, gradientDirection,         &
       nConstraint
  integer(kind = MPI_OFFSET_KIND) :: gradientFileOffset
  real(WP), allocatable :: controlMollifier(:,:), controlForcing(:,:),                       &
       gradientBuffer(:,:,:), instantaneousSensitivity(:), currentSensitivity(:),            &
       controlGradient(:), baselineValue(:)
  logical :: spaceTimeGradient
  logical, dimension(:), allocatable :: equalityConstraint
  character(len = str_medium) :: gradientFilename
  character(len = str_medium), allocatable :: sensitivityParameter(:)

contains

  ! Verify the controller patch is correct
  ! --------------------------------------
  subroutine verify_controller_patch(patch)

    ! External modules
    use string
    use simulation_flags
    use geometry

    implicit none

    ! Arguments
    type(t_Patch), intent(in) :: patch

    ! Local variables
    integer :: i, n, extent(6)
    character(len = str_long) :: message

    if (predictionOnly) call die('verify_controller_patch: controller requires adjoint run')

    extent(1) = patch%iMin; extent(2) = patch%iMax
    extent(3) = patch%jMin; extent(4) = patch%jMax
    extent(5) = patch%kMin; extent(6) = patch%kMax

    n = nDimensions
    do i = 1, nDimensions
       if (extent((i-1)*2+1) .lt. 0 .or. extent((i-1)*2+2) .gt. globalGridSize(i) .or.       &
            extent((i-1)*2+1) .gt. extent((i-1)*2+2))                                        &
            call die("verify_controller_patch: Invalid extent on '" //                       &
            trim(patch%name) // "'!")
       if (extent((i-1)*2+1) .eq. extent((i-1)*2+2)) n = n - 1
    end do


    if (n .ne. nDimensions) then
       write(message, '(2(A,I0.0),3A)') "verify_controller_patch: Expected a ", nDimensions, &
            "D patch, but extent represents a ", n, "D patch on '", trim(patch%name), "'!"
       call die(trim(message))
    end if

    return
  end subroutine verify_controller_patch


  ! Normalize the control mollifier
  ! -------------------------------
  subroutine normalize_control_mollifier

    ! External modules
    use parallel

    implicit none

    ! Local variables
    real(WP) :: mollifierNorm
    logical :: hasNegativeMollifier

    hasNegativeMollifier = any(controlMollifier(:,1) .lt. 0.0_WP)
    call parallel_lor(hasNegativeMollifier)
    if (hasNegativeMollifier)                                                                &
         call die('Control mollifying support is not non-negative everywhere!')

    call patch_quadrature(controllerPatch, controlMollifier(:,1), mollifierNorm)

    if (mollifierNorm .gt. 0.0_WP) controlMollifier = controlMollifier / mollifierNorm

    return
  end subroutine normalize_control_mollifier


  ! Load the sensitivity gradient on the controller patch
  ! -----------------------------------------------------
  subroutine load_sensitivity_gradient

    ! External modules
    use parallel

    implicit none

    ! Local variables
    integer :: extent(6), arrayOfSizes(5), arrayOfSubsizes(5), arrayOfStarts(5),             &
         mpiRealSubarrayType, mpiFileHandle, dataSize, ierror
    integer(kind = MPI_OFFSET_KIND) :: nBytesToRead

    if (controllerPatch%comm .eq. MPI_COMM_NULL) return

    extent(1) = controllerPatch%iMin; extent(2) = controllerPatch%iMax
    extent(3) = controllerPatch%jMin; extent(4) = controllerPatch%jMax
    extent(5) = controllerPatch%kMin; extent(6) = controllerPatch%kMax

    arrayOfSizes(1:3) = controllerPatch%globalSize
    arrayOfSizes(4) = size(gradientBuffer, 2)
    arrayOfSizes(5) = iGradientBuffer
    arrayOfSubsizes(1:3) = controllerPatch%localSize
    arrayOfSubsizes(4) = size(gradientBuffer, 2)
    arrayOfSubsizes(5) = iGradientBuffer
    arrayOfStarts(1:3) = controllerPatch%offset - extent(1::2) + 1
    arrayOfStarts(4:5) = 0
    call MPI_Type_create_subarray(5, arrayOfSizes, arrayOfSubsizes, arrayOfStarts,           &
         MPI_ORDER_FORTRAN, MPI_REAL_WP, mpiRealSubarrayType, ierror)
    call MPI_Type_commit(mpiRealSubarrayType, ierror)

    call MPI_File_open(controllerPatch%comm, trim(mpiiofs)//trim(gradientFilename)//char(0), &
         MPI_MODE_RDONLY, mpiInfo, mpiFileHandle, ierror)

    nBytesToRead = int(WP, MPI_Offset_kind) * product(int(controllerPatch%globalSize,        &
         MPI_OFFSET_KIND)) * size(gradientBuffer, 2) * iGradientBuffer
    if (gradientFileOffset - nBytesToRead .le. int(0, MPI_OFFSET_KIND)) then
       iGradientBuffer = int(gradientFileOffset / (int(WP, MPI_Offset_kind) *                &
            product(int(controllerPatch%globalSize, MPI_OFFSET_KIND)) *                      &
            size(gradientBuffer, 2)))
       gradientFileOffset = 0
    else
       gradientFileOffset = gradientFileOffset - nBytesToRead
    end if

    call MPI_File_set_view(mpiFileHandle, gradientFileOffset, MPI_REAL_WP,                   &
         mpiRealSubarrayType, "native", mpiInfo, ierror)

    dataSize = controllerPatch%nPatchPoints * size(gradientBuffer, 2) * iGradientBuffer
    call MPI_File_read_all(mpiFileHandle, gradientBuffer, dataSize, MPI_REAL_WP,             &
         MPI_STATUS_IGNORE, ierror)

    call MPI_File_close(mpiFileHandle, ierror)

    call MPI_Type_free(mpiRealSubarrayType, ierror)

    return
  end subroutine load_sensitivity_gradient


  ! Save the sensitivity gradient on the controller patch
  ! -----------------------------------------------------
  subroutine save_sensitivity_gradient

    ! External modules
    use parallel

    implicit none

    ! Local variables
    integer :: extent(6), arrayOfSizes(5), arrayOfSubsizes(5), arrayOfStarts(5),             &
         mpiRealSubarrayType, mpiFileHandle, dataSize, ierror

    if (controllerPatch%comm .eq. MPI_COMM_NULL .or. iGradientBuffer .eq. 0) return

    extent(1) = controllerPatch%iMin; extent(2) = controllerPatch%iMax
    extent(3) = controllerPatch%jMin; extent(4) = controllerPatch%jMax
    extent(5) = controllerPatch%kMin; extent(6) = controllerPatch%kMax

    arrayOfSizes(1:3) = controllerPatch%globalSize
    arrayOfSizes(4) = size(gradientBuffer, 2)
    arrayOfSizes(5) = iGradientBuffer
    arrayOfSubsizes(1:3) = controllerPatch%localSize
    arrayOfSubsizes(4) = size(gradientBuffer, 2)
    arrayOfSubsizes(5) = iGradientBuffer
    arrayOfStarts(1:3) = controllerPatch%offset - extent(1::2) + 1
    arrayOfStarts(4:5) = 0
    call MPI_Type_create_subarray(5, arrayOfSizes, arrayOfSubsizes, arrayOfStarts,           &
         MPI_ORDER_FORTRAN, MPI_REAL_WP, mpiRealSubarrayType, ierror)
    call MPI_Type_commit(mpiRealSubarrayType, ierror)

    call MPI_File_open(controllerPatch%comm, trim(mpiiofs)//trim(gradientFilename)//char(0), &
         MPI_MODE_WRONLY, mpiInfo, mpiFileHandle, ierror)

    call MPI_File_set_view(mpiFileHandle, gradientFileOffset, MPI_REAL_WP,                   &
         mpiRealSubarrayType, "native", mpiInfo, ierror)

    dataSize = controllerPatch%nPatchPoints * size(gradientBuffer, 2) * iGradientBuffer
    call MPI_File_write_all(mpiFileHandle, gradientBuffer, dataSize, MPI_REAL_WP,            &
         MPI_STATUS_IGNORE, ierror)

    gradientFileOffset = gradientFileOffset + int(WP, MPI_Offset_kind) *                     &
         product(int(controllerPatch%globalSize, MPI_OFFSET_KIND)) *                         &
         size(gradientBuffer, 2) * iGradientBuffer

    call MPI_File_close(mpiFileHandle, ierror)

    call MPI_Type_free(mpiRealSubarrayType, ierror)

    return
  end subroutine save_sensitivity_gradient


  ! Load individual cost sensitivities
  ! ----------------------------------
  subroutine load_cost_sensitivity

    ! External modules
    use parallel
    use fileio
    use time_info

    implicit none

    ! Local variables
    integer :: i, fileUnit, iostat, iBuffer
    real(WP), dimension(nControlParameters) :: rBuffer
    character(len = str_medium), dimension(nControlParameters) :: cBuffer
    character(len = str_medium) :: filename

    ! This is for non space-time gradients only
    if (spaceTimeGradient) return

    write(filename, '(A)') "cost_sensitivity.txt"

    ! Only the root process reads/writes
    if (iRank .eq. iRoot) then
       fileUnit = iOpen()
       open(unit = fileUnit, file = trim(filename), action = 'read', status = 'unknown',     &
            iostat = iostat)
    end if

    call parallel_bc(iostat)
    if (iostat .ne. 0)                                                                       &
         call die('load_cost_sensitivity: Failed to open cost sensitivity file!')

    ! Read the file
    if (iRank .eq. iRoot) then
       read(fileUnit, *, iostat = iostat)
       do i = 1, nControlParameters
          read(fileUnit, *, iostat = iostat) iBuffer, cBuffer(i), rBuffer(i),                &
               currentSensitivity(i)
       end do
    end if
    call parallel_bc(iostat)
    if (iostat .ne. 0) call die(trim(filename) //                                           &
         ': Failed to read cost sensitivity from file!')
    call parallel_bc(currentSensitivity)

    ! Close the file
    if (iRank .eq. iRoot) close(iclose(fileUnit))

    return
  end subroutine load_cost_sensitivity


  ! Save individual cost sensitivities
  ! ------------------------------------
  subroutine save_cost_sensitivity

    ! External modules
    use parallel
    use fileio
    use time_info

    implicit none

    ! Local variables
    integer :: i, fileUnit, iostat
    character(len = str_medium) :: filename

    ! This is for non space-time gradients only
    if (spaceTimeGradient) return

    ! Only the root process reads/writes
    if (iRank .eq. iRoot) then
       write(filename, '(A)') "cost_sensitivity.txt"
       fileUnit = iOpen()
       open(unit = fileUnit, file = trim(filename), action = 'write', status = 'unknown',    &
            iostat = iostat)
    end if

    call parallel_bc(iostat)
    if (iostat .ne. 0)                                                                       &
         call die('save_cost_sensitivity: Failed to open cost sensitivity file!')

    if (iRank .eq. iRoot) then
       ! Write the file
       write(fileUnit, '(A6,3A24)') 'Number', 'Name', 'Parameter', 'Sensitivity'
       do i = 1, nControlParameters
          write(fileUnit, '(I6,1X,A23,2(1X,SP,(SP,ES23.15E3)))') i,                          &
               trim(sensitivityParameter(i)), baselineValue(i), currentSensitivity(i)
       end do

       ! Close the file
       close(iclose(fileUnit))
    end if

    return
  end subroutine save_cost_sensitivity

end module controller


! ==================== !
! Setup the controller !
! ==================== !
subroutine controller_setup

  ! Internal modules
  use controller

  ! External modules
  use string
  use parser
  use parallel
  use simulation_flags
  use solver_options

  implicit none

  ! Local varianbles
  integer :: i, j, k, gridIndex, gradientBufferSize
  real(WP) :: steepness, fraction
  character(len = str_medium) :: type, key, param
  character(len = str_short) :: xyz(3)

  ! Number of control parameters
  call parser_read('number of control parameters', nControlParameters, 1)
  if (nControlParameters .le. 0)                                                             &
       call die('controller_setup: number of control parameters must be positive!')

  ! By default, there is no constaint
  nConstraint = 0

  ! Initialize the sensitivity arrays
  allocate(instantaneousSensitivity(nControlParameters))
  allocate(currentSensitivity(nControlParameters))
  allocate(controlGradient(nControlParameters))
  instantaneousSensitivity = 0.0_WP
  currentSensitivity = 0.0_WP
  controlGradient = 0.0_WP

  ! Setup the sensitivity parameter names
  allocate(sensitivityParameter(nControlParameters))
  do i = 1, nControlParameters
     write(key, '(A,I1.1)') "sensitivity parameter ", i
     write(param, '(A,I1.1)') "param_", i
     call parser_read(trim(key), sensitivityParameter(i), trim(param))
  end do

  ! By default the sensitivity gradient is not space-time dependent
  spaceTimeGradient = .false.

  ! By default the controller attempts to minimize the cost functional
  gradientDirection = -1

  ! Connect the controller patch
  allocate(controllerPatch)
  do i = 1, nPatches
     if (patches(i)%patchType .eq. ACTUATOR) exit
  end do
  if (i .gt. nPatches) return
  controllerPatch => patches(i)

  ! Ensure only one controller is specified
  do j = i + 1, nPatches
     if (patches(j)%patchType .eq. ACTUATOR)                                                 &
          call die('controller_setup: only one actuator patch can be specified!')
  end do

  ! Verify the controller patch
  call verify_controller_patch(controllerPatch)

  ! Setup the control mollifier
  allocate(controlMollifier(nGridPoints, 1))
  controlMollifier = 0.0_WP
  do k = controllerPatch%iStart(3), controllerPatch%iEnd(3)
     do j = controllerPatch%iStart(2), controllerPatch%iEnd(2)
        do i = controllerPatch%iStart(1), controllerPatch%iEnd(1)
           gridIndex = i - gridOffset(1) + localGridSize(1) *                                &
                (j - 1 - gridOffset(2) + localGridSize(2) *                                  &
                (k - 1 - gridOffset(3)))
           controlMollifier(gridIndex, 1) = 1.0_WP
        end do
     end do
  end do
  call parser_read('control mollifier steepness', steepness, 20.0_WP)
  call parser_read('control mollifier fraction', fraction, 0.1_WP)
  xyz(1) = 'x'; xyz(2) = 'y'; xyz(3) = 'z'
  do i = 1, nDimensions
     call parser_read('control support type in ' // trim(xyz(i)), type, 'none')
     select case (trim(type))
     case ('tanh')
        call patch_tanh_support(controllerPatch, i, controlMollifier(:,1),                   &
             steepness, fraction)

     case ('cubic')
        call patch_cubic_bspline_support(controllerPatch, i, controlMollifier(:,1))

     case ('none')

     case default
        call die("Unknown control support type '" // trim(type) // "'")

     end select
  end do
  call normalize_control_mollifier

  ! Get the controller type and setup
  call parser_read('controller type', type, '')

  select case (trim(type))

  case ('thermal actuator')
     controllerType = THERMAL_ACTUATOR
     call thermal_controller_setup

  case ('momentum actuator')
     controllerType = MOMENTUM_ACTUATOR
     call momentum_controller_setup

  case ('fuel actuator')
     controllerType = FUEL_ACTUATOR
     call fuel_controller_setup

  case ('ignition actuator')
     controllerType = IGNITION_ACTUATOR
     call ignition_controller_setup

  case ('chemical actuator')
     controllerType = CHEMICAL_ACTUATOR
     call chemical_controller_setup

  case ('perturbation actuator')
     controllerType = PERTURBATION_ACTUATOR
     call perturbation_controller_setup

  case ('inital condition actuator')
     controllerType = GENERIC_IC_ACTUATOR
     call ic_controller_setup

  case default

     controllerType = NULL_ACTUATOR

  end select

  ! Initialize parameters for sensitivity I/O
  iGradientBuffer = 0
  gradientFileOffset = int(0, MPI_OFFSET_KIND)

  if (spaceTimeGradient) then ! ... setup the control forcing for space-time gradient

     if (nControlParameters .ne. 1) call die('controller_setup: number of control parameters &
          &must be 1 to output 3D gradient file!')

     write(gradientFilename, '(3A)') "gradient_", trim(controllerPatch%name), ".dat"

     call parser_read('gradient buffer size', gradientBufferSize, 20)
     if (controllerPatch%nPatchPoints .gt. 0) then
        allocate(gradientBuffer(controllerPatch%nPatchPoints, 1, gradientBufferSize))
        gradientBuffer = 0.0_WP

        allocate(controlForcing(controllerPatch%nPatchPoints, nUnknowns))
        controlForcing = 0.0_WP
     end if

  else ! ... setup the control forcing for a vector of parameters

     if (.not. predictionOnly) then
        allocate(controlForcing(nControlParameters, 1))
        controlForcing = 0.0_WP
     end if
     
  end if

  return
end subroutine controller_setup


! ====================== !
! Cleanup the controller !
! ====================== !
subroutine controller_cleanup

  ! Internal modules
  use controller

  implicit none

  if (allocated(baselineValue)) deallocate(baselineValue)
  if (allocated(controlGradient)) deallocate(controlGradient)
  if (allocated(sensitivityParameter)) deallocate(sensitivityParameter)
  if (allocated(controlForcing)) deallocate(controlForcing)
  if (allocated(gradientBuffer)) deallocate(gradientBuffer)
  if (allocated(instantaneousSensitivity)) deallocate(instantaneousSensitivity)
  if (allocated(currentSensitivity)) deallocate(currentSensitivity)
  if (allocated(controlMollifier)) deallocate(controlMollifier)

  if (associated(controllerPatch)) nullify(controllerPatch)

  gradientDirection = -1
  iGradientBuffer = 0
  gradientFileOffset = int(0, MPI_OFFSET_KIND)

  return
end subroutine controller_cleanup


! ===================================== !
! Update control parameters using the   !
! adjoint gradient and acutation amount !
! ===================================== !
subroutine controller_update_parameters(actuationAmount)

  ! Internal modules
  use controller

  ! External modules
  use simulation_flags

  implicit none

  ! Arguments
  real(WP), intent(in) :: actuationAmount

  if (predictionOnly .or. .not. associated(controllerPatch)) return

  ! Start the controller timer
  call timing_start('controller')

  select case (controllerType)

  case (THERMAL_ACTUATOR)
     call thermal_controller_update_parameters(actuationAmount)

  case (MOMENTUM_ACTUATOR)
     call momentum_controller_update_parameters(actuationAmount)

  case (FUEL_ACTUATOR)
     call fuel_controller_update_parameters(actuationAmount)

  case (IGNITION_ACTUATOR)
     call ignition_controller_update_parameters(actuationAmount)

  case (CHEMICAL_ACTUATOR)
     call chemical_controller_update_parameters(actuationAmount)

  case (PERTURBATION_ACTUATOR)
     call perturbation_controller_update_parameters(actuationAmount)

  case (GENERIC_IC_ACTUATOR)
     call ic_controller_update_parameters(actuationAmount)

  end select

  ! Stop the controller timer
  call timing_stop('controller')

  return
end subroutine controller_update_parameters


! ===================================== !
! Update initial conditions using the   !       
! adjoint gradient and acutation amount !
! ===================================== !
subroutine controller_update_ic(actuationAmount, controlIteration)

  ! Internal modules
  use controller

  ! External modules
  use simulation_flags

  implicit none

  ! Arguments
  integer, intent(in)  :: controlIteration
  real(WP), intent(in) :: actuationAmount

  if (predictionOnly .or. .not. associated(controllerPatch)) return

  select case (controllerType)

  case (PERTURBATION_ACTUATOR)
     ! Update parameters before updating initial conditions
     call perturbation_controller_update_parameters(actuationAmount)
     call perturbation_controller_update_ic(controlIteration)

  case (GENERIC_IC_ACTUATOR)
     call ic_controller_update_ic(actuationAmount)

  end select

  return
end subroutine controller_update_ic


! ========================================== !
! Add the control forcing source term to the ! 
! forward state if necessary                 !
! ========================================== !
subroutine controller_add_source(mode, source)

  ! Internal modules
  use controller

  ! External modules
  use simulation_flags
  use solver_options
  use geometry

  implicit none

  ! Arguments
  integer, intent(in) :: mode
  real(WP), dimension(nGridPoints, nUnknowns), intent(inout) :: source

  ! Local variables
  integer :: i, j, k, l, gridIndex, patchIndex

  if (.not. allocated(controlForcing) .or. .not. spaceTimeGradient .or.                      &
       mode .eq. ADJOINT .or. (mode .eq. FORWARD .and. predictionOnly)) return

  ! Start the controller timer
  call timing_start('controller')

  do l = 1, nUnknowns
     do k = controllerPatch%iStart(3), controllerPatch%iEnd(3)
        do j = controllerPatch%iStart(2), controllerPatch%iEnd(2)
           do i = controllerPatch%iStart(1), controllerPatch%iEnd(1)
              gridIndex = i - gridOffset(1) + localGridSize(1) *                             &
                   (j - 1 - gridOffset(2) + localGridSize(2) *                               &
                   (k - 1 - gridOffset(3)))
              patchIndex = i - controllerPatch%offset(1) + controllerPatch%localSize(1) *    &
                   (j - 1 - controllerPatch%offset(2) + controllerPatch%localSize(2) *       &
                   (k - 1 - controllerPatch%offset(3)))

              source(gridIndex, l) = source(gridIndex, l) +                                  &
                   controlMollifier(gridIndex, 1) * controlForcing(patchIndex, l)
           end do
        end do
     end do
  end do

  ! Stop the controller timer
  call timing_stop('controller')

  return
end subroutine controller_add_source


! ======================================== !
! Compute the adjoint sensitivity gradient !
! ======================================== !
subroutine controller_compute_sensitivity

  ! Internal modules
  use controller

  ! External modules
  use simulation_flags

  implicit none

  if (predictionOnly .or. .not. associated(controllerPatch)) return

  ! Start the controller timer
  call timing_start('controller')

  select case (controllerType)

  case (THERMAL_ACTUATOR)
     call thermal_controller_compute_sensitivity

  case (MOMENTUM_ACTUATOR)
     call momentum_controller_compute_sensitivity

  case (FUEL_ACTUATOR)
     call fuel_controller_compute_sensitivity

  case (IGNITION_ACTUATOR)
     call ignition_controller_compute_sensitivity

  case (CHEMICAL_ACTUATOR)
     call chemical_controller_compute_sensitivity

  case (PERTURBATION_ACTUATOR)
     call perturbation_controller_compute_sensitivity

  case (GENERIC_IC_ACTUATOR)
     call ic_controller_compute_sensitivity

  end select

  ! Save the gradient of necessary
  if (spaceTimeGradient .and. iGradientBuffer .eq. size(gradientBuffer, 3)) then
     call save_sensitivity_gradient
     iGradientBuffer = 0
  end if

  ! Stop the controller timer
  call timing_stop('controller')

  return
end subroutine controller_compute_sensitivity


! ======================================================== !
! Compute the adjoint sensitivity of the initial condition !
! ======================================================== !
subroutine controller_compute_sensitivity_ic

  ! Internal modules
  use controller

  ! External modules
  use simulation_flags

  implicit none

  if (predictionOnly .or. .not. associated(controllerPatch)) return

  ! Start the controller timer
  call timing_start('controller')

  instantaneousSensitivity = 0.0_WP

  select case (controllerType)

  case (PERTURBATION_ACTUATOR)
     call perturbation_controller_compute_sensitivity_ic

  case (GENERIC_IC_ACTUATOR)
     call ic_controller_compute_sensitivity_ic

  end select

  ! Stop the controller timer
  call timing_stop('controller')

  return
end subroutine controller_compute_sensitivity_ic


! ================================================= !
! Call controller hooks before time marching starts !
! ================================================= !
subroutine controller_hook_before_timemarch(mode)

  ! Internal modules
  use controller

  ! External modules
  use parallel
  use fileio
  use simulation_flags
  use solver_options

  implicit none

  ! Arguments
  integer, intent(in) :: mode

  ! Local variables
  integer :: stat, fileUnit, mpiFileHandle, ierror
  logical :: fileExists

  if (predictionOnly .or. .not. spaceTimeGradient .or.                                       &
       controllerPatch%comm .eq. MPI_COMM_NULL) return

  select case (mode)

  case (FORWARD)
     if (iRank .eq. controllerPatch%masterRank) then
        inquire(file = trim(gradientFilename), exist = fileExists)
        if (.not. fileExists) call die("No gradient information available on patch '" //     &
             trim(controllerPatch%name) // "'")
     end if
     iGradientBuffer = size(gradientBuffer, 3) + 1
     call MPI_File_open(controllerPatch%comm, trim(mpiiofs)//trim(gradientFilename)//char(0),&
          MPI_MODE_WRONLY, mpiInfo, mpiFileHandle, ierror)
     call MPI_File_get_size(mpiFileHandle, gradientFileOffset, ierror)
     call MPI_File_close(mpiFileHandle, ierror)

  case (ADJOINT)
     if (iRank .eq. controllerPatch%masterRank) then
        fileUnit = iOpen()
        open(unit = fileUnit, file = trim(gradientFilename), iostat = stat, status = 'old')
        if (stat .eq. 0) close(iClose(fileUnit), status = 'delete')
        fileUnit = iOpen()
        open(unit = fileUnit, file = trim(gradientFilename), action = 'write',               &
             status = 'unknown')
        close(iClose(fileUnit))
     end if
     call MPI_Barrier(controllerPatch%comm, ierror)
     iGradientBuffer = 0
     gradientFileOffset = int(0, MPI_OFFSET_KIND)

  end select

  return
end subroutine controller_hook_before_timemarch


! ================================================ !
! Call controller hooks after time marching starts !
! ================================================ !
subroutine controller_hook_after_timemarch(mode)

  ! Internal modules
  use controller

  ! External modules
  use simulation_flags
  use solver_options

  implicit none

  ! Arguments
  integer, intent(in) :: mode

  ! Save individual sensitivities if not space-time gradient
  if (.not. spaceTimeGradient .and. mode.eq.ADJOINT)  call save_cost_sensitivity

  if (predictionOnly .or. .not. spaceTimeGradient .or.                                       &
       controllerPatch%comm .eq. MPI_COMM_NULL) return

  select case (mode)

  case (ADJOINT)
     call save_sensitivity_gradient

  end select

  return
end subroutine controller_hook_after_timemarch


! ========================= !
! Set the parameters bounds !
! ========================= !
subroutine controller_parameter_bound(lowerBound, upperBound, scale)

  ! Internal modules
  use controller

  implicit none

  ! Arguments
  real(WP), dimension(nControlParameters), intent(out) :: lowerBound, upperBound, scale

  ! Give default value
  lowerBound = -huge(1.0_WP)
  upperBound =  huge(1.0_WP)
  scale      =  1.0_WP

  select case (controllerType)

  case (PERTURBATION_ACTUATOR)
     call perturbation_controller_parameter_bound(lowerBound, upperBound, scale)

  case (CHEMICAL_ACTUATOR)
     call chemical_controller_parameter_bound(lowerBound, upperBound, scale)

  end select

  return
end subroutine controller_parameter_bound


! ============================ !
! Compute constraint functions !
! ============================ !
subroutine controller_constraint(controlParameters, functionValue)

  ! Internal modules
  use controller

  implicit none

  ! Arguments
  real(WP), dimension(nControlParameters), intent(in) :: controlParameters
  real(WP), dimension(nConstraint), intent(out) :: functionValue

  if (nConstraint.le.0) return

  select case (controllerType)

  case (PERTURBATION_ACTUATOR)
     call perturbation_controller_constraint(controlParameters, functionValue)

  end select

  return
end subroutine controller_constraint


! ======================================== !
! Compute gradient of constraint functions !
! ======================================== !
subroutine controller_constraint_gradient(controlParameters, gradientValue)

  ! Internal modules
  use controller

  implicit none

  ! Arguments
  real(WP), dimension(nControlParameters), intent(in) :: controlParameters
  real(WP), dimension(nConstraint,nControlParameters), intent(out) :: gradientValue

  if (nConstraint.le.0) return

  select case (controllerType)

  case (PERTURBATION_ACTUATOR)
     call perturbation_controller_constraint_gradient(controlParameters, gradientValue)

  end select

  return
end subroutine controller_constraint_gradient
