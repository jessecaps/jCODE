module excitation

  ! External modules
  use precision
  use grid_patch
  use state, only : targetState

  implicit none

  ! Excitation types
  integer, parameter ::                                                                      &
       EXCITATION_UNIFORM           = 0,                                                     &
       EXCITATION_AXIAL_FORCING     = 1,                                                     &
       EXCITATION_HELICAL_FORCING   = 2,                                                     &
       EXCITATION_HELICAL_FORCING_2 = 3,                                                     &
       EXCITATION_VORTEX_RING       = 4

  integer :: nExcitations
  real(WP) :: defaultExcitationAmplitude, defaultExcitationVelocity, defaultExcitationRadius,&
       defaultExcitationFrequency, defaultExcitationMomentumThickness
  real(WP), allocatable :: baseTargetMomentum(:,:), baseTargetEnergy(:)

  type, private :: t_Excitation
     integer :: type, nModes
     real(WP) :: frequency
     real(WP), allocatable :: strength(:,:), angularCoordinate(:), amplitudes(:), phases(:)
  end type t_Excitation

  type(t_Excitation), allocatable :: excitationData(:)
  type(t_Patch), pointer :: excitationPatch(:)

contains

  subroutine setup_excitation_patch(patch, bc)

    ! External modules
    use parser
    use math, only: pi
    use geometry

    implicit none

    ! Arguments
    type(t_Patch), intent(inout) :: patch
    type(t_Excitation), intent(inout) :: bc

    ! Local variables
    integer :: i, j, k, gridIndex, patchIndex, ierror
    real (WP) :: Uj, Rj, x0, x, y, z, r, ln2, Delta2, temp, maxVel, theta0, RjI, slope,      &
         amplitude
    character(len = str_medium) :: key

    ! Verify the excitation patch
    call verify_excitation_patch(patch)

    ! Read excitation type
    call parser_read(trim(patch%name) // ' excitation type', key)
    select case (key)

    case ('uniform')
       ! Uniformly distributed random fluctuations
       bc%type = EXCITATION_UNIFORM

    case ('axial')
       ! Axial forcing from Tyliszczak and Geurts 2014 (Flow Turbulence Combust)
       bc%type = EXCITATION_AXIAL_FORCING
       if (nDimensions.ne.3) call die('setup_excitation_patch: &
            &"axial" excitation only works in 3D for now')

    case ('helical')
       ! Helical forcing from Tyliszczak 2015 (International Journal of Heat and Fluid Flow)
       bc%type = EXCITATION_HELICAl_FORCING
       if (nDimensions.ne.3) call die('setup_excitation_patch: &
            &"helical" excitation only works in 3D for now')

    case ('helical 2')
       ! Helical forcing similar to Tyliszczak 2015 (Combustion and Flame)
       bc%type = EXCITATION_HELICAl_FORCING_2
       if (nDimensions.ne.3) call die('setup_excitation_patch: &
            &"helical 2" excitation only works in 3D for now')

    case ('vortex ring')
       ! Inflow forcing based on the unit vortex ring, for more information, look at:
       ! Bogey et al. 2003 (Theoret. Comput. Fluid Dynamics), Bogey and Bailly 2003 (AIAA)
       bc%type = EXCITATION_VORTEX_RING
       if (nDimensions.ne.3) call die('setup_excitation_patch: &
            &"vortex ring" excitation only works in 3D for now')

    case default
       call die('setup_excitation_patch: Unknown excitation type: ' // key)

    end select

    ! Read excitation parameters
    call parser_read(trim(patch%name) // ' excitation amplitude ratio', amplitude,           &
         defaultExcitationAmplitude)
    call parser_read(trim(patch%name) // ' excitation velocity', Uj,defaultExcitationVelocity)
    call parser_read(trim(patch%name) // ' excitation radius', Rj, defaultExcitationRadius)

    amplitude = amplitude * Uj

    select case (bc%type)

    case (EXCITATION_UNIFORM)

       allocate(bc%strength(patch%nPatchPoints, nDimensions)); bc%strength = 0.0_WP
       allocate(bc%amplitudes(1))
       bc%amplitudes(1) = amplitude

    case (EXCITATION_AXIAL_FORCING)

       allocate(bc%strength(patch%nPatchPoints, 1)); bc%strength = 0.0_WP

       call parser_read(trim(patch%name) // ' excitation momentum thickness', theta0,        &
            defaultExcitationMomentumThickness)
       call parser_read(trim(patch%name) // ' excitation Strouhal number', bc%frequency,     &
            defaultExcitationFrequency)

       slope = Rj / (4.0_WP * theta0)
       RjI = 1.0_WP / Rj
       bc%frequency = 2.0_WP * pi * bc%frequency * Uj / (2.0_WP * Rj)

       do k = patch%iStart(3), patch%iEnd(3)
          do j = patch%iStart(2), patch%iEnd(2)
             do i = patch%iStart(1), patch%iEnd(1)
                gridIndex = i - gridOffset(1) + localGridSize(1) *                           &
                     (j - 1 - gridOffset(2) + localGridSize(2) *                             &
                     (k - 1 - gridOffset(3)))
                patchIndex = i - patch%offset(1) + patch%localSize(1) *                      &
                     (j - 1 - patch%offset(2) + patch%localSize(2) *                         &
                     (k - 1 - patch%offset(3)))

                y = coordinates(gridIndex, 2)
                z = coordinates(gridIndex, 3)

                r = sqrt(y**2 + z**2)
                bc%strength(patchIndex, 1) = amplitude * 0.5_WP * (1.0_WP +                  &
                     tanh(slope * (Rj/r - RjI*r)))
             end do
          end do
       end do

    case (EXCITATION_HELICAL_FORCING)

       allocate(bc%strength(patch%nPatchPoints, 1)); bc%strength = 0.0_WP
       allocate(bc%angularCoordinate(patch%nPatchPoints)); bc%angularCoordinate = 0.0_WP

       call parser_read(trim(patch%name) // ' excitation momentum thickness', theta0,        &
            defaultExcitationMomentumThickness)
       call parser_read(trim(patch%name) // ' excitation Strouhal number', bc%frequency,     &
            defaultExcitationFrequency)

       slope = Rj / (4.0_WP * theta0)
       RjI = 1.0_WP / Rj
       bc%frequency = 2.0_WP * pi * bc%frequency * Uj / (2.0_WP * Rj)

       do k = patch%iStart(3), patch%iEnd(3)
          do j = patch%iStart(2), patch%iEnd(2)
             do i = patch%iStart(1), patch%iEnd(1)
                gridIndex = i - gridOffset(1) + localGridSize(1) *                           &
                     (j - 1 - gridOffset(2) + localGridSize(2) *                             &
                     (k - 1 - gridOffset(3)))
                patchIndex = i - patch%offset(1) + patch%localSize(1) *                      &
                     (j - 1 - patch%offset(2) + patch%localSize(2) *                         &
                     (k - 1 - patch%offset(3)))

                y = coordinates(gridIndex, 2)
                z = coordinates(gridIndex, 3)

                r = sqrt(y**2 + z**2)
                if (r.gt.0.0_WP) then
                   bc%angularCoordinate(patchIndex) = atan2(z,y)
                else
                   bc%angularCoordinate(patchIndex) = 0.0_WP
                end if
                bc%strength(patchIndex, 1) = amplitude * sin(pi * r * (0.5_WP*RjI)) *        &
                     0.5_WP * (1.0_WP + tanh(slope * (Rj/r - RjI*r)))
             end do
          end do
       end do

    case (EXCITATION_HELICAL_FORCING_2)

       allocate(bc%strength(patch%nPatchPoints, 1)); bc%strength = 0.0_WP
       allocate(bc%angularCoordinate(patch%nPatchPoints)); bc%angularCoordinate = 0.0_WP

       call parser_read(trim(patch%name) // ' excitation momentum thickness', theta0,        &
            defaultExcitationMomentumThickness)
       call parser_read(trim(patch%name) // ' excitation Strouhal number', bc%frequency,     &
            defaultExcitationFrequency)

       slope = Rj / (4.0_WP * theta0)
       RjI = 1.0_WP / Rj
       bc%frequency = 2.0_WP * pi * bc%frequency * Uj / (2.0_WP * Rj)

       do k = patch%iStart(3), patch%iEnd(3)
          do j = patch%iStart(2), patch%iEnd(2)
             do i = patch%iStart(1), patch%iEnd(1)
                gridIndex = i - gridOffset(1) + localGridSize(1) *                           &
                     (j - 1 - gridOffset(2) + localGridSize(2) *                             &
                     (k - 1 - gridOffset(3)))
                patchIndex = i - patch%offset(1) + patch%localSize(1) *                      &
                     (j - 1 - patch%offset(2) + patch%localSize(2) *                         &
                     (k - 1 - patch%offset(3)))

                y = coordinates(gridIndex, 2)
                z = coordinates(gridIndex, 3)

                r = sqrt(y**2 + z**2)
                bc%angularCoordinate(patchIndex) = 0.0_WP
                bc%strength(patchIndex, 1) = amplitude * sin(pi * y * (0.5_WP*RjI)) *        &
                     0.5_WP * (1.0_WP + tanh(slope * (Rj/r - RjI*r)))
             end do
          end do
       end do

    case (EXCITATION_VORTEX_RING)

       allocate(bc%strength(patch%nPatchPoints, 1)); bc%strength = 0.0_WP
       allocate(bc%angularCoordinate(patch%nPatchPoints)); bc%angularCoordinate = 0.0_WP

       bc%nModes = 8
       allocate(bc%amplitudes(0:bc%nModes)); bc%amplitudes = 0.0_WP
       allocate(bc%phases(0:bc%nModes)); bc%phases = 0.0_WP

       x0 = 0.8_WP * Rj
       ln2 = -log(2.0_WP)

       maxVel = 0.0_WP !... to store maximum magnitude of velocity

       do k = patch%iStart(3), patch%iEnd(3)
          do j = patch%iStart(2), patch%iEnd(2)
             do i = patch%iStart(1), patch%iEnd(1)
                gridIndex = i - gridOffset(1) + localGridSize(1) *                             &
                     (j - 1 - gridOffset(2) + localGridSize(2) *                               &
                     (k - 1 - gridOffset(3)))
                patchIndex = i - patch%offset(1) + patch%localSize(1) *                        &
                     (j - 1 - patch%offset(2) + patch%localSize(2) *                           &
                     (k - 1 - patch%offset(3)))

                x = coordinates(gridIndex, 1)
                y = coordinates(gridIndex, 2)
                z = coordinates(gridIndex, 3)

                r = sqrt(y**2 + z**2)
                if (r.gt.0.0_WP) then
                   bc%angularCoordinate(patchIndex) = atan2(z,y)
                   temp = 1.0_WP / (gridSpacing(gridIndex,2)**2 + gridSpacing(gridIndex,3)**2)
                   Delta2 = (x - x0)**2 + (r-Rj)**2
                   temp = 2.0_WP * Rj * sqrt(temp) * exp(ln2 * Delta2 * temp) / r
                else
                   bc%angularCoordinate(patchIndex) = 0.0_WP
                   temp = 0.0_WP
                end if

                bc%strength(patchIndex, 1) = temp * (r - Rj)
                bc%strength(patchIndex, 2) = temp * (x0 - x) *                               &
                     cos(bc%angularCoordinate(patchIndex))
                bc%strength(patchIndex, 3) = temp * (x0 - x) *                               &
                     sin(bc%angularCoordinate(patchIndex))

                maxVel = max(maxVel, sum(bc%strength(patchIndex, :))**2)
             end do
          end do
       end do

       call MPI_ALLREDUCE(MPI_IN_PLACE, maxVel, 1, MPI_REAL_WP, MPI_MAX, patch%comm, ierror)

       ! Normalize excitation velocity
       maxVel = 1.0_WP / sqrt(maxVel)
       bc%strength = bc%strength * maxVel

       ! Include excitation amplitude
       bc%strength = bc%strength * amplitude

    end select

    return
  end subroutine setup_excitation_patch


  subroutine cleanup_excitation_patch(patch, bc)

    implicit none

    ! Arguments
    type(t_Patch), intent(inout) :: patch
    type(t_Excitation), intent(inout) :: bc

    if (allocated(bc%strength)) deallocate(bc%strength)
    if (allocated(bc%angularCoordinate)) deallocate(bc%angularCoordinate)
    if (allocated(bc%amplitudes)) deallocate(bc%amplitudes)
    if (allocated(bc%phases)) deallocate(bc%phases)

    return
  end subroutine cleanup_excitation_patch


  ! Verify the excitation source patch is correct
  ! ---------------------------------------------
  subroutine verify_excitation_patch(patch)

    ! External modules
    use string
    use simulation_flags
    use geometry

    implicit none

    ! Arguments
    type(t_Patch), intent(in) :: patch

    ! Local variables
    integer :: i, n, extent(6)

    extent(1) = patch%iMin; extent(2) = patch%iMax
    extent(3) = patch%jMin; extent(4) = patch%jMax
    extent(5) = patch%kMin; extent(6) = patch%kMax

    n = nDimensions
    do i = 1, nDimensions
       if (extent((i-1)*2+1) .lt. 0 .or. extent((i-1)*2+2) .gt. globalGridSize(i) .or.       &
            extent((i-1)*2+1) .gt. extent((i-1)*2+2))                                        &
            call die("verify_excitation_patch: Invalid extent on '" //                       &
            trim(patch%name) // "'!")
       if (extent((i-1)*2+1) .eq. extent((i-1)*2+2)) n = n - 1
    end do

    return
  end subroutine verify_excitation_patch


  subroutine excitation_patch_forward(patch, bc, time, source)

    ! External modules
    use geometry
    use math, only : pi

    implicit none

    ! Arguments
    type(t_Patch), intent(inout) :: patch
    type(t_Excitation), intent(inout) :: bc
    real(WP), intent(in) :: time
    real(WP), dimension(:,:), intent(inout) :: source

    ! Local variables
    integer :: i, j, k, l, gridIndex, patchIndex
    real(WP) :: temp, oldV, newV

    ! Generate random numbers
    select case (bc%type)
    case (EXCITATION_UNIFORM)
       call random_number(bc%strength)
       bc%strength = 2.0_WP * bc%amplitudes(1) * (bc%strength - 0.5_WP)
    case (EXCITATION_VORTEX_RING)
       call random_number(bc%amplitudes); bc%amplitudes = 2.0_WP * (bc%amplitudes - 0.5_WP)
       call random_number(bc%phases); bc%phases = 2.0_WP * pi * bc%phases
    end select

    do k = patch%iStart(3), patch%iEnd(3)
       do j = patch%iStart(2), patch%iEnd(2)
          do i = patch%iStart(1), patch%iEnd(1)
             gridIndex = i - gridOffset(1) + localGridSize(1) *                                &
                  (j - 1 - gridOffset(2) + localGridSize(2) *                                  &
                  (k - 1 - gridOffset(3)))
             patchIndex = i - patch%offset(1) + patch%localSize(1) *                           &
                  (j - 1 - patch%offset(2) + patch%localSize(2) *                              &
                  (k - 1 - patch%offset(3)))

             ! Add excitation source to momentum equations
             select case (bc%type)

             case (EXCITATION_UNIFORM)
                do l = 1, nDimensions
                   source(gridIndex,l+1) = source(gridIndex,l+1) + bc%strength(patchIndex,l)
                end do

             case (EXCITATION_AXIAL_FORCING)
                oldV = targetState(gridIndex, 2) / targetState(gridIndex, 1)
                newV = oldV + bc%strength(patchIndex, 1) * sin(bc%frequency * time)

                ! Update the target momentum
                targetState(gridIndex, 2) = targetState(gridIndex, 1) * newV

                ! Correct the target energy to avoid changing the target pressure
                targetState(gridIndex, nDimensions+2) = targetState(gridIndex, nDimensions+2)&
                     + 0.5_WP * targetState(gridIndex, 1) * (newV**2 - oldV**2)

             case (EXCITATION_HELICAL_FORCING, EXCITATION_HELICAL_FORCING_2)
                oldV = targetState(gridIndex, 2) / targetState(gridIndex, 1)
                newV = oldV + bc%strength(patchIndex, 1) *                                   &
                     sin(bc%frequency * time + bc%angularCoordinate(patchIndex))

                ! Update the target momentum
                targetState(gridIndex, 2) = targetState(gridIndex, 1) * newV

                ! Correct the target energy to avoid changing the target pressure
                targetState(gridIndex, nDimensions+2) = targetState(gridIndex, nDimensions+2)&
                     + 0.5_WP * targetState(gridIndex, 1) * (newV**2 - oldV**2)

             case (EXCITATION_VORTEX_RING)
                temp = 0.0_WP
                do l = 0, bc%nModes
                   temp = temp + bc%amplitudes(l) * cos(l * bc%angularCoordinate(patchIndex) &
                        + bc%phases(l))
                end do
                do l = 1, nDimensions
                   source(gridIndex, l+1) = source(gridIndex, l+1) + temp *                  &
                        bc%strength(patchIndex, l)
                end do

             end select
          end do
       end do
    end do

    return
  end subroutine excitation_patch_forward

end module excitation


! ==================== !
! Setup the excitation !
! ==================== !
subroutine excitation_setup

  ! Internal modules
  use excitation

  ! External modules
  use parser
  use simulation_flags

  implicit none

  ! Local variables
  integer :: i, j, n, seed
  integer, allocatable :: seed_(:)

  ! Find the number of excitations
  nExcitations = 0; j = 0
  do i = 1, nPatches
     if (patches(i)%patchType .eq. EXCITATION_PATCH) then
        nExcitations = nExcitations + 1
        if (j .eq. 0) j = i
     end if
  end do
  if (nExcitations .eq. 0) return

  ! Allocate the excitation type
  allocate(excitationData(nExcitations))

  ! Read in default parameters
  call parser_read('default excitation amplitude ratio', defaultExcitationAmplitude, 0.1_WP)
  call parser_read('default excitation velocity', defaultExcitationVelocity, 0.1_WP)
  call parser_read('default excitation radius', defaultExcitationRadius, 0.5_WP)
  call parser_read('default excitation Strouhal Number', defaultExcitationFrequency, 1.0_WP)
  call parser_read('default excitation momentum thickness',                                  &
       defaultExcitationMomentumThickness, 0.05_WP * defaultExcitationRadius)

  ! Excitation random seed
  if (irank.eq.iroot) then
     call parser_read('excitation random seed', seed, -1)
     call random_seed(size = n)
     allocate(seed_(n))
     seed_ = seed
     call random_seed(put = seed_)
     deallocate(seed_)
  end if

  ! Connect the boundary patch
  excitationPatch => patches(j:j+nExcitations-1)

  ! Setup the boundary conditions
  do i = 1, nExcitations
     call setup_excitation_patch(excitationPatch(i), excitationData(i))
  end do

  do i = 1, nExcitations
     if (excitationData(i)%type .eq. EXCITATION_AXIAL_FORCING .or.                           &
          excitationData(i)%type .eq. EXCITATION_HELICAL_FORCING .or.                        &
          excitationData(i)%type .eq. EXCITATION_HELICAL_FORCING_2) then
        ! A target state is required!
        if (.not. useTargetState)                                                            &
             call die('excitation_setup: No target state available for excitation patch:' // &
             trim(excitationPatch(i)%name) // '!')

        allocate(baseTargetMomentum(nGridPoints, nDimensions))
        baseTargetMomentum(:,1:nDimensions) = targetState(:,2:nDimensions+1)

        allocate(baseTargetEnergy(nGridPoints))
        baseTargetEnergy = targetState(:,nDimensions+2)

        exit
     end if
  end do

  return
end subroutine excitation_setup


! ====================== !
! Cleanup the excitation !
! ====================== !
subroutine excitation_cleanup

  ! Internal modules
  use excitation

  implicit none

  ! Local variables
  integer :: i

  if (nExcitations .gt. 0) then
     do i = 1, nExcitations
        call cleanup_excitation_patch(excitationPatch(i), excitationData(i))
     end do
     deallocate(excitationData)
     nullify(excitationPatch)
  end if

  if (allocated(baseTargetMomentum)) deallocate(baseTargetMomentum)
  if (allocated(baseTargetEnergy)) deallocate(baseTargetEnergy)

  nExcitations = 0

  return
end subroutine excitation_cleanup


! ========================================= !
! Add the excitation during the forward run !
! ========================================= !
subroutine excitation_forward(time, source)

  ! Internal modules
  use excitation

  ! External modules
  use geometry
  use solver_options

  implicit none

  ! Arguments
  real(WP), intent(inout) :: time
  real(WP), dimension(nGridPoints, nUnknowns), intent(inout) :: source

  ! Local variables
  integer :: i

  if (allocated(baseTargetMomentum)) then
     targetState(:,2:nDimensions+1) = baseTargetMomentum(:,1:nDimensions)
     targetState(:,nDimensions+2) = baseTargetEnergy
  end if

  do i = 1, nExcitations
     call excitation_patch_forward(excitationPatch(i), excitationData(i), time, source)
  end do

  return
end subroutine excitation_forward
