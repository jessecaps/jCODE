module slip

  ! External modules
  use boundary
  use impenetrable, only : t_ImpenetrableWall
  use precision

  implicit none

  integer :: nSlipWalls

  type, private :: t_SlipWall
     type(t_ImpenetrableWall) :: impenetrableWallData
     real(WP) :: viscousPenaltyAmount
     real(WP), allocatable, dimension(:,:,:) :: viscousFluxes
  end type t_SlipWall

  type(t_SlipWall), allocatable :: slipWallData(:)
  type(t_Patch), pointer :: slipWallPatch(:)

contains

  subroutine setup_slip_patch(patch, bc)

    ! External modules
    use parser
    use simulation_flags
    use solver_options
    use first_derivative
    use impenetrable, only : setup_impenetrable_patch 

    implicit none

    ! Arguments
    type(t_Patch), intent(inout) :: patch
    type(t_SlipWall), intent(inout) :: bc

    ! Verify the patch type
    call verify_slip_patch(patch)

    ! Slip wall requires impenetrable wall data
    call setup_impenetrable_patch(patch, bc%impenetrableWallData)

    ! Viscous penalty amount
    if (useViscosity) then
       bc%viscousPenaltyAmount = 1.0_WP /                                                    &
            firstDerivative(abs(patch%normalDirection))%normBoundary(1)
       bc%viscousPenaltyAmount = sign(bc%viscousPenaltyAmount, real(patch%normalDirection,WP))

       ! Setup memory for storing viscous fluxes on the patch
       if (patch%nPatchPoints .gt. 0) then
          allocate(bc%viscousFluxes(patch%nPatchPoints, nUnknowns, nDimensions))
       end if
    else
       bc%viscousPenaltyAmount = 0.0_WP
    end if

    return
  end subroutine setup_slip_patch


  subroutine cleanup_slip_patch(patch, bc)

    ! External modules
    use impenetrable, only : cleanup_impenetrable_patch

    implicit none

    ! Arguments
    type(t_Patch), intent(inout) :: patch
    type(t_SlipWall), intent(inout) :: bc

    if (allocated(bc%viscousFluxes)) deallocate(bc%viscousFluxes)
    call cleanup_impenetrable_patch(patch, bc%impenetrableWallData)

    return
  end subroutine cleanup_slip_patch


  subroutine verify_slip_patch(patch)

    ! External modules
    use simulation_flags
    use geometry

    implicit none

    ! Arguments
    type(t_Patch), intent(in) :: patch

    if (patch%patchType .ne. SAT_SLIP_WALL)                                                  &
         call die('verify_slip_patch: patch type mismatch for SAT_SLIP_WALL')

    return
  end subroutine verify_slip_patch


  subroutine slip_patch_forward(patch, bc, source)

    ! External modules
    use solver_options
    use geometry
    use grid
    use grid_patch
    use state
    use impenetrable, only : impenetrable_patch_forward

    implicit none

    ! Arguments
    type(t_Patch), intent(inout) :: patch
    type(t_SlipWall), intent(inout) :: bc
    real(WP), dimension(:,:), intent(inout) :: source

    ! Local variables
    integer :: i, j, k, direction, gridIndex, patchIndex
    real(WP), allocatable :: metricsAlongNormalDirection(:), unitNormal(:), tangent(:)

    ! Enforce the no-penetration condition
    call impenetrable_patch_forward(patch, bc%impenetrableWallData, source)

    if (.not. useViscosity) return

    direction = abs(patch%normalDirection)

    allocate(metricsAlongNormalDirection(nDimensions))
    allocate(unitNormal(nDimensions))
    allocate(tangent(nUnknowns)); tangent = 1.0_WP
    
    do k = patch%iStart(3), patch%iEnd(3)
       do j = patch%iStart(2), patch%iEnd(2)
          do i = patch%iStart(1), patch%iEnd(1)
             gridIndex = i - gridOffset(1) + localGridSize(1) *                              &
                  (j - 1 - gridOffset(2) + localGridSize(2) *                                &
                  (k - 1 - gridOffset(3)))
             patchIndex = i - patch%offset(1) + patch%localSize(1) *                         &
                  (j - 1 - patch%offset(2) + patch%localSize(2) *                            &
                  (k - 1 - patch%offset(3)))
             if (allocated(patch%hole)) then
                if (patch%hole(patchIndex) .eq. 1) cycle
             end if

             metricsAlongNormalDirection =                                                   &
                  metrics(gridIndex,1+nDimensions*(direction-1):nDimensions*direction)
             unitNormal = metricsAlongNormalDirection /                                      &
                  sqrt(sum(metricsAlongNormalDirection ** 2))
             tangent(2:nDimensions+1) = 1.0_WP - unitNormal
             if (twoWayCoupling) tangent = tangent * volumeFraction(gridIndex, 1)

             ! Zero flux condition (penalize tangent directions)
             source(gridIndex,:) = source(gridIndex,:) + bc%viscousPenaltyAmount *           &
                  jacobian(gridIndex, 1) * matmul(bc%viscousFluxes(patchIndex,:,:),          &
                  metricsAlongNormalDirection) * tangent

          end do
       end do
    end do

    deallocate(metricsAlongNormalDirection, unitNormal, tangent)

    return
  end subroutine slip_patch_forward


  subroutine slip_patch_adjoint(patch, bc, source)

    ! External modules
    use simulation_flags
    use solver_options
    use geometry
    use grid
    use state
    use state_jacobian
    use impenetrable, only : impenetrable_patch_adjoint
    use first_derivative

    implicit none

    ! Arguments
    type(t_Patch), intent(inout) :: patch
    type(t_SlipWall), intent(inout) :: bc
    real(WP), dimension(:,:), intent(inout) :: source

    ! Local variables
    integer :: i, j, k, l, direction, gridIndex, patchIndex
    real(WP), allocatable :: metricsAlongNormalDirection(:), metricsAlongDirection2(:),      &
         localViscousFluxJacobian(:,:), localViscousFluxJacobian2(:,:), localMassFraction(:),&
         localEnthalpyFlux(:), localSpeciesFlux(:,:), temp1(:,:,:), temp2(:,:),              &
         unitNormal(:), tangent(:)

    ! Apply no-penetration conditions
    call impenetrable_patch_adjoint(patch, bc%impenetrableWallData, source)

    if (.not. useViscosity) return

    direction = abs(patch%normalDirection)

    allocate(metricsAlongNormalDirection(nDimensions))
    allocate(metricsAlongDirection2(nDimensions))
    allocate(localViscousFluxJacobian(nUnknowns, nUnknowns))
    allocate(localViscousFluxJacobian2(nUnknowns-1, nUnknowns-1))
    if (nSpecies .gt. 0) then
       allocate(localMassFraction(nSpecies))
       allocate(localSpeciesFlux(nSpecies,nDimensions))
       if (equationOfState .eq. IDEAL_GAS_MIXTURE) allocate(localEnthalpyFlux(nDimensions))
    end if
    allocate(temp1(nGridPoints, nUnknowns - 1, nDimensions)); temp1 = 0.0_WP
    allocate(unitNormal(nDimensions))
    allocate(tangent(nUnknowns)); tangent = 1.0_WP

    do k = patch%iStart(3), patch%iEnd(3)
       do j = patch%iStart(2), patch%iEnd(2)
          do i = patch%iStart(1), patch%iEnd(1)
             gridIndex = i - gridOffset(1) + localGridSize(1) *                              &
                  (j - 1 - gridOffset(2) + localGridSize(2) *                                &
                  (k - 1 - gridOffset(3)))
             patchIndex = i - patch%offset(1) + patch%localSize(1) *                         &
                  (j - 1 - patch%offset(2) + patch%localSize(2) *                            &
                  (k - 1 - patch%offset(3)))
             if (allocated(patch%hole)) then
                if (patch%hole(patchIndex) .eq. 1) cycle
             end if

             metricsAlongNormalDirection =                                                   &
                  metrics(gridIndex,1+nDimensions*(direction-1):nDimensions*direction)
             unitNormal = metricsAlongNormalDirection /                                      &
                  sqrt(sum(metricsAlongNormalDirection ** 2))
             tangent(2:nDimensions+1) = 1.0_WP - unitNormal

             if (nSpecies .gt. 0) then
                localMassFraction = massFraction(gridIndex,:)
                localSpeciesFlux = speciesFlux(gridIndex,:,:)
                if (equationOfState .eq. IDEAL_GAS_MIXTURE)                               &
                     localEnthalpyFlux = enthalpyFlux(gridIndex,:)
             end if
             call compute_first_partial_viscous_jacobian(conservedVariables(gridIndex,:), &
                  metricsAlongNormalDirection, stressTensor(gridIndex,:),                 &
                  heatFlux(gridIndex,:), localEnthalpyFlux, localSpeciesFlux,             &
                  localViscousFluxJacobian, specificVolume(gridIndex,1),                  &
                  velocity(gridIndex,:), temperature(gridIndex,1), localMassFraction)
             source(gridIndex,:) = source(gridIndex,:) - bc%viscousPenaltyAmount *        &
                  jacobian(gridIndex, 1) * matmul(transpose(localViscousFluxJacobian),    &
                  adjointVariables(gridIndex,:)) * tangent

             if (useContinuousAdjoint) cycle

             do l = 1, nDimensions
                metricsAlongDirection2 =                                                  &
                     metrics(gridIndex,1+nDimensions*(l-1):nDimensions*l)
                call compute_second_partial_viscous_jacobian(velocity(gridIndex,:),       &
                     dynamicViscosity(gridIndex,1),                                       &
                     secondCoefficientOfViscosity(gridIndex,1),                           &
                     thermalDiffusivity(gridIndex,1), massDiffusivity(gridIndex,:),       &
                     temperature(gridIndex,1), jacobian(gridIndex,1),                     &
                     metricsAlongNormalDirection, metricsAlongDirection2,                 &
                     localViscousFluxJacobian2)
                temp1(gridIndex,:,l) = - matmul(transpose(localViscousFluxJacobian2),     &
                     adjointVariables(gridIndex,2:nUnknowns)) * tangent
             end do

          end do
       end do
    end do

    deallocate(metricsAlongNormalDirection, tangent, unitNormal)
    if (allocated(metricsAlongDirection2)) deallocate(metricsAlongDirection2)
    if (allocated(localViscousFluxJacobian)) deallocate(localViscousFluxJacobian)
    if (allocated(localViscousFluxJacobian2)) deallocate(localViscousFluxJacobian2)
    if (allocated(localMassFraction)) deallocate(localMassFraction)
    if (allocated(localSpeciesFlux)) deallocate(localSpeciesFlux)
    if (allocated(localEnthalpyFlux)) deallocate(localEnthalpyFlux)

    if (useContinuousAdjoint .or. .not. useViscosity) return

    allocate(temp2(nGridPoints, nUnknowns - 1)); temp2 = 0.0_WP

    do i = 1, nDimensions
       call adjoint_first_derivative_apply(i, temp1(:,:,i))
    end do
    temp2 = sum(temp1, dim = 3) !... divergence of the adjoint flux

    deallocate(temp1)  !... no longer needed

    select case (equationOfState)

    case (IDEAL_GAS)

       temp2(:,nDimensions+1) = ratioOfSpecificHeats * specificVolume(:,1) *                 &
            temp2(:,nDimensions+1)
       do i = 1, nDimensions
          temp2(:,i) = specificVolume(:,1) * temp2(:,i) - velocity(:,i) *                    &
               temp2(:,nDimensions+1)
       end do
       do k = 1, nSpecies
          temp2(:,nDimensions+1+k) = specificVolume(:,1) * temp2(:,nDimensions+1+k)
       end do

    case (IDEAL_GAS_MIXTURE)

       temp2(:,nDimensions+1) = ratioOfSpecificHeats * specificVolume(:,1) *                 &
            mixtureMolecularWeight(:,1) * temp2(:,nDimensions+1)
       do i = 1, nDimensions
          temp2(:,i) = specificVolume(:,1) * temp2(:,i) - velocity(:,i) *                    &
               temp2(:,nDimensions+1)
       end do
       do k = 1, nSpecies
          temp2(:,nDimensions+1+k) = specificVolume(:,1) * temp2(:,nDimensions+1+k) +        &
               temperature(:,1) * (molecularWeightInverse(nSpecies+1) -                      &
               molecularWeightInverse(k)) * temp2(:,nDimensions+1) /                         &
               ratioOfSpecificHeats
       end do

    end select

    do i = 1, size(temp2, 2)
       temp2(:,i) = bc%viscousPenaltyAmount * jacobian(:,1) * temp2(:,i)
    end do

    source(:,2:nUnknowns) = source(:,2:nUnknowns) + temp2
    source(:,1) = source(:,1) - specificVolume(:,1) * conservedVariables(:,nDimensions+2)    &
         * temp2(:,nDimensions+1) - sum(velocity * temp2(:,1:nDimensions), dim = 2)
    if (nSpecies .gt. 0) source(:,1) = source(:,1) -                                         &
         sum(massFraction * temp2(:,nDimensions+2:nUnknowns-1), dim = 2)

    deallocate(temp2)

    return
  end subroutine slip_patch_adjoint

end module slip


! =============================== !
! Setup the slip boundary patches !
! =============================== !
subroutine slip_setup

  ! Internal modules
  use slip

  ! External modules
  use parser

  implicit none

  ! Local variables
  integer :: i, j

  ! Find the number of boundaries of this type
  nSlipWalls = 0; j = 0
  do i = 1, nPatches
     if (patches(i)%patchType .eq. SAT_SLIP_WALL) then
        nSlipWalls = nSlipWalls + 1
        if (j .eq. 0) j = i
     end if
  end do
  if (nSlipWalls .eq. 0) return

  ! Allocate the slip wall type
  allocate(slipWallData(nSlipWalls))

  ! Connect the boundary patch
  slipWallPatch => patches(j:j+nSlipWalls-1)

  ! Setup the boundary conditions
  do i = 1, nSlipWalls
     call setup_slip_patch(slipWallPatch(i), slipWallData(i))
  end do

  return
end subroutine slip_setup


! ================================= !
! Cleanup the slip boundary patches !
! ================================= !
subroutine slip_cleanup

  ! Internal modules
  use slip

  ! Local variables
  integer :: i

  if (nSlipWalls .gt. 0) then
     do i = 1, nSlipWalls
        call cleanup_slip_patch(slipWallPatch(i), slipWallData(i))
     end do
     deallocate(slipWallData)
     nullify(slipWallPatch)
  end if

  nSlipWalls = 0

  return
end subroutine slip_cleanup


! ========================================== !
! Add the slip source during the forward run !
! ========================================== !
subroutine slip_forward(source)

  ! Internal modules
  use slip

  ! External modules
  use solver_options
  use geometry

  implicit none

  ! Arguments
  real(WP), dimension(nGridPoints,nUnknowns), intent(inout) :: source

  ! Local variables
  integer :: i

  do i = 1, nSlipWalls
     call slip_patch_forward(slipWallPatch(i), slipWallData(i), source)
  end do

  return
end subroutine slip_forward


! ========================================== !
! Add the slip source during the adjoint run !
! ========================================== !
subroutine slip_adjoint(source)

  ! Internal modules
  use slip

  ! External modules
  use solver_options
  use geometry

  implicit none

  ! Arguments
  real(WP), dimension(nGridPoints, nUnknowns), intent(inout) :: source

  ! Local variables
  integer :: i

  do i = 1, nSlipWalls
     call slip_patch_adjoint(slipWallPatch(i), slipWallData(i), source)
  end do

  return
end subroutine slip_adjoint


! ================================================= !
! Store the viscous fluxes on the slip wall patches !
! ================================================= !
subroutine slip_store_viscous_fluxes(fluxes)

  ! Internal modules
  use slip

  ! External modules
  use solver_options
  use geometry

  implicit none

  ! Arguments
  real(WP), dimension(nGridPoints, nUnknowns, nDimensions), intent(in) :: fluxes

  ! Local variables
  integer :: i

  do i = 1, nSlipWalls
     if (allocated(slipWallData(i)%viscousFluxes))                                           &
          call patch_collect(slipWallPatch(i), fluxes, slipWallData(i)%viscousFluxes)
  end do

  return
end subroutine slip_store_viscous_fluxes
