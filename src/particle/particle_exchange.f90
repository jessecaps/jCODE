module particle_exchange

  ! External modules
  use particle
  use precision
  use simulation_flags
  use geometry
  use grid

  implicit none

  ! Global variables used during particle exchange
  integer :: nInterp, nOverlap(3,2), pStart(3), pEnd(3)
  integer, dimension(:,:,:), allocatable :: nPartInCell
  integer, dimension(:,:,:,:), allocatable :: partInCell
  integer :: vaporIndex
  real(WP) :: diffusionAmount, filterSize, volumeFactor
  real(WP), dimension(:,:), allocatable :: cartesianCoordinates, particleVelocity,           &
       granularTemperature
  real(WP), dimension(:,:,:,:), allocatable :: fluidDensity, fluidTemperature,               &
       fluidViscosity, fluidVelocity, fluidVorticity, fluidStress, fluidVolumeFraction,      &
       momentumSource, energySource, massSource, fluidMassFraction, fluidPressure,           &
       cartScalar, cartVector, fluidGradRho, fluidDivRhoU, fluidLevelset, fluidLevelsetNormal
  logical :: collisionsOn, useGranularTemperature, filterFluid

contains

  ! Fluid => particle interpolation routine
  ! ---------------------------------------
  subroutine interpolate_fluid_to_particle(particleIndex, particlePosition, density,         &
       temperature, viscosity, velocity, stress, volumeFraction, vorticity, massFraction,    &
       pressure, PTKE, gradRho, divRhoU, levelset, levelsetNormal)

    ! External modules
    use math

    implicit none

    ! Arguments
    integer, intent(in) :: particleIndex(3)
    real(WP), intent(in) :: particlePosition(3)
    real(WP), intent(out), optional :: density, temperature, viscosity, velocity(3),         &
         stress(3), volumeFraction, vorticity(3), massFraction, pressure, PTKE,              &
         gradRho(3), divRhoU, levelset, levelsetNormal(3)

    ! Local variables
    integer :: i, j, k, i1, j1, k1, i2, j2, k2, n, nInterp_(nDimensions)
    real(WP) :: interpCoeff(nInterp, nInterp, nInterp), weight(nInterp, nDimensions)
    
    select case (nDimensions)

    case (1)

       ! Get the interpolation points
       n = ceiling(0.5_WP * real(nInterp, WP))
       i1 = particleIndex(1) - n + 1
       if (mod(nInterp, 2) .gt. 0) then
          if (cartesianCoordinates(particleIndex(1)+1, 1) - particlePosition(1) .gt.         &
               particlePosition(1) - cartesianCoordinates(particleIndex(1),1)) i1 = i1 + 1
       end if
       i2 = i1 + nInterp - 1
       if (nOverlap(1,1) .eq. 0) i1 = max(i1, iStart(1))
       if (nOverlap(1,2) .eq. 0) i2 = min(i2, iEnd(1))
       nInterp_(1) = i2-i1+1

       ! Compute the interpolation coefficients
       call lagrange_basis(nInterp_(1), cartesianCoordinates(i1:i2,1), particlePosition(1),  &
            weight(1:nInterp_(1),1))

       ! Combine the interpolation coefficients to form a linear interpolation
       interpCoeff = 0.0_WP
       interpCoeff(1:nInterp_(1),1,1) = weight(1:nInterp_(1),1)

       ! Interpolate the fluid variables
       if (present(density)) density =                                                       &
            sum(interpCoeff(1:nInterp_(1),1,1) * fluidDensity(i1:i2,1,1,1))

       if (present(temperature)) temperature =                                               &
            sum(interpCoeff(1:nInterp_(1),1,1) * fluidTemperature(i1:i2,1,1,1))

       if (present(viscosity)) viscosity =                                                   &
            sum(interpCoeff(1:nInterp_(1),1,1) * fluidViscosity(i1:i2,1,1,1))

       if (present(massFraction)) massFraction =                                             &
            sum(interpCoeff(1:nInterp_(1),1,1) * fluidMassFraction(i1:i2,1,1,1))

       if (present(pressure)) pressure =                                                     &
            sum(interpCoeff(1:nInterp_(1),1,1) * fluidPressure(i1:i2,1,1,1))

       if (present(velocity)) then
          velocity = 0.0_WP
          velocity(1) = sum(interpCoeff(1:nInterp_(1),1,1) * fluidVelocity(i1:i2,1,1,1))
       end if

       if (present(gradRho)) then
          gradRho = 0.0_WP
          gradRho(1) = sum(interpCoeff(1:nInterp_(1),1,1) * fluidGradRho(i1:i2,1,1,1))
       end if

       if (present(divRhoU)) divRhoU =                                                       &
            sum(interpCoeff(1:nInterp_(1),1,1) * fluidDivRhoU(i1:i2,1,1,1))

       if (present(vorticity)) then
          vorticity = 0.0_WP
          do i = 1, 3
             vorticity(i) = sum(interpCoeff(1:nInterp_(1),1,1) * fluidVorticity(i1:i2,1,1,i))
          end do
       end if

       if (present(stress)) then
          stress = 0.0_WP
          stress(1) = sum(interpCoeff(1:nInterp_(1),1,1) * fluidStress(i1:i2,1,1,1))
       end if

       if (present(volumeFraction)) volumeFraction =                                         &
            sum(interpCoeff(1:nInterp_(1),1,1) * fluidVolumeFraction(i1:i2,1,1,1))

       if (present(levelset)) levelset =                                                     &
            sum(interpCoeff(1:nInterp_(1),1,1) * fluidLevelset(i1:i2,1,1,1))

       if (present(levelsetNormal)) then
          levelsetNormal = 0.0_WP
          levelsetNormal(1) = sum(interpCoeff(1:nInterp_(1),1,1) *                           &
               fluidLevelsetNormal(i1:i2,1,1,1))
       end if

    case (2)
       
       ! Get the interpolation points
       n = ceiling(0.5_WP * real(nInterp, WP))
       i1 = particleIndex(1) - n + 1
       j1 = particleIndex(2) - n + 1
       if (mod(nInterp, 2) .gt. 0) then
          if (cartesianCoordinates(particleIndex(1)+1, 1) - particlePosition(1) .gt.         &
               particlePosition(1) - cartesianCoordinates(particleIndex(1),1)) i1 = i1 + 1
          if (cartesianCoordinates(particleIndex(2)+1, 2) - particlePosition(2) .gt.         &
               particlePosition(2) - cartesianCoordinates(particleIndex(2),2)) j1 = j1 + 1
       end if
       i2 = i1 + nInterp - 1
       j2 = j1 + nInterp - 1
       if (nOverlap(1,1) .eq. 0) i1 = max(i1, iStart(1))
       if (nOverlap(1,2) .eq. 0) i2 = min(i2, iEnd(1))
       if (nOverlap(2,1) .eq. 0) j1 = max(j1, iStart(2))
       if (nOverlap(2,2) .eq. 0) j2 = min(j2, iEnd(2))
       nInterp_(1) = i2-i1+1
       nInterp_(2) = j2-j1+1

       ! Compute the interpolation coefficients
       call lagrange_basis(nInterp_(1), cartesianCoordinates(i1:i2,1), particlePosition(1),  &
            weight(1:nInterp_(1),1))
       call lagrange_basis(nInterp_(2), cartesianCoordinates(j1:j2,2), particlePosition(2),  &
            weight(1:nInterp_(2),2))

       ! Combine the interpolation coefficients
       interpCoeff = 0.0_WP
       do j = 1, nInterp_(2)
          do i = 1, nInterp_(1)
             interpCoeff(i,j,1) = weight(i,1) * weight(j,2)
          end do
       end do

       ! Interpolate the fluid variables
       if (present(density)) density = sum(                                                  &
            interpCoeff(1:nInterp_(1),1:nInterp_(2),1) * fluidDensity(i1:i2,j1:j2,1,1) )

       if (present(temperature)) temperature = sum(                                          &
            interpCoeff(1:nInterp_(1),1:nInterp_(2),1) * fluidTemperature(i1:i2,j1:j2,1,1) )

       if (present(viscosity)) viscosity = sum(                                              &
            interpCoeff(1:nInterp_(1),1:nInterp_(2),1) * fluidViscosity(i1:i2,j1:j2,1,1) )

       if (present(massFraction)) massFraction = sum(                                        &
            interpCoeff(1:nInterp_(1),1:nInterp_(2),1) * fluidMassFraction(i1:i2,j1:j2,1,1) )

       if (present(pressure)) pressure = sum(                                                &
            interpCoeff(1:nInterp_(1),1:nInterp_(2),1) * fluidPressure(i1:i2,j1:j2,1,1) )

       if (present(velocity)) then
          velocity = 0.0_WP
          do i = 1, nDimensions
             velocity(i) = sum(interpCoeff(1:nInterp_(1),1:nInterp_(2),1) *                  &
                  fluidVelocity(i1:i2,j1:j2,1,i))
          end do
       end if

       if (present(divRhoU)) divRhoU = sum(                                                  &
            interpCoeff(1:nInterp_(1),1:nInterp_(2),1) * fluidDivRhoU(i1:i2,j1:j2,1,1) )

       if (present(gradRho)) then
          gradRho = 0.0_WP
          do i = 1, nDimensions
             gradRho(i) = sum(interpCoeff(1:nInterp_(1),1:nInterp_(2),1) *                   &
                  fluidGradRho(i1:i2,j1:j2,1,i) )
          end do
       end if

       if (present(vorticity)) then
          vorticity = 0.0_WP
          do i = 1, 3
             vorticity(i) = sum(interpCoeff(1:nInterp_(1),1:nInterp_(2),1) *                 &
                  fluidVorticity(i1:i2,j1:j2,1,i))
          end do
       end if

       if (present(stress)) then
          stress = 0.0_WP
          do i = 1, nDimensions
             stress(i) = sum(interpCoeff(1:nInterp_(1),1:nInterp_(2),1) *                    &
                  fluidStress(i1:i2,j1:j2,1,i))
          end do
       end if

       if (present(volumeFraction)) volumeFraction = sum(                                    &
            interpCoeff(1:nInterp_(1),1:nInterp_(2),1) * fluidVolumeFraction(i1:i2,j1:j2,1,1))

       if (present(levelset)) levelset = sum(                                                &
            interpCoeff(1:nInterp_(1),1:nInterp_(2),1) * fluidLevelset(i1:i2,j1:j2,1,1))

       if (present(levelsetNormal)) then
          levelsetNormal = 0.0_WP
          do i = 1, nDimensions
             levelsetNormal(i) = sum(interpCoeff(1:nInterp_(1),1:nInterp_(2),1) *            &
                  fluidLevelsetNormal(i1:i2,j1:j2,1,i))
          end do
       end if

    case (3)

       ! Get the interpolation points
       n = ceiling(0.5_WP * real(nInterp, WP))
       i1 = particleIndex(1) - n + 1
       j1 = particleIndex(2) - n + 1
       k1 = particleIndex(3) - n + 1
       if (mod(nInterp, 2) .gt. 0) then
          if (cartesianCoordinates(particleIndex(1)+1, 1) - particlePosition(1) .gt.         &
               particlePosition(1) - cartesianCoordinates(particleIndex(1),1)) i1 = i1 + 1
          if (cartesianCoordinates(particleIndex(2)+1, 2) - particlePosition(2) .gt.         &
               particlePosition(2) - cartesianCoordinates(particleIndex(2),2)) j1 = j1 + 1
          if (cartesianCoordinates(particleIndex(3)+1, 3) - particlePosition(3) .gt.         &
               particlePosition(3) - cartesianCoordinates(particleIndex(3),3)) k1 = k1 + 1
       end if
       i2 = i1 + nInterp - 1
       j2 = j1 + nInterp - 1
       k2 = k1 + nInterp - 1
       if (nOverlap(1,1) .eq. 0) i1 = max(i1, iStart(1))
       if (nOverlap(1,2) .eq. 0) i2 = min(i2, iEnd(1))
       if (nOverlap(2,1) .eq. 0) j1 = max(j1, iStart(2))
       if (nOverlap(2,2) .eq. 0) j2 = min(j2, iEnd(2))
       if (nOverlap(3,1) .eq. 0) k1 = max(k1, iStart(3))
       if (nOverlap(3,2) .eq. 0) k2 = min(k2, iEnd(3))
       nInterp_(1) = i2-i1+1
       nInterp_(2) = j2-j1+1
       nInterp_(3) = k2-k1+1

       ! Compute the interpolation coefficients
       call lagrange_basis(nInterp_(1), cartesianCoordinates(i1:i2,1), particlePosition(1),  &
            weight(1:nInterp_(1),1))
       call lagrange_basis(nInterp_(2), cartesianCoordinates(j1:j2,2), particlePosition(2),  &
            weight(1:nInterp_(2),2))
       call lagrange_basis(nInterp_(3), cartesianCoordinates(k1:k2,3), particlePosition(3),  &
            weight(1:nInterp_(3),3))

       ! Combine the interpolation coefficients
       interpCoeff = 0.0_WP
       do k = 1, nInterp_(1)
          do j = 1, nInterp_(2)
             do i = 1, nInterp_(3)
                interpCoeff(i,j,k) = weight(i,1) * weight(j,2) * weight(k,3)
             end do
          end do
       end do

       ! Interpolate the fluid variables
       if (present(density)) density =                                                       &
            sum(interpCoeff(1:nInterp_(1),1:nInterp_(2),1:nInterp_(3)) *                     &
            fluidDensity(i1:i2,j1:j2,k1:k2,1))

       if (present(temperature)) temperature =                                               &
            sum(interpCoeff(1:nInterp_(1),1:nInterp_(2),1:nInterp_(3)) *                     &
            fluidTemperature(i1:i2,j1:j2,k1:k2,1))

       if (present(viscosity)) viscosity =                                                   &
            sum(interpCoeff(1:nInterp_(1),1:nInterp_(2),1:nInterp_(3)) *                     &
            fluidViscosity(i1:i2,j1:j2,k1:k2,1))

       if (present(massFraction)) massFraction =                                             &
            sum(interpCoeff(1:nInterp_(1),1:nInterp_(2),1:nInterp_(3)) *                     &
            fluidMassFraction(i1:i2,j1:j2,k1:k2,1))

       if (present(pressure)) pressure =                                                     &
            sum(interpCoeff(1:nInterp_(1),1:nInterp_(2),1:nInterp_(3)) *                     &
            fluidPressure(i1:i2,j1:j2,k1:k2,1))

       if (present(velocity)) then
          do i = 1, nDimensions
             velocity(i) = sum(interpCoeff(1:nInterp_(1),1:nInterp_(2),1:nInterp_(3)) *      &
                  fluidVelocity(i1:i2,j1:j2,k1:k2,i))
          end do
       end if

       if (present(divRhoU)) divRhoU =                                                       &
            sum(interpCoeff(1:nInterp_(1),1:nInterp_(2),1:nInterp_(3)) *                     &
            fluidDivRhoU(i1:i2,j1:j2,k1:k2,1))

       if (present(gradRho)) then
          do i = 1, nDimensions
             gradRho(i) = sum(interpCoeff(1:nInterp_(1),1:nInterp_(2),1:nInterp_(3)) *       &
                  fluidGradRho(i1:i2,j1:j2,k1:k2,i))
          end do
       end if
       
       if (present(vorticity)) then
          do i = 1, 3
             vorticity(i) = sum(interpCoeff(1:nInterp_(1),1:nInterp_(2),1:nInterp_(3)) *     &
                  fluidVorticity(i1:i2,j1:j2,k1:k2,i))
          end do
       end if

       if (present(stress)) then
          do i = 1, nDimensions
             stress(i) = sum(interpCoeff(1:nInterp_(1),1:nInterp_(2),1:nInterp_(3)) *        &
                  fluidStress(i1:i2,j1:j2,k1:k2,i))
          end do
       end if

       if (present(volumeFraction)) volumeFraction =                                         &
            sum(interpCoeff(1:nInterp_(1),1:nInterp_(2),1:nInterp_(3)) *                     &
            fluidVolumeFraction(i1:i2,j1:j2,k1:k2,1))

       if (present(levelset)) levelset =                                                     &
            sum(interpCoeff(1:nInterp_(1),1:nInterp_(2),1:nInterp_(3)) *                     &
            fluidLevelset(i1:i2,j1:j2,k1:k2,1))

       if (present(levelsetNormal)) then
          do i = 1, nDimensions
             levelsetNormal(i) = sum(interpCoeff(1:nInterp_(1),1:nInterp_(2),1:nInterp_(3)) *&
                  fluidLevelsetNormal(i1:i2,j1:j2,k1:k2,i))
          end do
       end if

    end select

    return
  end subroutine interpolate_fluid_to_particle


  ! Particle => grid extrapolation routine
  ! Second-order tri-linear extrapolation
  ! --------------------------------------
  subroutine extrapolate_particle_to_grid(particleIndex, particlePosition,                   &
       particleQuantity, source)

    ! External modules
    use math

    implicit none

    ! Arguments
    integer, intent(in) :: particleIndex(3)
    real(WP), intent(in) :: particlePosition(3), particleQuantity
    real(WP), intent(inout) :: source(pStart(1):pEnd(1), pStart(2):pEnd(2), pStart(3):pEnd(3))

    ! Local variables
    integer :: i, j, k, n, i1, j1, k1, i2, j2, k2
    integer, parameter :: nInterp_ = 2
    real(WP) :: interpCoeff(nInterp_, nInterp_, nInterp_), weight(nInterp_, nDimensions), norm

    select case (nDimensions)

    case (1)

       ! Get the interpolation points
       n = ceiling(0.5_WP * real(nInterp_, WP))
       i1 = particleIndex(1) - n + 1; i2 = particleIndex(1) + n
       if (nOverlap(1,1) .eq. 0) i1 = max(i1, iStart(1))
       if (nOverlap(1,2) .eq. 0) i2 = min(i2, iEnd(1))

       ! Compute the interpolation coefficients
       call lagrange_basis(i2-i1+1, cartesianCoordinates(i1:i2, 1), particlePosition(1),     &
            weight(:,1))

       ! Combine the interpolation coefficients to form a linear interpolation
       interpCoeff = 0.0_WP
       interpCoeff(:,1,1) = weight(:,1)

       ! Correct for walls with Neumann condition
       norm = sum(interpCoeff) + epsilon(1.0_WP)
       interpCoeff = interpCoeff / norm

       ! Perform extrapolation on source
       source(i1:i2, 1, 1) = source(i1:i2, 1, 1) + interpCoeff(:, 1, 1) * particleQuantity

    case (2)

       ! Get the interpolation points
       n = ceiling(0.5_WP * real(nInterp_, WP))
       i1 = particleIndex(1) - n + 1; i2 = particleIndex(1) + n
       j1 = particleIndex(2) - n + 1; j2 = particleIndex(2) + n
       if (nOverlap(1,1) .eq. 0) i1 = max(i1, iStart(1))
       if (nOverlap(1,2) .eq. 0) i2 = min(i2, iEnd(1))
       if (nOverlap(2,1) .eq. 0) j1 = max(j1, iStart(2))
       if (nOverlap(2,2) .eq. 0) j2 = min(j2, iEnd(2))

       ! Compute the interpolation coefficients
       call lagrange_basis(i2-i1+1, cartesianCoordinates(i1:i2, 1), particlePosition(1),     &
            weight(:,1))
       call lagrange_basis(j2-j1+1, cartesianCoordinates(j1:j2, 2), particlePosition(2),     &
            weight(:,2))

       ! Combine the interpolation coefficients
       interpCoeff = 0.0_WP
       do j = 1, nInterp_
          do i = 1, nInterp_
             interpCoeff(i,j,1) = weight(i,1) * weight(j,2)
          end do
       end do

       ! Correct for walls with Neumann condition
       norm = sum(interpCoeff) + epsilon(1.0_WP)
       interpCoeff = interpCoeff / norm

       ! Perform extrapolation on source
       source(i1:i2, j1:j2, 1) = source(i1:i2, j1:j2, 1) + interpCoeff(:, :, 1) *            &
            particleQuantity

    case (3)

       ! Get the interpolation points
       n = ceiling(0.5_WP * real(nInterp_, WP))
       i1 = particleIndex(1) - n + 1; i2 = particleIndex(1) + n
       j1 = particleIndex(2) - n + 1; j2 = particleIndex(2) + n
       k1 = particleIndex(3) - n + 1; k2 = particleIndex(3) + n
       if (nOverlap(1,1) .eq. 0) i1 = max(i1, iStart(1))
       if (nOverlap(1,2) .eq. 0) i2 = min(i2, iEnd(1))
       if (nOverlap(2,1) .eq. 0) j1 = max(j1, iStart(2))
       if (nOverlap(2,2) .eq. 0) j2 = min(j2, iEnd(2))
       if (nOverlap(3,1) .eq. 0) k1 = max(k1, iStart(3))
       if (nOverlap(3,2) .eq. 0) k2 = min(k2, iEnd(3))

       ! Compute the interpolation coefficients
       call lagrange_basis(i2-i1+1, cartesianCoordinates(i1:i2, 1), particlePosition(1),     &
            weight(:,1))
       call lagrange_basis(j2-j1+1, cartesianCoordinates(j1:j2, 2), particlePosition(2),     &
            weight(:,2))
       call lagrange_basis(k2-k1+1, cartesianCoordinates(k1:k2, 3), particlePosition(3),     &
            weight(:,3))

       ! Combine the interpolation coefficients
       do k = 1, nInterp_
          do j = 1, nInterp_
             do i = 1, nInterp_
                interpCoeff(i,j,k) = weight(i,1) * weight(j,2) * weight(k,3)
             end do
          end do
       end do

       ! Correct for walls with Neumann condition
       norm = sum(interpCoeff) + epsilon(1.0_WP)
       interpCoeff = interpCoeff / norm

       ! Perform extrapolation on source
       source(i1:i2, j1:j2, k1:k2) = source(i1:i2, j1:j2, k1:k2) + interpCoeff *             &
            particleQuantity

    end select

    return
  end subroutine extrapolate_particle_to_grid


  ! Mollification function
  ! -----------------------
  subroutine mollify(particlePosition, ic, jc, kc, fs)

    ! External modules
    use grid_functions

    implicit none

    ! Arguments
    integer, intent(in) :: ic, jc, kc
    real(WP), intent(in) :: particlePosition(3)
    real(WP), intent(out) :: fs
    real(WP) :: gridPosition(3)

    ! Local variables
    real(WP) :: s

    ! Get normalized Euclidean distance
    select case (nDimensions)
    case (1)
       gridPosition(1) = cartesianCoordinates(ic,1)
    case (2)
       gridPosition(1) = cartesianCoordinates(ic,1)
       gridPosition(2) = cartesianCoordinates(jc,2)
    case (3)
       gridPosition(1) = cartesianCoordinates(ic,1)
       gridPosition(2) = cartesianCoordinates(jc,2)
       gridPosition(3) = cartesianCoordinates(kc,3)
    end select

    s = sqrt(sum((particlePosition(1:nDimensions)-gridPosition(1:nDimensions))**2)) /        &
         filterSize * 2.0_WP

    ! Integrate the quartic spline kernel (Pepiot and Desjardins, Powder Technology, 2012)
    if (s.le.0.5_WP) then
       fs = 0.25*s**4 - 5.0_WP/8.0_WP*s**2 + 115.0_WP/192.0_WP
    else if (s.le.1.5_WP) then
       fs = -1.0_WP/6.0_WP*s**4 + 5.0_WP/6.0_WP*s**3 - 5.0_WP/4.0_WP*s**2 +                  &
            5.0_WP/24.0_WP*s + 55.0_WP/96.0_WP
    else if (s.le.2.5_WP) then
       fs = (2.5_WP-s)**4/24.0_WP
    else
       fs = 0.0_WP
    end if

    ! Integrate
    !fs = fs / jacobian(grid_index(ic,jc,kc), 1)

    return
  end subroutine mollify


  ! Filter particle data after initial extrapolation (2nd order)
  ! ------------------------------------------------------------
  subroutine filter_extrapolated_field(f)

    ! External modules
    use grid_functions, only : implicit_diffusion

    implicit none

    ! Arguments
    real(WP), intent(inout) :: f(:,:)

    if (diffusionAmount .gt. 0.0_WP) call implicit_diffusion(f, diffusionAmount)

    return
  end subroutine filter_extrapolated_field
  

  ! Recompute viscous & pressure stresses for interphase exchange
  ! -------------------------------------------------------------
  subroutine recompute_fluid_stress(fluidStress)

    ! External modules
    use first_derivative
    use grid_functions, only : transform_fluxes, laplacian
    use state

    implicit none

    ! Arguments
    real(WP), intent(out) :: fluidStress(nGridPoints, nUnknowns)

    ! Local variables
    integer :: i, k
    real(WP), dimension(nGridPoints, nUnknowns, nDimensions) :: temp, temp2

    ! Zero out the array
    temp = 0.0_WP

    ! Compute inviscid & viscous fluxes
    select case (nDimensions)

    case (1)
       temp(:,2,1) = -pressure(:,1)
       temp(:,3,1) = -velocity(:,1) * pressure(:,1)
       if (useViscosity) then
          temp(:,2,1) = temp(:,2,1) + stressTensor(:,1)
          temp(:,3,1) = temp(:,3,1) + velocity(:,1) * stressTensor(:,1) - heatFlux(:,1)
          if (equationOfState .eq. IDEAL_GAS_MIXTURE) then
             temp(:,3,1) = temp(:,3,1) - enthalpyFlux(:,1)
          end if
          do k = 1, nSpecies
             temp(:,k+3,1) = temp(:,k+3,1) - speciesFlux(:,k,1)
          end do
       end if

    case (2)
       temp(:,2,1) = -pressure(:,1)
       temp(:,3,2) = -pressure(:,1)
       temp(:,4,1) = -velocity(:,1) * pressure(:,1)
       temp(:,4,2) = -velocity(:,2) * pressure(:,1)
       if (useViscosity) then
          temp(:,2,1) = temp(:,2,1) + stressTensor(:,1)
          temp(:,2,2) = temp(:,2,2) + stressTensor(:,2)
          temp(:,3,1) = temp(:,3,1) + stressTensor(:,3)
          temp(:,3,2) = temp(:,3,2) + stressTensor(:,4)
          temp(:,4,1) = temp(:,4,1) + velocity(:,1) * stressTensor(:,1) +                    &
               velocity(:,2) * stressTensor(:,3) - heatFlux(:,1)
          temp(:,4,2) = temp(:,4,2) - velocity(:,1) * stressTensor(:,2) +                    &
               velocity(:,2) * stressTensor(:,4) - heatFlux(:,2)
          if (equationOfState .eq. IDEAL_GAS_MIXTURE) then
             temp(:,4,1) = temp(:,4,1) - enthalpyFlux(:,1)
             temp(:,4,2) = temp(:,4,2) - enthalpyFlux(:,2)
          end if
          do k = 1, nSpecies
             temp(:,k+4,1) = temp(:,k+4,1) - speciesFlux(:,k,1)
             temp(:,k+4,2) = temp(:,k+4,2) - speciesFlux(:,k,2)
          end do
       end if

    case (3)
       temp(:,2,1) = -pressure(:,1)
       temp(:,3,2) = -pressure(:,1)
       temp(:,4,3) = -pressure(:,1)
       temp(:,5,1) = -velocity(:,1) * pressure(:,1)
       temp(:,5,2) = -velocity(:,2) * pressure(:,1)
       temp(:,5,3) = -velocity(:,3) * pressure(:,1)
       if (useViscosity) then
          temp(:,2,1) = temp(:,2,1) + stressTensor(:,1)
          temp(:,2,2) = temp(:,2,2) + stressTensor(:,2)
          temp(:,2,3) = temp(:,2,3) + stressTensor(:,3)
          temp(:,3,1) = temp(:,3,1) + stressTensor(:,4)
          temp(:,3,2) = temp(:,3,2) + stressTensor(:,5)
          temp(:,3,3) = temp(:,3,3) + stressTensor(:,6)
          temp(:,4,1) = temp(:,4,1) + stressTensor(:,7)
          temp(:,4,2) = temp(:,4,2) + stressTensor(:,8)
          temp(:,4,3) = temp(:,4,3) + stressTensor(:,9)
          temp(:,5,1) = temp(:,5,1) + velocity(:,1) * stressTensor(:,1) +                    &
               velocity(:,2) * stressTensor(:,4) +                                           &
               velocity(:,3) * stressTensor(:,7) - heatFlux(:,1)
          temp(:,5,2) = temp(:,5,2) - velocity(:,1) * stressTensor(:,2) +                    &
               velocity(:,2) * stressTensor(:,5) +                                           &
               velocity(:,3) * stressTensor(:,8) - heatFlux(:,2)
          temp(:,5,3) = temp(:,5,3) - velocity(:,1) * stressTensor(:,3) +                    &
               velocity(:,2) * stressTensor(:,6) +                                           &
               velocity(:,3) * stressTensor(:,9) - heatFlux(:,3)
          if (equationOfState .eq. IDEAL_GAS_MIXTURE) then
             temp(:,5,1) = temp(:,5,1) - enthalpyFlux(:,1)
             temp(:,5,2) = temp(:,5,2) - enthalpyFlux(:,2)
             temp(:,5,3) = temp(:,5,3) - enthalpyFlux(:,3)
          end if
          do k = 1, nSpecies
             temp(:,k+5,1) = temp(:,k+5,1) - speciesFlux(:,k,1)
             temp(:,k+5,2) = temp(:,k+5,2) - speciesFlux(:,k,2)
             temp(:,k+5,3) = temp(:,k+5,3) - speciesFlux(:,k,3)
          end do
       end if

    end select

    ! Transform stress tensor from Cartesian to contravariant form
    call transform_fluxes(temp, temp2)

    ! Take derivatives of stress tensor
    do i = 1, nDimensions
       call first_derivative_apply(i, temp2(:,:,i))
    end do
    fluidStress = sum(temp2, dim = 3)

    ! Multiply by Jacobian
    do i = 1, nUnknowns
       fluidStress(:,i) = sum(temp2(:,i,:), dim = 2) * jacobian(:,1)
    end do

    return
  end subroutine recompute_fluid_stress


  ! Compute the pDV work term
  ! -------------------------
  subroutine compute_pdv_work(source)

    ! External modules
    use parallel
    use grid_functions, only : gradient
    use state

    implicit none

    ! Arguments
    real(WP), intent(inout) :: source(nGridPoints, nUnknowns)

    ! Local variables
    integer :: i
    real(WP), dimension(nGridPoints, nDimensions**2) :: gradUp
    real(WP), dimension(nGridPoints, nDimensions) :: alphaUp

    ! Compute the gradient
    do i = 1, nDimensions
       alphaUp(:,i) = particleVelocity(:,i) * (1.0_WP - volumeFraction(:,1))
    end do
    call gradient(alphaUp, gradUp)

    ! Add to right-hand-side
    select case (nDimensions)
    case (1)
       source(:,3) = source(:,3) + (stressTensor(:,1) - pressure(:,1)) * gradUp(:,1)
    case (2)
       source(:,4) = source(:,4) + (stressTensor(:,1) - pressure(:,1)) * gradUp(:,1) +       &
            stressTensor(:,2) * gradUp(:,2) + stressTensor(:,3) * gradUp(:,3) +              &
            (stressTensor(:,4) - pressure(:,1)) * gradUp(:,4)
    case (3)
       source(:,5) = source(:,5) + (stressTensor(:,1) - pressure(:,1)) * gradUp(:,1) +       &
            stressTensor(:,2) * gradUp(:,2) + stressTensor(:,3) * gradUp(:,3) +              &
            stressTensor(:,4) * gradUp(:,4) +                                                &
            (stressTensor(:,5) - pressure(:,1)) * gradUp(:,5) +                              &
            stressTensor(:,6) * gradUp(:,6) + stressTensor(:,7) * gradUp(:,7) +              &
            stressTensor(:,8) * gradUp(:,8) +                                                &
            (stressTensor(:,9) - pressure(:,1)) * gradUp(:,9)
    end select

    return
  end subroutine compute_pdv_work

end module particle_exchange


! =================================== !
! Setup the particle exchange routine !
! =================================== !
subroutine particle_exchange_setup

  ! Internal modules
  use particle_exchange

  ! External modules
  use parallel
  use parser
  use math, only : pi

  ! Local variables
  integer ::gridIndex, i, j, n, ijk(3), ierror

  if (.not. useParticles) return

  ! Particle volume factor (treat particles as cylinders in 2D, spheres in 3D)
  select case(nDimensions)
  case (1)
     volumeFactor = 1.0_WP
  case (2)
     volumeFactor = 0.25_WP * pi
  case (3)
     volumeFactor = pi / 6.0_WP
  end select

  ! Get the interpolation order
  call parser_read('particle interp order', nInterp, 2)
  if (nInterp .lt. 2) call die('particle_exchange_setup: intrp order must be >1!')
  nOverlap = ceiling(0.5_WP * real(nInterp, WP))

  ! Include an additional point for particle motion during RK4
  nOverlap = nOverlap + 1

  ! Adjust the ghost points at the boundaries
  do i = 1, 3
     if (.not. isPeriodic(i)) then
        if (procCoords(i) .eq. 0) nOverlap(i,1) = 0
        if (procCoords(i) .eq. nProcsDir(i) - 1) nOverlap(i,2) = 0
     end if
  end do

  ! Prepare Cartesian vectors for particle localization
  allocate(cartesianCoordinates(1 - maxval(nOverlap) : maxval(globalGridSize) +              &
       maxval(nOverlap), nDimensions))
  cartesianCoordinates = -huge(1.0_WP)
  do j = 1, nDimensions
     do i = iStart(j), iEnd(j)
        ijk = iStart; ijk(j) = i
        gridIndex = ijk(1) - gridOffset(1) + localGridSize(1) *                              &
             (ijk(2) - 1 - gridOffset(2) + localGridSize(2) *                                &
             (ijk(3) - 1 - gridOffset(3)))
        cartesianCoordinates(i, j) = coordinates(gridIndex, j)
     end do
     !call MPI_Allgather(cartesianCoordinates(iStart(j):iEnd(j), j), localGridSize(j),        &
     !     MPI_REAL_WP, cartesianCoordinates(1:globalGridSize(3), j), globalGridSize(j),      &
     !     MPI_REAL_WP, commDir(j), ierror)
     call MPI_Allreduce(MPI_IN_PLACE, cartesianCoordinates(:, j),                            &
          size(cartesianCoordinates, 1), MPI_REAL_WP, MPI_MAX,                               &
          commDir(j), ierror)
     if (isPeriodic(j) .and. procCoords(j) .eq. 0) then
        do n = 1, nOverlap(j, 1)
           cartesianCoordinates(1 - n, j) = - periodicLength(j) +                            &
                cartesianCoordinates(globalGridSize(j) - n + 1, j)
        end do
     end if
     if (isPeriodic(j) .and. procCoords(j) .eq. nProcsDir(j) - 1) then
        do n = 1, nOverlap(j, 2)
           cartesianCoordinates(globalGridSize(j) + n, j) = periodicLength(j) +              &
                cartesianCoordinates(n, j)
        end do
     end if
     call MPI_Allreduce(MPI_IN_PLACE, cartesianCoordinates(:, j),                            &
          size(cartesianCoordinates, 1), MPI_REAL_WP, MPI_MAX, commDir(j), ierror)
  end do

  ! Get local indices of Cartesian arrays with ghost points
  pStart = 1; pEnd = 1
  select case (nDimensions)
  case (1)
     pStart(1) = gridOffset(1) + 1 - nOverlap(1,1)
     pEnd(1) = gridOffset(1) + localGridSize(1) + nOverlap(1,2)
  case (2)
     pStart(1) = gridOffset(1) + 1 - nOverlap(1,1)
     pEnd(1) = gridOffset(1) + localGridSize(1) + nOverlap(1,2)
     pStart(2) = gridOffset(2) + 1 - nOverlap(2,1)
     pEnd(2) = gridOffset(2) + localGridSize(2) + nOverlap(2,2)
  case (3)
     pStart(1) = gridOffset(1) + 1 - nOverlap(1,1)
     pEnd(1) = gridOffset(1) + localGridSize(1) + nOverlap(1,2)
     pStart(2) = gridOffset(2) + 1 - nOverlap(2,1)
     pEnd(2) = gridOffset(2) + localGridSize(2) + nOverlap(2,2)
     pStart(3) = gridOffset(3) + 1 - nOverlap(3,1)
     pEnd(3) = gridOffset(3) + localGridSize(3) + nOverlap(3,2)
  end select

  ! Allocate fluid arrays with extra layer of ghost cells
  allocate(fluidDensity(pStart(1) : pEnd(1), pStart(2) : pEnd(2), pStart(3) : pEnd(3), 1))
  allocate(fluidTemperature(pStart(1) : pEnd(1), pStart(2) : pEnd(2),                        &
       pStart(3) : pEnd(3), 1))
  allocate(fluidViscosity(pStart(1) : pEnd(1), pStart(2) : pEnd(2), pStart(3) : pEnd(3), 1))
  allocate(fluidVelocity(pStart(1) : pEnd(1), pStart(2) : pEnd(2),                           &
       pStart(3) : pEnd(3), nDimensions))
  allocate(fluidStress(pStart(1) : pEnd(1), pStart(2) : pEnd(2), pStart(3) : pEnd(3),        &
       nDimensions))
  allocate(fluidVolumeFraction(pStart(1) : pEnd(1), pStart(2) : pEnd(2),                     &
       pStart(3) : pEnd(3), 1))
  if (useSaffmanLift .or. useFriction) allocate(fluidVorticity(pStart(1) : pEnd(1),          &
       pStart(2) : pEnd(2), pStart(3) : pEnd(3), 3))
  if (useAddedMass) then
     allocate(fluidGradRho(pStart(1) : pEnd(1), pStart(2) : pEnd(2), pStart(3) : pEnd(3),    &
          nDimensions))
     allocate(fluidDivRhoU(pStart(1) : pEnd(1), pStart(2) : pEnd(2), pStart(3) : pEnd(3), 1))
  end if
  if (usePhaseChange) then
     allocate(fluidMassFraction(pStart(1) : pEnd(1), pStart(2) : pEnd(2),                    &
          pStart(3) : pEnd(3), 1))
     allocate(fluidPressure(pStart(1) : pEnd(1), pStart(2) : pEnd(2),                        &
          pStart(3) : pEnd(3), 1))
  end if
  if (useIBM .and. collisionsOn) then
     allocate(fluidLevelset(pStart(1) : pEnd(1), pStart(2) : pEnd(2),                        &
          pStart(3) : pEnd(3), 1))
     allocate(fluidLevelsetNormal(pStart(1) : pEnd(1), pStart(2) : pEnd(2),                  &
          pStart(3) : pEnd(3), nDimensions))
  end if

  ! Allocate the fluid source terms with extra layer of ghost cells
  if (twoWayCoupling) then
     allocate(momentumSource(pStart(1) : pEnd(1), pStart(2) : pEnd(2), pStart(3) : pEnd(3),  &
          nDimensions))
     allocate(energySource(pStart(1) : pEnd(1), pStart(2) : pEnd(2), pStart(3) : pEnd(3), 1))
     momentumSource = 0.0_WP
     energySource = 0.0_WP
     if (usePhaseChange .or. usePTKE) then
        allocate(massSource(pStart(1) : pEnd(1), pStart(2) : pEnd(2), pStart(3) : pEnd(3), 1))
        massSource = 0.0_WP
     end if
  end if

  ! Allocate the temporary Cartesian arrays
  allocate(cartScalar(pStart(1) : pEnd(1), pStart(2) : pEnd(2), pStart(3) : pEnd(3), 1))
  allocate(cartVector(pStart(1) : pEnd(1), pStart(2) : pEnd(2), pStart(3) : pEnd(3), nDimensions))

  ! Allocate the particle localization arrays
  if (collisionsOn) then
     allocate(nPartInCell(pStart(1)-1 : pEnd(1)+1, pStart(2)-1 : pEnd(2)+1,                  &
          pStart(3)-1 : pEnd(3)+1))
     allocate(partInCell(pStart(1)-1 : pEnd(1)+1, pStart(2)-1 : pEnd(2)+1,                   &
          pStart(3)-1 : pEnd(3)+1, 400))
  end if

  ! Allocate particle dependent variables
  allocate(particleVelocity(nGridPoints, nDimensions)); particleVelocity = 0.0_WP
  allocate(granularTemperature(nGridPoints, 1))

  return
end subroutine particle_exchange_setup


! ===================================== !
! Cleanup the particle exchange routine !
! ===================================== !
subroutine particle_exchange_cleanup

  ! Internal modules
  use particle_exchange

  implicit none

  if (allocated(cartesianCoordinates)) deallocate(cartesianCoordinates)
  if (allocated(fluidLevelset)) deallocate(fluidLevelset)
  if (allocated(fluidLevelsetNormal)) deallocate(fluidLevelsetNormal)
  if (allocated(fluidGradRho)) deallocate(fluidGradRho)
  if (allocated(fluidDivRhoU)) deallocate(fluidDivRhoU)
  if (allocated(fluidDensity)) deallocate(fluidDensity)
  if (allocated(fluidTemperature)) deallocate(fluidTemperature)
  if (allocated(fluidViscosity)) deallocate(fluidViscosity)
  if (allocated(fluidMassFraction)) deallocate(fluidMassFraction)
  if (allocated(fluidPressure)) deallocate(fluidPressure)
  if (allocated(fluidVelocity)) deallocate(fluidVelocity)
  if (allocated(fluidVorticity)) deallocate(fluidVorticity)
  if (allocated(fluidStress)) deallocate(fluidStress)
  if (allocated(fluidVolumeFraction)) deallocate(fluidVolumeFraction)
  if (allocated(momentumSource)) deallocate(momentumSource)
  if (allocated(energySource)) deallocate(energySource)
  if (allocated(massSource)) deallocate(massSource)
  if (allocated(nPartInCell)) deallocate(nPartInCell)
  if (allocated(partInCell)) deallocate(partInCell)
  if (allocated(particleVelocity)) deallocate(particleVelocity)
  if (allocated(granularTemperature)) deallocate(granularTemperature)
  if (allocated(cartScalar)) deallocate(cartScalar)
  if (allocated(cartVector)) deallocate(cartVector)

  return
end subroutine particle_exchange_cleanup


! ==================================================== !
! Prepare the Cartesian arrays for interphase exchange !
! ==================================================== !
subroutine prepare_cartesian_arrays

  ! Internal modules
  use particle_exchange

  ! External modules
  use parallel
  use state
  use state_functions
  use grid_levelset
  use grid_functions, only : gradient, divergence

  implicit none

  ! Local variables
  integer :: gridIndex, i, j, k, n, nVariables
  real(WP), dimension(nGridPoints, nUnknowns) :: stress
  real(WP), allocatable :: vorticity(:,:), divRhoU(:), gradRho(:,:), temp(:,:),              &
       gridBuffer(:,:,:,:)

  if (.not. useParticles) return

  ! Zero-out the source terms
  if (twoWayCoupling) then
     momentumSource = 0.0_WP
     energySource   = 0.0_WP
     if (usePhaseChange .or. usePTKE) massSource = 0.0_WP
  end if

  ! Zero-out the temporary arrays
  cartVector = 0.0_WP
  cartScalar = 0.0_WP

  ! Recompute the fluid stresses
  call recompute_fluid_stress(stress)

  ! Particle solver needs density, temperature, viscosity, volume fraction,
  ! velocity, and stress + any optional variables
  nVariables = 4 + 2 * nDimensions
  if (useSaffmanLift .or. useFriction) nVariables = nVariables + 3
  if (useAddedMass) nVariables = nVariables + 1 + nDimensions
  if (usePhaseChange) nVariables = nVariables + 2
  if (useIBM .and. collisionsOn) nVariables = nVariables + 1 + nDimensions
  allocate(temp(nGridPoints, nVariables))
  allocate(gridBuffer(pStart(1) : pEnd(1), pStart(2) : pEnd(2), pStart(3) : pEnd(3),         &
       nVariables))

  ! Pack data belonging to local processor in temporary array
  temp(:,1) = volumeFraction(:,1)
  temp(:,2) = conservedVariables(:,1) / volumeFraction(:,1)
  temp(:,3) = temperature(:,1)
  temp(:,4) = dynamicViscosity(:,1)
  if (useShockCapturing) temp(:,4) = temp(:,4) - artificialShearViscosity(:,1)
  if (useLES) temp(:,4) = temp(:,4) - turbulentViscosity(:,1)
  temp(:,5:4+nDimensions) = velocity(:,1:nDimensions)
  temp(:,5+nDimensions:2*nDimensions+4) = stress(:,2:nDimensions+1)
  n = 4 + 2 * nDimensions
  
  ! Include vorticity if using lift
  if (useSaffmanLift .or. useFriction) then
     allocate(vorticity(nGridPoints, 3))
     call compute_vorticity(velocityGradient, vorticity)
     temp(:,n+1:n+3) = vorticity(:,1:3)
     n = n + 3
  end if

  ! Include divergence of rho*u and gradient of density if using added mass
  if (useAddedMass) then
     allocate(divRhoU(nGridPoints))
     allocate(gradRho(nGridPoints, 3))
     call gradient(conservedVariables(:,1), gradRho)
     call divergence(conservedVariables(:,2:nDimensions+1), divRhoU)
     temp(:,n+1) = divRhoU
     temp(:,n+2:n+1+nDimensions) = gradRho(:,1:nDimensions)
     n = n + 1 + nDimensions
  end if
  
  if (usePhaseChange) then
     temp(:,n+1) = massFraction(:,vaporIndex)
     temp(:,n+2) = pressure(:,1)
     n = n + 2
  end if
  
  if (useIBM .and. collisionsOn) then
     temp(:,n+1) = levelset(:, 1)
     temp(:,n+2:n+1+nDimensions) = levelsetNormal(:, 1:nDimensions)
     n = n + 1 + nDimensions
  end if

  ! Filter the fluid quantities (except stress)
  if (filterFluid) then
     call filter_extrapolated_field(temp)
     !temp(:,5+nDimensions:2*nDimensions+4) = stress(:,2:nDimensions+1)
  end if

  ! Send to Cartesian array
  do k = iStart(3), iEnd(3)
     do j = iStart(2), iEnd(2)
        do i = iStart(1), iEnd(1)
           gridIndex = i - gridOffset(1) + localGridSize(1) *                                &
                (j - 1 - gridOffset(2) + localGridSize(2) *                                  &
                (k - 1 - gridOffset(3)))
           gridBuffer(i,j,k,:) = temp(gridIndex,:)
        end do
     end do
  end do

  ! Exhange data in direction `1`
  if (nDimensions.ge.1) call fill_ghost_points(gridBuffer(:, iStart(2) : iEnd(2),            &
       iStart(3) : iEnd(3), 1:nVariables), 1, (/nOverlap(1,1), nOverlap(1,2) /))

  ! Exhange data in direction `2`
  if (nDimensions.ge.2) call fill_ghost_points(gridBuffer(iStart(1) - nOverlap(1,1) :        &
       iEnd(1) + nOverlap(1,2), :, iStart(3):iEnd(3), 1:nVariables), 2,                      &
       (/nOverlap(2,1), nOverlap(2,2) /))

  ! Exhange data in direction `3`
  if (nDimensions.ge.3) call fill_ghost_points(gridBuffer(iStart(1) - nOverlap(1,1) :        &
       iEnd(1) + nOverlap(1,2), iStart(2) - nOverlap(2,1) : iEnd(2) + nOverlap(2,2), :,      &
       1:nVariables), 3, (/nOverlap(3,1), nOverlap(3,2) /))

  ! Unpack the temporary array
  fluidVolumeFraction(:,:,:,1) = gridBuffer(:,:,:,1)
  fluidDensity(:,:,:,1) = gridBuffer(:,:,:,2)
  fluidTemperature(:,:,:,1) = gridBuffer(:,:,:,3)
  fluidViscosity(:,:,:,1) = gridBuffer(:,:,:,4)
  fluidVelocity(:,:,:,1:nDimensions) = gridBuffer(:,:,:,5:4+nDimensions)
  fluidStress(:,:,:,1:nDimensions) = gridBuffer(:,:,:,5+nDimensions:2*nDimensions+4)
  n = 4 + 2 * nDimensions
  if (useSaffmanLift .or. useFriction) then
     fluidVorticity(:,:,:,1:3) = gridBuffer(:,:,:,n+1:n+3)
     n = n + 3
  end if
  if (useAddedMass) then
     fluidDivRhoU(:,:,:,1) = gridBuffer(:,:,:,n+1)
     fluidGradRho(:,:,:,1:nDimensions) = gridBuffer(:,:,:,n+2:n+1+nDimensions)
     n = n + 1 + nDimensions
  end if
  if (usePhaseChange) then
     fluidMassFraction(:,:,:,1) = gridBuffer(:,:,:,n+1)
     fluidPressure(:,:,:,1) = gridBuffer(:,:,:,n+2)
     n = n + 2
  end if
  if (useIBM .and. collisionsOn) then
     fluidLevelset(:,:,:,1) = gridBuffer(:,:,:,n+1)
     fluidLevelsetNormal(:,:,:,1:nDimensions) = gridBuffer(:,:,:,n+2:n+1+nDimensions)
     n = n + 1 + nDimensions
  end if
  
  ! Cleanup
  deallocate(gridBuffer, temp)
  if (allocated(vorticity)) deallocate(vorticity)
  if (allocated(gradRho)) deallocate(gradRho)
  if (allocated(divRhoU)) deallocate(divRhoU)

  return
end subroutine prepare_cartesian_arrays


! ====================================== !
! Transfer particle sources to the fluid !
! ====================================== !
subroutine transfer_interphase_exchange(part, dudt, dTdt, dmdt, dwdt, Cpd, Lv)

  ! Internal modules
  use particle_exchange

  ! External modules
  use solver_options

  implicit none

  ! Arguments
  type(t_Particle), intent(in) :: part
  real(WP), dimension(3), intent(in) :: dudt, dwdt
  real(WP), intent(in) :: dTdt, dmdt, Cpd, Lv

  ! Local variables
  integer :: i
  real(WP) :: mass, volume, massSource_, momentumSource_(3), energySource_

  if (.not. twoWayCoupling .or. part%stop.eq.1) return

  ! Get particle volume
  volume = volumeFactor * part%diameter**nDimensions

  ! Create the source terms
  mass = volume * particleDensity
  momentumSource_ = mass * dudt
  energySource_ = mass * Cpd * dTdt
  if (usePhaseChange) then
     massSource_ = dmdt
     momentumSource_(1:nDimensions) = momentumSource_(1:nDimensions) +                       &
          part%velocity(1:nDimensions) * dmdt
     energySource_ = energySource_ + dmdt*part%temperature - dmdt*Lv
  end if
  if (useFriction) then
     energySource_ = energySource_ + sum(dwdt * part%angularVelocity) * mass
  end if

  ! Extrapolate to the grid
  do i = 1, nDimensions

     ! Send momentum to the fluid
     call extrapolate_particle_to_grid(part%gridIndex, part%position, momentumSource_(i),    &
          momentumSource(:,:,:,i))

     ! Add work due to drag to the energy exchange term
     energySource_ = energySource_ + momentumSource_(i) * part%velocity(i)

  end do

  ! Send energy to the fluid
  call extrapolate_particle_to_grid(part%gridIndex, part%position, energySource_,            &
       energySource(:,:,:,1))

  ! Send mass to the fluid
  if (usePhaseChange .or. usePTKE) then
     call extrapolate_particle_to_grid(part%gridIndex, part%position, massSource_,           &
          massSource(:,:,:,1))
  end if

  return
end subroutine transfer_interphase_exchange


! ======================================================== !
! Add particle sources to the fluid during the forward run !
! ======================================================== !
subroutine particle_source_forward(source)

  ! Internal modules
  use particle_exchange

  ! External modules
  use parallel
  use first_derivative
  use grid_functions, only : transform_fluxes, gradient
  use state

  implicit none

  ! Arguments
  real(WP), intent(inout) :: source(nGridPoints, nUnknowns)

  ! Local variables
  integer :: i, j, k, gridIndex
  real(WP), dimension(nGridPoints, nUnknowns) :: temp
  real(WP), dimension(nGridPoints, nDimensions) :: gradVolumeFraction

  if (.not.useParticles .or. .not.twoWayCoupling) return
  
  ! Start the particle timer
  call timing_start('particles')

  ! Communicate at the borders
  do i = 1, nDimensions
     call border_summation(momentumSource, i, (/nOverlap(i,1), nOverlap(i,2) /))
     call border_summation(energySource, i, (/nOverlap(i,1), nOverlap(i,2) /))
     if (usePhaseChange .or. usePTKE) then
        call border_summation(massSource, i, (/nOverlap(i,1), nOverlap(i,2) /))
     end if
  end do

  ! Normalize the source terms
  temp = 0.0_WP
  do k = iStart(3), iEnd(3)
     do j = iStart(2), iEnd(2)
        do i = iStart(1), iEnd(1)
           gridIndex = i - gridOffset(1) + localGridSize(1) *                                &
                (j - 1 - gridOffset(2) + localGridSize(2) *                                  &
                (k - 1 - gridOffset(3)))
           temp(gridIndex, 2:nDimensions+1) = momentumSource(i,j,k,1:nDimensions) *          &
                jacobian(gridIndex, 1)
           temp(gridIndex, nDimensions+2) = energySource(i,j,k,1) * jacobian(gridIndex, 1)
           if (usePhaseChange) temp(gridIndex, nDimensions+3) =                              &
                massSource(i,j,k,1) * jacobian(gridIndex, 1)
        end do
     end do
  end do

  ! Filter the source terms
  i = nDimensions + 2
  if (usePhaseChange) i = i + 1
  call filter_extrapolated_field(temp(:,2:i))

  ! Add PTKE source
  if (usePTKE) temp(:, nDimensions+3) = temp(:, nDimensions+2) -                             &
       sum(temp(:, 2:nDimensions+1) * velocity(:, 1:nDimensions), dim = 2)

  ! Add non-conservative terms
  call gradient(volumeFraction(:,1), gradVolumeFraction)
  select case (nDimensions)
  case (1)
     temp(:,2) = temp(:,2) + (pressure(:,1) - stressTensor(:,1)) * gradVolumeFraction(:,1)
     temp(:,3) = temp(:,3) + heatFlux(:,1) * gradVolumeFraction(:,1)
     do k = 1, nSpecies
        temp(:,3+k) = temp(:,3+k) + speciesFlux(:,k,1) * gradVolumeFraction(:,1)
     end do
  case (2)
     temp(:,2) = temp(:,2) + (pressure(:,1) - stressTensor(:,1)) * gradVolumeFraction(:,1) - &
          stressTensor(:,2) * gradVolumeFraction(:,2)
     temp(:,3) = temp(:,3) - stressTensor(:,3) * gradVolumeFraction(:,1) +                   &
          (pressure(:,1) - stressTensor(:,4)) * gradVolumeFraction(:,2)
     temp(:,4) = temp(:,4) + heatFlux(:,1) * gradVolumeFraction(:,1)
     temp(:,4) = temp(:,4) + heatFlux(:,2) * gradVolumeFraction(:,2)
     do k = 1, nSpecies
        temp(:,4+k) = temp(:,4+k) + speciesFlux(:,k,1) * gradVolumeFraction(:,1) +           &
             speciesFlux(:,k,2) * gradVolumeFraction(:,2)
     end do
     if (usePTKE) then

     end if
  case (3)
     temp(:,2) = temp(:,2) + (pressure(:,1) - stressTensor(:,1)) * gradVolumeFraction(:,1) - &
          stressTensor(:,4) * gradVolumeFraction(:,2) -                                      &
          stressTensor(:,7) * gradVolumeFraction(:,3)
     temp(:,3) = temp(:,3) - stressTensor(:,2) * gradVolumeFraction(:,1) +                   &
          (pressure(:,1) - stressTensor(:,5)) * gradVolumeFraction(:,2) -                    &
          stressTensor(:,8) * gradVolumeFraction(:,3)
     temp(:,4) = temp(:,4) - stressTensor(:,3) * gradVolumeFraction(:,1) -                   &
          stressTensor(:,6) * gradVolumeFraction(:,2) +                                      &
          (pressure(:,1) - stressTensor(:,9)) * gradVolumeFraction(:,3)
     temp(:,5) = temp(:,5) + heatFlux(:,1) * gradVolumeFraction(:,1)
     temp(:,5) = temp(:,5) + heatFlux(:,2) * gradVolumeFraction(:,2)
     temp(:,5) = temp(:,5) + heatFlux(:,3) * gradVolumeFraction(:,3)
     do k = 1, nSpecies
        temp(:,5+k) = temp(:,5+k) + speciesFlux(:,k,1) * gradVolumeFraction(:,1) +           &
             speciesFlux(:,k,2) * gradVolumeFraction(:,2) +                                  &
             speciesFlux(:,k,3) * gradVolumeFraction(:,3)
     end do
     if (usePTKE) then
        
     end if
  end select

  ! Compute pDV work term
  call compute_pdv_work(temp)
 
  ! Add them to the right-hand side
  do i = 1, nGridPoints
     source(i,:) = source(i,:) + temp(i,:)
  end do

  ! Stop the particle timer
  call timing_stop('particles')

  return
end subroutine particle_source_forward


! ======================================================== !
! Add particle sources to the fluid during the adjoint run !
! ======================================================== !
subroutine particle_source_adjoint(source)

  ! Internal modules
  use particle_exchange

  ! External modules
  use parallel
  use state

  implicit none

  ! Arguments
  real(WP), intent(inout) :: source(nGridPoints, nUnknowns)

  if (.not. twoWayCoupling) return

  ! Start the particle timer
  call timing_start('particles')

  ! Not yet implemented!

  ! Stop the particle timer
  call timing_stop('particles')

  return
end subroutine particle_source_adjoint
