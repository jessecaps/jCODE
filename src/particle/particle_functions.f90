module particle_functions

  implicit none

  ! Empty module

end module particle_functions


! =========================== !
! Compute the volume fraction !
! =========================== !
subroutine compute_volume_fraction

  ! Internal modules
  use particle_exchange

  ! External modules
  use parallel
  use math, only : pi
  use state

  implicit none

  ! Local variables
  integer :: i, j, k, gridIndex
  real(WP) :: volume
  real(WP), parameter :: vfClip = 0.2_WP

  if (.not. useParticles) return

  ! Get the particle concentration
  cartScalar = 0.0_WP
  do i = 1, nParticles
     if (particles(i)%stop .eq. 1) cycle
     volume = volumeFactor * particles(i)%diameter**nDimensions
     call extrapolate_particle_to_grid(particles(i)%gridIndex, particles(i)%position,        &
          volume, cartScalar(:,:,:,1))
  end do

  ! Communicate at the borders
  do i = 1, nDimensions
     call border_summation(cartScalar, i, (/nOverlap(i,1), nOverlap(i,2) /))
  end do

  ! Send back the fluid volume fraction
  do k = iStart(3), iEnd(3)
     do j = iStart(2), iEnd(2)
        do i =  iStart(1), iEnd(1)
           gridIndex = i - gridOffset(1) + localGridSize(1) *                                &
                (j - 1 - gridOffset(2) + localGridSize(2) *                                  &
                (k - 1 - gridOffset(3)))
           volumeFraction(gridIndex, 1) = cartScalar(i,j,k,1) * jacobian(gridIndex, 1)
        end do
     end do
  end do

  ! Filter
  call filter_extrapolated_field(volumeFraction(:,1:1))

  ! Store the fluid volume fraction
  volumeFraction(:,1) = 1.0_WP - volumeFraction(:,1)

  ! Clip volume fraction
  where (volumeFraction(:,1) .lt. vfClip)
     volumeFraction(:,1) = vfClip
  end where

  return
end subroutine compute_volume_fraction


! ========================================== !
! Compute particle-phase dependent variables !
! Currently computes volume fraction,        !
! Eulerian particle velocity, and the        !
! granular temperature                       !
! ========================================== !
subroutine compute_particle_dependent_variables

  ! Internal modules
  use particle_exchange

  ! External modules
  use parallel
  use state

  implicit none

  ! Local variables
  integer :: i, j, k, gridIndex
  real(WP) :: volume, ep, up(3)
  real(WP), parameter :: vfClip = 0.2_WP

  if (.not. useParticles) return

  ! Get the particle concentration and velocity
  cartScalar = 0.0_WP
  cartVector = 0.0_WP
  do i = 1, nParticles
     if (particles(i)%stop .eq. 1) cycle
     volume = volumeFactor * particles(i)%diameter**nDimensions
     call extrapolate_particle_to_grid(particles(i)%gridIndex, particles(i)%position,        &
          volume, cartScalar(:,:,:,1))
     do j = 1, nDimensions
        call extrapolate_particle_to_grid(particles(i)%gridIndex, particles(i)%position,     &
             particles(i)%velocity(j) * volume, cartVector(:,:,:,j))
       end do
  end do

  ! Communicate at the borders
  do i = 1, nDimensions
     call border_summation(cartScalar, i, (/nOverlap(i,1), nOverlap(i,2) /))
     call border_summation(cartVector, i, (/nOverlap(i,1), nOverlap(i,2) /))
  end do

  ! Send back the fluid volume fraction
  do k = iStart(3), iEnd(3)
     do j = iStart(2), iEnd(2)
        do i =  iStart(1), iEnd(1)
           gridIndex = i - gridOffset(1) + localGridSize(1) *                                &
                (j - 1 - gridOffset(2) + localGridSize(2) *                                  &
                (k - 1 - gridOffset(3)))
           volumeFraction(gridIndex, 1) = cartScalar(i,j,k,1) * jacobian(gridIndex, 1)
           particleVelocity(gridIndex, 1:nDimensions) = cartVector(i,j,k,1:nDimensions) *    &
                jacobian(gridIndex, 1)
        end do
     end do
  end do

  ! Filter and remove volume fraction
  call filter_extrapolated_field(volumeFraction(:,1:1))
  call filter_extrapolated_field(particleVelocity(:,1:nDimensions))
  do i = 1, nDimensions
     particleVelocity(:,i) = particleVelocity(:,i) / (volumeFraction(:,1) + epsilon(1.0_WP))
  end do
  
  ! Store the fluid volume fraction
  volumeFraction(:,1) = 1.0_WP - volumeFraction(:,1)

  ! Compute granular temperature
  if (useGranularTemperature) then
     
     ! Get data belonging to local processor
     do k = iStart(3), iEnd(3)
        do j = iStart(2), iEnd(2)
           do i = iStart(1), iEnd(1)
              gridIndex = i - gridOffset(1) + localGridSize(1) *                             &
                   (j - 1 - gridOffset(2) + localGridSize(2) *                               &
                   (k - 1 - gridOffset(3)))
              fluidVelocity(i,j,k,1:nDimensions) = particleVelocity(gridIndex,1:nDimensions)
              fluidVolumeFraction(i,j,k,1) = volumeFraction(gridIndex, 1)
           end do
        end do
     end do

     ! Get data from the neighboring processors
     if (nDimensions .ge. 1) then
        ! Exhange data in direction `1`
        call fill_ghost_points(fluidVelocity(:, iStart(2) : iEnd(2), iStart(3) : iEnd(3), :),&
             1, (/nOverlap(1,1), nOverlap(1,2) /))
        call fill_ghost_points(fluidVolumeFraction(:, iStart(2) : iEnd(2),                   &
             iStart(3) : iEnd(3), :), 1, (/nOverlap(1,1), nOverlap(1,2) /))
     end if
     if (nDimensions .ge. 2) then
        ! Exhange data in direction `2`
        call fill_ghost_points(fluidVelocity(iStart(1) - nOverlap(1,1) : iEnd(1) +           &
             nOverlap(1,2), :, iStart(3):iEnd(3), :) , 2, (/nOverlap(2,1), nOverlap(2,2) /))
        call fill_ghost_points(fluidVolumeFraction(iStart(1) - nOverlap(1,1) : iEnd(1) +     &
             nOverlap(1,2), :, iStart(3):iEnd(3), :) , 2, (/nOverlap(2,1), nOverlap(2,2) /))
     end if
     if (nDimensions .ge. 3) then
        ! Exhange data in direction `3`
        call fill_ghost_points(fluidVelocity(iStart(1) - nOverlap(1,1) : iEnd(1) +           &
             nOverlap(1,2), iStart(2) - nOverlap(2,1) : iEnd(2) + nOverlap(2,2), :, :) , 3,  &
             (/nOverlap(3,1), nOverlap(3,2) /))
        call fill_ghost_points(fluidVolumeFraction(iStart(1) - nOverlap(1,1) : iEnd(1) +     &
             nOverlap(1,2), iStart(2) - nOverlap(2,1) : iEnd(2) + nOverlap(2,2), :, :) , 3,  &
             (/nOverlap(3,1), nOverlap(3,2) /))
     end if

     ! Send granular temperature to the grid
     cartScalar = 0.0_WP
     do i = 1, nParticles
        if (particles(i)%stop .eq. 1) cycle
        ! Interpolate
        call interpolate_fluid_to_particle(particles(i)%gridIndex, particles(i)%position,    &
             velocity = up, volumeFraction = ep)
        ep = 1.0_WP - ep
        volume = volumeFactor * particles(i)%diameter**nDimensions
        do j = 1, nDimensions
           call extrapolate_particle_to_grid(particles(i)%gridIndex, particles(i)%position,  &
                (particles(i)%velocity(j) - up(j))**2 * volume / ep, cartScalar(:,:,:,1))
        end do
     end do

     ! Communicate at the borders
     do i = 1, nDimensions
        call border_summation(cartScalar, i, (/nOverlap(i,1), nOverlap(i,2) /))
     end do

     ! Send back granular temperature and normalize
     do k = iStart(3), iEnd(3)
        do j = iStart(2), iEnd(2)
           do i =  iStart(1), iEnd(1)
              gridIndex = i - gridOffset(1) + localGridSize(1) *                             &
                   (j - 1 - gridOffset(2) + localGridSize(2) *                               &
                   (k - 1 - gridOffset(3)))
              granularTemperature(gridIndex, 1) = cartScalar(i,j,k,1) /                      &
                   real(nDimensions, WP) * jacobian(gridIndex, 1) 
           end do
        end do
     end do

     ! Filter
     call filter_extrapolated_field(granularTemperature(:,1:1))

  end if

  ! Clip volume fraction
  where (volumeFraction(:,1) .lt. vfClip)
     volumeFraction(:,1) = vfClip
  end where

  return
end subroutine compute_particle_dependent_variables


! ==================================================================== !
! Compute pseudo turbulent Reynolds stress using the model from        !
! M. Mehrabadi et al. Pseudo-turbulent gas-phase velocity fluctuations !
! in homogeneous gas–solid flow: fixed particle assemblies and freely  !
! evolving suspensions, JFM (2015)                                     !
! ==================================================================== !
subroutine compute_pseudo_reynolds_stress

  ! Internal modules
  use particle_exchange

  ! External modules
  use state

  implicit none

  ! Local variables
  integer :: i, dir
  real(WP) :: rhoK, bPar, bPerp, buf
  real(WP), dimension(3) :: u1, u2, u3, u1dot, u2dot, u3dot
  real(WP), dimension(3,3) :: Q, temp, bij
  real(WP) :: particleReynoldsNumber(nGridPoints)
  real(WP), parameter :: oneThird = 1.0_WP / 3.0_WP
  real(WP), parameter :: a = 0.523_WP
  real(WP), parameter :: b = 0.305_WP
  real(WP), parameter :: c = 0.114_WP
  real(WP), parameter :: d = 3.511_WP
  real(WP), parameter :: e = 1.801_WP
  real(WP), parameter :: f = 0.005_WP

  if (.not.usePTKE) return

  ! Assume we transport pseudo-TKE as the first scalar after energy
  
  select case(nDimensions)

  case (2)
     
     !  Assume the Reynolds stresses are isotropic in 2D
     do i = 1, nGridPoints
        
        rhoK = conservedVariables(i,nDimensions+3) / volumeFraction(i,1)
        reynoldsStress(i,1) = rhoK
        reynoldsStress(i,4) = rhoK

     end do

  case (3)

     do i = 1, nGridPoints

        ! Compute the anisotropy tensor components (smoothly decay to 0 away from particles)
        buf = 0.5_WP * ((tanh(10.0_WP * ((1.0_WP - volumeFraction(i,1)) - 0.01_WP))) + 1.0_WP)
        bPar = buf * a / (1.0_WP + b * exp(-c * particleReynoldsNumber(i))) *                &
             exp(-d * (1.0_WP - volumeFraction(i,1)) / (1.0_WP + e *                         &
             exp(-f * particleReynoldsNumber(i))))
        bPerp = -0.5_WP * bPar

        ! Add trace
        bPar = bPar + oneThird
        bPerp = bPerp + oneThird

        ! Get the local PTKE
        rhoK = conservedVariables(i,nDimensions+3) / volumeFraction(i,1)

        ! Compute the reference axes
        u1 = particleVelocity(i,:) - velocity(i,:)
        u1dot = dot_product(u1, u1) + epsilon(1.0_WP)

        ! Generate an orthonormal set based on max slip component
        dir = maxloc(abs(u1), dim = 1)
        select case (dir)
        case (1)
           ! Max slip in direction `1`
           u2    = -(u1(2) / u1dot) * u1
           u2(2) = 1.0_WP + u2(2)
           u2dot = dot_product(u2, u2) + epsilon(1.0_WP)
        case(2)
           ! Max slip in direction `2`
           u2    = -(u1(1) / u1dot) * u1
           u2(1) = 1.0_WP + u2(1)
           u2dot = dot_product(u2, u2) + epsilon(1.0_WP)
        case(3)
           ! Max slip in direction `23`
           u2    = -(u1(1) / u1dot) * u1
           u2(1) = 1.0_WP + u2(1)
           u2dot = dot_product(u2, u2) + epsilon(1.0_WP)
        end select

        ! Right-hand coordinate system for U3 (cross-product)
        u3(1)  = u1(2) * u2(3) - u1(3) * u2(2)
        u3(2)  = u1(3) * u2(1) - u1(1) * u2(3)
        u3(3)  = u1(1) * u2(2) - u1(2) * u2(1)
        u3dot = dot_product(u3, u3) + epsilon(1.0_WP)

        ! Normalize basis vectors
        u1 = u1 / sqrt(u1dot)
        u2 = u2 / sqrt(u2dot)
        u3 = u3 / sqrt(u3dot)

        ! Construct rotation matrices
        Q(:,1) = u1
        Q(:,2) = u2
        Q(:,3) = u3

        ! Direction `1` parallel to slip velocity (bPar)
        ! Multiply diagonal b^dag by Q^T
        temp(1,1) = bPar *Q (1,1)
        temp(1,2) = bPar *Q (2,1)
        temp(1,3) = bPar *Q (3,1)

        temp(2,1) = bPerp * Q(1,2)
        temp(2,2) = bPerp * Q(2,2)
        temp(2,3) = bPerp * Q(3,2)

        temp(3,1) = bPerp * Q(1,3)
        temp(3,2) = bPerp * Q(2,3)
        temp(3,3) = bPerp * Q(3,3)

        ! Multiply Q by b^dag*Q^T (temp) to get bij tensor
        bij(1,1) = Q(1,1) * temp(1,1) + Q(1,2) * temp(2,1) + Q(1,3) * temp(3,1)
        bij(1,2) = Q(1,1) * temp(1,2) + Q(1,2) * temp(2,2) + Q(1,3) * temp(3,2)
        bij(1,3) = Q(1,1) * temp(1,3) + Q(1,2) * temp(2,3) + Q(1,3) * temp(3,3)

        bij(2,1) = Q(2,1) * temp(1,1) + Q(2,2) * temp(2,1) + Q(2,3) * temp(3,1)
        bij(2,2) = Q(2,1) * temp(1,2) + Q(2,2) * temp(2,2) + Q(2,3) * temp(3,2)
        bij(2,3) = Q(2,1) * temp(1,3) + Q(2,2) * temp(2,3) + Q(2,3) * temp(3,3)

        bij(3,1) = Q(3,1) * temp(1,1) + Q(3,2) * temp(2,1) + Q(3,3) * temp(3,1)
        bij(3,2) = Q(3,1) * temp(1,2) + Q(3,2) * temp(2,2) + Q(3,3) * temp(3,2)
        bij(3,3) = Q(3,1) * temp(1,3) + Q(3,2) * temp(2,3) + Q(3,3) * temp(3,3)

        ! Update the Reynolds stresses
        reynoldsStress(i,1) = 2.0_WP * rhoK * bij(1,1)
        reynoldsStress(i,2) = 2.0_WP * rhoK * bij(1,2)
        reynoldsStress(i,3) = 2.0_WP * rhoK * bij(1,3)
        reynoldsStress(i,4) = reynoldsStress(i,2)
        reynoldsStress(i,5) = 2.0_WP * rhoK * bij(2,2)
        reynoldsStress(i,6) = 2.0_WP * rhoK * bij(2,3)
        reynoldsStress(i,7) = reynoldsStress(i,3)
        reynoldsStress(i,8) = reynoldsStress(i,6)
        reynoldsStress(i,9) = 2.0_WP * rhoK * bij(3,3)

     end do

  end select

  return
end subroutine compute_pseudo_reynolds_stress


! ================================================================== !
! Compute the correction to Stokes drag for an individual particle.  !
! This routine contains a library of drag models that returns        !
! dragCorrection = Cd * Rep / 24 in order to update acceleration:    !
! du/dt = dragCorrection * alpha * (uf - up) / taup in particle_rhs. !
! ================================================================== !
subroutine compute_drag_correction(part, density, temperature, viscosity, velocity, vf,     &
     dragCorrection)

  ! Internal modules
  use particle_exchange

  ! External modules
  use parallel
  use math, only : pi
  use solver_options
  use particle_solver
  use particle

  implicit none

  ! Arguments
  type(t_Particle), intent(in) :: part
  real(WP), intent(in) :: density, temperature, viscosity, velocity(3), vf
  real(WP), intent(out) :: dragCorrection

  ! Non-dimensional numbers
  real(WP) :: Kn, Ma, MaInf, Rep, ReInf, vp, CD

  ! Dimensional parameters
  real(WP) :: relativeVelocity, soundSpeed

  ! Model parameters
  real(WP) :: Mcr, CDMcr, zsub, fsub(3), Csub(3), fsup(3), Csup(3), zsup, CDM1,              &
       CDM175, CDstd, CD1, CD2, CM, GM, HM, JM, fKn, Knp, sDrag, b1, b2, b3, CD3, Rem,       &
       fsupM1(3), fsupM175(3), fsubM1(3), fsubMcr(3),                                        &
       delt_o, Ao, alpha_o, alpha_hoc, u_rat, T_s_rat, Ma_s, alpha, C2, alpha_inf,           &
       Theta_singh, Re_tilde, Cd_c, Wr, f_low_Kn, s, Cd_FM, a_singh, Br_singh, omega_exp

  ! Particle volume fraction
  vp = 1.0_WP - vf

  ! Relative superficial velocity
  relativeVelocity = sqrt(sum((part%velocity(1:nDimensions)-velocity(1:nDimensions))**2))

  ! Particle Reynolds number
  Rep = density * part%diameter * relativeVelocity / viscosity + epsilon(1.0_WP)

  ! Particle Mach number
  soundSpeed = sqrt((ratioOfSpecificHeats-1.0_WP) * temperature)
  Ma = relativeVelocity / soundSpeed + epsilon(1.0_WP)

  ! Particle Knudsen number
  Kn = sqrt(0.5_WP*pi*ratioOfSpecificHeats) * Ma / Rep

  ! Compute the drag correction
  select case (dragModel)
  case (NO_DRAG)

     ! No drag
     dragCorrection = 0.0_WP

  case (STOKES)

     ! Stokes drag
     dragCorrection = 1.0_WP

  case (SCHILLER_NAUMANN)

     ! Clift, Grace & Weber (1978)
     if (Rep.lt.1000.0_WP) then
        dragCorrection = 1.0_WP + 0.15_WP * Rep**(0.687_WP)
     else
        dragCorrection = 0.44_WP * Rep / 24.0_WP
     end if

  case (GIDASPOW)

     ! Gidaspow (1994)
     ! This model is usually written as d(rho u)/dt = ... + beta (up-uf)
     ! Here, dragCorrection = beta / (rhop * vp) * taup
     Rep = Rep * vf
     if (vf.ge.0.8_WP) then
        ! Wen & Yu
        if (Rep.lt.1000.0_WP) then
           dragCorrection = (1.0_WP + 0.15_WP * Rep**(0.687_WP)) * vf**(-2.65_WP)
        else
           dragCorrection = 0.44_WP * Rep / 24.0_WP * vf**(-2.65_WP)
        end if
     else
        ! Ergun
        dragCorrection = (150.0_WP * vp / vf + 1.75_WP * Rep / vf) / 18.0_WP
     end if

  case (TENNETI)

     ! Tenneti & Subramaniam (2011)
     Rep = Rep * vf
     b1 = 5.81_WP*vp/vf**3+0.48_WP*vp**(1.0_WP/3.0_WP)/vf**4
     b2 = vp**3*Rep*(0.95_WP+0.61_WP*vp**3/vf**2)

     dragCorrection = vf * ((1.0_WP+0.15_WP*Rep**(0.687_WP))/vf**3+b1+b2)

  case (BASSET)

     ! Account for compressibility & rarefaction effects
     ! Cunningham correction factor (Crowe et al. Multiphase Flow)
     dragCorrection = 1.0_WP + Kn * (2.49_WP + 0.84_WP * exp(-1.74_WP/Kn))
     dragCorrection = 1.0_WP  / dragCorrection

  case (PARMAR)

     ! Based on Mach-dependent drag coefficient proposed Parmar (2010)
     ! Unsteady forces on a particle in compressible flows. PhD Thesis, U. Florida.
     ! Includes volume fraction correction by Sangani. Valid for Map <= 1.75.
     ! ---------------------------------------------------------------------------

     ! Rep-dependent drag coefficients based on Mach number
     Ma = min(Ma, 1.75_WP)
     Mcr = 0.6_WP
     CDMcr = (1.0_WP + 0.15_WP * Rep**(0.684_WP)) + 0.513_WP /                               &
          (1.0_WP + 483.0_WP / Rep**(0.669_WP))
     CDM1 = (1.0_WP + 0.118_WP * Rep**(0.813_WP)) + 0.69_WP /                                &
          (1.0_WP + 3550.0_WP / Rep**(0.793_WP))
     CDM175 = (1.0_WP + 0.107_WP * Rep**(0.867_WP)) + 0.646_WP /                             &
          (1.0_WP + 861.0_WP / Rep**(0.634_WP))

     ! Compute the drag coefficient
     if (Ma .le. Mcr) then
        ! Subcritical regime (standard drag of Clift and Gauvin)
        CDstd = (1.0_WP + 0.15_WP * Rep**(0.687_WP)) + 0.42_WP /                             &
             (1.0_WP + 42500.0_WP / Rep**(1.16_WP))
        dragCorrection = CDstd + (CDMcr - CDstd) * Ma / Mcr

     else if (Ma.gt.Mcr .and. Ma.le.1.0_WP) then
        ! Intermediate regime (supercritical and subsonic)
        Csub(1) = 6.48_WP; Csub(2) = 9.28_WP; Csub(3) = 12.21_WP
        fsub(1) = -0.087_WP + 2.92_WP*Ma - 4.75_WP*Ma**2 + 2.83_WP*Ma**3
        fsub(2) = -0.12_WP  + 2.66_WP*Ma - 4.36_WP*Ma**2 + 2.53_WP*Ma**3
        fsub(3) =  1.84_WP - 5.13_WP*Ma + 6.05_WP*Ma**2 - 1.91_WP*Ma**3
        fsubM1(1) = 0.913_WP; fsubM1(2) = 0.71_WP; fsubM1(3) = 0.85_WP
        fsubMcr(1) = 0.56628_WP; fsubMcr(2) = 0.45288_WP; fsubMcr(3) = 0.52744_WP
        zsub =                                                                               &
             (fsub(1) - fsubMcr(1)) / (fsubM1(1) - fsubMcr(1)) *                             &
             (log(Rep) - Csub(2))/(Csub(1)-Csub(2)) * (log(Rep)-Csub(3))/(Csub(1)-Csub(3)) + &
             (fsub(2) - fsubMcr(2)) / (fsubM1(2) - fsubMcr(2)) *                             &
             (log(Rep) - Csub(1))/(Csub(2)-Csub(1)) * (log(Rep)-Csub(3))/(Csub(2)-Csub(3)) + &
             (fsub(3) - fsubMcr(3)) / (fsubM1(3) - fsubMcr(3)) *                             &
             (log(Rep) - Csub(1))/(Csub(3)-Csub(1)) * (log(Rep)-Csub(2))/(Csub(3)-Csub(2))
        dragCorrection = CDMcr + (CDM1 - CDMcr) * zsub

     else
        ! Supersonic regime
        Csup(1) = 6.48_WP; Csup(2) = 8.93_WP; Csup(3) = 12.21_WP
        fsup(1) = 0.126_WP + 1.15_WP*Ma - 0.306_WP*Ma**2 - 0.007_WP*Ma**3 -                  &
             0.061_WP*exp((1.0_WP-Ma)/0.011_WP)
        fsup(2) = -0.901_WP + 2.93_WP*Ma - 1.573_WP*Ma**2 + 0.286_WP*Ma**3 -                 &
             0.042_WP*exp((1.0_WP-Ma)/0.01_WP)
        fsup(3) = 0.13_WP + 1.42_WP*Ma - 0.818_WP*Ma**2 + 0.161_WP*Ma**3 -                   &
             0.043_WP*exp((1.0_WP-Ma)/0.012_WP)
        fsupM1(1) = 0.902_WP; fsupM1(2) = 0.7_WP; fsupM1(3) = 0.85_WP
        fsupM175(1) = 1.163859375_WP; fsupM175(2) = 0.94196875_WP; fsupM175(3) = 0.972734375_WP
        zsup =                                                                               &
             (fsup(1) - fsupM1(1)) / (fsupM175(1) - fsupM1(1)) *                             &
             (log(Rep) - Csup(2))/(Csup(1)-Csup(2)) * (log(Rep)-Csup(3))/(Csup(1)-Csup(3)) + &
             (fsup(2) - fsupM1(2)) / (fsupM175(2) - fsupM1(2)) *                             &
             (log(Rep) - Csup(1))/(Csup(2)-Csup(1)) * (log(Rep)-Csup(3))/(Csup(2)-Csup(3)) + &
             (fsup(3) - fsupM1(3)) / (fsupM175(3) - fsupM1(3)) *                             &
             (log(Rep) - Csup(1))/(Csup(3)-Csup(1)) * (log(Rep)-Csup(2))/(Csup(3)-Csup(2))
        dragCorrection = CDM1 + (CDM175 - CDM1) * zsup
     end if

     ! Correct for volume fraction (Sangani et al., 1991)
     dragCorrection = dragCorrection * (1.0_WP + 2.0_WP * vp) / (1.0_WP - vp)**2

  case (THEO)

     CD = -6.9791_WP*Ma**6 + 22.967_WP*Ma**5 - 29.792_WP*Ma**4 + 19.563_WP*Ma**3 -           &
          6.2377_WP*Ma**2 + 0.9322_WP*Ma + 0.3597_WP

     dragCorrection = CD * Rep / 24.0_WP


!!$  case (LOTH)
!!$
!!$     ! Old drag model proposed by Loth (2008)
!!$     ! Loth, E. (2008). Compressibility and rarefaction effects on drag of a
!!$     ! spherical particle. AIAA journal, 46(9), 2219-2228.
!!$     ! ---------------------------------------------------------------------
!!$
!!$     ! Rarefraction-dominant regime
!!$     Knp = sqrt(pi*ratioOfSpecificHeats/2.0_WP)*(Ma/Rep)
!!$     fKn = 1.0_WP/(1.0_WP + Knp*(2.514_WP + 0.8_WP*exp(-0.55_WP/Knp)))
!!$     CD1 = 24.0_WP/Rep * (1.0_WP + 0.15_WP*Rep**(0.687_WP)) * fKn
!!$     sDrag = Ma*sqrt(ratioOfSpecificHeats/2.0_WP)
!!$
!!$     CD2prime = (1.0_WP + 2.0_WP*(sDrag)**2) * exp(-1.0_WP*(sDrag)**2) /                     &
!!$          ((sDrag)**3*sqrt(pi) + epsilon(1.0_WP)) + (4.0_WP*(sDrag)**4 +                     &
!!$          4.0_WP*(sDrag)**2 - 1.0_WP) * erf(sDrag)/(2.0_WP*(sDrag)**4 + epsilon(1.0_WP)) 
!!$
!!$     CD2primeT = CD2prime +                                                                  &
!!$          2.0_WP/(3.0_WP*(sDrag) + epsilon(1.0_WP)) * sqrt(pi*part%temperature/temperature)
!!$
!!$     CD2 = CD2primeT/(1.0_WP + (CD2prime/(1.63_WP) - 1.0_WP)*sqrt(Rep/45.0_WP))
!!$
!!$     ! Compression-dominant regime
!!$     if (Ma .le. 1.45_WP) then
!!$        CM = 5.0_WP/3.0_WP + 2.0_WP/3.0_WP*tanh(3.0_WP*log(Ma + 0.1_WP))
!!$     else
!!$        CM = 2.044_WP + 0.2_WP*exp(-1.8_WP*(log(Ma/1.5_WP))**2)
!!$     end if
!!$
!!$     if (Ma .lt. 0.89_WP) then
!!$        GM = 1.0_WP - 1.525_WP*Ma**4
!!$     else
!!$        GM = 0.0002_WP + 0.0008_WP*tanh(12.77_WP*(Ma - 2.02_WP))
!!$     end if
!!$
!!$     HM = 1.0_WP - (0.258_WP*CM)/(1.0_WP + 514_WP*GM)
!!$
!!$     if (Rep .gt. 45.0_WP) then
!!$        CD = 24/Rep*(1.0_WP + 0.15_WP * Rep**(0.687_WP))*HM + 0.42_WP*CM/                    &
!!$             (1+(42500.0_WP*GM/Rep**(1.16_WP)))
!!$     else
!!$        CD = CD1/(1.0_WP + Ma**4) + (Ma**4*CD2)/(1.0_WP + Ma**4) 
!!$     end if
!!$
!!$     dragCorrection = CD * Rep / 24.0_WP

  case (LOTH)

     ! New drag model proposed by Loth (2021)
     ! Loth, E., Tyler Daspit, J., Jeong, M., Nagata, T., & Nonomura, T. (2021).
     ! Supersonic and Hypersonic Drag Coefficients for a Sphere. AIAA Journal, 1-14
     ! -----------------------------------------------------------------------------

     if (Rep .le. 45.0_WP) then
        ! Rarefraction-dominated regime
        Knp = sqrt(0.5_WP * pi * ratioOfSpecificHeats) * Ma / Rep
        fKn = 1.0_WP / (1.0_WP + Knp*(2.514_WP + 0.8_WP*exp(-0.55_WP/Knp)))
        CD1 = (1.0_WP + 0.15_WP*Rep**(0.687_WP)) * fKn
        sDrag = Ma * sqrt(0.5_WP * ratioOfSpecificHeats)
        if (Ma .le. 1.0_WP) then
           JM = 2.26_WP - 0.1_WP/Ma + 0.14_WP/Ma**3
        else
           JM = 1.6_WP + 0.25_WP/Ma + 0.11_WP/Ma**2 + 0.44_WP/Ma**3
        end if
        CD2 = (1.0_WP + 2.0_WP*sDrag**2) * exp(-sDrag**2) / (sDrag**3*sqrt(pi) +             &
             epsilon(1.0_WP)) + (4.0_WP*sDrag**4 + 4.0_WP*sDrag**2 - 1.0_WP) *               &
             erf(sDrag) / (2.0_WP*sDrag**4 + epsilon(1.0_WP)) +                              &
             2.0_WP / (3.0_WP * sDrag + epsilon(1.0_WP)) * sqrt(pi)
        CD2 = CD2 / (1.0_WP + (CD2/JM - 1.0_WP) * sqrt(Rep/45.0_WP))
        dragCorrection = CD1 / (1.0_WP + Ma**4) +                                            &
             Rep / 24.0_WP * Ma**4 * CD2 / (1.0_WP + Ma**4)

     else
        ! Compression-dominated regime
        if (Ma .lt. 1.5_WP) then
           CM = 1.65_WP + 0.65_WP * tanh(4.0_WP*Ma - 3.4_WP)
        else
           CM = 2.18_WP - 0.13_WP * tanh(0.9_WP*Ma - 2.7_WP)
        end if
        if (Ma .lt. 0.8_WP) then
           GM = 166.0_WP*Ma**3 + 3.29_WP*Ma**2 - 10.9_WP*Ma + 20.0_WP
        else
           GM = 5.0_WP + 40.0_WP*Ma**(-3.0_WP)
        end if

        if (Ma .lt. 1.0_WP) then
           HM = 0.0239_WP*Ma**3 + 0.212_WP*Ma**2 - 0.074_WP*Ma + 1.0_WP
        else
           HM = 0.93_WP + 1.0_WP / (3.5_WP + Ma**5)
        end if
        dragCorrection = (1.0_WP + 0.15_WP * Rep**(0.687_WP))*HM + Rep / 24.0_WP *           &
             0.42_WP*CM / (1.0_WP+42500.0_WP/Rep**(1.16_WP*CM) + GM/sqrt(Rep))
     end if

  case (HENDERSON)

     ! Non-relative dimensionless quantities
     MaInf = sqrt(sum((velocity(1:nDimensions))**2))
     ReInf = density * part%diameter * MaInf / viscosity
     ReInf = ReInf + epsilon(1.0_WP)
     MaInf = MaInf / soundSpeed

     if (Ma .le. 1.0_WP) then

        CD = 24.0_WP / (Rep + Ma*sqrt(0.5_WP*ratioOfSpecificHeats)*(4.33_WP                  &
             + (3.65_WP - 1.53_WP*(part%temperature/temperature))                            &
             / (1.0_WP + 0.353_WP*(part%temperature/temperature)))                           &
             * exp(-0.247_WP*Rep/(Ma*sqrt(0.5_WP*ratioOfSpecificHeats))))                    &
             + exp(-0.5_WP*Ma/sqrt(Rep))*((4.5_WP + 0.38_WP                                  &
             * (0.03_WP*Rep + 0.48_WP*sqrt(Rep)))                                            &
             / (1.0_WP + 0.03_WP*Rep + 0.48_WP*sqrt(Rep)) + 0.1_WP*Ma**2 + 0.2_WP*Ma**8)     &
             + (1.0_WP - exp(-1.0_WP*Ma/Rep)) * 0.6_WP * Ma*sqrt(0.5_WP*ratioOfSpecificHeats)

     else if (Ma .gt. 1.0 .and. Ma .lt. 1.75_WP) then

        CD = (24.0_WP / (Rep + 1.0_WP*sqrt(0.5_WP*ratioOfSpecificHeats)*(4.33_WP             &
             + (3.65_WP - 1.53_WP*(part%temperature/temperature))                            &
             / (1.0_WP + 0.353_WP*(part%temperature/temperature)))                           &
             * exp(-0.247_WP*Rep/(sqrt(0.5_WP*ratioOfSpecificHeats))))                &
             + exp(-0.5_WP*1.0_WP/sqrt(Rep))*((4.5_WP + 0.38_WP                              &
             * (0.03_WP*Rep + 0.48_WP*sqrt(Rep)))                                            &
             / (1.0_WP + 0.03_WP*Rep + 0.48_WP*sqrt(Rep))                                    &
             + 0.1_WP*1.0_WP**2 + 0.2_WP*1.0_WP**8) + (1.0_WP - exp(-1.0_WP*1.0_WP/Rep))     &
             * 0.6_WP * 1.0_WP*sqrt(0.5_WP*ratioOfSpecificHeats))                            &
             + 4.0_WP/3.0_WP * (MaInf - 1.0_WP)                                              &
             * ((0.9_WP + 0.34_WP/1.75_WP**2 + 1.86_WP * sqrt(1.75_WP/ReInf)                 &
             * (2.0_WP + 2.0_WP/(1.75_WP*sqrt(0.5_WP*ratioOfSpecificHeats))**2               &
             + 1.058_WP/(1.75_WP*sqrt(0.5_WP*ratioOfSpecificHeats))                          &
             * sqrt(part%temperature/temperature)                                            &
             - 1.0_WP/(1.75_WP*sqrt(0.5_WP*ratioOfSpecificHeats))**4))                       &
             / (1.0_WP + 1.86_WP * sqrt(1.75_WP/ReInf))                                      &
             - (24.0_WP / (Rep + 1.0_WP*sqrt(0.5_WP*ratioOfSpecificHeats)*(4.33_WP           &
             + (3.65_WP - 1.53_WP*(part%temperature/temperature))                            &
             / (1.0_WP + 0.353_WP*(part%temperature/temperature)))                           &
             * exp(-0.247_WP*Rep/(1.0_WP*sqrt(0.5_WP*ratioOfSpecificHeats))))                &
             + exp(-0.5_WP*1.0_WP/sqrt(Rep))*((4.5_WP + 0.38_WP                              &
             * (0.03_WP*Rep + 0.48_WP*sqrt(Rep)))                                            &
             / (1.0_WP + 0.03_WP*Rep + 0.48_WP*sqrt(Rep))                                    &
             + 0.1_WP*1.0_WP**2 + 0.2_WP*1.0_WP**8) + (1.0_WP - exp(-1.0_WP/Rep))            &
             * 0.6_WP * 1.0_WP*sqrt(0.5_WP*ratioOfSpecificHeats)))

     else

        CD = (0.9_WP + 0.34_WP/MaInf**2 + 1.86_WP * sqrt(MaInf/ReInf)                        &
             * (2.0_WP + 2.0_WP/(MaInf*sqrt(0.5_WP*ratioOfSpecificHeats))**2                 &
             + 1.058_WP/(MaInf*sqrt(0.5_WP*ratioOfSpecificHeats))                            &
             * sqrt(part%temperature/temperature)                                            &
             - 1.0_WP/(MaInf*sqrt(0.5_WP*ratioOfSpecificHeats))**4))                         &
             / (1.0_WP + 1.86_WP * sqrt(MaInf/ReInf))

     end if

     dragCorrection = CD * Rep / 24.0_WP

  case(SINGH)

     ! Singh, N., Kroells, M., Li, C., Ching, E., Ihme, M., Hogan, C. J., &
     ! Schwartzentruber, T. E. (2022). General drag coefficient for flow over
     ! spherical particles. AIAA journal, 60(2), 587-597.
     ! ----------------------------------------------------------------------

     if (Rep .ge. 10000.0_WP)                                                                &
          call die('compute_drag_correction: Singh model is not intended for Rep > 10k')

     omega_exp = 0.74_WP
     delt_o = 9.4_WP
     Ao = 0.271616115889543_WP
     alpha_o = 0.3555_WP
     alpha_hoc = 0.9071_WP
     if (Ma .lt. 1.0_WP) then
        u_rat = 1.0_WP
        T_s_rat = 1.0_WP
        Ma_s = Ma
        alpha = 1.0_WP
        C2 = 1.0_WP
     else
        u_rat = (2.0_WP + (ratioOfSpecificHeats - 1.0_WP) *Ma**2) /                          &
             ((ratioOfSpecificHeats+1.0_WP)*Ma**2)   
        T_s_rat = ((ratioOfSpecificHeats - 1.0_WP) *Ma**2 + 2.0_WP) *                        &
             (2.0_WP * ratioOfSpecificHeats * Ma**2 - (ratioOfSpecificHeats - 1.0_WP)) /     &
             ((ratioOfSpecificHeats + 1.0_WP)**2 * Ma**2)
        Ma_s = sqrt(((ratioOfSpecificHeats - 1.0_WP) * Ma**2 +2.0_WP) /                      &
             (2.0_WP * ratioOfSpecificHeats *Ma**2 - (ratioOfSpecificHeats - 1.0_WP)))
        alpha = 1.0_WP / (alpha_o * Ma + 1.0_WP - alpha_o)
        alpha_inf = 1.0_WP / alpha_o/ Ma
        C2 = (0.9_WP - Ao * (1 + (ratioOfSpecificHeats - 1.0_WP)**2 / 4.0_WP /               &
             ratioOfSpecificHeats)**(ratioOfSpecificHeats / (ratioOfSpecificHeats - 1.0_WP)))&
             / (1.0_WP - alpha_inf * (ratioOfSpecificHeats - 1.0_WP) /                       &
             (ratioOfSpecificHeats + 1.0_WP))
     end if
     theta_singh = (1.0_WP + (ratioOfSpecificHeats - 1.0_WP) * Ma_s**2 / 2.0_WP)**           &
          (ratioOfSpecificHeats / (ratioOfSpecificHeats - 1.0_WP))
     Re_tilde = Rep * (1.0_WP / alpha**2 / T_s_rat)**omega_exp *                             &
          theta_singh**((ratioOfSpecificHeats + 1.0_WP) / 2.0_WP / ratioOfSpecificHeats-     &
          (ratioOfSpecificHeats - 1.0_WP) * omega_exp / ratioOfSpecificHeats)
     Cd_c = C2 * (1.0_WP - alpha * u_rat) + (Ao * theta_singh) * (1.0_WP + delt_o /          &
          sqrt(Re_tilde))**2
     Wr = Ma**(2.0_WP * omega_exp) / Rep

     if (Ma .gt. 1.0_WP) then
        Wr = Wr * (1.0_WP + 1.0_WP / T_s_rat)**omega_exp
     end if

     f_low_Kn = 1.0_WP / (1.0_WP + Kn * (2.514_WP + 0.8_WP * exp(-0.55_WP / Kn))) /          &
          (1.0_WP + alpha_hoc * Wr)
     s = sqrt(0.5_WP * ratioOfSpecificHeats) * Ma
     Cd_FM = (1.0_WP + 2.0_WP * s**2) * exp(-s**2) / (s**3 * sqrt(pi)) +                     &
          (4.0_WP * s**4 + 4.0_WP*s**2 - 1.0_WP) * erf(s) / (2.0_WP * s**4) +                &
          2.0_WP / (3.0_WP*s) * sqrt(pi)

     a_singh = 1.8_WP

     Br_singh = Wr * (Ma**(2.0_WP * omega_exp - 1.0_WP) + 1.0_WP) /                          &
          (Ma**(2.0_WP * omega_exp - 1.0_WP))
     CD = Cd_FM * Br_singh**a_singh / (1.0_WP + Br_singh**a_singh) + Cd_c * f_low_Kn /       &
          (1.0_WP + Br_singh**a_singh)

     dragCorrection = CD * Rep / 24.0_WP

  case(KC)

     ! Volume fraction and Mach number dependent drag by Khalloufi and Capecelatro
     ! Augments single-particle drag coefficient of Singh (2022) with volume fraction
     ! correction of Sangani et al. (1991). Also replaced Mach number with a
     ! volume fraction-dependent effective Mach
     ! ----------------------------------------------------------------------

     if (Rep .ge. 10000.0_WP)                                                                &
          call die('compute_drag_correction: KC model is not intended for Rep > 10k')

     ! Effective Mach number
     Ma = Ma * (1.0_WP + 0.33_WP * (1.0_WP - exp(-34.0_WP * vp)))

     omega_exp = 0.74_WP
     delt_o = 9.4_WP
     Ao = 0.271616115889543_WP
     alpha_o = 0.3555_WP
     alpha_hoc = 0.9071_WP
     if (Ma .lt. 1.0_WP) then
        u_rat = 1.0_WP
        T_s_rat = 1.0_WP
        Ma_s = Ma
        alpha = 1.0_WP
        C2 = 1.0_WP
     else
        u_rat = (2.0_WP + (ratioOfSpecificHeats - 1.0_WP) *Ma**2) /                          &
             ((ratioOfSpecificHeats+1.0_WP)*Ma**2)   
        T_s_rat = ((ratioOfSpecificHeats - 1.0_WP) *Ma**2 + 2.0_WP) *                        &
             (2.0_WP * ratioOfSpecificHeats * Ma**2 - (ratioOfSpecificHeats - 1.0_WP)) /     &
             ((ratioOfSpecificHeats + 1.0_WP)**2 * Ma**2)
        Ma_s = sqrt(((ratioOfSpecificHeats - 1.0_WP) * Ma**2 +2.0_WP) /                      &
             (2.0_WP * ratioOfSpecificHeats *Ma**2 - (ratioOfSpecificHeats - 1.0_WP)))
        alpha = 1.0_WP / (alpha_o * Ma + 1.0_WP - alpha_o)
        alpha_inf = 1.0_WP / alpha_o/ Ma
        C2 = (0.9_WP - Ao * (1 + (ratioOfSpecificHeats - 1.0_WP)**2 / 4.0_WP /               &
             ratioOfSpecificHeats)**(ratioOfSpecificHeats / (ratioOfSpecificHeats - 1.0_WP)))&
             / (1.0_WP - alpha_inf * (ratioOfSpecificHeats - 1.0_WP) /                       &
             (ratioOfSpecificHeats + 1.0_WP))
     end if
     theta_singh = (1.0_WP + (ratioOfSpecificHeats - 1.0_WP) * Ma_s**2 / 2.0_WP)**           &
          (ratioOfSpecificHeats / (ratioOfSpecificHeats - 1.0_WP))
     Re_tilde = Rep * (1.0_WP / alpha**2 / T_s_rat)**omega_exp *                             &
          theta_singh**((ratioOfSpecificHeats + 1.0_WP) / 2.0_WP / ratioOfSpecificHeats-     &
          (ratioOfSpecificHeats - 1.0_WP) * omega_exp / ratioOfSpecificHeats)
     Cd_c = C2 * (1.0_WP - alpha * u_rat) + (Ao * theta_singh) * (1.0_WP + delt_o /          &
          sqrt(Re_tilde))**2
     Wr = Ma**(2.0_WP * omega_exp) / Rep

     if (Ma .gt. 1.0_WP) then
        Wr = Wr * (1.0_WP + 1.0_WP / T_s_rat)**omega_exp
     end if

     f_low_Kn = 1.0_WP / (1.0_WP + Kn * (2.514_WP + 0.8_WP * exp(-0.55_WP / Kn))) /          &
          (1.0_WP + alpha_hoc * Wr)
     s = sqrt(0.5_WP * ratioOfSpecificHeats) * Ma
     Cd_FM = (1.0_WP + 2.0_WP * s**2) * exp(-s**2) / (s**3 * sqrt(pi)) +                     &
          (4.0_WP * s**4 + 4.0_WP*s**2 - 1.0_WP) * erf(s) / (2.0_WP * s**4) +                &
          2.0_WP / (3.0_WP*s) * sqrt(pi)

     a_singh = 1.8_WP

     Br_singh = Wr * (Ma**(2.0_WP * omega_exp - 1.0_WP) + 1.0_WP) /                          &
          (Ma**(2.0_WP * omega_exp - 1.0_WP))
     CD = Cd_FM * Br_singh**a_singh / (1.0_WP + Br_singh**a_singh) + Cd_c * f_low_Kn /       &
          (1.0_WP + Br_singh**a_singh)

     ! Convert to drag correction and augment by Sangani
     dragCorrection = CD * Rep / 24.0_WP * (1.0_WP + 2.0_WP * vp) / (1.0_WP - vp)**2

  case(OSNES)

     ! Osnes, A., M. Vartdal, M.  Khalloufi, J. Capecelatro, and S. Balachandar. 
     ! Comprehensive quasi-steady force correlations for 
     ! compressible flow through random particle suspensions (2023) 
     ! International Journal of Multiphase Flow 
     ! -------------------------------------------------------------------------

     if (Rep .le. 45.0_WP) then
        ! Rarefraction-dominated regime
        Knp = sqrt(0.5_WP * pi * ratioOfSpecificHeats) * Ma / Rep
        fKn = 1.0_WP / (1.0_WP + Knp*(2.514_WP + 0.8_WP*exp(-0.55_WP/Knp)))
        CD1 = (1.0_WP + 0.15_WP*Rep**(0.687_WP)) * fKn
        sDrag = Ma * sqrt(0.5_WP * ratioOfSpecificHeats)
        if (Ma .le. 1.0_WP) then
           JM = 2.26_WP - 0.1_WP/Ma + 0.14_WP/Ma**3
        else
           JM = 1.6_WP + 0.25_WP/Ma + 0.11_WP/Ma**2 + 0.44_WP/Ma**3
        end if
        CD2 = (1.0_WP + 2.0_WP*sDrag**2) * exp(-sDrag**2) / (sDrag**3*sqrt(pi) +             &
             epsilon(1.0_WP)) + (4.0_WP*sDrag**4 + 4.0_WP*sDrag**2 - 1.0_WP) *               &
             erf(sDrag) / (2.0_WP*sDrag**4 + epsilon(1.0_WP)) +                              &
             2.0_WP / (3.0_WP * sDrag + epsilon(1.0_WP)) * sqrt(pi)
        CD2 = CD2 / (1.0_WP + (CD2/JM - 1.0_WP) * sqrt(Rep/45.0_WP))
        dragCorrection = CD1 / (1.0_WP + Ma**4) +                                            &
             Rep / 24.0_WP * Ma**4 * CD2 / (1.0_WP + Ma**4)
        CD3 = dragCorrection * 24.0_WP/Rep
     else
        ! Compression-dominated regime
        if (Ma .lt. 1.5_WP) then
           CM = 1.65_WP + 0.65_WP * tanh(4.0_WP*Ma - 3.4_WP)
        else
           CM = 2.18_WP - 0.13_WP * tanh(0.9_WP*Ma - 2.7_WP)
        end if
        if (Ma .lt. 0.8_WP) then
           GM = 166.0_WP*Ma**3 + 3.29_WP*Ma**2 - 10.9_WP*Ma + 20.0_WP
        else
           GM = 5.0_WP + 40.0_WP*Ma**(-3.0_WP)
        end if

        if (Ma .lt. 1.0_WP) then
           HM = 0.0239_WP*Ma**3 + 0.212_WP*Ma**2 - 0.074_WP*Ma + 1.0_WP
        else
           HM = 0.93_WP + 1.0_WP / (3.5_WP + Ma**5)
        end if
        CD3 = 24.0_WP/Rep*(1.0_WP + 0.15_WP * Rep**(0.687_WP))*HM + Rep / 24.0_WP *           &
             0.42_WP*CM / (1.0_WP+42500.0_WP/Rep**(1.16_WP*CM) + GM/sqrt(Rep))
     end if
     vp = min(0.4_WP, vp)
     CD3 = CD3/(1.0_WP-vp)
     Rem = (1.0_WP-vp)*Rep
     b1 = (5.81_WP*vp/(1.0_WP-vp)**2.0_WP+0.48_WP*vp**(1.0_WP/3.0_WP)/(1.0_WP-vp)**3.0_WP)
     b2 = (1.0_WP-vp)*(vp**3.0_WP*Rem*(0.95_WP+0.61_WP*vp**3.0_WP/(1.0_WP-vp)**2.0_WP))

     b1 = b1*(1.0_WP-vp)
     b2 = b2*(1.0_WP-vp)

     b3 = min((Ma/0.05_WP)**0.5_WP, 1.0_WP)*                                                  &
          (5.65_WP*vp-22.0_WP*vp**2.0_WP+23.4_WP*vp**3.0_WP)                                   &
          *(1.0_WP+tanh((Ma-(0.65_WP-0.24_WP*vp))/0.35_WP))
     dragCorrection = CD3*Rep/24.0_WP + b1 + b2 + b3*Rep/24.0_WP

  end select

  return

end subroutine compute_drag_correction


! ===================================================== !
! Overwrite particle velocity with local fluid velocity !
! ===================================================== !
subroutine overwrite_particle_velocity

  ! External modules
  use particle_exchange
  use state, only : update_state

  implicit none

  ! Local variables
  integer :: i

  if (.not. useParticles) return

  ! Prepare the fluid arrays
  call update_state
  call prepare_cartesian_arrays

  ! Interpolate fluid velocity to particle location
  do i = 1, nParticles
     call interpolate_fluid_to_particle(particles(i)%gridIndex, particles(i)%position,       &
          velocity = particles(i)%velocity)
  end do

  return
end subroutine overwrite_particle_velocity
