module monitor_homogeneous

  ! External modules
  use monitor

  implicit none

end module monitor_homogeneous


! =================================================== !
! Setup the routine to monitor homogeneous turbulence !
! =================================================== !
subroutine monitor_homogeneous_setup

  ! Internal modules
  use monitor_homogeneous

  ! External modules
  use parser
  use simulation_flags
  use solver_options
  use geometry

  implicit none

  ! Local variables
  integer :: n

  if (.not. all(isPeriodic(1:nDimensions))) return
  if (.not. useViscosity) return

  ! Set the monitor names
  n = 16
  if (useIBM) n = n + 2
  call monitor_create('homogeneous', n)
  call monitor_set_header(1 , 'TKE', 'r')
  call monitor_set_header(2 , 'urms', 'r')
  call monitor_set_header(3 , 'R11', 'r')
  call monitor_set_header(4 , 'R22', 'r')
  call monitor_set_header(5 , 'R33', 'r')
  call monitor_set_header(6 , 'R12', 'r')
  call monitor_set_header(7 , 'R13', 'r')
  call monitor_set_header(8 , 'R23', 'r')
  call monitor_set_header(9 , 'eta', 'r')
  call monitor_set_header(10, 'Re_lambda', 'r')
  call monitor_set_header(11, 'Ma_t', 'r')
  call monitor_set_header(12 , 'enstrophy', 'r')
  call monitor_set_header(13, "<u'div(rhouu)>", 'r')
  call monitor_set_header(14, "<p*div(u')>", 'r')
  call monitor_set_header(15, "dissipation", 'r')
  call monitor_set_header(16, 'num-diss', 'r')
  if (useIBM) then
     call monitor_set_header(17, 'drag-prod-P', 'r')
     call monitor_set_header(18, 'drag-prod-V', 'r')
  end if

  return
end subroutine monitor_homogeneous_setup


! ========================================= !
! Compute homogeneous turbulence statistics !
! ========================================= !
subroutine monitor_homogeneous_timestep

  ! Internal modules
  use monitor_homogeneous

  ! External modules
  use parser
  use parallel
  use solver_options
  use geometry
  use state
  use state_functions
  use time_info
  use grid
  use grid_functions
  use first_derivative
  use dissipation

  implicit none

  ! Local variables
  integer :: i, j, k
  real(WP) :: volume, meanDensity, favreViscosity, favreVelocity(3), enstrophy, TKE, urms,   &
       vorticity(nGridPoints, 3), Rij(3, 3), uPrime(3), divU, eta, lambda, pressureStrain,   &
       viscousDissipation, meanSoundSpeed, numericalDissipation, dragP, dragVisc, uPrimeConv
  real(WP), dimension(nGridPoints, nDimensions, nDimensions) :: fluxes1, fluxes2
  real(WP), dimension(nGridPoints, nDimensions) :: convectiveFlux, fluidStress
  real(WP), dimension(nGridPoints) :: divUprimeTau, divUprimeP

  if (.not. all(isPeriodic(1:nDimensions))) return
  if (.not. useViscosity) return

  ! Compute vorticity
  call compute_vorticity(velocityGradient, vorticity)

  ! Compute mean statistics
  volume = 0.0_WP
  meanDensity = 0.0_WP
  meanSoundSpeed = 0.0_WP
  favreViscosity = 0.0_WP
  favreVelocity = 0.0_WP
  enstrophy = 0.0_WP
  do i = 1, nGridPoints
     volume = volume + gridNorm(i, 1)
     meanDensity = meanDensity + conservedVariables(i, 1) * gridNorm(i, 1)
     meanSoundSpeed = meanSoundSpeed + sqrt(ratioOfSpecificHeats * pressure(i, 1) *          &
          specificVolume(i, 1)) * gridNorm(i, 1)
     favreViscosity = favreViscosity + dynamicViscosity(i, 1) * gridNorm(i, 1)
     favreVelocity(1:nDimensions) = favreVelocity(1:nDimensions) +                           &
          conservedVariables(i,2:nDimensions+1) * gridNorm(i, 1)
     enstrophy = enstrophy + 0.5_WP * conservedVariables(i, 1) * sum(vorticity(i,:)**2) *    &
          gridNorm(i, 1)
  end do

  ! Domain & Favre average
  call parallel_sum(volume)
  call parallel_sum(meanDensity); meanDensity = meanDensity / volume
  call parallel_sum(meanSoundSpeed); meanSoundSpeed = meanSoundSpeed / volume
  call parallel_sum(favreViscosity); favreViscosity = favreViscosity / meanDensity / volume
  call parallel_sum(enstrophy); enstrophy = enstrophy / meanDensity / volume
  call parallel_sum(favreVelocity); favreVelocity = favreVelocity / meanDensity / volume

  ! Compute convetive term
  select case (nDimensions)
  case (1)
     fluxes1(:,1,1) = conservedVariables(:,2) * velocity(:,1)
  case (2)
     fluxes1(:,1,1) = conservedVariables(:,2) * velocity(:,1)
     fluxes1(:,2,1) = conservedVariables(:,2) * velocity(:,2)
     fluxes1(:,1,2) = conservedVariables(:,2) * velocity(:,2)
     fluxes1(:,2,2) = conservedVariables(:,3) * velocity(:,2)
  case (3)
     fluxes1(:,1,1) = conservedVariables(:,2) * velocity(:,1)
     fluxes1(:,2,1) = conservedVariables(:,2) * velocity(:,2)
     fluxes1(:,3,1) = conservedVariables(:,2) * velocity(:,3)
     fluxes1(:,1,2) = conservedVariables(:,2) * velocity(:,2)
     fluxes1(:,2,2) = conservedVariables(:,3) * velocity(:,2)
     fluxes1(:,3,2) = conservedVariables(:,3) * velocity(:,3)
     fluxes1(:,1,3) = conservedVariables(:,2) * velocity(:,3)
     fluxes1(:,2,3) = conservedVariables(:,3) * velocity(:,3)
     fluxes1(:,3,3) = conservedVariables(:,4) * velocity(:,3)
  end select
  call transform_fluxes(fluxes1, fluxes2)
  do i = 1, nDimensions
     call first_derivative_apply(i, fluxes2(:,:,i))
  end do
  do i = 1, nDimensions
     convectiveFlux(:,i) = -sum(fluxes2(:,i,:), dim = 2) * jacobian(:, 1)
  end do

  ! Compute stress tensor for drag production
  if (useIBM) then
     select case (nDimensions)
     case (1)
        fluxes1(:,1,1) = -pressure(:,1)
        if (useViscosity) then
           fluxes2(:,1,1) = stressTensor(:,1)
        end if
     case (2)
        fluxes1(:,1,1) = -pressure(:,1)
        fluxes1(:,2,2) = -pressure(:,1)
        if (useViscosity) then
           fluxes2(:,1,1) = stressTensor(:,1)
           fluxes2(:,1,2) = stressTensor(:,2)
           fluxes2(:,2,1) = stressTensor(:,3)
           fluxes2(:,2,2) = stressTensor(:,4)
        end if
     case (3)
        fluxes1(:,1,1) = -pressure(:,1)
        fluxes1(:,2,2) = -pressure(:,1)
        fluxes1(:,3,3) = -pressure(:,1)
        if (useViscosity) then
           fluxes2(:,1,1) = stressTensor(:,1)
           fluxes2(:,1,2) = stressTensor(:,2)
           fluxes2(:,1,3) = stressTensor(:,3)
           fluxes2(:,2,1) = stressTensor(:,4)
           fluxes2(:,2,2) = stressTensor(:,5)
           fluxes2(:,2,3) = stressTensor(:,6)
           fluxes2(:,3,1) = stressTensor(:,7)
           fluxes2(:,3,2) = stressTensor(:,8)
           fluxes2(:,3,3) = stressTensor(:,9)
        end if
     end select

     ! Compute div(u'P)
     do i = 1, nGridPoints
        uPrime(1:nDimensions) = velocity(i,1:nDimensions) - favreVelocity(1:nDimensions)
        do j = 1, nDimensions
           fluidStress(i,j) = sum(fluxes1(i,j,:) * uPrime(1:nDimensions))
        end do
     end do
     call divergence(fluidStress, divUprimeP)

     ! Compute div(u'tau)
     do i = 1, nGridPoints
        uPrime(1:nDimensions) = velocity(i,1:nDimensions) - favreVelocity(1:nDimensions)
        do j = 1, nDimensions
           fluidStress(i,j) = sum(fluxes2(i,j,:) * uPrime(1:nDimensions))
        end do
     end do
     call divergence(fluidStress, divUprimeTau)
  end if

  ! Compute higher-order statistics
  Rij = 0.0_WP
  pressureStrain = 0.0_WP
  viscousDissipation = 0.0_WP
  lambda = 0.0_WP
  numericalDissipation = 0.0_WP
  uPrimeConv = 0.0_WP
  dragP = 0.0_WP
  dragVisc = 0.0_WP
  do i = 1, nGridPoints

     ! Compute fluctating velocity (about Favre average)
     uPrime(1:nDimensions) = velocity(i,1:nDimensions) - favreVelocity(1:nDimensions)

     ! Compute Reynolds stresses
     do j = 1, nDimensions
        do k = 1, nDimensions
           Rij(j,k) = Rij(j,k) + conservedVariables(i,1) * uPrime(j) * uPrime(k) * gridNorm(i, 1)
        end do
     end do
     
     ! Compute local divergence
     select case (nDimensions)
     case (1)
        divU = velocityGradient(i,1)
        lambda = lambda + velocityGradient(i,1)**2 * gridNorm(i, 1)
     case (2)
        divU = velocityGradient(i,1) + velocityGradient(i,4)
        lambda = lambda + (velocityGradient(i,1)**2 + velocityGradient(i,4)**2) *            &
             gridNorm(i, 1)
     case (3)
        divU = velocityGradient(i,1) + velocityGradient(i,5) + velocityGradient(i,9)
        lambda = lambda + (velocityGradient(i,1)**2 + velocityGradient(i,5)**2 +             &
             velocityGradient(i,9)**2) * gridNorm(i, 1)
     end select

     ! Pressure-strain correlation
     pressureStrain = pressureStrain + pressure(i,1) * divU * gridNorm(i, 1)

     ! Viscous dissipation
     if (useViscosity) then
        viscousDissipation = viscousDissipation - sum(stressTensor(i,:) *                    &
             velocityGradient(i,:)) * gridNorm(i,1)
     end if

     ! Numerical dissipation
     if (useDissipation) then
        numericalDissipation = numericalDissipation + sum( uPrime(1:nDimensions) *           &
             dissipationSource(i, 2:nDimensions+1) * gridNorm(i, 1) )
     end if

     ! Convective flux
     uPrimeConv = uPrimeConv + sum(uPrime(1:nDimensions) * convectiveFlux(i,1:nDimensions)) *&
          gridNorm(i,1)

     ! Drag production
     if (useIBM) then
        dragP = dragP + divUprimeP(i) * gridNorm(i, 1)
        dragVisc = dragVisc + divUprimeTau(i) * gridNorm(i, 1)
     end if

  end do

  ! Domain & Favre average
  call parallel_sum(lambda); lambda = lambda / volume
  call parallel_sum(pressureStrain); pressureStrain = pressureStrain / volume
  call parallel_sum(viscousDissipation); viscousDissipation = viscousDissipation / volume
  call parallel_sum(numericalDissipation); numericalDissipation = numericalDissipation / volume
  call parallel_sum(uPrimeConv); uPrimeConv = uPrimeConv / volume
  do i = 1, nDimensions
     do j = 1, nDimensions
        call parallel_sum(Rij(i,j)); Rij(i,j) = Rij(i,j) / volume / meanDensity
     end do
  end do
  if (useIBM) then
     call parallel_sum(dragP); dragP = dragP / volume
     call parallel_sum(dragVisc); dragVisc = dragVisc / volume
  end if

  ! Turbulent kinetic energy and root-mean-square velocity
  TKE = 0.5_WP * (Rij(1,1) + Rij(2,2) + Rij(3,3))
  urms = sqrt(2.0_WP * TKE / real(nDimensions, WP))

  ! Compute Taylor microscale
  lambda = sqrt(2.0_WP * TKE) / sqrt(lambda)

  ! Kolmogorov length scale
  eta =  (favreViscosity**3 / abs(viscousDissipation))**(0.25_WP)

  ! Set the homogeneous turbulence parameters
  call monitor_select('homogeneous')
  call monitor_set_single_value(1 , TKE)
  call monitor_set_single_value(2 , urms)
  call monitor_set_single_value(3 , Rij(1,1))
  call monitor_set_single_value(4 , Rij(2,2))
  call monitor_set_single_value(5 , Rij(3,3))
  call monitor_set_single_value(6 , Rij(1,2))
  call monitor_set_single_value(7 , Rij(1,3))
  call monitor_set_single_value(8 , Rij(2,3))
  call monitor_set_single_value(9 , eta)
  call monitor_set_single_value(10, urms * lambda / favreViscosity)
  call monitor_set_single_value(11, sqrt(2.0_WP * TKE) / meanSoundSpeed)
  call monitor_set_single_value(12, enstrophy)
  call monitor_set_single_value(13, uPrimeConv)
  call monitor_set_single_value(14, pressureStrain)
  call monitor_set_single_value(15, viscousDissipation)
  call monitor_set_single_value(16, numericalDissipation)
  if (useIBM) then
     call monitor_set_single_value(17, dragP)
     call monitor_set_single_value(18, dragVisc)
  end if

  return
end subroutine monitor_homogeneous_timestep
