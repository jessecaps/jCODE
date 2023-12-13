module rayleigh_taylor

  ! External modules
  use precision
  use string
  use parser
  use simulation_flags
  use solver_options
  use geometry
  use grid
  use grid_functions

  implicit none

  ! Global variables
  integer :: nx, ny, nz
  real(WP) :: Lx, Ly, Lz, dx, dy, dz

end module rayleigh_taylor

subroutine rayleigh_taylor_grid

  ! Internal modules
  use rayleigh_taylor

  ! External modules
  use parallel

  implicit none

  ! Local variables
  integer :: i, j, k
  real(WP) :: sigma, b, c, minMeshsize, maxMeshsize, y1,  y2
  real(WP), allocatable, dimension(:) :: s, g
  logical :: stretchGrid

  ! Simplify
  nx = globalGridSize(1)
  ny = globalGridSize(2)
  nz = globalGridSize(3)

  ! Read in the grid size
  call parser_read('Lx', Lx, 0.0_WP)
  call parser_read('Ly', Ly, 0.0_WP)
  call parser_read('Lz', Lz, 0.0_WP)

  ! Compute the grid spacing
  if (periodicityType(1) .eq. PLANE) then
     dx = Lx / real(nx, WP)
  else
     dx = Lx / real(nx-1, WP)
  end if
  if (periodicityType(2) .eq. PLANE) then
     dy = Ly / real(ny, WP)
  else
     dy = Ly / real(ny-1, WP)
  end if
  if (periodicityType(3) .eq. PLANE) then
     dz = Lz / real(nz, WP)
  else
     dz = Lz / real(nz-1, WP)
  end if

  ! Should we stretch the mesh?
  call parser_read('grid stretching', stretchGrid, .false.)

  ! Generate the grid
  do k = iStart(3), iEnd(3)
     do j = iStart(2), iEnd(2)
        do i = iStart(1), iEnd(1)

           ! Create X
           if (nx .gt. 1) then
              if (periodicityType(1) .eq. PLANE) then
                 coordinates(grid_index(i,j,k), 1) = (Lx - dx) *  real(i - 1, WP) /          &
                      real(nx - 1, WP) - 0.5_WP * Lx
              else
                 coordinates(grid_index(i,j,k), 1) = Lx * real(i - 1, WP) /                  &
                      real(nx - 1, WP) - 0.5_WP * Lx
              end if
           end if

           ! Create Y
           if (ny .gt. 1 .and. .not. stretchGrid) then
              coordinates(grid_index(i,j,k), 2) = Ly * real(j - 1, WP) /                     &
                   real(ny - 1, WP) - 0.5_WP * Ly
           end if

           ! Create Z
           if (nz .gt. 1) then
              if (periodicityType(3) .eq. PLANE) then
                 coordinates(grid_index(i,j,k), 3) = (Lz - dz) *  real(k - 1, WP) /          &
                      real(nz - 1, WP) - 0.5_WP * Lz
              else
                 coordinates(grid_index(i,j,k), 3) = Lz * real(k - 1, WP) /                  &
                      real(nz - 1, WP) - 0.5_WP * Lz
              end if
           end if

        end do
     end do
  end do

  ! Grid stretching
  if (stretchGrid) then
     ! Grid stretching parameters
     sigma = 0.21_WP
     b = 12.0_WP
     c = 0.6_WP

     ! Create uniform spacing
     allocate(s(ny))
     do j = 1, ny
        s(j) = real(j - 1, WP) / real(ny - 1, WP)
     end do

     ! Compute mapping g(s)
     allocate(g(ny))
     call mapping_function(s, b, c, sigma, g)

     ! Find min/max spacing
     minMeshsize =  huge(1.0_WP)
     maxMeshsize = -huge(1.0_WP)

     do k = iStart(3), iEnd(3)
        do j = iStart(2), iEnd(2)
           do i = iStart(1), iEnd(1)
              ! Create y
              coordinates(grid_index(i,j,k), 2) = 0.5_WP * Ly * (1.0_WP + g(j)) - 0.5_WP * Ly

              ! Find min/max spacing
              if (j .gt. iStart(2) + 1) then
                 y1 = coordinates(grid_index(i,j  ,k), 2)
                 y2 = coordinates(grid_index(i,j-1,k), 2)
                 minMeshsize = min(minMeshsize, abs(y2 - y1))
                 maxMeshsize = max(maxMeshsize, abs(y2 - y1))
              end if
           end do
        end do
     end do
     call parallel_max(maxMeshsize)
     call parallel_min(minMeshsize)
     if (iRank .eq. iRoot) then
        print *
        print *, 'min/max y-spacing:', minMeshsize, maxMeshsize
        print *        
     end if

     deallocate(s)
     deallocate(g)
  end if

  ! Setup the grid metrics
  call grid_metrics_setup

contains

  ! Mapping functional for grid stretching
  ! --------------------------------------
  subroutine mapping_function(s, b, c, sigma, g)

    ! External modules
    use math, only : pi

    implicit none

    real(WP), intent(in) :: s(:), b, c, sigma
    real(WP), intent(out) :: g(size(s))

    g = ((s - 0.5_WP) * (1.0_WP + 2.0_WP * b) - b * sigma *                                  &
         (exp(- ((s - 0.5_WP + c) / sigma) ** 2) / sqrt(pi) +                                &
         ((s - 0.5_WP + c) / sigma) * erf((s - 0.5_WP + c) / sigma) -                        &
         exp(- ((s - 0.5_WP - c) / sigma) ** 2) / sqrt(pi) -                                 &
         ((s - 0.5_WP - c) / sigma) * erf((s - 0.5_WP - c) / sigma))) /                      &
         (0.5_WP + b - b * sigma * (exp(- ((0.5_WP + c) / sigma) ** 2) /                     &
         sqrt(pi) + ((0.5_WP + c) / sigma) * erf((0.5_WP + c) / sigma) -                     &
         exp(- ((0.5_WP - c) / sigma) ** 2) / sqrt(pi) - ((0.5_WP - c) / sigma) *            &
         erf((0.5_WP - c) / sigma)))

    return
  end subroutine mapping_function

end subroutine rayleigh_taylor_grid

subroutine rayleigh_taylor_data

  ! Internal modules
  use rayleigh_taylor

  ! External modules
  use math, only : pi
  use parallel
  use parser
  use state, only : conservedVariables
  use equation_of_state
  use random
  use string

  implicit none

  ! Local variables
  integer  :: i, j, m, nModes, ierror
  real(WP) :: rand, T0(1), P0, Rho0, x, y, z, massDiffusivity(1,1), Pmin, massFraction1,        &
       molecularWeightOfMixture, initialThicknessInverse, minAmplitude, maxAmplitude, eta,      &
       constant1, constant2, constant3, coefficient1, coefficient2                              
  real(WP), allocatable :: Wi(:), densityGradient(:,:), pressure(:),                            &
       waveNumber(:), amplitude(:), phase(:)
  logical :: isothermalInitialConditions, randomPerturbation, velocityCorrection
  character(len = str_medium) :: key

  ! Set the number of conserved variables and allocate the array
  nSpecies = 1
  nUnknowns = nDimensions + 2 + nSpecies
  allocate(conservedVariables(nGridPoints, nUnknowns))

  ! Initialize solver options
  call solver_options_setup
  
  ! Get molecular weights
  if (equationOfState .eq. IDEAL_GAS_MIXTURE) then
     allocate(Wi(nSpecies+1))
     Wi = molecularWeightInverse
  else
     call die('raylaigh_taylor_data: Initializationr requires ideal gas mixture')
  end if
  
  ! Initialize the random number generator
  call random_setup

  ! Check in isothermal initial conditions
  call parser_read('isothermal initial conditions', isothermalInitialConditions, .true.)

  ! Read in reference pressure
  call parser_read('reference pressure', P0, 1.0_WP / ratioOfSpecificHeats)
 
  if (isothermalInitialConditions) then
     ! Read in reference temperature
     call parser_read('reference temperature', T0(1), 1.0_WP / (ratioOfSpecificHeats - 1.0_WP))
  
     ! Calculate the reference density: fluid interface was chosen as the reference
     Rho0 = ratioOfSpecificHeats * P0 / ((ratioOfSpecificHeats - 1.0_WP) * T0(1))                &
          / (0.5_WP*Wi(1) + 0.5_WP*Wi(2)) 
  else
     ! Calculate the reference density and temperature, density equals to the molecular weight
     ! and fluid interface was chosen as the reference
     Rho0 = 1.0_WP / (0.5_WP*Wi(1) + 0.5_WP*Wi(2)) 
     T0 = ratioOfSpecificHeats * P0 / (ratioOfSpecificHeats - 1.0_WP)

     ! Estimate the minimum pressure value which is around the top wall
     Pmin = P0 - 0.5_WP * Ly * froudeNumberInverse / (1.0_WP*Wi(1) + 0.0_WP*Wi(2))
     write(key,*) Pmin
     if (Pmin .lt. 0.0_WP)       call die("rayleigh_taylor_data: minimum pressure value would be &
          &probably negavitve: " // trim(key) //  "!")

  end if

  ! Read in initial diffusion thickness
  call parser_read('initial diffusion thickness', initialThicknessInverse)
  if (initialThicknessInverse .le. 1.0e-9_WP) call die("rayleigh_taylor_data: initial diffusion & 
       & thickness is too small!")
  InitialThicknessInverse = 1.0_WP / initialThicknessInverse

  ! Read in perturbation properties 
  call parser_read('number of modes', nModes, 1)
  if (nModes .le. 0) call die("rayleigh_taylor_data: number of modes must be positive!") 
  
  allocate(waveNumber(nModes), amplitude(nModes), phase(nModes))

  randomPerturbation = .false.
  if (nModes .gt. 1) call parser_read('random perturbation', randomPerturbation, .false.)
  
  do i = 1, nModes
     write(key, '(A,I1.1)') 'perturbation wave number ', i
     call parser_read(trim(key), waveNumber(i), real(i, WP))
  end do

  call parser_read('maximum perturbation amplitude', maxAmplitude, 0.0_WP)
  if (maxAmplitude .lt. 0.0_WP)                                                                  &
       call die("rayleigh_taylor_data: maximum perturbation amplitude must NOT be negative!")

  if (randomPerturbation) then

     call parser_read('min perturbation amplitude', minAmplitude, 0.0_WP)
     if (minAmplitude .lt. 0.0_WP) call die("Minimum perturbation amplitude must NOT be negative!")
     
     if (iRank .eq. iRoot) then     
        do i=1, nModes
           call random_number(rand)
           amplitude(i) = minAmplitude + rand * (maxAmplitude - minAmplitude)
           call random_number(rand)
           phase(i) = rand * 2.0_WP * pi
        end do
     end if
     call MPI_BCAST(amplitude, nModes, MPI_REAL_WP, iRoot, comm, ierror)
     call MPI_BCAST(phase, nModes, MPI_REAL_WP, iRoot, comm, ierror)

  else 
     do i=1, nModes
        amplitude(i) = maxAmplitude
        phase(i) = 0.0_WP
     end do
  end if
  
  ! Compute constants to be used later
  constant1 = ( Wi(1) + Wi(2) ) / 2.0_WP
  constant2 = ( Wi(1) - Wi(2) ) / 2.0_WP
  constant3 = constant2 * initialThicknessInverse * 2.0_WP / sqrt(pi) 

  ! Compute mass fraction, density, and pressure
  allocate(pressure(nGridPoints))
  do i = 1, nGridPoints

     ! Coordinates
     x = coordinates(i, 1)
     y = coordinates(i, 2)
     z = 0.0_WP
     if (nz .gt. 1) z = coordinates(i, 3)

     ! Compute the perturbated fluid interface
     eta = 0.0_WP
     do m=1, nModes
        eta = eta + amplitude(m) * cos( waveNumber(m) * x + phase(m) )
     end do

     ! Compute the mass fraction distribution
     massFraction1 = 0.5_WP * (1.0_WP + erf( (y - eta) * initialThicknessInverse ))

     ! Compute the molecular weight of the mixture and its y derivative
     molecularWeightOfMixture = 1.0_WP / ( constant1 + constant2 *                              &
          erf( (y - eta) * initialThicknessInverse)) 

     ! Compute coefficient1
     coefficient1 = constant3 * molecularWeightOfMixture *                                      &
          exp( - (y - eta)**2 * initialThicknessInverse**2 )                                ! B     

     ! Compute the pressure and set the density
     if (isothermalInitialConditions) then     
        coefficient2 = ratioOfSpecificHeats * molecularWeightOfMixture /                        &
             ( (ratioOfSpecificHeats - 1.0_WP) * T0(1) )                                    ! G

        ! Compute the pressure
        pressure(i) = P0 * exp( - coefficient2 * froudeNumberInverse *                          &
             ( y - eta + (0.5_WP / initialThicknessInverse + eta) *                             &
             coefficient1  / initialThicknessInverse ) ) 

        ! Set the density
        conservedVariables(i,1) = coefficient2 * pressure(i)

     else

        ! Set the density
        conservedVariables(i,1) = molecularWeightOfMixture

        ! Compute the pressure 
        pressure(i) = P0 - molecularWeightOfMixture * froudeNumberInverse *                     &
             ( y - eta + (0.5_WP / initialThicknessInverse + eta) *                             &
             coefficient1  / initialThicknessInverse )
         
     end if

     ! Set the mass fraction
     conservedVariables(i,nDimensions+2+1) = conservedVariables(i,1) * massFraction1

  end do
  deallocate(waveNumber, amplitude, phase)
    
  ! Initialize the velocity fiels
  conservedVariables(:,2:nDimensions+1) = 0.0_WP

  ! Correct the velocity to prevent unphysical waves at the interface (Pooya Movahed, 2014)
  call parser_read('rayleigh taylor velocity correction', velocityCorrection, .true.)

  if (velocityCorrection) then
     if (.not. useViscosity)                                                                  &
          call die("rayleigh_taylor_data: velocity correction cannot be used with '           &
          & 'inviscid flow!")
       
     ! Compute the density gradient
     allocate(densityGradient(nGridPoints, nDimensions))
     call gradient(conservedVariables(:,1), densityGradient)

     ! Get the mass diffusivity
     call compute_transport_variables(T0, massDiffusivity = massDiffusivity)
      
     ! Set the momentum
     do i = 1, nGridPoints
        do j = 1, nDimensions
           conservedVariables(i,j+1) = - massDiffusivity(1,1) * densityGradient(i,j)
        end do
     end do
     deallocate(densityGradient)
  end if

  ! Set the energy   
  do i = 1, nGridPoints
     conservedVariables(i, nDimensions+2) = pressure(i) /                                    &
          (ratioOfSpecificHeats - 1.0_WP) + 0.5_WP *                                         &
          sum(conservedVariables(i,2:nDimensions+1)**2) /                                    &
          conservedVariables(i,1)
  end do
  deallocate(Wi, pressure)

  return
end subroutine rayleigh_taylor_data
