program test_adjoint

  ! External modules
  use precision
  use string
  use parser
  use parallel
  use geometry
  use grid
  use random
  use state
  use state_functions, only : compute_cartesian_viscous_fluxes
  use grid_functions
  use fileio
  use first_derivative

  implicit none

  ! Test mode
  integer, parameter ::                                                                      &
       NULL          = 0,                                                                    &
       DERIVATIVE    = 1,                                                                    &
       !DERIVATIVE_BC = 2,                                                                    &
       FLUX          = 3,                                                                    &
       INVISCID_FLUX = 4,                                                                    &
       VISCOUS_FLUX  = 5,                                                                    &
       DISSIPATION   = 6,                                                                    &
       COMBUSTION_   = 7,                                                                    &
       BOUNDARY      = 8

  ! Local variables
  integer :: i, k, nIterations, restartIteration, fileUnit, iostat, testMode
  real(WP) :: geometricGrowthFactor, initialActuationAmount, scalar1, scalar2, error,        &
       actuationAmount
  real(WP), allocatable :: rand(:),deltaConservedVariables(:,:),deltaPrimitiveVariables(:,:),&
       forwardSource(:,:), adjointSource(:,:), conservedVariables0(:,:), forwardSource0(:,:),&
       viscousFlux(:,:,:), temp(:,:), targetViscousFluxes(:,:,:)
  logical :: readData
  character(len = str_medium) :: input, gridFile, val, filename
  character(len = 1) :: xyz(3)

  ! Initialize parallel environment and parse the input file
  call parallel_init
  call parallel_get_inputname(input)
  call parser_init
  call parser_parsefile(input)
  disableManualDecomp = .true.

  if (iRank .eq. iRoot) then
     write (*,*)
     write (*,*) '========================'
     write (*,*) '| jCODE - Test Adjoint |'
     write (*,*) '========================'
     write (*,*)
  end if

  ! Initialize the random number generator
  call random_setup

  ! Initialize simulation flags
  call simulation_flags_setup

  ! Check if grid and solution files are available
  call parser_read('read data', readData, .false.)

  if (.not. readData) then
     ! Get the grid dimensions from input
     xyz(1) = 'x'; xyz(2) = 'y'; xyz(3) = 'z'
     nDimensions = 0
     do i = 1, 3
        call parser_read('n'//xyz(i), globalGridSize(i), 1)
        if (globalGridSize(i) .gt. 1) nDimensions = nDimensions + 1
     end do

     ! Get periodicity information from input
     periodicityType = NONE
     periodicLength = 0.0_WP
     do i = 1, nDimensions
        call parser_read('periodicity type in ' // xyz(i), val, '')
        periodicityType(i) = NONE
        if (trim(val) .eq. 'plane') then
           periodicityType(i) = PLANE
        else if (trim(val) .eq. 'overlap') then
           periodicityType(i) = OVERLAP
        else if (trim(val) .eq. 'polar') then
           periodicityType(i) = POLAR
        end if
        if (periodicityType(i) == PLANE) call parser_read('L' // xyz(i), periodicLength(i))
     end do
     isPeriodic = (periodicityType .ne. NONE)

     ! Is the domain curvilinear?
     call parser_read('curvilinear domain', isDomainCurvilinear, .false.)  
  else
     call get_dimensions(gridFile)
  end if

  ! Initialize the geometric parameters
  call geometry_setup

  ! Get/Set the total number of conserved variables
  if (.not. readData) then
     ! Get the number of species
     call parser_read('number of species', nSpecies, 0)

     ! Set the number of conserved variables
     nUnknowns = nDimensions + 2 + nSpecies
  else
     call get_nUnknowns
  end if

  ! Initialize solver options
  call solver_options_setup

  ! Setup the stecil operators
  call operator_setup

  ! Setup the grid
  call grid_setup

  ! Generate/Read the grid
  if (.not. readData) then
     call generate_grid
     ! Write the grid
     call parser_read('grid file', gridFile)
     call simulation_write(IO_GRID, gridFile)
  else
     call simulation_read(IO_GRID, trim(gridFile))
  end if

  ! Setup patches
  call grid_patch_setup
  
  ! Compute the grid metrics
  call grid_metrics_setup

  ! Setup the state
  call state_setup

  ! Setup the particles
  call particle_setup

  ! Get the conserved variables filename
  call parser_read('solution file to read', filename)

  if (.not.readData) allocate(rand(nGridPoints))

  ! Randomize/Read conserved variables
  if (.not. readData) then
     call random_number(rand)
     conservedVariables(:,1) = rand * (10.0_WP - 0.01_WP) + 0.01_WP ! ... density
     do i = 1, nDimensions
        call random_number(rand)
        rand = rand * (10.0_WP - (-10.0_WP)) - 10.0_WP ! ... velocity in i-direction
        conservedVariables(:,i+1) = conservedVariables(:,1) * rand
     end do
     call random_number(rand)
     rand = rand * (10.0_WP - 0.01_WP) + 0.01_WP ! ... pressure
     conservedVariables(:,nDimensions+2) = rand / (ratioOfSpecificHeats - 1.0_WP) + 0.5_WP * &
          sum(conservedVariables(:,2:nDimensions+1)**2, dim = 2) / conservedVariables(:,1)
     do k = 1, nSpecies
        call random_number(rand) ! ... mass fraction for k species
        conservedVariables(:,nDimensions+2+k) = conservedVariables(:,1) * rand
     end do

     call simulation_write(IO_FORWARD_STATE, trim(filename))
  else
     call simulation_read(IO_FORWARD_STATE, trim(filename))
  end if

  ! Update the state
  call update_state

  ! Get the delta conserved variables filename
  call parser_read('delta solution file to read', filename, 'data.delta')
  allocate(temp(nGridPoints, nUnknowns)); temp = conservedVariables ! ... tempelate for I/O

  ! Randomize/Read delta conserved variables
  if (.not. readData) then
     ! Randomize delta primitive variables: density, velocity, pressure, mass fraction
     allocate(deltaPrimitiveVariables(nGridPoints,nUnknowns))
     call random_number(deltaPrimitiveVariables)
     deltaPrimitiveVariables = deltaPrimitiveVariables * (1.0_WP - (-1.0_WP)) - 1.0_WP

     ! Randomize delta conserved variables
     allocate(deltaConservedVariables(nGridPoints,nUnknowns))  
     deltaConservedVariables(:,1) = deltaPrimitiveVariables(:,1) ! ... delta density
     do i = 1, nDimensions
        deltaConservedVariables(:,i+1) =                                                     &
             conservedVariables(:,i+1)/conservedVariables(:,1) * deltaPrimitiveVariables(:,1)&
             + conservedVariables(:,1) * deltaPrimitiveVariables(:,i+1)
     end do
     deltaConservedVariables(:,nDimensions+2) = deltaPrimitiveVariables(:,nDimensions+2) /   &
          (ratioOfSpecificHeats - 1.0_WP) + 0.5_WP * deltaPrimitiveVariables(:,1) *          &
          sum(conservedVariables(:,2:nDimensions+1)**2, dim = 2) / conservedVariables(:,1) + &
          sum(conservedVariables(:,2:nDimensions+1) *                                        &
          deltaPrimitiveVariables(:,2:nDimensions+1), dim = 2) * conservedVariables(:,1)
     do k = 1, nSpecies
        deltaConservedVariables(:,nDimensions+2+k) = conservedVariables(:,nDimensions+2+k) / &
             conservedVariables(:,1) * deltaPrimitiveVariables(:,1) +                        &
             conservedVariables(:,1) * deltaPrimitiveVariables(:,nDimensions+2+k)
        conservedVariables(:,nDimensions+2+k) = conservedVariables(:,1) * rand
     end do
     deallocate(deltaPrimitiveVariables)

     conservedVariables = deltaConservedVariables
     call simulation_write(IO_FORWARD_STATE, trim(filename))
  else
     call simulation_read(IO_FORWARD_STATE, trim(filename))
     deltaConservedVariables = conservedVariables
  end if

  ! Return back the conserved variables
  conservedVariables = temp
  deallocate(temp)

  ! Get the target state filename
  call parser_read('target state file', filename)

  ! Randomize/read target state variables
  if (.not. readData) then
     call random_number(rand)
     targetState(:,1) = rand * (10.0_WP - 0.01_WP) + 0.01_WP ! ... density
     do i = 1, nDimensions
        call random_number(rand)
        rand = rand * (10.0_WP - (-10.0_WP)) - 10.0_WP ! ... velocity in i-direction
        targetState(:,i+1) = targetState(:,1) * rand
     end do
     call random_number(rand)
     rand = rand * (10.0_WP - 0.01_WP) + 0.01_WP ! ... pressure
     targetState(:,nDimensions+2) = rand / (ratioOfSpecificHeats - 1.0_WP) + 0.5_WP *        &
          sum(targetState(:,2:nDimensions+1)**2, dim = 2) / targetState(:,1)
     do k = 1, nSpecies
        call random_number(rand) ! ... mass fraction for k species
        targetState(:,nDimensions+2+k) = targetState(:,1) * rand
     end do
     call simulation_write(IO_TARGET_STATE, trim(filename))
  else
     call simulation_read(IO_TARGET_STATE, trim(filename))
  end if

  if (allocated(rand)) deallocate(rand)

  ! Setup the boundary conditions after setting the target values
  call boundary_setup

  ! Compute the target viscous fluxes
  if (useViscosity) then
     allocate(targetViscousFluxes(nGridPoints, nUnknowns, nDimensions))
     call update_state(targetState)
     if (equationOfState .eq. IDEAL_GAS) then
        call compute_cartesian_viscous_fluxes(velocity, stressTensor, heatFlux,              &
             targetViscousFluxes, speciesFlux = speciesFlux)
     else
        call compute_cartesian_viscous_fluxes(velocity, stressTensor, heatFlux,              &
             targetViscousFluxes, speciesFlux, enthalpyFlux)
     end if
     ! Recompute state variables using the conserved variables
     call update_state
  end if

  ! Setup the control mollifier for the flame thickening case
  call controller_setup

  ! Recompute state variables using the conserved variables
  call update_state

  ! Get the adjoint filename
  call parser_read('adjoint file to write', filename)  

  ! Randomize/Read adjoint variables
  if (.not. readData) then
     call random_number(adjointVariables)
     call simulation_write(IO_ADJOINT_STATE, trim(filename))
  else
     call simulation_read(IO_ADJOINT_STATE, trim(filename))
  end if

  ! Randomly change the flame thickening factor and the flame efficiency factor
  !call random_number(flameThickening)
  !call random_number(flameEfficiency)
  call update_state

  call parser_read('number of control iterations', nIterations)
  if (nIterations.lt.0) call die('Number of control iterations must be a non-negative number')

  call parser_read('restart control iteration', restartIteration, 0)
  restartIteration = max(restartIteration, 0)

  filename = 'gradient_error.txt'
  if (iRank .eq. iRoot) then
     if (restartIteration .eq. 0) then
        fileUnit = iOpen()
        open(unit = fileUnit, file = trim(filename), action = 'write', status = 'unknown',   &
             iostat = iostat)
     else
        fileUnit = iOpen()
        open(unit = fileUnit, file = trim(filename), action = 'readwrite', status = 'old',   &
             position = 'rewind', iostat = iostat)
     end if
  end if

  call parallel_bc(iostat)
  if (iostat .ne. 0) call die('Failed to open file for writing!')

  if (nIterations .gt. 0) then
     call parser_read('initial actuation amount', initialActuationAmount)
     call parser_read('actuation amount geometric growth', geometricGrowthFactor,            &
          1.0_WP / (10.0_WP ** 0.25_WP)) ! ... default yields 4 iterations per decade
  end if

  ! Set the test mode
  call parser_read('test mode', val)
  select case(trim(val))
  case('')
     testMode = NULL
  case('derivative')
     testMode = DERIVATIVE
  case('flux')
     testMode = FLUX
  case('inviscid flux')
     testMode = INVISCID_FLUX
  case('viscous flux')
     testMode = VISCOUS_FLUX
  case('dissipation')
     testMode = DISSIPATION
  case('combustion')
     testMode = COMBUSTION_
  case('boundary')
     testMode = BOUNDARY
     if (useViscosity) allocate(viscousFlux(nGridPoints, nUnknowns, nDimensions))
  case default
     call die("unknown test model '" // trim(val) //"'")
  end select

  ! Allocate and initilize forward source term
  allocate(forwardSource(nGridPoints, nUnknowns))
  forwardSource = 0.0_WP

  ! Apply the forward operator
  select case(testMode)
  case(NULL)
     ! Nothing to do
  case(DERIVATIVE)
     call compute_divergence(FORWARD, forwardSource)
  case(FLUX)
     call add_fluxes_forward(forwardSource)
  case(INVISCID_FLUX)
     call add_inviscid_flux(FORWARD, forwardSource)
  case(VISCOUS_FLUX)
     call add_viscous_flux(FORWARD, forwardSource)
  case(DISSIPATION)
     call add_dissipation(FORWARD, forwardSource)
  case(COMBUSTION_)
     call combustion_source(FORWARD, forwardSource)
  case(BOUNDARY)
     ! Compute and send viscous fluxes to appropriate patches if required
     if (useViscosity) then
        call compute_cartesian_viscous_fluxes(velocity, stressTensor, heatFlux, viscousFlux, &
             speciesFlux, enthalpyFlux)
        call farfield_store_viscous_fluxes(viscousFlux)
     end if
     call boundary_sources(FORWARD, forwardSource)
  end select

  ! Allocate and initilize adjoint source term
  allocate(adjointSource(nGridPoints, nUnknowns))
  adjointSource = 0.0_WP

  ! Apply the adjoint operator
  select case(testMode)
  case(NULL)
     ! Nothing to do
  case(DERIVATIVE)
     call compute_divergence(ADJOINT, adjointSource)
  case(FLUX)
     call add_fluxes_adjoint(adjointSource)
  case(INVISCID_FLUX)
     call add_inviscid_flux(ADJOINT, adjointSource)
  case(VISCOUS_FLUX)
     call add_viscous_flux(ADJOINT, adjointSource)
  case(DISSIPATION)
     call add_dissipation(ADJOINT, adjointSource)
  case(COMBUSTION_)
     call combustion_source(ADJOINT, adjointSource)
  case(BOUNDARY)
     call boundary_sources(ADJOINT, adjointSource)
  end select

  ! Compute < R^dagger[Q,Q^dagger], delta Q >
  scalar1 = inner_product(adjointSource, deltaConservedVariables)
  deallocate(adjointSource)

  if (iRank .eq. iRoot) then
     write(fileUnit, '(A4,4A24)') 'i', 'alpha', 'alpha<Q^dagger,delta R>',                   &
          '<R^dagger,delta Q>', 'Error'
     write(fileUnit, '(I4,4(1X,SP,(SP,ES23.15E3)))') 0, 0.0_WP, 0.0_WP, scalar1, 0.0_WP
     flush(fileUnit)
  end if

  if (restartIteration .eq. 0) restartIteration = restartIteration + 1

  do i = 1, restartIteration - 1
     if (iRank .eq. iRoot)                                                                   &
          read(fileUnit, *, iostat = iostat) k, actuationAmount, scalar2, scalar1, error
     call parallel_bc(iostat)
     if (iostat .ne. 0) call die(trim(filename) //                                           &
          ': history is too short for the specified restart iteration!')
  end do

  ! Store the baseline conserved vairables and forward source
  allocate(conservedVariables0(nGridPoints, nUnknowns)); conservedVariables0 = conservedVariables
  allocate(forwardSource0(nGridPoints, nUnknowns)); forwardSource0 = forwardSource
  do i = restartIteration, restartIteration + nIterations - 1
     actuationAmount = initialActuationAmount * geometricGrowthFactor ** real(i - 1, WP)

     ! Finite difference on conserved variables
     conservedVariables = conservedVariables0 + actuationAmount * deltaConservedVariables

     ! Update state values
     call update_state

     ! Compute the forward source for the perturbed conserved variables
     forwardSource = 0.0_WP
     select case(testMode)
     case(NULL)
        ! Nothing to do
     case(DERIVATIVE)
        call compute_divergence(FORWARD, forwardSource)
     case(FLUX)
        call add_fluxes_forward(forwardSource)
     case(INVISCID_FLUX)
        call add_inviscid_flux(FORWARD, forwardSource)
     case(VISCOUS_FLUX)
        call add_viscous_flux(FORWARD, forwardSource)
     case(DISSIPATION)
        call add_dissipation(FORWARD, forwardSource)
     case(COMBUSTION_)
        call combustion_source(FORWARD, forwardSource)
     case(BOUNDARY)
        ! Compute and send viscous fluxes to appropriate patches if required
        if (useViscosity) then
           call compute_cartesian_viscous_fluxes(velocity, stressTensor, heatFlux,           &
                viscousFlux, speciesFlux, enthalpyFlux)
           call farfield_store_viscous_fluxes(viscousFlux)
        end if
        call boundary_sources(FORWARD, forwardSource)
     end select

     ! Compute < Q^dagger, delta R[Q,delta Q] >
     scalar2 = inner_product(adjointVariables, forwardSource - forwardSource0)

     ! Compute the error
     error = abs(scalar2 / actuationAmount + scalar1)

     if (iRank .eq. iRoot) then
        write(fileUnit, '(I4,4(1X,SP,(SP,ES23.15E3)))') i, actuationAmount, scalar2,         &
             scalar2 / actuationAmount, error
        flush(fileUnit)
     end if
  end do
  deallocate(conservedVariables0, forwardSource0)

  if (iRank .eq. iRoot) close(iclose(fileUnit))

  ! Clean up
  deallocate(deltaConservedVariables, forwardSource)
  if (allocated(viscousFlux)) deallocate(viscousFlux)
  if (allocated(targetViscousFluxes)) deallocate(targetViscousFluxes)

  ! Finalize the parallel environment
  call parallel_finalize

end program test_adjoint


! A subroutine for generating grid
! --------------------------------
subroutine generate_grid

  ! External modules
  use precision
  use parallel
  use string
  use parser
  use geometry
  use grid
  use grid_functions

  implicit none

  ! Local variables
  integer :: nx, ny, nz, i, j, k
  real(WP) :: Lx, Ly, Lz, dx, dy, dz, sigma, b, c, minMeshsize, maxMeshsize, x1,  x2
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
           if (.not. stretchGrid) then
              if (periodicityType(1) .eq. PLANE) then
                 coordinates(grid_index(i,j,k), 1) = (Lx - dx) *  real(i - 1, WP) /          &
                      real(nx - 1, WP)
              else
                 coordinates(grid_index(i,j,k), 1) = Lx * real(i - 1, WP) /                  &
                      real(nx - 1, WP)
              end if
           end if

           ! Create Y
           if (ny .gt. 1) then
              if (periodicityType(2) .eq. PLANE) then
                 coordinates(grid_index(i,j,k), 2) = (Ly - dy) *  real(j - 1, WP) /          &
                      real(ny - 1, WP)
              else
                 coordinates(grid_index(i,j,k), 2) = Ly * real(j - 1, WP) /                  &
                      real(ny - 1, WP)
              end if
           end if

           ! Create Z
           if (nz .gt. 1) then
              if (periodicityType(3) .eq. PLANE) then
                 coordinates(grid_index(i,j,k), 3) = (Lz - dz) *  real(k - 1, WP) /          &
                      real(nz - 1, WP)
              else
                 coordinates(grid_index(i,j,k), 3) = Lz * real(k - 1, WP) /                  &
                      real(nz - 1, WP)
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
     allocate(s(nx))
     do i = 1, nx
        s(i) = real(i - 1, WP) / real(nx - 1, WP)
     end do

     ! Compute mapping g(s)
     allocate(g(nx))
     call mapping_function(s, b, c, sigma, g)

     ! Find min/max spacing
     minMeshsize =  huge(1.0_WP)
     maxMeshsize = -huge(1.0_WP)

     do k = iStart(3), iEnd(3)
        do j = iStart(2), iEnd(2)
           do i = iStart(1), iEnd(1)
              ! Create y
              coordinates(grid_index(i,j,k), 1) = 0.5_WP * Lx * (1.0_WP + g(j))

              ! Find min/max spacing
              if (i .gt. iStart(1) + 1) then
                 x1 = coordinates(grid_index(i  ,j,k), 1)
                 x2 = coordinates(grid_index(i-1,j,k), 1)
                 minMeshsize = min(minMeshsize, abs(x2 - x1))
                 maxMeshsize = max(maxMeshsize, abs(x2 - x1))
              end if
           end do
        end do
     end do
     call parallel_max(maxMeshsize)
     call parallel_min(minMeshsize)
     if (iRank .eq. iRoot) then
        print *
        print *, 'min/max x-spacing:', minMeshsize, maxMeshsize
        print *
     end if

     deallocate(s)
     deallocate(g)
  end if

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

end subroutine generate_grid



! Compute divergence of the conserved variables
! ---------------------------------------------
subroutine compute_divergence(mode, source)

  ! External modules
  use state_jacobian
  use grid
  use first_derivative
  use state_functions
  use state

  implicit none

  ! Arguments
  integer, intent(in) :: mode
  real(WP), dimension(nGridPoints, nUnknowns), intent(inout) :: source

  ! Local variables
  integer :: i
  real(WP), allocatable :: temp1(:,:,:), temp2(:,:)

  select case (mode)

  case (FORWARD)

     allocate(temp1(nGridPoints, nUnknowns, nDimensions))
     allocate(temp2(nGridPoints, nUnknowns))

     ! Partial derivatives of conserved variables w.r.t. *computational* coordinates
     do i = 1, nDimensions
        temp1(:,:,i) = conservedVariables
        call first_derivative_apply(i, temp1(:,:,i))
     end do
     temp2 = sum(temp1, dim = 3) !... divergence of the conserved variables flux

     deallocate(temp1) !... no longer needed

     do i = 1, nUnknowns
        source(:,i) = source(:,i) - temp2(:,i)
     end do

     deallocate(temp2) !... no longer needed

  case (ADJOINT)

     allocate(temp1(nGridPoints, nUnknowns, nDimensions))
     allocate(temp2(nGridPoints, nUnknowns))

     ! Partial derivatives of adjoint variables w.r.t. *computational* coordinates
     do i = 1, nDimensions
        temp1(:,:,i) = adjointVariables
        call adjoint_first_derivative_apply(i, temp1(:,:,i))
     end do
     temp2 = sum(temp1, dim = 3) !... divergence of the adjoint flux

     deallocate(temp1) ! ... no longer needed

     do i = 1, nUnknowns
        source(:,i) = source(:,i) + temp2(:,i)
     end do

     deallocate(temp2) ! ... no longer needed

  end select

  return
end subroutine compute_divergence




! Add the inviscid fluxes
! -----------------------
subroutine add_inviscid_flux(mode, source)

  ! External modules
  use state_jacobian
  use grid
  use first_derivative
  use grid_functions, only : transform_fluxes
  use state_functions
  use state

  implicit none

  ! Arguments
  integer, intent(in) :: mode
  real(WP), dimension(nGridPoints, nUnknowns), intent(inout) :: source

  ! Local variables
  integer :: i, j
  real(WP), allocatable :: fluxes1(:,:,:), fluxes2(:,:,:), temp1(:,:,:),                     &
       localFluxJacobian1(:,:), localConservedVariables(:), localVelocity(:),                &
       localMassFraction(:), localMetricsAlongDirection1(:)

  select case (mode)

  case (FORWARD)

     allocate(fluxes1(nGridPoints, nUnknowns, nDimensions))
     allocate(fluxes2(nGridPoints, nUnknowns, nDimensions))

     ! Compute Cartesian form of inviscid fluxes
     call compute_cartesian_inviscid_fluxes(conservedVariables, velocity, pressure(:,1),     &
          fluxes1)

     ! Transform fluxes from Cartesian to contravariant form: `fluxes1` has the Cartesian
     ! form of total fluxes... upon return, `fluxes2` has the contravariant form
     call transform_fluxes(fluxes1, fluxes2)

     deallocate(fluxes1) !... no longer needed

     ! Take derivatives of fluxes
     do i = 1, nDimensions
        call first_derivative_apply(i, fluxes2(:,:,i))
     end do
     source = source - sum(fluxes2, dim = 3)

     deallocate(fluxes2) !... no longer needed

  case (ADJOINT)

     ! Partial derivatives of adjoint variables w.r.t. *computational* coordinates
     allocate(temp1(nGridPoints, nUnknowns, nDimensions))
     do i = 1, nDimensions
        temp1(:,:,i) = adjointVariables
        call adjoint_first_derivative_apply(i, temp1(:,:,i))
     end do

     allocate(localFluxJacobian1(nUnknowns, nUnknowns))
     allocate(localConservedVariables(nUnknowns))
     allocate(localVelocity(nDimensions))
     if (nSpecies .gt. 0) allocate(localMassFraction(nSpecies))
     allocate(localMetricsAlongDirection1(nDimensions))

     do j = 1, nGridPoints

        localConservedVariables = conservedVariables(j,:)
        localVelocity = velocity(j,:)
        if (nSpecies .gt. 0) localMassFraction = massFraction(j,:)

        do i = 1, nDimensions

           localMetricsAlongDirection1 = metrics(j,1+nDimensions*(i-1):nDimensions*i)

           call compute_jacobian_of_inviscid_flux(localConservedVariables,                   &
                localMetricsAlongDirection1, localFluxJacobian1,                             &
                specificVolume = specificVolume(j,1), velocity = localVelocity,              &
                pressure = pressure(j,1), massFraction = localMassFraction)

           source(j,:) = source(j,:) + matmul(transpose(localFluxJacobian1), temp1(j,:,i))

        end do

     end do

     if (allocated(localConservedVariables)) deallocate(localConservedVariables)
     if (allocated(localFluxJacobian1)) deallocate(localFluxJacobian1)
     if (allocated(localMassFraction)) deallocate(localMassFraction)
     deallocate(localMetricsAlongDirection1)
     deallocate(localVelocity)
     deallocate(temp1)

  end select

  do i = 1, nUnknowns
     source(:,i) = source(:,i) * jacobian(:,1)
  end do

  return
end subroutine add_inviscid_flux



! Add the viscous flux
! --------------------
subroutine add_viscous_flux(mode, source)

  ! External modules
  use state_jacobian
  use grid
  use first_derivative
  use grid_functions, only : transform_fluxes, laplacian
  use state_functions
  use state

  implicit none

  ! Arguments
  integer, intent(in) :: mode
  real(WP), dimension(nGridPoints, nUnknowns), intent(inout) :: source

  ! Local variables
  integer :: i, j, k
  real(WP), allocatable :: fluxes1(:,:,:), fluxes2(:,:,:), temp1(:,:,:), temp2(:,:),         &
       localFluxJacobian1(:,:), localFluxJacobian2(:,:), localConservedVariables(:),         &
       localVelocity(:), localMassFraction(:), localMetricsAlongDirection1(:),               &
       localMetricsAlongDirection2(:), localStressTensor(:), localHeatFlux(:),               &
       localEnthalpyFlux(:), localSpeciesFlux(:,:), localAdjointDiffusion(:,:)

  select case (mode)

  case (FORWARD)

     allocate(fluxes1(nGridPoints, nUnknowns, nDimensions))
     allocate(fluxes2(nGridPoints, nUnknowns, nDimensions))

     fluxes1 = 0.0_WP

     ! Compute Cartesian form of viscous fluxes if viscous terms are included
     if (useViscosity) then
        call compute_cartesian_viscous_fluxes(velocity, stressTensor, heatFlux,              &
             fluxes2, speciesFlux, enthalpyFlux)

        ! Update Cartesian form of total fluxes
        if (.not. useSplitViscosity) fluxes1 = fluxes1 - fluxes2

        ! Send viscous fluxes to appropriate patches
        call farfield_store_viscous_fluxes(fluxes2)
     end if

     ! Transform fluxes from Cartesian to contravariant form: `fluxes1` has the Cartesian
     ! form of total fluxes... upon return, `fluxes2` has the contravariant form
     call transform_fluxes(fluxes1, fluxes2)

     deallocate(fluxes1) !... no longer needed

     ! Take derivatives of fluxes
     do i = 1, nDimensions
        call first_derivative_apply(i, fluxes2(:,:,i))
     end do
     source = source - sum(fluxes2, dim = 3)

     deallocate(fluxes2) !... no longer needed

  case (ADJOINT)

     ! Partial derivatives of adjoint variables w.r.t. *computational* coordinates
     allocate(temp1(nGridPoints, nUnknowns, nDimensions))
     do i = 1, nDimensions
        temp1(:,:,i) = adjointVariables
        call adjoint_first_derivative_apply(i, temp1(:,:,i))
     end do

     allocate(localFluxJacobian1(nUnknowns, nUnknowns))
     allocate(localConservedVariables(nUnknowns))
     allocate(localVelocity(nDimensions))
     if (nSpecies .gt. 0) allocate(localMassFraction(nSpecies))
     allocate(localMetricsAlongDirection1(nDimensions))

     if (useViscosity) then
        allocate(localFluxJacobian2(nUnknowns, nUnknowns))
        allocate(localStressTensor(nDimensions ** 2))
        allocate(localHeatFlux(nDimensions))
        if (nSpecies .gt. 0) then
           allocate(localSpeciesFlux(nSpecies,nDimensions))
           if (equationOfState .eq. IDEAL_GAS_MIXTURE)                                       &
                allocate(localEnthalpyFlux(nDimensions))
        end if
     end if

     do j = 1, nGridPoints

        localConservedVariables = conservedVariables(j,:)
        localVelocity = velocity(j,:)
        if (nSpecies .gt. 0) localMassFraction = massFraction(j,:)
        if (useViscosity) then
           localStressTensor = stressTensor(j,:)
           localHeatFlux = heatFlux(j,:)
           if (nSpecies .gt. 0) then
              localSpeciesFlux = speciesFlux(j,:,:)
              if (equationOfState .eq. IDEAL_GAS_MIXTURE)                                    &
                   localEnthalpyFlux = enthalpyFlux(j,:)
           end if
        end if

        do i = 1, nDimensions

           localMetricsAlongDirection1 = metrics(j,1+nDimensions*(i-1):nDimensions*i)

           localFluxJacobian1 = 0.0_WP

           if (useViscosity) then
              call compute_first_partial_viscous_jacobian(localConservedVariables,           &
                   localMetricsAlongDirection1, localStressTensor, localHeatFlux,            &
                   localEnthalpyFlux, localSpeciesFlux, localFluxJacobian2,                  &
                   specificVolume(j,1), localVelocity, temperature(j,1), localMassFraction)
              localFluxJacobian1 = localFluxJacobian1 - localFluxJacobian2
           end if

           source(j,:) = source(j,:) + matmul(transpose(localFluxJacobian1), temp1(j,:,i))

        end do

     end do

     if (allocated(localConservedVariables)) deallocate(localConservedVariables)
     if (allocated(localFluxJacobian1)) deallocate(localFluxJacobian1)
     if (allocated(localHeatFlux)) deallocate(localHeatFlux)
     if (allocated(localEnthalpyFlux)) deallocate(localEnthalpyFlux)
     if (allocated(localStressTensor)) deallocate(localStressTensor)
     if (allocated(localMassFraction)) deallocate(localMassFraction)
     if (allocated(localSpeciesFlux)) deallocate(localSpeciesFlux)
     if (allocated(localFluxJacobian2)) deallocate(localFluxJacobian2)

     if (useViscosity) then

        allocate(temp2(nGridPoints, nUnknowns - 1)); temp2 = 0.0_WP

        allocate(localMetricsAlongDirection2(nDimensions))
        allocate(localFluxJacobian2(nUnknowns - 1, nUnknowns - 1))
        allocate(localAdjointDiffusion(nUnknowns - 1, nDimensions))

        do k = 1, nGridPoints

           localVelocity = velocity(k,:)
           localAdjointDiffusion = 0.0_WP

           do j = 1, nDimensions

              localMetricsAlongDirection2 = metrics(k,1+nDimensions*(j-1):nDimensions*j)

              do i = 1, nDimensions

                 localMetricsAlongDirection1 = metrics(k,1+nDimensions*(i-1):nDimensions*i)

                 call compute_second_partial_viscous_jacobian(localVelocity,                 &
                      dynamicViscosity(k,1), secondCoefficientOfViscosity(k,1),              &
                      thermalDiffusivity(k,1), massDiffusivity(k,:), temperature(k,1),       &
                      jacobian(k,1), localMetricsAlongDirection1,                            &
                      localMetricsAlongDirection2, localFluxJacobian2)
                 localAdjointDiffusion(:,j) = localAdjointDiffusion(:,j) +                   &
                      matmul(transpose(localFluxJacobian2), temp1(k,2:nUnknowns,i))

              end do

           end do

           do j = 1, nDimensions
              temp1(k,2:nUnknowns,j) = localAdjointDiffusion(:,j)
           end do

        end do

        do j = 1, nDimensions
           call adjoint_first_derivative_apply(j, temp1(:,2:nUnknowns,j))
        end do
        temp2 = sum(temp1(:,2:nUnknowns,:), dim = 3) !... divergence of the adjoint flux

        select case (equationOfState)

        case (IDEAL_GAS)

           temp2(:,nDimensions+1) = ratioOfSpecificHeats * specificVolume(:,1) *             &
                temp2(:,nDimensions+1)
           do i = 1, nDimensions
              temp2(:,i) = specificVolume(:,1) * temp2(:,i) - velocity(:,i) *                &
                   temp2(:,nDimensions+1)
           end do
           do k = 1, nSpecies
              temp2(:,nDimensions+1+k) = specificVolume(:,1) * temp2(:,nDimensions+1+k)
           end do

        case (IDEAL_GAS_MIXTURE)

           temp2(:,nDimensions+1) = ratioOfSpecificHeats * specificVolume(:,1) *             &
                mixtureMolecularWeight(:,1) * temp2(:,nDimensions+1)
           do i = 1, nDimensions
              temp2(:,i) = specificVolume(:,1) * temp2(:,i) - velocity(:,i) *                &
                   temp2(:,nDimensions+1)
           end do
           do k = 1, nSpecies
              temp2(:,nDimensions+1+k) = specificVolume(:,1) * temp2(:,nDimensions+1+k) +    &
                   temperature(:,1) * (molecularWeightInverse(nSpecies+1) -                  &
                   molecularWeightInverse(k)) * temp2(:,nDimensions+1) /                     &
                   ratioOfSpecificHeats
           end do

        end select

        source(:,2:nUnknowns) = source(:,2:nUnknowns) - temp2
        source(:,1) = source(:,1) + specificVolume(:,1) *                                    &
             conservedVariables(:,nDimensions+2) * temp2(:,nDimensions+1) +                  &
             sum(velocity * temp2(:,1:nDimensions), dim = 2)
        if (nSpecies .gt. 0) source(:,1) = source(:,1) +                                     &
             sum(massFraction * temp2(:,nDimensions+2:nUnknowns-1), dim = 2)

        deallocate(temp2)

        deallocate(localAdjointDiffusion)
        deallocate(localFluxJacobian2)
        deallocate(localMetricsAlongDirection2)

     end if

     deallocate(localMetricsAlongDirection1)
     deallocate(localVelocity)
     deallocate(temp1)

  end select

  ! Multiply by Jacobian
  do i = 1, nUnknowns
     source(:,i) = source(:,i) * jacobian(:,1)
  end do

  return
end subroutine add_viscous_flux



! -------------------------------------------
subroutine die(errorText)

  ! External modules
  use parallel

  implicit none

  ! Arguments
  character(len = *), intent(in) :: errorText

  call parallel_kill(errorText)

  return
end subroutine die
