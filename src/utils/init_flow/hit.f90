module hit

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

end module hit

subroutine hit_grid

  ! Internal modules
  use hit

  implicit none

  ! Local variables
  integer :: i, j, k

  ! Check for periodicity
  do i = 1, nDimensions
     if (periodicityType(i) .ne. PLANE) call die('HIT must be periodic in direction ' //     &
          num2str(i))
  end do

  ! Simplify
  nx = globalGridSize(1)
  ny = globalGridSize(2)
  nz = globalGridSize(3)

  ! Read in the grid size
  call parser_read('Lx', Lx)
  Ly = Lx
  call parser_read('Lz', Lz, 0.0_WP)

  ! Compute the grid spacing (assume periodic)
  dx = Lx / real(nx, WP)
  dy = Ly / real(ny, WP)
  dz = Lz / real(nz, WP)

  ! Generate the grid
  do k = iStart(3), iEnd(3)
     do j = iStart(2), iEnd(2)
        do i = iStart(1), iEnd(1)
           if (nx .gt. 1) coordinates(grid_index(i,j,k), 1) = (Lx - dx) *  real(i - 1, WP) / &
                real(nx - 1, WP)
           if (ny .gt. 1) coordinates(grid_index(i,j,k), 2) = (Ly - dy) *  real(j - 1, WP) / &
                real(ny - 1, WP)
           if (nz .gt. 1) coordinates(grid_index(i,j,k), 3) = (Lz - dz) *  real(k - 1, WP) / &
                real(nz - 1, WP)
        end do
     end do
  end do

  ! Setup the grid metrics
  call grid_metrics_setup

  return
end subroutine hit_grid


subroutine hit_data

  ! Internal modules
  use hit

  ! External modules
  use math
  use random
  use parallel
  use state, only : conservedVariables

  implicit none

  ! Libraries
#ifdef USE_FFTW
  include 'fftw3.f'

  ! Local variables
  integer :: i, j, k, dim
  real(WP) :: gamma, P0, rho0, Ut, le
  real(WP), dimension(:,:), allocatable :: velocity

  ! Spectrum computation
  integer  :: nk
  real(WP) :: psr, ps1, ps2
  real(WP) :: ke, dk, kc, kk, kx, ky, kz, kk2
  real(WP) :: energy_spec, spec_amp, eps, amp_disc, rand
  complex(WP), dimension(:,:,:), pointer :: ak, bk

  ! Fourier coefficients
  integer(KIND=8) :: plan_r2c, plan_c2r
  complex(WP) :: ii = (0.0_WP, 1.0_WP)
  complex(WP), dimension(:,:,:), pointer :: Uk, Vk, Wk
  complex(WP), dimension(:,:,:), pointer :: Cbuf
  real(WP), dimension(:,:,:), pointer :: Rbuf

  if (nProcs .gt. 1) call die('HIT init requires nProcs=1')

  ! Set the number of conserved variables and allocate the array
  nUnknowns = nDimensions + 2
  allocate(conservedVariables(nGridPoints, nUnknowns))
  allocate(velocity(nGridPoints, nDimensions))

  ! Initialize homogeneous isotropic turbulence using spectrum proposed by
  ! Passot, T., & Pouquet, A. (1987). Numerical simulation of compressible
  ! homogeneous flows in the turbulent regime. Journal of Fluid Mechanics.

  call parser_read('fluctuations', Ut)
  call parser_read('energetic scale', le, Lx/twoPi)

  ke = 2.0_WP * pi / le
  dk = 2.0_WP * pi / Lx
  kc = real(nx / 2, WP) * dk
  eps = ke / 1000000.0_WP
  spec_amp = 16.0_WP * sqrt(2.0_WP / pi) * Ut**2 / ke
  amp_disc = sqrt(dk)**3
  
  ! Compute spectrum
  nk = nx / 2 + 1
  allocate(ak(nk,ny,nz), bk(nk,ny,nz))
  do k=1,nz
     do j=1,ny
        do i=1,nk
           ! Random numbers
           call random_number(rand)
           psr = twoPi * (rand - 0.5_WP)
           call random_number(rand)
           ps1 = twoPi * (rand - 0.5_WP)
           call random_number(rand)
           ps2 = twoPi * (rand - 0.5_WP)
           ! Wavenumbers
           kx = real(i - 1, WP) * dk
           ky = real(j - 1, WP) * dk
           if (j.gt.nk) ky=-real(nx+1-j,WP) * dk
           kz = real(k - 1, WP) * dk
           if (k.gt.nk) kz=-real(nx+1-k, WP) * dk
           kk = sqrt(kx**2 + ky**2 + kz**2)
           ! Spectrum
           energy_spec=spec_amp*(kk/ke)**4*exp(-2.0_WP*(kk/ke)**2)
           ! Coeff
           ak(i,j,k) = 0.0_WP
           bk(i,j,k) = 0.0_WP
           if ((kk.gt.eps) .and. (kk.le.kc)) then
              ak(i,j,k)=amp_disc*sqrt(energy_spec/(2.0_WP*pi*kk**2))*exp(ii*ps1)*cos(psr)
              bk(i,j,k)=amp_disc*sqrt(energy_spec/(2.0_WP*pi*kk**2))*exp(ii*ps2)*sin(psr)
           end if
        end do
     end do
  end do
  
  ! Compute 3D velocity field
  allocate(Uk(nk,ny,nz))
  allocate(Vk(nk,ny,nz))
  allocate(Wk(nk,ny,nz))
  Uk=(0.0_WP,0.0_WP)
  Vk=(0.0_WP,0.0_WP)
  Wk=(0.0_WP,0.0_WP)
  
  do dim = 1, 3
     
     ! Compute the Fourier coefficients
     do k=1,nz
        do j=1,ny
           do i=1,nk
              
              ! Wavenumbers
              kx=real(i-1,WP)*dk
              ky=real(j-1,WP)*dk
              if (j.gt.nk) ky=-real(nx+1-j,WP)*dk
              kz=real(k-1,WP)*dk
              if (k.gt.nk) kz=-real(nx+1-k,WP)*dk
              kk =sqrt(kx**2+ky**2+kz**2)
              kk2=sqrt(kx**2+ky**2)
              
              ! Compute the Fourier coefficients
              if ((kk.gt.eps).and.(kk.le.kc)) then
                 if (dim.eq.1) then
                    if (kk2.lt.eps) then
                       Uk(i,j,k)=(ak(i,j,k)+bk(i,j,k))/sqrt(2.0_WP)
                    else
                       Uk(i,j,k)=(ak(i,j,k)*kk*ky+bk(i,j,k)*kx*kz)/(kk*kk2)
                    end if
                 end if
                 if (dim.eq.2) then
                    if (kk2.lt.eps) then
                       Vk(i,j,k)=(bk(i,j,k)-ak(i,j,k))/sqrt(2.0_WP)
                    else
                       Vk(i,j,k)=(bk(i,j,k)*ky*kz-ak(i,j,k)*kk*kx)/(kk*kk2)
                    end if
                 end if
                 if (dim.eq.3) Wk(i,j,k)=-bk(i,j,k)*kk2/kk
              end if
              
           end do
        end do
     end do
     
  end do
  
  ! Oddball
  do k=2,nz
     do j=nk+1,ny
        Uk(1,j,k)=conjg(Uk(1,ny+2-j,nz+2-k))
        Vk(1,j,k)=conjg(Vk(1,ny+2-j,nz+2-k))
        Wk(1,j,k)=conjg(Wk(1,ny+2-j,nz+2-k))
     end do
  end do
  do k=nk+1,nz
     Uk(1,1,k)=conjg(Uk(1,1,nz+2-k))
     Vk(1,1,k)=conjg(Vk(1,1,nz+2-k))
     Wk(1,1,k)=conjg(Wk(1,1,nz+2-k))
  end do
  
  ! Inverse Fourier transform
  allocate(Cbuf(nk,ny,nz))
  allocate(Rbuf(nx,ny,nz))
  call dfftw_plan_dft_c2r_3d(plan_c2r,nx,ny,nz,Cbuf,Rbuf,FFTW_ESTIMATE)
  call dfftw_plan_dft_r2c_3d(plan_r2c,nx,ny,nz,Rbuf,Cbuf,FFTW_ESTIMATE)
  
  ! Execute the plans
  do k=1,nz
     do j=1,ny
        do i=1,nk
           Cbuf(i,j,k) = Uk(i,j,k)
        end do
     end do
  end do
  call dfftw_execute(plan_c2r)
  do k=1,nz
     do j=1,ny
        do i=1,nx
           velocity(grid_index(i,j,k), 1) = Rbuf(i,j,k)
        end do
     end do
  end do
  do k=1,nz
     do j=1,ny
        do i=1,nk
           Cbuf(i,j,k) = Vk(i,j,k)
        end do
     end do
  end do
  call dfftw_execute(plan_c2r)
  do k=1,nz
     do j=1,ny
        do i=1,nx
           velocity(grid_index(i,j,k), 2) = Rbuf(i,j,k)
        end do
     end do
  end do
  if (nDimensions .gt. 2) then
     do k=1,nz
        do j=1,ny
           do i=1,nk
              Cbuf(i,j,k) = Wk(i,j,k)
           end do
        end do
     end do
     call dfftw_execute(plan_c2r)
     do k=1,nz
        do j=1,ny
           do i=1,nx
              velocity(grid_index(i,j,k), 3) = Rbuf(i,j,k)
           end do
        end do
     end do
  end if

  ! Ambient pressure
  call parser_read('ratio of specific heats', gamma, 1.4_WP)
  P0 = 1.0_WP / gamma

  ! Constant density
  rho0 = 1.0_WP

  ! Set the density
  conservedVariables(:, 1) = rho0

  ! Set the momentum
  do i = 1, nDimensions
     conservedVariables(:, i+1) = rho0 * velocity(:,i)
  end do

  ! Set the energy
  conservedVariables(:, nDimensions+2) = P0 / (gamma - 1.0_WP) + 0.5_WP * rho0 *             &
       sum(velocity(:,1:nDimensions)**2, dim = 2)

  ! Clean up
  deallocate(velocity, Uk, Vk, Wk, ak, bk, Cbuf, Rbuf)

#endif
  
  return
end subroutine hit_data


subroutine hit_ibm

  ! Internal modules
  use hit

  ! External modules
  use random
  use math
  use parallel
  use ibm

  implicit none

  ! Local variables
  integer :: i, j, k, npartx, nparty, npartz, ix, iy, iz, ierr
  real(WP) :: volume, volumeFraction, particleVolume, volumeFactor, sumParticleVolume, Lp,   &
       Lpy, Lpz, rand, r, distance(3), dMean, dStd, dMin, dMax, dShift, particleSpacing
  character(len = str_medium) :: filename, particleDistribution
  logical :: success, preventOverlap, setSeed
  integer :: mySeed(33)=123456789

  ! Particle type
  type :: t_Particle
     real(WP) :: diameter
     real(WP), dimension(3) :: position
     real(WP), dimension(3) :: velocity
  end type t_Particle
  type(t_Particle), dimension(:), allocatable :: particles
  integer :: nParticles

  ! Return if not writing an ibm file
  call parser_read('init ibm file', filename, '')
  if (len_trim(filename) .eq. 0) return

  ! Read particle diameter
  call parser_read('particle mean diameter', dMean)
  call parser_read('particle std diameter', dStd, 0.0_WP)
  call parser_read('particle min diameter', dMin, 0.0_WP)
  call parser_read('particle max diameter', dMax, 0.0_WP)
  call parser_read('particle diameter shift', dShift, 0.0_WP)

  ! Get the particle distribution type
  call parser_read('particle distribution', particleDistribution, 'random')
  call parser_read('prevent particle overlap', preventOverlap, .true.)

  ! Initialize the random number generator
  call parser_read('set seed', setSeed, .false.)
  if (.not. setSeed) then
     call random_setup
  else
     call random_seed(put=mySeed)
  end if

  ! Initial Volume that particles occupy
  if (nz .gt. 1) then
     volume = Lx * Ly * Lz
  else
     volume = Lx * Ly
  end if

  ! Get particle volume
  if (nDimensions .eq. 2) then
     volumeFactor = 1.0_WP / 4.0_WP
  else if (nDimensions .eq. 3) then
     volumeFactor = 1.0_WP / 6.0_WP
  end if
  particleVolume = pi * volumeFactor * dMean ** nDimensions

  select case (trim(particleDistribution))

  case ('random')

     ! Get the number of particles based on the volume fraction and distribute to processors
     call parser_read('particle volume fraction', volumeFraction)
     nParticles = int(volumeFraction * volume / particleVolume)

  case ('uniform')

     ! Mean interparticle distance
     call parser_read('particle volume fraction', volumeFraction)
     npartx = 1; nparty = 1; npartz = 1
     Lp = (particleVolume / volumeFraction)**(1.0_WP / real(nDimensions, WP))
     npartx = int(Lx / Lp)
     Lp = Lx / real(npartx, WP)
     if (ny .gt. 1) nparty = int(Ly / Lp)
     Lpy = Ly / real(nparty, WP)
     if (nz .gt. 1) npartz = int(Lz / Lp)
     Lpz = Lz / real(npartz, WP)
     nParticles = npartx * nparty * npartz

  case default

     call die("Unknown particle distribution '" // trim(particleDistribution) // "'")

  end select

  ! Allocate the particle vector
  allocate(particles(nParticles))

  ! Only root process generates IBM particle data
  if (iRank .eq. iRoot) then

     ! Initialize the particle parameters
     sumParticleVolume = 0.0_WP
     do i = 1, nParticles

        ! Set particle size distribution (compact support lognormal)
        particles(i)%diameter = random_lognormal(m=dMean,sd=dStd) + dShift
        if (dStd.gt.0.0_WP) then
           do while (particles(i)%diameter.gt.dMax+epsilon(1.0_WP) .or.                      &
                particles(i)%diameter.lt.dMin-epsilon(1.0_WP))
              particles(i)%diameter=random_lognormal(m=dMean,sd=dStd) + dShift
           end do
        else
           particles(i)%diameter = dMean
        end if

        ! Position (this will be reset later)
        particles(i)%position = 0.0_WP

        ! Particle velocity
        particles(i)%velocity = 0.0_WP

        ! Sum particle volume
        sumParticleVolume = sumParticleVolume + pi * volumeFactor *                          &
             particles(i)%diameter ** nDimensions
     end do

     ! Compute effective volume fraction
     volumeFraction = sumParticleVolume / volume

     ! Distribute the particles
     select case (trim(particleDistribution))

     case ('random')
        i=1
        do while (i.le.nParticles)

           ! Give the particle a random position
           if (nx .gt. 1) then
              call random_number(rand)
              particles(i)%position(1) = Lx * rand
           end if
           if (ny .gt. 1) then
              call random_number(rand)
              particles(i)%position(2) = Ly * rand
           end if
           if (nz .gt. 1) then
              call random_number(rand)
              particles(i)%position(3) = Lz * rand
           end if

           ! Prevent particle overlap
           success = .true.
           if (preventOverlap) then
              part2: do j = 1, i-1
                 ! Compute separation distance (account for periodicity)
                 distance = 0.0_WP
                 do k = 1, nDimensions
                    distance(k) = abs(particles(i)%position(k) - particles(j)%position(k))
                    if (isPeriodic(k)) then
                       distance(k) = min(distance(k), periodicLength(k) -                    &
                            abs(particles(i)%position(k) - particles(j)%position(k)))
                    end if
                 end do
                 r = sqrt(sum(distance**2))
                 particleSpacing = 0.5_WP*(particles(i)%diameter+particles(j)%diameter)      &
                      + 3.0_WP * minGridSpacing
                 if (setSeed) then
                    particleSpacing = 0.5375_WP*(particles(i)%diameter+particles(j)%diameter)
                 end if
                 if (r.le.particleSpacing) then
                    i=i-1
                    success = .false.
                    exit part2
                 end if
              end do part2
           end if
           i=i+1

           if (success .and. modulo(i,1000).eq.0)                                            &
                print *, real(i,SP)/real(nParticles,SP)*100.0_SP,'%'

        end do

     case ('uniform')
        do i = 1, nParticles
           ix = (i - 1) / (nparty * npartz)
           iy = (i - 1 - nparty * npartz * ix) / npartz
           iz = i - 1 - nparty * npartz * ix - npartz * iy
           particles(i)%position(1) = (ix + 0.5_WP) * Lp 
           particles(i)%position(2) = (iy + 0.5_WP) * Lpy
           particles(i)%position(3) = (iz + 0.5_WP) * Lpz
        end do

     end select

     ! Output stuff to the screen
     print *
     print *, 'Number of particles: ', nParticles
     print *, 'Volume fraction: ', volumeFraction
     print *
  end if

  ! Set IBM object data for i/o
  nObjects = nParticles
  allocate(object(nObjects))
  if (iRank .eq. iRoot) then
     do i = 1, nObjects
        select case (nDimensions)
        case (2)
           object(i)%volume = 0.25_WP * pi * particles(i)%diameter**2
        case (3)
           object(i)%volume = pi * particles(i)%diameter**3 / 6.0_WP
        end select
        object(i)%position = particles(i)%position
        object(i)%velocity = particles(i)%velocity
        object(i)%angularVelocity = 0.0_WP
     end do
  end if
  call MPI_BCAST(object, nObjects, MPI_OBJECT, iRoot, MPI_COMM_WORLD, ierr)

  return
end subroutine hit_ibm
