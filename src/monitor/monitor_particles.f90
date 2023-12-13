module monitor_particles

  ! External modules
  use monitor
  use particle

  implicit none

end module monitor_particles


! ========================================== !
! Setup the routine to monitor the particles !
! ========================================== !
subroutine monitor_particles_setup(mode, controlIteration)

  ! Internal modules
  use monitor_particles

  ! external modules
  use simulation_flags
  use solver_options
  use geometry

  implicit none

  ! Arguments
  integer, intent(in) :: mode, controlIteration

  ! Local variables
  integer :: i
  character(len = 1) :: vel(3), xyz(3)

  if (.not. useParticles) return

  select case (mode)

  case (FORWARD)

     if (controlIteration .eq. 0) then

        ! ----------------------------------------
        ! Monitor the baseline particle prediction

        nParticlesPrevious = nParticlesGlobal
        
        ! Basic particle information
        call monitor_create('particles', 11)
        call monitor_set_header(1 , 'npart', 'i')
        call monitor_set_header(2 , 'npart-in', 'i')
        call monitor_set_header(3 , 'npart-out', 'i')
        call monitor_set_header(4 , 'ncol', 'i')
        call monitor_set_header(5 , 'min(Rep)', 'r')
        call monitor_set_header(6 , 'max(Rep)', 'r')
        call monitor_set_header(7 , 'min(Map)', 'r')
        call monitor_set_header(8 , 'max(Map)', 'r')
        call monitor_set_header(9 , 'min(Knp)', 'r')
        call monitor_set_header(10, 'max(Knp)', 'r')
        call monitor_set_header(11, 'minDragTime', 'r')

        ! Monitor the particle velocity
        call monitor_create('particle_vel', nDimensions * 4)

        ! Velocity
        vel = (/ 'U', 'V', 'W' /)
        do i = 1, nDimensions
           call monitor_set_header((i-1) * 4 + 1, 'min'//vel(i), 'r')
           call monitor_set_header((i-1) * 4 + 2, 'max'//vel(i), 'r')
           call monitor_set_header((i-1) * 4 + 3, 'mean'//vel(i), 'r')
           call monitor_set_header((i-1) * 4 + 4, 'var'//vel(i), 'r')
        end do

        ! Monitor the particle temperature
        call monitor_create('particle_heat', 3)
        call monitor_set_header(1, 'minT', 'r')
        call monitor_set_header(2, 'maxT', 'r')
        call monitor_set_header(3, 'meanT', 'r')


        ! Monitor the particle acceleration
        call monitor_create('particle_force', nDimensions * 6)
        xyz = (/ 'x', 'y', 'z' /)
        do i = 1, nDimensions
           call monitor_set_header((i-1) * 6 + 1, 'drag_'//xyz(i), 'r')
           call monitor_set_header((i-1) * 6 + 2, 'stress_'//xyz(i), 'r')
           call monitor_set_header((i-1) * 6 + 3, 'gravity_'//xyz(i), 'r')
           call monitor_set_header((i-1) * 6 + 4, 'collision_'//xyz(i), 'r')
           call monitor_set_header((i-1) * 6 + 5, 'lift_'//xyz(i), 'r')
           call monitor_set_header((i-1) * 6 + 6, 'added_mass_'//xyz(i), 'r')
        end do

        ! Reset values
        nParticleCollisions = 0
        do i = 1, nDimensions
           monitorParticle(i)%drag      = 0.0_WP
           monitorParticle(i)%stress    = 0.0_WP
           monitorParticle(1)%gravity   = 0.0_WP
           monitorParticle(i)%collision = 0.0_WP
           monitorParticle(i)%lift      = 0.0_WP
           monitorParticle(i)%addedMass = 0.0_WP
        end do

        ! Monitor pseudo-turbulent Reynolds stresses
        if (usePTKE) then
           select case (nDimensions)
           case (2)
              call monitor_create('ptke', 4)
              call monitor_set_header(1,'<k>','r')
              call monitor_set_header(2,'<R11>','r')
              call monitor_set_header(3,'<R22>','r')
              call monitor_set_header(4,'<R12>','r')
           case (3)
              call monitor_create('ptke', 7)
              call monitor_set_header(1,'<k>','r')
              call monitor_set_header(2,'<R11>','r')
              call monitor_set_header(3,'<R22>','r')
              call monitor_set_header(4,'<R33>','r')
              call monitor_set_header(5,'<R12>','r')
              call monitor_set_header(6,'<R13>','r')
              call monitor_set_header(7,'<R23>','r')
           end select
        end if

     else

        ! -------------------------------------
        ! Monitor the forward control iteration

     end if

  case (ADJOINT)

     if (controlIteration .eq. 0) then

        ! -------------------------------------
        ! Monitor the baseline adjoint solution

     else

        ! -------------------------------------
        ! Monitor the adjoint control iteration

     end if

  end select

  return
end subroutine monitor_particles_setup


! ================================================================ !
! Compute various statistics of the particles during a forward run !
! ================================================================ !
subroutine monitor_particles_forward(controlIteration)

  ! Internal modules
  use monitor_particles

  ! External modules
  use parallel
  use simulation_flags
  use solver_options
  use geometry
  use grid, only : gridNorm
  use state
  use time_integrator, only : nStages
  use particle_solver, only : minDragTime, minRep, maxRep, minMap, maxMap, minKnp, maxKnp

  implicit none

  ! Arguments
  integer, intent(in) :: controlIteration

  ! Local variables
  integer :: i, j, nParticlesOut
  real(WP) ::  minU(3), maxU(3), meanU(3), varU(3), minT, maxT, meanT, meanPTKE, meanR(6),   &
       volumeInverse, buf, normalization

  if (.not. useParticles) return

  if (controlIteration .eq. 0) then

     ! -------------------------------
     ! Monitor the baseline prediction

     ! Determine the total number of particles to leave since we last monitored
     nParticlesOut = nParticlesPrevious - nParticlesGlobal + nParticlesIn

     ! Set basic particle info
     call monitor_select('particles')
     call monitor_set_single_value(1 , real(nParticlesGlobal, WP))
     call monitor_set_single_value(2 , real(nParticlesIn, WP))
     call monitor_set_single_value(3 , real(nParticlesOut, WP))
     call monitor_set_single_value(4 , real(nParticleCollisions, WP))
     call monitor_set_single_value(5 , minRep)
     call monitor_set_single_value(6 , maxRep)
     call monitor_set_single_value(7 , minMap)
     call monitor_set_single_value(8 , maxMap)
     call monitor_set_single_value(9 , minKnp)
     call monitor_set_single_value(10 ,maxKnp)
     call monitor_set_single_value(11, minDragTime)

     ! Reset
     nParticlesPrevious = nParticlesGlobal

     ! Particle stats
     minU = huge(1.0_WP); maxU = -huge(1.0_WP); meanU = 0.0_WP; varU = 0.0_WP
     minT = huge(1.0_WP); maxT = -huge(1.0_WP); meanT = 0.0_WP
     do i = 1, nParticles
        do j = 1, nDimensions
           minU(j) = min(minU(j), particles(i)%velocity(j))
           maxU(j) = max(maxU(j), particles(i)%velocity(j))
           meanU(j) = meanU(j) + particles(i)%velocity(j)
        end do
        minT = min(minT, particles(i)%temperature)
        maxT = max(maxT, particles(i)%temperature)
        meanT = meanT + particles(i)%temperature
     end do
     do i = 1, nDimensions
        call parallel_min(minU(i))
        call parallel_max(maxU(i))
        call parallel_sum(meanU(i))
        meanU(i) = meanU(i) / real(nParticlesGlobal, WP)
     end do
     call parallel_min(minT)
     call parallel_max(maxT)
     call parallel_sum(meanT)
     meanT = meanT / real(nParticlesGlobal, WP)

     do i = 1, nParticles
        do j = 1, nDimensions
           varU(j) = varU(j) + (particles(i)%velocity(j) - meanU(j))**2
        end do
     end do
     do i = 1, nDimensions
        call parallel_sum(varU(i))
        varU(i) = varU(i) / real(nParticlesGlobal, WP)
     end do

     ! Set velocity values to monitor
     call monitor_select('particle_vel')
     do i = 1, nDimensions
        call monitor_set_single_value((i-1) * 4 + 1, minU(i))
        call monitor_set_single_value((i-1) * 4 + 2, maxU(i))
        call monitor_set_single_value((i-1) * 4 + 3, meanU(i))
        call monitor_set_single_value((i-1) * 4 + 4, varU(i))
     end do

     ! Set temperature values to monitor
     call monitor_select('particle_heat')
     call monitor_set_single_value(1, minT)
     call monitor_set_single_value(2, maxT)
     call monitor_set_single_value(3, meanT)

     ! Set particle acceleration to monitor
     normalization = 1.0_WP / real(nParticlesGlobal, WP) / real(reportInterval, WP) /        &
          real(nStages, WP)
     call monitor_select('particle_force')
     do i = 1, nDimensions
        ! Drag
        call parallel_sum(monitorParticle(i)%drag, buf); buf = buf * normalization
        call monitor_set_single_value((i-1) * 6 + 1, buf)

        ! Stress
        call parallel_sum(monitorParticle(i)%stress, buf); buf = buf * normalization
        call monitor_set_single_value((i-1) * 6 + 2, buf)

        ! Gravity
        call parallel_sum(monitorParticle(i)%gravity, buf); buf = buf * normalization
        call monitor_set_single_value((i-1) * 6 + 3, buf)

        ! Collision
        call parallel_sum(monitorParticle(i)%collision, buf); buf = buf * normalization
        call monitor_set_single_value((i-1) * 6 + 4, buf)

        ! Lift
        call parallel_sum(monitorParticle(i)%lift, buf); buf = buf * normalization
        call monitor_set_single_value((i-1) * 6 + 5, buf)

        ! Added mass
        call parallel_sum(monitorParticle(i)%addedMass, buf); buf = buf * normalization
        call monitor_set_single_value((i-1) * 6 + 6, buf)

        ! Reset values
        monitorParticle(i)%drag      = 0.0_WP
        monitorParticle(i)%stress    = 0.0_WP
        monitorParticle(i)%gravity   = 0.0_WP
        monitorParticle(i)%collision = 0.0_WP
        monitorParticle(i)%lift      = 0.0_WP
        monitorParticle(i)%addedMass = 0.0_WP
     end do

     ! Pesudo-turbulent Reynolds stresses
     if (usePTKE) then
        call monitor_select('ptke')
        meanPTKE = 0.0_WP
        meanR = 0.0_WP
        volumeInverse = 0.0_WP
        select case (nDimensions)
        case (2)
           do i = 1, nGridPoints
              buf = gridNorm(i,1) / conservedVariables(i,1)
              volumeInverse = volumeInverse + gridNorm(i,1)
              meanPTKE = meanPTKE + conservedVariables(i,nDimensions+3) * buf
              meanR(1) = meanR(1) + reynoldsStress(i,1) * buf
              meanR(2) = meanR(2) + reynoldsStress(i,4) * buf
              meanR(3) = meanR(3) + reynoldsStress(i,2) * buf
           end do
           call parallel_sum(volumeInverse); volumeInverse = 1.0_WP / volumeInverse
           call parallel_sum(meanPTKE); meanPTKE = meanPTKE * volumeInverse
           call parallel_sum(meanR(1)); meanR(1) = meanR(1) * volumeInverse
           call parallel_sum(meanR(2)); meanR(2) = meanR(2) * volumeInverse
           call parallel_sum(meanR(3)); meanR(3) = meanR(3) * volumeInverse
           call monitor_set_single_value(1, meanPTKE)
           call monitor_set_single_value(2, meanR(1))
           call monitor_set_single_value(3, meanR(2))
           call monitor_set_single_value(4, meanR(3))
        case (3)
           do i = 1, nGridPoints
              buf = gridNorm(i,1) / conservedVariables(i,1)
              volumeInverse = volumeInverse + gridNorm(i,1)
              meanPTKE = meanPTKE + conservedVariables(i,nDimensions+3) * buf
              meanR(1) = meanR(1) + reynoldsStress(i,1) * buf
              meanR(2) = meanR(2) + reynoldsStress(i,5) * buf
              meanR(3) = meanR(3) + reynoldsStress(i,9) * buf
              meanR(4) = meanR(4) + reynoldsStress(i,2) * buf
              meanR(5) = meanR(5) + reynoldsStress(i,3) * buf
              meanR(6) = meanR(6) + reynoldsStress(i,6) * buf
           end do
           call parallel_sum(volumeInverse); volumeInverse = 1.0_WP / volumeInverse
           call parallel_sum(meanPTKE); meanPTKE = meanPTKE * volumeInverse
           call parallel_sum(meanR(1)); meanR(1) = meanR(1) * volumeInverse
           call parallel_sum(meanR(2)); meanR(2) = meanR(2) * volumeInverse
           call parallel_sum(meanR(3)); meanR(3) = meanR(3) * volumeInverse
           call parallel_sum(meanR(4)); meanR(4) = meanR(4) * volumeInverse
           call parallel_sum(meanR(5)); meanR(5) = meanR(5) * volumeInverse
           call parallel_sum(meanR(6)); meanR(6) = meanR(6) * volumeInverse
           call monitor_set_single_value(1, meanPTKE)
           call monitor_set_single_value(2, meanR(1))
           call monitor_set_single_value(3, meanR(2))
           call monitor_set_single_value(4, meanR(3))
           call monitor_set_single_value(5, meanR(4))
           call monitor_set_single_value(6, meanR(5))
           call monitor_set_single_value(7, meanR(6))
        end select
     end if

  else

     ! --------------------------------------
     ! Monitor the forward control iteration


  end if

  return
end subroutine monitor_particles_forward
