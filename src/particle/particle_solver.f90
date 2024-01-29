module particle_solver

  ! External modules
  use particle
  use particle_exchange

  implicit none

  ! Drag parameters
  integer :: dragModel
  integer, parameter ::                                                                      &
       NO_DRAG          = 0,                                                                 &
       STOKES           = 1,                                                                 &
       SCHILLER_NAUMANN = 2,                                                                 &
       GIDASPOW         = 3,                                                                 &
       TENNETI          = 4,                                                                 &
       BASSET           = 5,                                                                 &
       PARMAR           = 6,                                                                 &
       THEO             = 7,                                                                 &
       LOTH             = 8,                                                                 &
       HENDERSON        = 9,                                                                 &
       SINGH            = 10,                                                                &
       KC               = 11,                                                                &
       OSNES            = 12
  ! Collision parameters
  real(WP) :: collisionTime, coefficientOfRestitution, coefficientOfFriction
  logical :: particleWall(3,2), interParticleCollisions

  ! Phase change
  real(WP) :: diameterCutoff

  ! Monitor min/max particle quantities
  real(WP) :: minDragTime, minRep, maxRep, minMap, MaxMap, minKnp, MaxKnp

end module particle_solver


! ========================= !
! Setup the particle solver !
! ========================= !
subroutine particle_solver_setup

  ! Internal modules
  use particle_solver

  ! External modules
  use parser
  use parallel
  use geometry
  use string
  use thermochem
  use solver_options

  implicit none

  ! Local variables
  integer :: i
  character(len = str_medium) :: val
  real(WP) :: Cpg, Cpd
  logical :: wall

  if (.not. useParticles) return

  ! Particles are not yet implemented for curvilinear grids
  if (isDomainCurvilinear)                                                                   &
       call die('Particle solver currently implemented for Cartesian grids only!')

  ! Particles require a viscous flow
  if (.not. useViscosity) call die('Particle solver requires a viscous flow!')

  ! Pseudo turbulence safety check
  if (usePTKE) then
     if (.not.twoWayCoupling) call die('Reynolds stress requires two-way couling!')
     if (nSpecies.lt.1) call die('PTKE requires at least one scalar!')
  end if

  ! Get the particle properties
  call parser_read('particle density', particleDensity)
  call parser_read('particle specific heat', particleSpecificHeat, 1.0_WP)

  ! Particle filter (Capecelatro & Desjardins, 2013)
  call parser_read('particle filter size', filterSize, -1.0_WP)
  if (twoWayCoupling .and. filterSize .gt. 0.0_WP) then
     diffusionAmount = max(filterSize**2-minGridSpacing**2, 0.0_WP) / (16.0_WP * log(2.0_WP))
     !diffusionAmount = filterSize**2 / (16.0_WP * log(2.0_WP))
     
     ! Decide if fluid variables should be filtered prior to interpolation
     call parser_read('filter fluid', filterFluid, .false.)
  else
     diffusionAmount = -1.0_WP
     filterFluid = .false.
  end if

  ! Initialize min/max quantities
  minDragTime = huge(1.0_WP)
  minRep = 0.0_WP; maxRep = 0.0_WP
  minMap = 0.0_WP; maxMap = 0.0_WP
  minKnp = 0.0_WP; maxKnp = 0.0_WP

  ! Include granular temperature
  call parser_read('use granular temperature', useGranularTemperature, .false.)

  ! Get the drag model
  call parser_read('drag model', val, 'stokes')
  select case (trim(val))

  case ('none', 'NONE', 'None')

     dragModel = NO_DRAG

  case ('stokes', 'Stokes')

     dragModel = STOKES

  case ('schiller naumann')

     dragModel = SCHILLER_NAUMANN

  case ('gidaspow', 'Gidaspow')

     dragModel = GIDASPOW

  case ('tenneti', 'Tenneti')

     dragModel = TENNETI

  case ('henderson', 'Henderson')

     dragModel = HENDERSON

  case ('loth', 'Loth')

     dragModel = LOTH

  case ('osnes', 'Osnes')

     dragModel = OSNES
     
     
  case ('basset', 'Basset')

     dragModel = BASSET

  case ('parmar', 'Parmar')

     dragModel = PARMAR
     
  case ('singh', 'SINGH')

     dragModel = SINGH

  case ('kc', 'KC')

     dragModel = KC

  case default

     call die("particle_setup: unknown drag model: '" // trim(val) // "'!")

  end select

  ! Lift model
  call parser_read('use saffman lift', useSaffmanLift, .false.)

  ! Added mass
  call parser_read('use added mass', useAddedMass, .false.)

  ! Interphase heat transfer
  call parser_read('use interphase heat transfer', useParticleHeat, .true.)

  ! Collisions
  call parser_read('particle collisions', collisionsOn, .false.)
  interParticleCollisions = .false.
  if (collisionsOn) then
     ! Get the collision parameters
     call parser_read('inter-particle collisions', interParticleCollisions, .true.)
     call parser_read('collision time', collisionTime)
     call parser_read('coefficient of restitution', coefficientOfRestitution, 0.85_WP)
     call parser_read('coefficient of friction', coefficientOfFriction, 0.0_WP)
     if (coefficientOfFriction .gt. 0.0_WP) then
        useFriction = .true.
     else
        useFriction = .false.
     end if

     ! Check for walls (by default add walls in non-periodic directions)
     particleWall = .false.
     do i = 1, nDimensions
        if (.not. isPeriodic(i)) then
           if (procCoords(i).eq.0) particleWall(i,1) = .true.
           if (procCoords(i).eq.nProcsDir(i)-1) particleWall(i,2) = .true.
        end if
     end do

     ! Option to remove walls from non-periodic directions 
     if (particleWall(1,1) .or. particleWall(1,2)) then
        call parser_read('particle wall in x', wall, .true.)
        if (.not.wall) particleWall(1,:) = .false.
     end if
     if (particleWall(2,1) .or. particleWall(2,2)) then
        call parser_read('particle wall in y', wall, .true.)
        if (.not.wall) particleWall(2,:) = .false.
     end if
     if (particleWall(3,1) .or. particleWall(3,2)) then
        call parser_read('particle wall in z', wall, .true.)
        if (.not.wall) particleWall(3,:) = .false.
     end if
  else
     useFriction = .false.
  end if

  ! Phase change
  call parser_read('use phase change', usePhaseChange, .false.)
  if (usePhaseChange) then
     ! Make sure we have at least one scalar
     if (nSpecies.le.0) call die('particle_solver_setup: phase change requires nSpecies > 0!')

     ! Get the scalar index corresponding to vapor mass fraction
     call get_species_index('water', vaporIndex)

     ! Overwrite particle specific heat
     call get_Cp('air', 300.0_WP, Cpg)
     call get_Cp('water', 300.0_WP, Cpd)
     particleSpecificHeat = Cpd / Cpg

     ! Determine the minimum diameter before removing the particle
     call parser_read('diameter cutoff', diameterCutoff, 0.0_WP)
  end if

  return
end subroutine particle_solver_setup


! =========================== !
! Cleanup the particle solver !
! =========================== !
subroutine particle_solver_cleanup

  ! Internal modules
  use particle_solver

  implicit none

  collisionsOn = .false.

  return
end subroutine particle_solver_cleanup


! =========================== !
! Compute particle collisions !
! =========================== !
subroutine compute_particle_collisions

  ! Internal modules
  use particle_solver

  ! External modules
  use parallel
  use math
  use geometry
  use time_info, only : timestepSize

  implicit none

  ! Local variables
  integer :: i, j, k, ii, jj, kk, i1, i2, j1, j2, k1, k2, ip, jp, n, ijk(3), gridIndex
  integer(WP) :: id1, id2
  real(WP), parameter :: oneSixth = 1.0_WP / 6.0_WP
  real(WP), parameter :: clipCol = 0.2_WP
  real(WP) :: d1, d2, mass1, mass2, particleSeparation, radiusOfInfluence, delta, rnv,       &
       effectiveMass, springForce, damper, rtv, buf
  real(WP), dimension(3) :: pos1, pos2, vel1, vel2, n12, v12, v12n, normalCollision,         &
       omega1, omega2, t12, tangentialCollision

  ! Reset particle collision counter
  nParticleCollisions = 0

  if (.not. collisionsOn) return

  ! Set collision time to be 20x the simulation timestep
  collisionTime = 20.0_WP * timestepSize

  ! Communicate ghost particles
  if (interParticleCollisions) then
     call communicate_ghost_particles

     ! Initialize arrays for particle localization
     nPartInCell = 0
     partInCell = 0

     do i = 1, nParticles

        ! Localize particles in cells
        ii = particles(i)%gridIndex(1)
        jj = particles(i)%gridIndex(2)
        kk = particles(i)%gridIndex(3)
        nPartInCell(ii,jj,kk) = nPartInCell(ii,jj,kk) + 1
        partInCell(ii,jj,kk,nPartInCell(ii,jj,kk)) = i

        ! Reset the collision force and torque
        particles(i)%collision = 0.0_WP
        particles(i)%torque = 0.0_WP
     end do

     ! Map ghost particle to cell
     do i = 1, nGhostParticles

        ! Localize particles in cells
        ii = ghostParticles(i)%gridIndex(1)
        jj = ghostParticles(i)%gridIndex(2)
        kk = ghostParticles(i)%gridIndex(3)
        nPartInCell(ii,jj,kk) = nPartInCell(ii,jj,kk) + 1
        partInCell(ii,jj,kk,nPartInCell(ii,jj,kk)) = -i

     end do

     ! Initialize index extents
     i1 = 1; i2 = 1
     j1 = 1; j2 = 1
     k1 = 1; k2 = 1
  else
     do i = 1, nParticles
        ! Reset the collision force and torque
        particles(i)%collision = 0.0_WP
        particles(i)%torque = 0.0_WP
     end do
  end if

  ! Implement the collision term
  do ip = 1, nParticles

     ! Cycle if negative id
     if (particles(ip)%id.le.0) cycle

     ! Localize particles in cells
     i = particles(ip)%gridIndex(1)
     j = particles(ip)%gridIndex(2)
     k = particles(ip)%gridIndex(3)

     ! Extract `particle 1` data
     id1    = particles(ip)%id
     pos1   = particles(ip)%position
     vel1   = particles(ip)%velocity
     omega1 = particles(ip)%angularVelocity
     d1     = particles(ip)%diameter
     mass1  = oneSixth * particleDensity * pi * d1**3

     ! Particle-wall collisions
     ! ------------------------

     ! Closest wall parameters
     vel2   = 0.0_WP
     omega2 = 0.0_WP
     d2     = 0.0_WP

     ! Particle-wall radius of influence
     radiusOfInfluence = 0.0_WP!0.2_WP * d1

     ! Loop through each wall
     do n = 1, nDimensions

        ! Check collisions with left wall
        if (particleWall(n,1)) then
           ijk(1) = i; ijk(2) = j; ijk(3) = k
           ijk(n) = iStart(n)
           gridIndex = ijk(1) - gridOffset(1) + localGridSize(1) *                           &
                (ijk(2) - 1 - gridOffset(2) + localGridSize(2) *                             &
                (ijk(3) - 1 - gridOffset(3)))
           pos2 = pos1
           pos2(n) = coordinates(gridIndex, n)

           ! Distance between particle and wall
           particleSeparation = sqrt(sum((pos1 - pos2) * (pos1 - pos2)))

           ! Distance of influence
           delta = 0.5_WP * d1 + radiusOfInfluence - particleSeparation
           delta = min(delta, clipCol * 0.5_WP * d1)

           if (delta .gt. 0.0_WP) then
              ! Recompute distance of influence
              n12 = (pos2 - pos1) / particleSeparation
              v12 = vel1 - vel2
              rnv = sum(v12 * n12)
              radiusOfInfluence = min(2.0_WP * abs(rnv) * timeStepSize, radiusOfInfluence)
              delta = 0.5_WP * d1 + radiusOfInfluence - particleSeparation
              delta = min(delta, clipCol * 0.5_WP * d1)
              if (delta .gt. 0.0_WP) then
                 ! Normal collision
                 v12n = rnv * n12
                 effectiveMass = mass1
                 springForce = effectiveMass / collisionTime**2 * (pi**2 +                   &
                      log(coefficientOfRestitution)**2)
                 damper = -2.0_WP * log(coefficientOfRestitution) *                          &
                      effectiveMass / collisionTime
                 normalCollision  = -springForce * delta * n12 - damper * v12n
                 ! Tangential collision
                 t12 = v12 - v12n + cross_product(0.5_WP * (d1 * omega1 + d2 * omega2), n12)
                 rtv = sqrt(sum(t12 * t12))
                 tangentialCollision = 0.0_WP
                 if (useFriction .and. rtv.gt.0.0_WP) tangentialCollision =                  &
                      -coefficientOfFriction * sqrt(sum(normalCollision * normalCollision)) *&
                      t12 / rtv
                 ! Calculate acceleration due to collisions
                 normalCollision = normalCollision / mass1
                 tangentialCollision = tangentialCollision / mass1
                 particles(ip)%collision(1:nDimensions) =                                    &
                      particles(ip)%collision(1:nDimensions) +                               &
                      normalCollision(1:nDimensions) + tangentialCollision(1:nDimensions)
                 ! Calculate collision torque
                 particles(ip)%torque(1) = particles(ip)%torque(1) + 0.5_WP * (d1 * n12(2) * &
                      tangentialCollision(3) - d1 * n12(3) * tangentialCollision(2))
                 particles(ip)%torque(2) = particles(ip)%torque(2) + 0.5_WP * (d1 * n12(3) * &
                      tangentialCollision(1) - d1 * n12(1) * tangentialCollision(3))
                 particles(ip)%torque(3) = particles(ip)%torque(3) + 0.5_WP * (d1 * n12(1) * &
                      tangentialCollision(2) - d1 * n12(2) * tangentialCollision(1))
              end if
           end if
        end if

        ! Check collisions with right wall
        if (particleWall(n,2)) then
           ijk(1) = i; ijk(2) = j; ijk(3) = k
           ijk(n) = iEnd(n)
           gridIndex = ijk(1) - gridOffset(1) + localGridSize(1) *                           &
                (ijk(2) - 1 - gridOffset(2) + localGridSize(2) *                             &
                (ijk(3) - 1 - gridOffset(3)))
           pos2 = pos1
           pos2(n) = coordinates(gridIndex, n)

           ! Distance between particle and wall
           particleSeparation = sqrt(sum((pos1 - pos2) * (pos1 - pos2)))

           ! Distance of influence
           delta = 0.5_WP * d1 + radiusOfInfluence - particleSeparation
           delta = min(delta, clipCol * 0.5_WP * d1)

           if (delta .gt. 0.0_WP) then
              ! Recompute distance of influence
              n12 = (pos2 - pos1) / particleSeparation
              v12 = vel1 - vel2
              rnv = sum(v12 * n12)
              radiusOfInfluence = min(2.0_WP * abs(rnv) * timestepSize, radiusOfInfluence)
              delta = 0.5_WP * d1 + radiusOfInfluence - particleSeparation
              delta = min(delta, clipCol * 0.5_WP * d1)
              if (delta .gt. 0.0_WP) then
                 ! Normal collision
                 v12n = rnv * n12
                 effectiveMass = mass1
                 springForce = effectiveMass / collisionTime**2 * (pi**2 +                   &
                      log(coefficientOfRestitution)**2)
                 damper = -2.0_WP * log(coefficientOfRestitution) *                          &
                      effectiveMass / collisionTime
                 normalCollision  = -springForce * delta * n12 - damper * v12n
                 ! Tangential collision
                 t12 = v12 - v12n + cross_product(0.5_WP * d1 * omega1, n12)
                 rtv = sqrt(sum(t12 * t12))
                 tangentialCollision = 0.0_WP
                 if (useFriction .and. rtv.gt.0.0_WP) tangentialCollision =                  &
                      -coefficientOfFriction * sqrt(sum(normalCollision * normalCollision)) *&
                      t12 / rtv
                 ! Calculate acceleration due to collisions
                 normalCollision = normalCollision / mass1
                 tangentialCollision = tangentialCollision / mass1
                 particles(ip)%collision(1:nDimensions) =                                    &
                      particles(ip)%collision(1:nDimensions) +                               &
                      normalCollision(1:nDimensions) + tangentialCollision(1:nDimensions)
                 ! Calculate collision torque
                 particles(ip)%torque(1) = particles(ip)%torque(1) + 0.5_WP * (d1 * n12(2) * &
                      tangentialCollision(3) - d1 * n12(3) * tangentialCollision(2))
                 particles(ip)%torque(2) = particles(ip)%torque(2) + 0.5_WP * (d1 * n12(3) * &
                      tangentialCollision(1) - d1 * n12(1) * tangentialCollision(3))
                 particles(ip)%torque(3) = particles(ip)%torque(3) + 0.5_WP * (d1 * n12(1) * &
                      tangentialCollision(2) - d1 * n12(2) * tangentialCollision(1))
              end if
           end if
        end if

     end do

     ! Particle-IBM collisions
     ! -----------------------
     if (useIBM) then
        ! IBM parameters
        vel2   = 0.0_WP
        omega2 = 0.0_WP
        d2     = 0.0_WP

        ! Particle-IBM radius of influence
        radiusOfInfluence = 0.0_WP!0.2_WP * d1

        ! Get distance to IB surface
        call interpolate_fluid_to_particle(particles(ip)%gridIndex, particles(ip)%position,  &
             levelset = particleSeparation)

        ! Distance of influence
        delta = 0.5_WP * d1 + radiusOfInfluence - particleSeparation
        delta = min(delta, clipCol * 0.5_WP * d1)

        if (delta .gt. 0.0_WP) then
           ! Normal vector
           call interpolate_fluid_to_particle(particles(ip)%gridIndex,particles(ip)%position,&
                levelsetNormal = n12)
           buf = sqrt(sum(n12*n12)) + epsilon(1.0_WP)
           n12 = -n12/buf

           ! Recompute distance of influence
           v12 = vel1 - vel2
           rnv = sum(v12 * n12)
           radiusOfInfluence = min(2.0_WP * abs(rnv) * timestepSize, radiusOfInfluence)
           delta = 0.5_WP * d1 + radiusOfInfluence - particleSeparation
           delta = min(delta, clipCol * 0.5_WP * d1)
           if (delta .gt. 0.0_WP) then
              ! Normal collision
              v12n = rnv * n12
              effectiveMass = mass1
              springForce = effectiveMass / collisionTime**2 * (pi**2 +                      &
                   log(coefficientOfRestitution)**2)
              damper = -2.0_WP * log(coefficientOfRestitution) *                             &
                   effectiveMass / collisionTime
              normalCollision  = -springForce * delta * n12 - damper * v12n
              ! Tangential collision
              t12 = v12 - v12n + cross_product(0.5_WP * d1 * omega1, n12)
              rtv = sqrt(sum(t12 * t12))
              tangentialCollision = 0.0_WP
              if (useFriction .and. rtv.gt.0.0_WP) tangentialCollision =                     &
                   -coefficientOfFriction * sqrt(sum(normalCollision * normalCollision)) *   &
                   t12 / rtv
              ! Calculate acceleration due to collisions
              normalCollision = normalCollision / mass1
              tangentialCollision = tangentialCollision / mass1
              particles(ip)%collision(1:nDimensions) =                                       &
                   particles(ip)%collision(1:nDimensions) +                                  &
                   normalCollision(1:nDimensions) + tangentialCollision(1:nDimensions)
              ! Calculate collision torque
              particles(ip)%torque(1) = particles(ip)%torque(1) + 0.5_WP * (d1 * n12(2) *    &
                   tangentialCollision(3) - d1 * n12(3) * tangentialCollision(2))
              particles(ip)%torque(2) = particles(ip)%torque(2) + 0.5_WP * (d1 * n12(3) *    &
                   tangentialCollision(1) - d1 * n12(1) * tangentialCollision(3))
              particles(ip)%torque(3) = particles(ip)%torque(3) + 0.5_WP * (d1 * n12(1) *    &
                   tangentialCollision(2) - d1 * n12(2) * tangentialCollision(1))
               ! Update collision counter
              nParticleCollisions = nParticleCollisions + 1
           end if
        end if
     end if
     
     ! Particle-particle collisions
     ! ----------------------------
     if (interParticleCollisions) then

        ! Define bounds to loop over neighboring cells
        select case (nDimensions)
        case (1)
           i1 = i - 1; i2 = i + 1
        case (2)
           i1 = i - 1; i2 = i + 1
           j1 = j - 1; j2 = j + 1
        case (3)
           i1 = i - 1; i2 = i + 1
           j1 = j - 1; j2 = j + 1
           k1 = k - 1; k2 = k + 1
        end select

        ! Loop over neighboring cells
        do kk = k1, k2
           do jj = j1, j2
              do ii = i1, i2

                 ! Loop over particles in cells
                 do n = 1, nPartInCell(ii,jj,kk)

                    ! Get particle index
                    jp = partInCell(ii,jj,kk,n)

                    if (jp .gt. 0) then
                       ! Extract `particle 2` data
                       id2    = particles(jp)%id
                       pos2   = particles(jp)%position
                       vel2   = particles(jp)%velocity
                       omega2 = particles(jp)%angularVelocity
                       d2     = particles(jp)%diameter
                       mass2  = oneSixth * particleDensity * pi * d2**3
                    else
                       ! Get data from ghost particles
                       jp     = -jp
                       id2    = ghostParticles(jp)%id
                       pos2   = ghostParticles(jp)%position
                       vel2   = ghostParticles(jp)%velocity
                       omega2 = ghostParticles(jp)%angularVelocity
                       d2     = ghostParticles(jp)%diameter
                       mass2  = oneSixth * particleDensity * pi * d2**3
                    end if

                    ! Distance between particles `1` and `2`
                    particleSeparation = sqrt(sum((pos1 - pos2) * (pos1 - pos2)))

                    ! Guess particle-particle radius of influence
                    radiusOfInfluence = 0.0_WP!0.1_WP * (d1 + d2)

                    ! Distance of influence
                    delta = 0.5_WP * (d1 + d2) + radiusOfInfluence - particleSeparation
                    delta = min(delta, clipCol * 0.5_WP * (d1 + d2))

                    if (delta .gt. 0.0_WP .and. id1 .ne. id2) then
                       ! Recompute distance of influence
                       n12 = (pos2 - pos1) / particleSeparation
                       v12 = vel1 - vel2
                       rnv = sum(v12 * n12)
                       radiusOfInfluence = min(abs(rnv) * timestepSize, radiusOfInfluence)
                       delta = 0.5_WP * (d1 + d2) + radiusOfInfluence - particleSeparation
                       delta = min(delta, clipCol * 0.5_WP * (d1 + d2))
                       if (delta .gt. 0.0_WP) then
                          ! Normal collision
                          v12n = rnv * n12
                          effectiveMass = mass1 * mass2 / (mass1 + mass2)
                          springForce = effectiveMass / collisionTime**2 * (pi**2 +          &
                               log(coefficientOfRestitution)**2)
                          damper = -2.0_WP * log(coefficientOfRestitution) *                 &
                               effectiveMass / collisionTime
                          normalCollision  = -springForce * delta * n12 - damper * v12n
                          ! Tangential collision
                          t12 = v12 - v12n + cross_product(0.5_WP * (d1*omega1 + d2*omega2), &
                               n12)
                          rtv = sqrt(sum(t12 * t12))
                          tangentialCollision = 0.0_WP
                          if (useFriction .and. rtv.gt.0.0_WP) tangentialCollision =         &
                               -coefficientOfFriction *                                      &
                               sqrt(sum(normalCollision * normalCollision)) * t12 / rtv
                          ! Calculate acceleration due to collisions
                          normalCollision = normalCollision / mass1
                          tangentialCollision = tangentialCollision / mass1
                          particles(ip)%collision(1:nDimensions) =                           &
                               particles(ip)%collision(1:nDimensions) +                      &
                               normalCollision(1:nDimensions) +                              &
                               tangentialCollision(1:nDimensions)
                          ! Calculate collision torque
                          particles(ip)%torque(1) = particles(ip)%torque(1) + 0.5_WP *       &
                               (d1 * n12(2) * tangentialCollision(3) -                       &
                               d1 * n12(3) * tangentialCollision(2))
                          particles(ip)%torque(2) = particles(ip)%torque(2) + 0.5_WP *       &
                               (d1 * n12(3) * tangentialCollision(1) -                       &
                               d1 * n12(1) * tangentialCollision(3))
                          particles(ip)%torque(3) = particles(ip)%torque(3) + 0.5_WP *       &
                               (d1 * n12(1) * tangentialCollision(2) -                       &
                               d1 * n12(2) * tangentialCollision(1))
                          ! Update collision counter
                          nParticleCollisions = nParticleCollisions + 1
                       end if
                    end if

                 end do

              end do
           end do
        end do
     end if

  end do

  ! Add up particle collision counter
  call parallel_sum(nParticleCollisions)

  return
end subroutine compute_particle_collisions


! ================================ !
! Compute particle right-hand side !
! ================================ !
subroutine particle_rhs(part, dxdt, dudt, dwdt, dTdt, dmdt)

  ! Internal modules
  use particle_solver

  ! External modules
  use math, only : pi
  use thermochem
  use gravity_source
  
  implicit none

  ! Arguments
  type(t_Particle), intent(inout) :: part
  real(WP), dimension(3), intent(out) :: dxdt, dudt, dwdt
  real(WP), intent(out) :: dTdt, dmdt

  ! Acceleration terms
  real(WP), dimension(3) :: dragForce, gravityForce, collisionForce, liftForce, addedMass

  ! Interphase exchange
  real(WP) :: heatExchange, massExchange
  real(WP), dimension(3) :: momentumExchange

  ! Fluid quantities
  real(WP) :: density, viscosity, temperature, velocity(3), stress(3), vorticity(3), Yv,     &
       pressure, gradRho(3), divRhoU, soundSpeed

  ! Non-dimensional numbers
  real(WP) :: Rep, Pr, Nu, Sc, Sh, vf, Ma, Kn

  ! Lift and drag parameters
  real(WP) :: responseTimeInverse, relativeVelocity, dragCorrection, omega, Cl

  ! Rotation
  real(WP) :: momentOfInertia
  real(WP), dimension(3) :: hTorque

  ! Phase change
  real(WP) :: mass, Pref, Td, Cpd, Cpg, Cpv, Cpmix, Lv, Yvs

  ! Added mass
  real(WP) :: Cm, eta, rhs(3)

  if (part%stop .eq. 1) return

  ! Zero out the forces
  stress = 0.0_WP
  dragForce = 0.0_WP
  gravityForce = 0.0_WP
  collisionForce = 0.0_WP
  liftForce = 0.0_WP
  addedMass = 0.0_WP

  ! Interpolate fluid quantities to particle location
  call interpolate_fluid_to_particle(part%gridIndex, part%position, density, temperature,    &
       viscosity, velocity, stress, vf)
  if (useSaffmanLift .or. useFriction) call interpolate_fluid_to_particle(part%gridIndex,    &
       part%position, vorticity = vorticity)
  if (useAddedMass) call interpolate_fluid_to_particle(part%gridIndex, part%position,        &
       gradRho = gradRho, divRhoU = divRhoU)
  if (usePhaseChange) then
     call interpolate_fluid_to_particle(part%gridIndex, part%position,                       &
          massFraction = Yv, pressure = pressure)
     Yv = max(Yv, 0.0_WP)
  end if

  ! Relative velocity
  relativeVelocity = sqrt(sum((part%velocity(1:nDimensions) - velocity(1:nDimensions))**2))

  ! Non-dimensional numbers
  Rep = density * part%diameter * relativeVelocity / viscosity
  Rep = Rep + epsilon(1.0_WP)
  Pr = 1.0_WP / prandtlNumberInverse
  soundSpeed = sqrt((ratioOfSpecificHeats-1.0_WP) * temperature)
  Ma = relativeVelocity / soundSpeed
  Kn = sqrt(0.5_WP * ratioOfSpecificHeats * pi) * Ma / Rep

  ! Particle response time
  responseTimeInverse = 18.0_WP * viscosity / (particleDensity * part%diameter**2)

  ! Drag contribution
  call compute_drag_correction(part, density, temperature, viscosity, velocity, vf,          &
       dragCorrection)
  dragForce(1:nDimensions) = dragCorrection * responseTimeInverse *                          &
       (velocity(1:nDimensions) - part%velocity(1:nDimensions))

  ! Store minimum drag time for stability (based on max eigenvalue)
  if (dragCorrection.gt.0.0_WP .and. vf.lt.1.0_WP)  minDragTime = min(minDragTime,           &
       density * part%diameter**2 / (18.0_WP * viscosity * dragCorrection * (1.0_WP - vf)))
  minRep = min(minRep, Rep)
  maxRep = max(maxRep, Rep)
  minMap = min(minMap, Ma)
  maxMap = max(maxMap, Ma)
  minKnp = min(minKnp, Kn)
  maxKnp = max(maxKnp, Kn)

  ! Gravity contribution
  if (allocated(gravity)) gravityForce(1:nDimensions) = gravity(1:nDimensions)

  ! Collision contribution
  if (collisionsOn) collisionForce(1:nDimensions) = part%collision(1:nDimensions)

  ! Angular velocity
  if (useFriction) then
     hTorque = 6.0_WP * viscosity * (0.5_WP * vorticity - part%angularVelocity) /            &
          particleDensity 
     momentOfInertia = 0.1_WP * part%diameter**2
     dwdt = (part%torque + hTorque) / momentOfInertia
  else
     hTorque = 0.0_WP
     dwdt = 0.0_WP
  end if

  ! Compute the lift force
  if (useSaffmanLift) then

     ! Vorticity magnitude
     omega = sqrt(sum(vorticity)**2)

     ! Saffman lift (1968)
     Cl = 9.69_WP / pi / part%diameter / particleDensity * sqrt(viscosity * density) *       &
          (omega+epsilon(1.0_WP))**(-0.5_WP)
     liftForce(1) = Cl * ((velocity(2) - part%velocity(2)) * vorticity(3) -                  &
          (velocity(3) - part%velocity(3)) * vorticity(2))
     liftForce(2) = Cl * ((velocity(3) - part%velocity(3)) * vorticity(1) -                  &
          (velocity(1) - part%velocity(1)) * vorticity(3))
     liftForce(3) = Cl * ((velocity(1) - part%velocity(1)) * vorticity(2) -                  &
          (velocity(2) - part%velocity(2)) * vorticity(1))

!!$        ! Corrected for high Re (Mei 1999)
!!$        beta = 0.5_WP * part%diameter * omega /                                              &
!!$             sqrt(sum((part%velocity(1:nDimensions) - velocity(1:nDimensions))**2))
!!$        b1 = (1.0_WP - 0.3314_WP * sqrt(beta)) * exp(-0.1_WP*Rep)
!!$        b2 = 0.3314_WP * sqrt(beta)
!!$        if (Rep.le.40.0_WP) then
!!$           liftForce(1) = liftForce(1) * b1 + b2
!!$           liftForce(2) = liftForce(2) * b1 + b2
!!$           liftForce(3) = liftForce(3) * b1 + b2
!!$        else
!!$           liftForce(1) = liftForce(1) * 0.0524 * sqrt(beta*Rep)
!!$           liftForce(2) = liftForce(2) * 0.0524 * sqrt(beta*Rep)
!!$           liftForce(3) = liftForce(3) * 0.0524 * sqrt(beta*Rep)
!!$        end if
  end if

  ! Position right-hand side
  dxdt = 0.0_WP
  dxdt(1:nDimensions) = part%velocity(1:nDimensions)

  ! Get the total particle acceleration
  dudt = dragForce + stress / particleDensity + gravityForce + collisionForce + liftForce
  
  ! Added mass contribution
  if (useAddedMass) then
     ! Mach/vol frac correction to the added mass coefficient from Koneru & Balachandar (2020)
     if (Ma .lt. 0.6_WP) then
        eta = 1.0_WP + 1.8_WP * Ma**2 + 7.6_WP * Ma**4
     else
        eta = 2.633_WP
     end if
     Cm = 0.5_WP * eta * (1.0_WP + 2.0_WP * (1.0_WP - vf)) / vf
     if (part%id .gt. 0) then
        rhs = dudt
        addedMass(1:nDimensions) = Cm / particleDensity * (stress(1:nDimensions) +              &
             part%velocity(1:nDimensions) * divRhoU - part%velocity(1:nDimensions) *            &
             sum(part%velocity(1:nDimensions) * gradRho(1:nDimensions)))
        dudt = (dudt + addedMass) / (1.0_WP + Cm * density / particleDensity)
        addedMass = dudt - rhs
     else
        addedMass = Cm * stress(1:nDimensions) / particleDensity
        dudt = dudt + addedMass
     end if
  end if

  ! Initialize mass/heat transfer
  dTdt = 0.0_WP
  dmdt = 0.0_WP
  Lv = 0.0_WP
  Cpd = particleSpecificHeat

  ! Compute mass/heat transfer sources
  if (usePhaseChange) then
     ! ---------------------------------------------------------------
     ! Taken from Russo et al. JFM (2014)
     ! Modeled from Bird et al. 1960, only valid for water/air mixture

     ! Non-dimensional numbers
     Sc = 1.0_WP / schmidtNumberInverse(vaporIndex)
     Nu = 2.0_WP+0.6_WP*sqrt(Rep)*(Pr)**(1.0_WP/3.0_WP)
     Sh = 2.0_WP+0.6_WP*sqrt(Rep)*(Sc)**(1.0_WP/3.0_WP)

     ! Thermal properties (use dimensional temperature!)
     Td =  part%temperature * (ratioOfSpecificHeats - 1.0_WP) * 293.15_WP
     call get_Cp('air', Td, Cpg)
     call get_Cp('water vapor', Td, Cpv)
     call get_Cp('water', Td, Cpd)
     Cpmix = Yv*Cpv+(1.0_WP-Yv)*Cpg
     call get_latent_heat('water', Td, Lv)

     ! Saturation vapour mass fraction using Antoine's relation
     Pref = 1.0_WP / ratioOfSpecificHeats
     Yvs = Pref / pressure * exp(11.6834_WP - 3816.44_WP / (Td - 46.28_WP))

     ! Mass transfer
     mass = particleDensity * pi * part%diameter**3 / 6.0_WP
     dmdt = -Sh * mass * log((1.0_WP - Yv) / (1.0_WP - Yvs)) * responseTimeInverse /         &
          (3.0_WP * Sc)

     ! Heat transfer
     dTdt = Nu * Cpmix * (temperature - part%temperature) * responseTimeInverse /            &
          (3.0_WP * Pr * Cpd) + dmdt / mass * Lv / Cpd /                                     &
          (ratioOfSpecificHeats - 1.0_WP) / 293.15_WP

     ! Non-dimensionalize thermal properties for interphase exchange
     Cpd = Cpd / Cpmix
     Lv = Lv / Cpmix / (ratioOfSpecificHeats - 1.0_WP) / 293.15_WP

  elseif (useParticleHeat) then
     !----------------------------
     ! Nusselt number (Gunn, 1978)
     Nu = (7.0_WP - 10.0_WP*vf + 5.0_WP*vf**2) * (1.0_WP + 0.7_WP*Rep**(0.2_WP) *            &
          Pr**(1.0_WP/3.0_WP)) + (1.33_WP - 2.4_WP*vf + 1.2_WP*vf**2) *                      &
          Rep**(0.7_WP)*Pr**(1.0_WP/3.0_WP)

     dTdt = Nu * responseTimeInverse * prandtlNumberInverse /                                &
          (3.0_WP * particleSpecificHeat) * (temperature - part%temperature)
  end if

!!$  ! Remove small droplets and transfer to fluid
!!$  if (usePhaseChange .and. part%diameter.le.diameterCutoff) then
!!$     part%stop = 1
!!$     !dmdt=max(-mass/dt_solve, 0.0_WP)
!!$  end if

  ! Transfer momentum to the fluid phase
  momentumExchange = -dragForce - liftForce - addedMass
  heatExchange = -dTdt
  massExchange = -dmdt
  call transfer_interphase_exchange(part, momentumExchange, heatExchange, massExchange,      &
       -hTorque, Cpd, Lv)

  ! Store particle forces for monitor
  monitorParticle(1:nDimensions)%drag = monitorParticle(1:nDimensions)%drag +                &
       dragForce(1:nDimensions)
  monitorParticle(1:nDimensions)%stress = monitorParticle(1:nDimensions)%stress +            &
       stress(1:nDimensions) / particleDensity
  monitorParticle(1:nDimensions)%gravity = monitorParticle(1:nDimensions)%gravity +          &
       gravityForce(1:nDimensions)
  monitorParticle(1:nDimensions)%collision = monitorParticle(1:nDimensions)%collision +      &
       collisionForce(1:nDimensions)
  monitorParticle(1:nDimensions)%lift = monitorParticle(1:nDimensions)%lift +                &
       liftForce(1:nDimensions)
  monitorParticle(1:nDimensions)%addedMass = monitorParticle(1:nDimensions)%addedMass +      &
       addedMass(1:nDimensions)

  return
end subroutine particle_rhs
