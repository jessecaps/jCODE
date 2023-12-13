module ibm_solver

  ! External modules
  use precision
  use ibm

  implicit none

  ! Global variables
  logical :: ibm_collisions, collisionWall(3,2)
  real(WP) :: ibmDensity, coefficientOfRestitution, coefficientOfFriction, collisionTime

end module ibm_solver


! ==================== !
! Setup the IBM solver !
! ==================== !
subroutine ibm_solver_setup

  ! Internal modules
  use ibm_solver

  ! External modules
  use parser
  use parallel
  use geometry
  
  implicit none

  ! Local variables
  logical :: wall

  ! Read in density
  if (ibm_move) call parser_read('ibm density', ibmDensity)

  ! Account for collisions
  call parser_read('ibm collisions', ibm_collisions, .false.)

  if (ibm_collisions) then
     ! Get the collision parameters
     call parser_read('coefficient of restitution', coefficientOfRestitution, 0.85_WP)
     call parser_read('coefficient of friction', coefficientOfFriction, 0.1_WP)
     call parser_read('collision time', collisionTime)

     ! Check for walls
     collisionWall = .false.
     call parser_read('ibm wall in x', wall, .false.)
     if (wall) collisionWall(1,:) = .true.
     call parser_read('ibm wall in y', wall, .false.)
     if (wall) collisionWall(2,:) = .true.
     call parser_read('ibm wall in z', wall, .false.)
     if (wall) collisionWall(3,:) = .true.
  else
     collisionWall = .false.
  end if

  return
end subroutine ibm_solver_setup


! ====================== !
! Cleanup the IBM solver !
! ====================== !
subroutine ibm_solver_cleanup

  ! Internal modules
  use ibm_solver
  
  implicit none

  return
end subroutine ibm_solver_cleanup


! ================================= !
! Compute IBM right-hand side terms !
! ================================= !
subroutine ibm_solver_rhs

  ! Internal modules
  use ibm_solver

  ! External modules
  use parallel
  use math
  use gravity_source

  implicit none

  ! Local variables
  integer :: i
  real(WP) :: mass, momentOfInertia

  if (.not. ibm_move) return

  ! Compute collision force on objects
  call compute_ibm_collisions

  ! Get forces acting on IBM objects
  call ibm_integrate_forces

  ! Update object acceleration
  do i = 1, nObjects

     ! Get object mass
     mass = ibmDensity * object(i)%volume

     ! Get moment of inertia
     select case (nDimensions)
     case (2)
        momentOfInertia = 0.5_WP * object(i)%volume / pi
     case (3)
        momentOfInertia = 0.1_WP * (6.0_WP * object(i)%volume / pi)**(2.0_WP / 3.0_WP)
     end select

     ! Fluid stresses
     object(i)%dudt(1:nDimensions) = (object(i)%pForce(1:nDimensions) +                      &
          object(i)%vForce(1:nDimensions)) / mass

     ! Collisions
     object(i)%dudt(1:nDimensions) = object(i)%dudt(1:nDimensions) +                         &
          object(i)%cForce(1:nDimensions) / mass

     ! Gravity
     if (allocated(gravity)) object(i)%dudt(1:nDimensions) = object(i)%dudt(1:nDimensions) + &
          gravity(1:nDimensions)

     ! Torque
     object(i)%dwdt(1:3) = (object(i)%hTorque(1:3) + object(i)%cTorque(1:3)) / momentOfInertia
     
  end do

  return
end subroutine ibm_solver_rhs


! ====================================== !
! Compute collision force on IBM objects !
! ====================================== !
subroutine compute_ibm_collisions

  ! Internal modules
  use ibm_solver

  ! External modules
  use simulation_flags
  use math
  use geometry
  use grid
  use grid_functions

  implicit none

  ! Local variables
  integer :: ip, jp, n, gridIndex
  real(WP), parameter :: oneSixth = 1.0_WP / 6.0_WP
  real(WP) :: d1, d2, mass1, mass2, objectSeparation, radiusOfInfluence, delta, rtv,         &
       effectiveMass, springForce, damper
  real(WP), dimension(3) :: dist, pos1, pos2, vel1, vel2, vel12, omega1, omega2,             &
       n12, v12n, t12, normalForce, tangentialForce

  ! Only implemented for cylinders and spheres for now
  if (.not. useIBM) return
  if (.not. ibm_collisions) return
  if (ibmType.ne.IBM_PARTICLE) return

  ! Zero out the forces
  do ip = 1, nObjects
     object(ip)%cForce = 0.0_WP
     object(ip)%cTorque = 0.0_WP
  end do

  ! Implement the collision term
  do ip = 1, nObjects

     ! Extract `particle 1` data
     pos1   = object(ip)%position
     vel1   = object(ip)%velocity
     omega1 = object(ip)%angularVelocity
     mass1 = ibmDensity * object(ip)%volume
     d1 = (2.0_WP * real(nDimensions, WP) * object(ip)%volume / pi)                          &
          ** (1.0_WP / real(nDimensions, WP))

     ! Particle-wall collisions
     ! ------------------------

     ! Closest wall parameters
     vel2   = 0.0_WP
     omega2 = 0.0_WP
     d2     = 0.0_WP

     ! Loop through each wall
     do n = 1, nDimensions

        ! Check collisions with left wall
        if (collisionWall(n,1)) then
           gridIndex = iStart(1) - gridOffset(1) + localGridSize(1) *                        &
                (iStart(2) - 1 - gridOffset(2) + localGridSize(2) *                          &
                (iStart(3) - 1 - gridOffset(3)))
           pos2 = pos1
           pos2(n) = coordinates(gridIndex, n)

           ! Distance between particle and wall
           dist = pos2 - pos1
           call correct_periodic_distance(dist)
           objectSeparation = sqrt(sum(dist**2))

           ! Particle-wall radius of influence
           radiusOfInfluence = gridSpacing(gridIndex, n)

           ! Distance of influence
           delta = 0.5_WP * d1 + radiusOfInfluence - objectSeparation
           delta = min(delta, 2.0_WP * minGridSpacing)

           if (delta .gt. 0.0_WP) then
              ! Collision parameters
              effectiveMass = mass1
              springForce = effectiveMass / collisionTime**2 * (pi**2 +                      &
                   log(coefficientOfRestitution)**2)
              damper = -2.0_WP * log(coefficientOfRestitution) *                             &
                   effectiveMass / collisionTime

              ! Normal force
              n12 = dist / objectSeparation
              vel12 = vel1 + cross_product(0.5_WP * d1 * omega1, n12) - (                    &
                   vel2 + cross_product(0.5_WP * d2 * omega2, -n12))
              v12n = sum(vel12 * n12) * n12
              normalForce = -springForce * delta * n12 - damper * v12n

              ! Tangential force
              t12 = vel12 - v12n
              rtv = sqrt(sum(t12 * t12))
              tangentialForce = 0.0_WP
              if (rtv.gt.0.0_WP) tangentialForce = -coefficientOfFriction *                  &
                   sqrt(sum(normalForce * normalForce)) * t12 / rtv

              ! Update the force
              object(ip)%cForce(1:nDimensions) = object(ip)%cForce(1:nDimensions) +          &
                   normalForce(1:nDimensions) + tangentialForce(1:nDimensions)
              
              ! Update the torque
              object(ip)%cTorque(1) = object(ip)%cTorque(1) + 0.5_WP * (d1 * n12(2) *        &
                   tangentialForce(3) - d1 * n12(3) * tangentialForce(2)) / mass1
              object(ip)%cTorque(2) = object(ip)%cTorque(2) + 0.5_WP * (d1 * n12(3) *        &
                   tangentialForce(1) - d1 * n12(1) * tangentialForce(3)) / mass1
              object(ip)%cTorque(3) = object(ip)%cTorque(3) + 0.5_WP * (d1 * n12(1) *        &
                   tangentialForce(2) - d1 * n12(2) * tangentialForce(1)) / mass1
              
           end if
        end if

        ! Check collisions with right wall
        if (collisionWall(n,2)) then
           gridIndex = iEnd(1) - gridOffset(1) + localGridSize(1) *                          &
                (iEnd(2) - 1 - gridOffset(2) + localGridSize(2) *                            &
                (iEnd(3) - 1 - gridOffset(3)))
           pos2 = pos1
           pos2(n) = coordinates(gridIndex, n)

           ! Distance between particle and wall
           dist = pos2 - pos1
           call correct_periodic_distance(dist)
           objectSeparation = sqrt(sum(dist**2))

           ! Particle-wall radius of influence
           radiusOfInfluence = gridSpacing(gridIndex, n)

           ! Distance of influence
           delta = 0.5_WP * d1 + radiusOfInfluence - objectSeparation
           delta = min(delta, 2.0_WP * minGridSpacing)

           if (delta .gt. 0.0_WP) then
              ! Collision parameters
              springForce = effectiveMass / collisionTime**2 * (pi**2 +                      &
                   log(coefficientOfRestitution)**2)
              damper = -2.0_WP * log(coefficientOfRestitution) *                             &
                   effectiveMass / collisionTime

              ! Normal force
              n12 = dist / objectSeparation
              vel12 = vel1 + cross_product(0.5_WP * d1 * omega1, n12) - (                    &
                   vel2 + cross_product(0.5_WP * d2 * omega2, -n12))
              v12n = sum(vel12 * n12) * n12
              normalForce = -springForce * delta * n12 - damper * v12n

              ! Tangential force
              t12 = vel12 - v12n
              rtv = sqrt(sum(t12 * t12))
              tangentialForce = 0.0_WP
              if (rtv.gt.0.0_WP) tangentialForce = -coefficientOfFriction *                  &
                   sqrt(sum(normalForce * normalForce)) * t12 / rtv

              ! Update the force
              object(ip)%cForce(1:nDimensions) = object(ip)%cForce(1:nDimensions) +          &
                   normalForce(1:nDimensions) + tangentialForce(1:nDimensions)

              ! Update the torque
              object(ip)%cTorque(1) = object(ip)%cTorque(1) + 0.5_WP * (d1 * n12(2) *        &
                   tangentialForce(3) - d1 * n12(3) * tangentialForce(2)) / mass1
              object(ip)%cTorque(2) = object(ip)%cTorque(2) + 0.5_WP * (d1 * n12(3) *        &
                   tangentialForce(1) - d1 * n12(1) * tangentialForce(3)) / mass1
              object(ip)%cTorque(3) = object(ip)%cTorque(3) + 0.5_WP * (d1 * n12(1) *        &
                   tangentialForce(2) - d1 * n12(2) * tangentialForce(1)) / mass1
           end if
        end if

     end do

     ! Particle-particle collisions
     ! ----------------------------

     ! Loop over all other particles
     do jp = ip + 1, nObjects

        ! Extract `particle 2` data
        pos2   = object(jp)%position
        vel2   = object(jp)%velocity
        omega2 = object(jp)%angularVelocity
        mass2  = ibmDensity * object(jp)%volume
        d2     = (2.0_WP * real(nDimensions, WP) * object(jp)%volume / pi)                   &
             ** (1.0_WP / real(nDimensions, WP))        

        ! Distance between particles `1` and `2`
        dist = pos2 - pos1
        call correct_periodic_distance(dist)
        objectSeparation = sqrt(sum(dist**2))

        ! Set radius of influence based on grid spacing
        radiusOfInfluence = sqrt(real(nDimensions, WP)) * minGridSpacing

        ! Distance of influence
        delta = 0.5_WP * (d1 + d2) + radiusOfInfluence - objectSeparation
        delta = min(delta, 2.0_WP * minGridSpacing)

        if (delta .gt. 0.0_WP) then
           ! Collision parameters
           effectiveMass = mass1 * mass2 / (mass1 + mass2)
           springForce = effectiveMass / collisionTime**2 * (pi**2 +                         &
                log(coefficientOfRestitution)**2)
           damper = -2.0_WP * log(coefficientOfRestitution) *                                &
                effectiveMass / collisionTime

           ! Normal force
           n12 = dist / objectSeparation
           vel12 = vel1 + cross_product(0.5_WP * d1 * omega1, n12) - (                       &
                vel2 + cross_product(0.5_WP * d2 * omega2, -n12))
           v12n = sum(vel12 * n12) * n12
           normalForce = -springForce * delta * n12 - damper * v12n

           ! Tangential force
           t12 = vel12 - v12n
           rtv = sqrt(sum(t12 * t12))
           tangentialForce = 0.0_WP
           if (rtv.gt.0.0_WP) tangentialForce = -coefficientOfFriction *                     &
                sqrt(sum(normalForce * normalForce)) * t12 / rtv

           ! Update the force on particles `i` and `j`
           object(ip)%cForce(1:nDimensions) = object(ip)%cForce(1:nDimensions) +             &
                normalForce(1:nDimensions) + tangentialForce(1:nDimensions)
           
           object(jp)%cForce(1:nDimensions) = object(jp)%cForce(1:nDimensions) -             &
                (normalForce(1:nDimensions) + tangentialForce(1:nDimensions))

           ! Update the torque on particles `i` and `j`
           object(ip)%cTorque(1) = object(ip)%cTorque(1) + 0.5_WP * (d1 * n12(2) *           &
                tangentialForce(3) - d1 * n12(3) * tangentialForce(2)) / mass1
           object(ip)%cTorque(2) = object(ip)%cTorque(2) + 0.5_WP * (d1 * n12(3) *           &
                tangentialForce(1) - d1 * n12(1) * tangentialForce(3)) / mass1
           object(ip)%cTorque(3) = object(ip)%cTorque(3) + 0.5_WP * (d1 * n12(1) *           &
                tangentialForce(2) - d1 * n12(2) * tangentialForce(1)) / mass1

           object(jp)%cTorque(1) = object(ip)%cTorque(1) - 0.5_WP * (d1 * n12(2) *           &
                tangentialForce(3) - d1 * n12(3) * tangentialForce(2)) / mass2
           object(jp)%cTorque(2) = object(ip)%cTorque(2) - 0.5_WP * (d1 * n12(3) *           &
                tangentialForce(1) - d1 * n12(1) * tangentialForce(3)) / mass2
           object(jp)%cTorque(3) = object(ip)%cTorque(3) - 0.5_WP * (d1 * n12(1) *           &
                tangentialForce(2) - d1 * n12(2) * tangentialForce(1)) / mass2
        end if
     end do

  end do
  
  return
  
end subroutine compute_ibm_collisions


! ========================================= !
! Correct a object position for periodicity !
! ========================================= !
subroutine correct_object_position(obj)

  ! Internal modules
  use ibm_solver

  ! External modules
  use geometry

  implicit none

  ! Arguments
  type(t_Object), intent(inout) :: obj

  ! Local variables
  integer :: i

  do i = 1, nDimensions
     if (isPeriodic(i)) obj%position(i) = domainExtent(i,1) +                              &
          modulo(obj%position(i) - domainExtent(i,1), periodicLength(i))
  end do

  return
end subroutine correct_object_position


! ====================================== !
! Check if an object has left the domain !
! ====================================== !
subroutine check_object_bounds(obj)

  ! Internal modules
  use ibm_solver

  ! External modules
  use math
  use geometry

  implicit none

  ! Arguments
  type(t_Object), intent(inout) :: obj

  ! Local variables
  integer :: i
  real(WP) :: R

  if (.not.ibm_move .or. obj%remove) return

  if (ibmType .eq. IBM_PARTICLE) then
     R = 0.5_WP * (2.0_WP * real(nDimensions, WP) * obj%volume / pi)**(1.0_WP /              &
          real(nDimensions, WP))
  else
     R = 0.0_WP
  end if

  do i = 1, nDimensions
     if (.not. isPeriodic(i)) then
        if (obj%position(i) + R .lt. domainExtent(i,1) .or.                                  &
             obj%position(i) - R .gt. domainExtent(i,2)) obj%remove = .true.
     end if
  end do

  return
end subroutine check_object_bounds
