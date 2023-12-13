module ibm_exchange

  ! External modules
  use precision
  use simulation_flags
  use geometry
  use grid
  
  implicit none

  ! Global variables
  integer :: nO, nOverlap(3,2), imin, imax, jmin, jmax, kmin, kmax, nBand
  integer, dimension(:,:,:), allocatable :: nMarkerInCell
  integer, dimension(:,:,:,:), allocatable :: markerInCell
  real(WP), dimension(:,:), allocatable :: cartesianCoordinates
  real(WP), dimension(:,:,:,:), allocatable :: gridDensity, gridPressure, gridVelocity,      &
       gridViscousStress, gridLevelset, gridDensityGradient, gridEnergyGradient,             &
       gridNormalVelocityGradient, gridBuffer, ibmSource

contains

  ! Compute regularized delta function
  ! ----------------------------------
  subroutine ibm_delta_function(delta,ic,jc,kc,xp,yp,zp,dxi)
    implicit none
    
    ! Arguments
    real(WP), intent(out) :: delta
    integer, intent(in) :: ic, jc, kc
    real(WP), intent(in) :: xp, yp, zp, dxi(3)

    ! Local variables
    real(WP) :: xc, yc, zc
    real(WP) :: deltax, deltay, deltaz, r

    ! Initialize delta function
    deltax = 1.0_WP; deltay = 1.0_WP; deltaz = 1.0_WP

    ! Compute in x
    xc = cartesianCoordinates(ic, 1)
    r = (xp - xc) * dxi(1) * 0.6_WP
    deltax = roma_kernel(r) * dxi(1)

    ! Compute in y
    if (nDimensions.gt.1) then
       yc = cartesianCoordinates(jc, 2)
       r = (yp - yc) * dxi(2) * 0.6_WP
       deltay = roma_kernel(r) * dxi(2)
    end if

    ! Compute in z
    if (nDimensions.gt.2) then
       zc = cartesianCoordinates(kc, 3)
       r = (zp - zc) * dxi(3) * 0.6_WP
       deltaz = roma_kernel(r) * dxi(3)
    end if

    ! Put it all together
    delta = deltax * deltay * deltaz

    return
  end subroutine ibm_delta_function
  

  ! Mollification kernels
  ! Roma A, Peskin C and Berger M 1999 J. Comput. Phys. 153 509–534
  ! ---------------------------------------------------------------
  function roma_kernel(d) result(phi)
    implicit none

    ! Arguments
    real(WP), intent(in) :: d
    real(WP)             :: phi

    ! Local arguments
    real(WP) :: r

    r = abs(d)
    if (r.le.0.5_WP) then
       phi = 1.0_WP/3.0_WP*(1.0_WP+sqrt(-3.0_WP*r**2+1.0_WP))
    else if (r.gt.0.5_WP .and. r.le.1.5_WP) then
       phi = 1.0_WP/6.0_WP*(5.0_WP-3.0_WP*r-sqrt(-3.0_WP*(1.0_WP-r)**2+1.0_WP))
    else
       phi = 0.0_WP
    end if

    return
  end function roma_kernel

  ! Juric D. Computations of phase change. University of Michigan; 1996. Ph.D. thesis
  ! ---------------------------------------------------------------------------------
  function juric_kernel(d) result(phi)
    implicit none
    
    ! Arguments
    real(WP), intent(in) :: d
    real(WP)             :: phi

    ! Local variables
    real(WP) :: r

    r = abs(d)
    if (r.le.1.0_WP) then
       phi = juric_function(r)
    else if (r.gt.1.0_WP .and. r.lt.2.0_WP) then
       phi = 0.5_WP - juric_function(2.0_WP - r)
    else
       phi = 0.0_WP
    end if

    return
  end function juric_kernel

  function juric_function(d) result(f)
    implicit none

    ! Arguments
    real(WP), intent(in) :: d
    real(WP)             :: f

    ! Local variables
    real(WP) :: r

    r = abs(d)
    f = (3.0_WP - 2.0_WP * r + sqrt(1.0_WP + 4.0_WP * r - 4.0_WP * r**2)) / 8.0_WP

    return
  end function juric_function
  
end module ibm_exchange


! ============================== !
! Setup the IBM exchange routine !
! ============================== !
subroutine ibm_exchange_setup(nOverlap_)

  ! Internal modules
  use ibm_exchange

  ! External modules
  use parallel
  use parser
  use solver_options
  use geometry

  ! Arguments
  integer, intent(inout) :: nOverlap_(3,2)
  
  ! Local variables
  integer :: gridIndex, i, j, n, ijk(3), ierror

  if (.not. useIBM) return

  ! Additional points for mollification
  nO = 3
  nOverlap = 3
  
  ! Adjust the ghost points at the boundaries
  do i = 1, 3
     if (.not. isPeriodic(i)) then
        if (procCoords(i) .eq. 0) nOverlap(i,1) = 0
        if (procCoords(i) .eq. nProcsDir(i) - 1) nOverlap(i,2) = 0
     end if
  end do

  ! Get local indices of grid arrays with overlap
  imin = gridOffset(1) + 1 - nOverlap(1,1)
  imax = gridOffset(1) + localGridSize(1) + nOverlap(1,2)
  if (nDimensions.gt.1) then
     jmin = gridOffset(2) + 1 - nOverlap(2,1)
     jmax = gridOffset(2) + localGridSize(2) + nOverlap(2,2)
  else
     jmin = 1
     jmax = 1
  end if
  if (nDimensions.gt.2) then
     kmin = gridOffset(3) + 1 - nOverlap(3,1)
     kmax = gridOffset(3) + localGridSize(3) + nOverlap(3,2)
  else
     kmin = 1
     kmax = 1
  end if

  ! Allocate grid arrays with extra layer of ghost cells
  allocate(gridLevelset(imin:imax, jmin:jmax, kmin:kmax, 1))
  allocate(gridDensity(imin:imax, jmin:jmax, kmin:kmax, 1))
  allocate(gridPressure(imin:imax, jmin:jmax, kmin:kmax, 1))
  allocate(gridVelocity(imin:imax, jmin:jmax, kmin:kmax, nDimensions))
  allocate(gridViscousStress(imin:imax, jmin:jmax, kmin:kmax, nDimensions**2))
  allocate(ibmSource(imin:imax, jmin:jmax, kmin:kmax, nUnknowns))
  allocate(gridBuffer(imin:imax, jmin:jmax, kmin:kmax, 3*nDimensions+5))
  allocate(gridDensityGradient(imin:imax, jmin:jmax, kmin:kmax, 1))
  allocate(gridEnergyGradient(imin:imax, jmin:jmax, kmin:kmax, 1))
  allocate(gridNormalVelocityGradient(imin:imax, jmin:jmax, kmin:kmax, nDimensions))
  
  ! Prepare Cartesian vectors for marker localization (use max overlap)
  nOverlap_ = max(nOverlap_, nOverlap)
  
  allocate(cartesianCoordinates(1 - maxval(nOverlap_) : maxval(globalGridSize) +             &
       maxval(nOverlap_), nDimensions))
  cartesianCoordinates = -huge(1.0_WP)
  do j = 1, nDimensions
     do i = iStart(j), iEnd(j)
        ijk = iStart; ijk(j) = i
        gridIndex = ijk(1) - gridOffset(1) + localGridSize(1) *                              &
             (ijk(2) - 1 - gridOffset(2) + localGridSize(2) *                                &
             (ijk(3) - 1 - gridOffset(3)))
        cartesianCoordinates(i, j) = coordinates(gridIndex, j)
     end do
     call MPI_Allreduce(MPI_IN_PLACE, cartesianCoordinates(:, j),                            &
          size(cartesianCoordinates, 1), MPI_REAL_WP, MPI_MAX,                               &
          commDir(j), ierror)
     if (isPeriodic(j) .and. procCoords(j) .eq. 0) then
        do n = 1, nOverlap_(j, 1)
           cartesianCoordinates(1 - n, j) = - periodicLength(j) +                            &
                cartesianCoordinates(globalGridSize(j) - n + 1, j)
        end do
     end if
     if (isPeriodic(j) .and. procCoords(j) .eq. nProcsDir(j) - 1) then
        do n = 1, nOverlap_(j, 2)
           cartesianCoordinates(globalGridSize(j) + n, j) = periodicLength(j) +              &
                cartesianCoordinates(n, j)
        end do
     end if
     call MPI_Allreduce(MPI_IN_PLACE, cartesianCoordinates(:, j),                            &
          size(cartesianCoordinates, 1), MPI_REAL_WP, MPI_MAX, commDir(j), ierror)
  end do
  
  ! Allocate the marker localization arrays (assume max 10 markers per cell)
  call parser_read('levelset band', nBand, 8)
  allocate(nMarkerInCell(imin:imax, jmin:jmax, kmin:kmax))
  allocate(markerInCell(imin:imax, jmin:jmax, kmin:kmax, 10))
     
  return
end subroutine ibm_exchange_setup


! ================================ !
! Cleanup the IBM exchange routine !
! ================================ !
subroutine ibm_exchange_cleanup

  ! Internal modules
  use ibm_exchange

  implicit none

  if (allocated(cartesianCoordinates)) deallocate(cartesianCoordinates)
  if (allocated(gridLevelset)) deallocate(gridLevelset)
  if (allocated(gridDensity)) deallocate(gridDensity)
  if (allocated(gridPressure)) deallocate(gridPressure)
  if (allocated(gridVelocity)) deallocate(gridVelocity)
  if (allocated(gridNormalVelocityGradient)) deallocate(gridNormalVelocityGradient)
  if (allocated(gridViscousStress)) deallocate(gridViscousStress)
  if (allocated(ibmSource)) deallocate(ibmSource)
  if (allocated(gridBuffer)) deallocate(gridBuffer)
  if (allocated(gridDensityGradient)) deallocate(gridDensityGradient)
  if (allocated(gridEnergyGradient)) deallocate(gridEnergyGradient)

  return
end subroutine ibm_exchange_cleanup


! ============================================== !
! Routine to send informaton from grid => marker !
! Takes in (xp,yp,zp) and (ip,jp,kp)             !
! as well as any 'big' array A                   !
! and returns Ap, interpolation of A @ p         !
! ============================================== !
subroutine ibm_interpolate(A,xp,yp,zp,ip,jp,kp,Ap,dir)

  ! Internal modules
  use ibm_exchange

  ! External modules
  use parallel
  use string

  implicit none

  ! Arguments
  real(WP), dimension(imin:imax, jmin:jmax, kmin:kmax), intent(in) :: A
  real(WP), intent(in) :: xp, yp, zp
  integer, intent(in) :: ip, jp, kp
  character(len=*), intent(in) :: dir
  real(WP), intent(out) :: Ap

  ! Local variables
  integer :: gridIndex
  integer :: di, dj, dk
  integer :: i1, i2, j1, j2, k1, k2
  real(WP) :: inverseGridSpacing(3), buf
  real(WP), dimension(-nO:+nO, -nO:+nO, -nO:+nO) :: delta

  ! Get inverse of local grid spacing
  gridIndex = ip - gridOffset(1) + localGridSize(1) *                                        &
       (jp - 1 - gridOffset(2) + localGridSize(2) *                                          &
       (kp - 1 - gridOffset(3)))
  inverseGridSpacing(1:nDimensions) = 1.0_WP / gridSpacing(gridIndex, 1:nDimensions)

  ! Get the neighboring points
  i1=ip-nO; i2=ip+nO
  j1=jp-nO; j2=jp+nO
  k1=kp-nO; k2=kp+nO

  ! Correct for walls
  select case (nDimensions)
  case (1)
     j1 = jp; j2 = jp
     k1 = kp; k2 = kp
     if (nOverlap(1,1) .eq. 0) i1 = max(i1, iStart(1))
     if (nOverlap(1,2) .eq. 0) i2 = min(i2, iEnd(1))
  case (2)
     k1 = kp; k2 = kp
     if (nOverlap(1,1) .eq. 0) i1 = max(i1, iStart(1))
     if (nOverlap(1,2) .eq. 0) i2 = min(i2, iEnd(1))
     if (nOverlap(2,1) .eq. 0) j1 = max(j1, iStart(2))
     if (nOverlap(2,2) .eq. 0) j2 = min(j2, iEnd(2))
  case (3)
     if (nOverlap(1,1) .eq. 0) i1 = max(i1, iStart(1))
     if (nOverlap(1,2) .eq. 0) i2 = min(i2, iEnd(1))
     if (nOverlap(2,1) .eq. 0) j1 = max(j1, iStart(2))
     if (nOverlap(2,2) .eq. 0) j2 = min(j2, iEnd(2))
     if (nOverlap(3,1) .eq. 0) k1 = max(k1, iStart(3))
     if (nOverlap(3,2) .eq. 0) k2 = min(k2, iEnd(3))
  end select

  ! Loop over neighboring cells and compute regularized delta function
  delta = 0.0_WP
  select case (trim(adjustl(dir)))
  case ('in')
     do dk = k1, k2
        do dj = j1, j2
           do di = i1, i2
              if (gridLevelset(di,dj,dk,1) .le. 0.0_WP) then
                 call ibm_delta_function(delta(di-ip,dj-jp,dk-kp), di, dj, dk, xp, yp, zp,   &
                      inverseGridSpacing)
              end if
           end do
        end do
     end do
  case ('out')
     do dk = k1, k2
        do dj = j1, j2
           do di = i1, i2
              if (gridLevelset(di,dj,dk,1) .gt. 0.0_WP) then
                 call ibm_delta_function(delta(di-ip,dj-jp,dk-kp), di, dj, dk, xp, yp, zp,   &
                      inverseGridSpacing)
              end if
           end do
        end do
     end do
  case ('inout')
     do dk = k1, k2
        do dj = j1, j2
           do di = i1, i2
              call ibm_delta_function(delta(di-ip,dj-jp,dk-kp), di, dj, dk, xp, yp, zp,      &
                   inverseGridSpacing)
           end do
        end do
     end do
  case default
     call die("ibm_interpolate: unknown direction '"//trim(dir)//"'")
  end select

  ! Normalize
  buf = sum(delta)
  if (buf.gt.0.0_WP) delta = delta / buf

  ! Perform the actual interpolation on Ap
  Ap = sum(delta(i1-ip:i2-ip, j1-jp:j2-jp, k1-kp:k2-kp) * A(i1:i2,j1:j2,k1:k2))! /            &
  !    jacobian(gridIndex, 1)

  return
end subroutine ibm_interpolate


! ============================================== !
! Routine to send informaton from marker => grid !
! ============================================== !
subroutine ibm_extrapolate(A,xp,yp,zp,ip,jp,kp,Ap)

  ! Internal modules
  use ibm_exchange

  ! External modules
  use parallel

  implicit none

  ! Arguments
  real(WP), dimension(imin:imax, jmin:jmax, kmin:kmax), intent(inout) :: A
  real(WP), intent(in) :: xp, yp, zp
  integer, intent(in) :: ip, jp, kp
  real(WP), intent(in) :: Ap

  ! Local variables
  integer :: i1, i2, j1, j2, k1, k2, di, dj, dk, gridIndex
  integer :: iw
  real(WP) :: inverseGridSpacing(3), buf, xpi
  real(WP), dimension(-nO:+nO, -nO:+nO, -nO:+nO) :: delta, deltaImage

  ! Get inverse of local grid spacing
  gridIndex = ip - gridOffset(1) + localGridSize(1) *                                        &
       (jp - 1 - gridOffset(2) + localGridSize(2) *                                          &
       (kp - 1 - gridOffset(3)))
  inverseGridSpacing(1:nDimensions) = 1.0_WP / gridSpacing(gridIndex, 1:nDimensions)

  ! Get the neighboring points
  i1=ip-nO; i2=ip+nO
  j1=jp-nO; j2=jp+nO
  k1=kp-nO; k2=kp+nO

  ! Correct for walls
  select case (nDimensions)
  case (1)
     j1 = jp; j2 = jp
     k1 = kp; k2 = kp
     if (nOverlap(1,1) .eq. 0) i1 = max(i1, iStart(1))
     if (nOverlap(1,2) .eq. 0) i2 = min(i2, iEnd(1))
  case (2)
     k1 = kp; k2 = kp
     if (nOverlap(1,1) .eq. 0) i1 = max(i1, iStart(1))
     if (nOverlap(1,2) .eq. 0) i2 = min(i2, iEnd(1))
     if (nOverlap(2,1) .eq. 0) j1 = max(j1, iStart(2))
     if (nOverlap(2,2) .eq. 0) j2 = min(j2, iEnd(2))
  case (3)
     if (nOverlap(1,1) .eq. 0) i1 = max(i1, iStart(1))
     if (nOverlap(1,2) .eq. 0) i2 = min(i2, iEnd(1))
     if (nOverlap(2,1) .eq. 0) j1 = max(j1, iStart(2))
     if (nOverlap(2,2) .eq. 0) j2 = min(j2, iEnd(2))
     if (nOverlap(3,1) .eq. 0) k1 = max(k1, iStart(3))
     if (nOverlap(3,2) .eq. 0) k2 = min(k2, iEnd(3))
  end select

  ! Loop over neighboring cells and compute regularized delta function
  delta = 0.0_WP
  do dk = k1, k2
     do dj = j1, j2
        do di = i1, i2
           if (gridLevelset(di,dj,dk,1) .le. 0.0_WP) then
              call ibm_delta_function(delta(di-ip,dj-jp,dk-kp), di, dj, dk, xp, yp, zp,         &
                   inverseGridSpacing)
           end if
        end do
     end do
  end do

  ! Get image of delta across non-periodic boundary
  deltaImage = 0.0_WP
  if (.not. isPeriodic(1)) then
     if (ip.le.3) then
        iw = iStart(1)
     elseif (ip.ge.globalGridSize(1)-1) then
        iw=iEnd(1)
     else
        iw = 0
     end if
     if (iw .gt. 0) then
        xpi = xp - 2.0_WP * (xp - cartesianCoordinates(iw, 1))
        do dk = k1, k2
           do dj = j1, j2
              do di = i1, i2
                 if (gridLevelset(di,dj,dk,1) .le. 0.0_WP) then
                    call ibm_delta_function(deltaImage(di-ip,dj-jp,dk-kp), di, dj, dk, xpi,  &
                         yp, zp, inverseGridSpacing)
                 end if
              end do
           end do
        end do
     end if
  end if
  delta = delta + deltaImage

  ! Normalize
  buf = sum(delta)
  if (buf.gt.0.0_WP) delta = delta / buf * jacobian(gridIndex, 1)

  ! Perform the actual extrapolation on A
  A(i1:i2, j1:j2, k1:k2) = A(i1:i2, j1:j2, k1:k2) +                                          &
       delta(i1-ip:i2-ip, j1-jp:j2-jp, k1-kp:k2-kp) * Ap

  return
end subroutine ibm_extrapolate


! =========================================================== !
! Same as ibm_extrapolate but acts on array of size nUnknowns !
! =========================================================== !
subroutine ibm_extrapolate_source(A,xp,yp,zp,ip,jp,kp,Ap)

  ! Internal modules
  use ibm_exchange

  ! External modules
  use parallel
  use solver_options

  implicit none

  ! Arguments
  real(WP), dimension(imin:imax, jmin:jmax, kmin:kmax, nUnknowns), intent(inout) :: A
  real(WP), intent(in) :: xp, yp, zp
  integer, intent(in) :: ip, jp, kp
  real(WP), dimension(nUnknowns), intent(in) :: Ap

  ! Local variables
  integer :: i1, i2, j1, j2, k1, k2, di, dj, dk, gridIndex
  integer :: n, iw
  real(WP) :: inverseGridSpacing(3), buf, xpi
  real(WP), dimension(-nO:+nO, -nO:+nO, -nO:+nO) :: delta, deltaImage

  ! Get inverse of local grid spacing
  gridIndex = ip - gridOffset(1) + localGridSize(1) *                                        &
       (jp - 1 - gridOffset(2) + localGridSize(2) *                                          &
       (kp - 1 - gridOffset(3)))
  inverseGridSpacing(1:nDimensions) = 1.0_WP / gridSpacing(gridIndex, 1:nDimensions)

  ! Get the neighboring points
  i1=ip-nO; i2=ip+nO
  j1=jp-nO; j2=jp+nO
  k1=kp-nO; k2=kp+nO

  ! Correct for walls
  select case (nDimensions)
  case (1)
     j1 = jp; j2 = jp
     k1 = kp; k2 = kp
     if (nOverlap(1,1) .eq. 0) i1 = max(i1, iStart(1))
     if (nOverlap(1,2) .eq. 0) i2 = min(i2, iEnd(1))
  case (2)
     k1 = kp; k2 = kp
     if (nOverlap(1,1) .eq. 0) i1 = max(i1, iStart(1))
     if (nOverlap(1,2) .eq. 0) i2 = min(i2, iEnd(1))
     if (nOverlap(2,1) .eq. 0) j1 = max(j1, iStart(2))
     if (nOverlap(2,2) .eq. 0) j2 = min(j2, iEnd(2))
  case (3)
     if (nOverlap(1,1) .eq. 0) i1 = max(i1, iStart(1))
     if (nOverlap(1,2) .eq. 0) i2 = min(i2, iEnd(1))
     if (nOverlap(2,1) .eq. 0) j1 = max(j1, iStart(2))
     if (nOverlap(2,2) .eq. 0) j2 = min(j2, iEnd(2))
     if (nOverlap(3,1) .eq. 0) k1 = max(k1, iStart(3))
     if (nOverlap(3,2) .eq. 0) k2 = min(k2, iEnd(3))
  end select

  ! Loop over neighboring cells and compute regularized delta function
  delta = 0.0_WP
  do dk = k1, k2
     do dj = j1, j2
        do di = i1, i2
          if (gridLevelset(di,dj,dk,1) .le. 0.0_WP) then
              call ibm_delta_function(delta(di-ip,dj-jp,dk-kp), di, dj, dk, xp, yp, zp,         &
                   inverseGridSpacing)
           end if
        end do
     end do
  end do

  ! Get image of delta across non-periodic boundary
  deltaImage = 0.0_WP
  if (.not. isPeriodic(1)) then
     if (ip.le.3) then
        iw = iStart(1)
     elseif (ip.ge.globalGridSize(1)-1) then
        iw=iEnd(1)
     else
        iw = 0
     end if
     if (iw .gt. 0) then
        xpi = xp - 2.0_WP * (xp - cartesianCoordinates(iw, 1))
        do dk = k1, k2
           do dj = j1, j2
              do di = i1, i2
                 if (gridLevelset(di,dj,dk,1) .le. 0.0_WP) then
                    call ibm_delta_function(deltaImage(di-ip,dj-jp,dk-kp), di, dj, dk, xpi,  &
                         yp, zp, inverseGridSpacing)
                 end if
              end do
           end do
        end do
     end if
  end if
  delta = delta + deltaImage

  ! Normalize
  buf = sum(delta)
  if (buf.gt.0.0_WP) delta = delta / buf * jacobian(gridIndex, 1)

  ! Perform the actual extrapolation on A
  do n = 1, nUnknowns
     A(i1:i2, j1:j2, k1:k2, n) = A(i1:i2, j1:j2, k1:k2, n) +                                 &
          delta(i1-ip:i2-ip, j1-jp:j2-jp, k1-kp:k2-kp) * Ap(n)
  end do

  return
end subroutine ibm_extrapolate_source


! ===================================== !
! Interpolate from grid to marker using !
! inverse distance function             !
! ===================================== !
subroutine ibm_inverse_distance_function(A,xp,yp,zp,ip,jp,kp,Ap)

  ! Internal modules
  use ibm_exchange

  ! External modules
  use parallel
  use geometry

  implicit none

  ! Arguments
  real(WP), dimension(imin:imax, jmin:jmax, kmin:kmax), intent(in) :: A
  real(WP), intent(in) :: xp, yp, zp
  integer, intent(in) :: ip, jp, kp
  real(WP), intent(out) :: Ap

    ! Local variables
    integer :: i1, j1, k1, i2, j2, k2
    real(WP) :: interpCoeff(2,2,2), dist(2,2,2), eta(2,2,2), alpha(2,2,2), buf
    real(WP), parameter :: eps=1.0e-10_WP

    select case (nDimensions)

    case (2)

       ! Get the interpolation points
       i1 = ip; i2 = i1 + 1
       j1 = jp; j2 = j1 + 1

       ! Get interpolation weights (Chaudhuri et al. 2011, JCP)
       dist(1,1,1) = sqrt(                                                                   &
            (cartesianCoordinates(i1, 1) - xp)**2 +                                          &
            (cartesianCoordinates(j1, 2) - yp)**2 )
       dist(2,1,1) = sqrt(                                                                   &
            (cartesianCoordinates(i2, 1) - xp)**2 +                                          &
            (cartesianCoordinates(j1, 2) - yp)**2 )
       dist(1,2,1) = sqrt(                                                                   &
            (cartesianCoordinates(i1, 1) - xp)**2 +                                          &
            (cartesianCoordinates(j2, 2) - yp)**2 )
       dist(2,2,1) = sqrt(                                                                   &
            (cartesianCoordinates(i2, 1) - xp)**2 +                                          &
            (cartesianCoordinates(j2, 2) - yp)**2 )
       interpCoeff = 0.0_WP
       if (dist(1,1,1) .le. eps * minGridSpacing) then
          interpCoeff(1,1,1) = 1.0_WP
       else if (dist(2,1,1) .le. eps * minGridSpacing) then
          interpCoeff(2,1,1) = 1.0_WP
       else if (dist(1,2,1) .le. eps * minGridSpacing) then
          interpCoeff(1,2,1) = 1.0_WP
       else if (dist(2,2,1) .le. eps * minGridSpacing) then
          interpCoeff(2,2,1) = 1.0_WP
       else
          eta(:,:,1) = 1.0_WP / dist(:,:,1)**2
          alpha = 1.0_WP
          if (gridLevelset(i1,j1,1,1).lt.0.0_WP) alpha(1,1,1) = 0.0_WP
          if (gridLevelset(i2,j1,1,1).lt.0.0_WP) alpha(2,1,1) = 0.0_WP
          if (gridLevelset(i1,j2,1,1).lt.0.0_WP) alpha(1,2,1) = 0.0_WP
          if (gridLevelset(i2,j2,1,1).lt.0.0_WP) alpha(2,2,1) = 0.0_WP
          buf = sum(alpha(:,:,1) * eta(:,:,1))
          interpCoeff(:,:,1) = alpha(:,:,1) * eta(:,:,1) / buf
       end if

       ! Interpolate neighboring points
       Ap = interpCoeff(1,1,1) * A(i1, j1, 1) + interpCoeff(2,1,1) * A(i2, j1, 1) +          &
            interpCoeff(1,2,1) * A(i1, j2, 1) + interpCoeff(2,2,1) * A(i2, j2, 1)

    case (3)

       ! Get the interpolation points.
       i1 = ip; i2 = i1 + 1
       j1 = jp; j2 = j1 + 1
       k1 = kp; k2 = k1 + 1

       ! Get interpolation weights (Chaudhuri et al. 2011, JCP)
       dist(1,1,1) = sqrt(                                                                   &
            (cartesianCoordinates(i1, 1) - xp)**2 +                                          &
            (cartesianCoordinates(j1, 2) - yp)**2 +                                          &
            (cartesianCoordinates(k1, 3) - zp)**2 )
       dist(2,1,1) = sqrt(                                                                   &
            (cartesianCoordinates(i2, 1) - xp)**2 +                                          &
            (cartesianCoordinates(j1, 2) - yp)**2 +                                          &
            (cartesianCoordinates(k1, 3) - zp)**2 )
       dist(1,2,1) = sqrt(                                                                   &
            (cartesianCoordinates(i1, 1) - xp)**2 +                                          &
            (cartesianCoordinates(j2, 2) - yp)**2 +                                          &
            (cartesianCoordinates(k1, 3) - zp)**2 )
       dist(2,2,1) = sqrt(                                                                   &
            (cartesianCoordinates(i2, 1) - xp)**2 +                                          &
            (cartesianCoordinates(j2, 2) - yp)**2 +                                          &
            (cartesianCoordinates(k1, 3) - zp)**2 )
       dist(1,1,2) = sqrt(                                                                   &
            (cartesianCoordinates(i1, 1) - xp)**2 +                                          &
            (cartesianCoordinates(j1, 2) - yp)**2 +                                          &
            (cartesianCoordinates(k2, 3) - zp)**2 )
       dist(2,1,2) = sqrt(                                                                   &
            (cartesianCoordinates(i2, 1) - xp)**2 +                                          &
            (cartesianCoordinates(j1, 2) - yp)**2 +                                          &
            (cartesianCoordinates(k2, 3) - zp)**2 )
       dist(1,2,2) = sqrt(                                                                   &
            (cartesianCoordinates(i1, 1) - xp)**2 +                                          &
            (cartesianCoordinates(j2, 2) - yp)**2 +                                          &
            (cartesianCoordinates(k2, 3) - zp)**2 )
       dist(2,2,2) = sqrt(                                                                   &
            (cartesianCoordinates(i2, 1) - xp)**2 +                                          &
            (cartesianCoordinates(j2, 2) - yp)**2 +                                          &
            (cartesianCoordinates(k2, 3) - zp)**2 )
       interpCoeff = 0.0_WP
       if (dist(1,1,1) .le. eps * minGridSpacing) then
          interpCoeff(1,1,1) = 1.0_WP
       else if (dist(2,1,1) .le. eps * minGridSpacing) then
          interpCoeff(2,1,1) = 1.0_WP
       else if (dist(1,2,1) .le. eps * minGridSpacing) then
          interpCoeff(1,2,1) = 1.0_WP
       else if (dist(2,2,1) .le. eps * minGridSpacing) then
          interpCoeff(2,2,1) = 1.0_WP
       else if (dist(1,1,2) .le. eps * minGridSpacing) then
          interpCoeff(1,1,2) = 1.0_WP
       else if (dist(2,1,2) .le. eps * minGridSpacing) then
          interpCoeff(2,1,2) = 1.0_WP
       else if (dist(1,2,2) .le. eps * minGridSpacing) then
          interpCoeff(1,2,2) = 1.0_WP
       else if (dist(2,2,2) .le. eps * minGridSpacing) then
          interpCoeff(2,2,2) = 1.0_WP
       else
          eta = 1.0_WP / dist**2
          alpha = 1.0_WP
          if (gridLevelset(i1,j1,k1,1).lt.0.0_WP) alpha(1,1,1) = 0.0_WP
          if (gridLevelset(i2,j1,k1,1).lt.0.0_WP) alpha(2,1,1) = 0.0_WP
          if (gridLevelset(i1,j2,k1,1).lt.0.0_WP) alpha(1,2,1) = 0.0_WP
          if (gridLevelset(i2,j2,k1,1).lt.0.0_WP) alpha(2,2,1) = 0.0_WP
          if (gridLevelset(i1,j1,k2,1).lt.0.0_WP) alpha(1,1,2) = 0.0_WP
          if (gridLevelset(i2,j1,k2,1).lt.0.0_WP) alpha(2,1,2) = 0.0_WP
          if (gridLevelset(i1,j2,k2,1).lt.0.0_WP) alpha(1,2,2) = 0.0_WP
          if (gridLevelset(i2,j2,k2,1).lt.0.0_WP) alpha(2,2,2) = 0.0_WP
          buf = sum(alpha * eta)
          interpCoeff = alpha * eta / buf
       end if

       ! Interpolate the fluid variables
       Ap = interpCoeff(1,1,1) * A(i1, j1, k1) + interpCoeff(2,1,1) * A(i2, j1, k1) +        &
            interpCoeff(1,2,1) * A(i1, j2, k1) + interpCoeff(2,2,1) * A(i2, j2, k1) +        &
            interpCoeff(1,1,2) * A(i1, j1, k2) + interpCoeff(2,1,2) * A(i2, j1, k2) +        &
            interpCoeff(1,2,2) * A(i1, j2, k2) + interpCoeff(2,2,2) * A(i2, j2, k2)

    end select

  return
end subroutine ibm_inverse_distance_function


! ============================================== !
! Prepare the arrays for grid => marker exchange !
! ============================================== !
subroutine prepare_grid_arrays(densityGradient, energyGradient, normalVelocityGradient)

  ! Internal modules
  use ibm_exchange

  ! External modules
  use parallel
  use grid_levelset
  use state

  implicit none

  ! Arguments
  real(WP), dimension(nGridPoints), intent(in) :: densityGradient, energyGradient
  real(WP), dimension(nGridPoints, nDimensions), intent(in) :: normalVelocityGradient
  
  ! Local variables
  integer :: gridIndex, i, j, k, nVariables

  if (.not. useIBM) return

  nVariables = 2 * nDimensions + 5

  ! Pack data belonging to local processor in temporary array
  do k = iStart(3), iEnd(3)
     do j = iStart(2), iEnd(2)
        do i = iStart(1), iEnd(1)
           gridIndex = i - gridOffset(1) + localGridSize(1) *                                &
                (j - 1 - gridOffset(2) + localGridSize(2) *                                  &
                (k - 1 - gridOffset(3)))
           ! Pack up the levelset
           gridBuffer(i,j,k,1) = levelset(gridIndex, 1)
           ! Pack up the density
           gridBuffer(i,j,k,2) = conservedVariables(gridIndex, 1)
           ! Pack up the pressure
           gridBuffer(i,j,k,3) = pressure(gridIndex, 1)
           ! Pack up the density gradient
           gridBuffer(i,j,k,4) = densityGradient(gridIndex)
           ! Pack up the energy gradient
           gridBuffer(i,j,k,5) = energyGradient(gridIndex)
           ! Pack up the velocity
           gridBuffer(i,j,k,6:nDimensions+5) = velocity(gridIndex,1:nDimensions)
           ! Pack the velocity gradient
           gridBuffer(i,j,k,nDimensions+6:nVariables) =                                      &
                normalVelocityGradient(gridIndex, 1:nDimensions)
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
  gridLevelset(:,:,:,1) = gridBuffer(:,:,:,1)
  gridDensity(:,:,:,1) = gridBuffer(:,:,:,2)
  gridPressure(:,:,:,1) = gridBuffer(:,:,:,3)
  gridDensityGradient(:,:,:,1) = gridBuffer(:,:,:,4)
  gridEnergyGradient(:,:,:,1) = gridBuffer(:,:,:,5)
  gridVelocity(:,:,:,1:nDimensions) = gridBuffer(:,:,:,6:nDimensions+5)
  gridNormalVelocityGradient(:,:,:,1:nDimensions) = gridBuffer(:,:,:,nDimensions+6:nVariables)

  return
end subroutine prepare_grid_arrays
