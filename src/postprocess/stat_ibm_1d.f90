module stat_ibm_1d

  ! External modules
  use stat
  use string
  
  implicit none

  ! Global variables
  integer :: nStatPoints
  real(WP) :: xmin, xmax, dx
  real(WP), dimension(:), allocatable :: grid1D, TKE_OLD, TKE_NEW
 
end module stat_ibm_1d


! =========================!
! Initialize 1D statistics !
! =========================!
subroutine stat_ibm_1d_setup

  ! Internal modules
  use stat_ibm_1d

  ! External modules
  use parser
  use parallel
  use simulation_flags
  use geometry
  use grid

  implicit none

  ! Local variables
  integer :: i
  
  ! Return if IBM is not used
  if (.not. useIBM) return

  ! Get stat size
  call parser_read('stat ibm bin size', nStatPoints, globalGridSize(1))

  ! Allocate arrays
  allocate(grid1D(nStatPoints))
  allocate(TKE_OLD(nStatPoints)); TKE_OLD = 0.0_WP
  allocate(TKE_NEW(nStatPoints)); TKE_NEW = 0.0_WP

  ! Generate the 1D grid
  xmin = minval(coordinates(:,1))
  call parallel_min(xmin)
  xmax = maxval(coordinates(:,1))
  call parallel_max(xmax)
  dx = (xmax - xmin) / real(nStatPoints-1, WP)
  do i = 1, nStatPoints
     grid1D(i) = xmin + real(i-1,WP) * dx
  end do
 
  return
end subroutine stat_ibm_1d_setup

! ======================!
! Cleanup 1D statistics !
! ======================!
subroutine stat_ibm_1d_cleanup

  ! Internal modules
  use stat_ibm_1d
  
  implicit none

  if (allocated(grid1D)) deallocate(grid1D)
  if (allocated(TKE_OLD)) deallocate(TKE_OLD)
  if (allocated(TKE_NEW)) deallocate(TKE_NEW)

  return
end subroutine stat_ibm_1d_cleanup

! ===================== !
! Compute 1D statistics !
! ===================== !
subroutine stat_ibm_1d_compute

  ! Internal modules
  use stat_ibm_1d

  ! External modules
  use parallel
  use simulation_flags
  use geometry
  use ibm
  use grid
  use state
  use grid_functions
  use grid_levelset
  use math

  implicit none

  ! Local variables
  integer :: i, j
  real(WP), dimension(nStatPoints, nDimensions) :: favreVelocity
  real(WP), dimension(nStatPoints) :: meanRho, totalVolume

  ! Return if IBM is not used
  if (.not. useIBM) return

  ! Update old value
  TKE_OLD = TKE_NEW

  ! Compute mean stats
  meanRho = 0.0_WP
  favreVelocity = 0.0_WP
  totalVolume = 0.0_WP
  do i = 1, nGridPoints
     j = 1 + nint((coordinates(i,1) - xmin) / dx)
     totalVolume(j) = totalVolume(j) + primitiveGridNorm(i,1)
     meanRho(j) = meanRho(j) + gridNorm(i,1) * conservedVariables(i,1)
     favreVelocity(j,1:nDimensions) = favreVelocity(j,1:nDimensions) + gridNorm(i, 1) *      &
          conservedVariables(i,2:nDimensions+1)
  end do
  call parallel_sum(meanRho)
  call parallel_sum(totalVolume)
  do i = 1, nDimensions
     call parallel_sum(favreVelocity(:,i))
     favreVelocity(:,i) = favreVelocity(:,i) / meanRho
  end do

  ! Update new <alpha> <rho> K
  TKE_NEW = 0.0_WP
  do i = 1, nGridPoints
     j = 1 + nint((coordinates(i,1) - xmin) / dx)
     TKE_NEW(j) = TKE_NEW(j) + 0.5_WP * conservedVariables(i, 1) * sum(                      &
          velocity(i,1:nDimensions) - favreVelocity(j,1:nDimensions))**2 * gridNorm(i, 1)
  end do
  call parallel_sum(TKE_NEW)
  TKE_NEW = TKE_NEW / totalVolume

  return
end subroutine stat_ibm_1d_compute


! =================== !
! Read ibm statistics !
! =================== !
subroutine stat_ibm_1d_read

  ! Internal modules
  use stat_ibm_1d

  ! External modules
  use parallel
  use ibm

  implicit none

  ! Nothing to do
  
  return
end subroutine stat_ibm_1d_read


! ================== !
! Read 1D statistics !
! ================== !
subroutine stat_ibm_1d_write

  ! Internal modules
  use stat_ibm_1d

  ! External modules
  use time_info
  use parallel
  use fileio
  use math
  use geometry
  use ibm
  use grid
  use state
  use grid_functions
  use grid_levelset
  use first_derivative

  implicit none

  ! Local variables
  character(len=str_medium) :: filename
  integer  :: i, j, n, iunit, ierror
  real(WP) :: D, F(3), pPrime
  real(WP), dimension(nDimensions) :: uPrime, uDoublePrime, objectVelocity
  real(WP), dimension(nDimensions**2) :: sigmaPrime
  real(WP), dimension(nStatPoints, nDimensions) :: meanVelocity, favreVelocity
    real(WP), dimension(nStatPoints, nDimensions**2) :: meanSigma
  real(WP), dimension(nStatPoints) :: meanRho, meanUU, meanVV, meanWW, meanUV, meanUW,       &
       meanVW, alpha, volume, totalVolume, rhoUUU, pU, uSigma, meanP,    &
       dragP, dragV, meanUprime, meandUdx, pDivU, sigmaDivU, dkdt
  real(WP), dimension(nGridPoints, nDimensions) :: gradI

  ! Return if IBM is not used
  if (.not. useIBM) return

  ! Compute gradient of the indicator function
  call gradient(indicatorFunction, gradI)
  
  ! Compute mean stats
  meanRho = 0.0_WP
  meanVelocity = 0.0_WP
  favreVelocity = 0.0_WP
  totalVolume = 0.0_WP
  volume = 0.0_WP
  meanSigma = 0.0_WP
  meanP = 0.0_WP
  meandUdx = 0.0_WP
  do i = 1, nGridPoints
     j = 1 + nint((coordinates(i,1) - xmin) / dx)

     ! Sum up grid volumes
     volume(j) = volume(j) + gridNorm(i,1)
     totalVolume(j) = totalVolume(j) + primitiveGridNorm(i,1)

     ! Mean stats
     meanRho(j) = meanRho(j) + gridNorm(i,1) * conservedVariables(i,1)
     meanVelocity(j,1:nDimensions) = meanVelocity(j,1:nDimensions) +                         &
          gridNorm(i,1) * velocity(i,1:nDimensions)
     favreVelocity(j,1:nDimensions) = favreVelocity(j,1:nDimensions) +                       &
          gridNorm(i,1) * conservedVariables(i,2:nDimensions+1)
     meanP(j) = meanP(j) + gridNorm(i,1) * pressure(i,1)
     meanSigma(j,:) = meanSigma(j,:) + gridNorm(i,1) * stressTensor(i,:)
     meandUdx(j) = meandUdx(j) + gridNorm(i, 1) * velocityGradient(i,1)
  end do

  ! Sum them over procs
  call parallel_sum(volume)
  call parallel_sum(totalVolume)
  call parallel_sum(meanRho)
  call parallel_sum(meanP)
  call parallel_sum(meandUdx)
  do i = 1, nDimensions
     call parallel_sum(meanVelocity(:,i))
     call parallel_sum(favreVelocity(:,i))
  end do
  do i = 1, nDimensions**2
     call parallel_sum(meanSigma(:,i))
  end do

  ! Normalize
  do i = 1, nStatPoints
     alpha(i) = volume(i) / totalVolume(i)
     meanRho(i) = meanRho(i) / volume(i)
     meanP(i) = meanP(i) / volume(i)
     meandUdx(i) = meandUdx(i) / volume(i)
     meanVelocity(i,:) = meanVelocity(i,:) / volume(i)
     favreVelocity(i,:) = favreVelocity(i,:) / volume(i) / meanRho(i)
     meanSigma(i,:) = meanSigma(i,:) / volume(i)
  end do
  
  ! Compute variance
  meanUU = 0.0_WP; meanVV = 0.0_WP; meanWW = 0.0_WP
  meanUV = 0.0_WP; meanUW = 0.0_WP; meanVW = 0.0_WP
  rhoUUU = 0.0_WP; pU = 0.0_WP; uSigma = 0.0_WP; meanUprime = 0.0_WP
  meanUPrime = 0.0_WP; sigmaDivU = 0.0_WP
  pDivU = 0.0_WP
  dragP = 0.0_WP; dragV = 0.0_WP
  do i = 1, nGridPoints
     j = 1 + nint((coordinates(i,1) - xmin) / dx)

     ! Get Reynolds & Favre fluctuations
     uPrime = velocity(i,1:nDimensions) - meanVelocity(j,1:nDimensions)
     uDoublePrime = velocity(i,1:nDimensions) - favreVelocity(j,1:nDimensions)
     pPrime = pressure(i,1) - meanP(j)
     sigmaPrime = stressTensor(i,:) - meanSigma(j,:)

     ! Reynolds stresses
     meanUU(j) = meanUU(j) + gridNorm(i,1) * conservedVariables(i,1) * uDoublePrime(1)**2
     if (nDimensions .gt. 1) then
        meanVV(j) = meanVV(j) + gridNorm(i,1) * conservedVariables(i,1) * uDoublePrime(2)**2
        meanUV(j) = meanUV(j) + gridNorm(i,1) * conservedVariables(i,1) * uDoublePrime(1) *  &
             uDoublePrime(2)
     end if
     if (nDimensions .gt. 2) then
        meanWW(j) = meanWW(j) + gridNorm(i,1) * conservedVariables(i,1) * uDoublePrime(3)**2
        meanUW(j) = meanUW(j) + gridNorm(i,1) * conservedVariables(i,1) * uDoublePrime(1) *  &
             uDoublePrime(3)
        meanVW(j) = meanVW(j) + gridNorm(i,1) * conservedVariables(i,1) * uDoublePrime(2) *  &
             uDoublePrime(3)
     end if

     ! Other terms
     meanUprime(j) = meanUprime(j) + uDoublePrime(1) * gridNorm(i,1)
     pU(j) = pU(j) + pPrime * uPrime(1) * gridNorm(i,1)
     uSigma(j) = uSigma(j) + sum(uPrime * sigmaPrime(1:nDimensions)) * gridNorm(i,1)
     rhoUUU(j) = rhoUUU(j) + conservedVariables(i,1) * sum(uDoublePrime) * uDoublePrime(1) * &
          gridNorm(i,1)
     pDivU(j) = pDivU(j) + pPrime * (velocityGradient(i,1) - meandUdx(j)) * gridNorm(i,1)
     sigmaDivU(j) = sigmaDivU(j) + (sum(sigmaPrime * velocityGradient(i,:)) -  &
          sigmaPrime(1) * meanDuDx(j)) * gridNorm(i,1)
     !sigmaDivU(l) = sigmaDivU(l) + ((sigmaPrime(1) *                                &
     !     velocityGradient(i,1)) - sigmaPrime(1) * meanDuDx(l)) *              &
     !     gridNorm(i,1) 
     select case (nDimensions)
     case (2)
        pDivU(j) = pDivU(j) + pPrime * gridNorm(i,1) * velocityGradient(i, 4)
        !sigmaDivU(l) = sigmaDivU(l) + (2.0_WP*(sigmaPrime(2) *                         &
        !  velocityGradient(i,2)) - sigmaPrime(1) * meanDuDx(l)) *            &
        !  gridNorm(i,1)

     case (3)
        pDivU(j) = pDivU(j) + pPrime * gridNorm(i,1) *                                       &
             (velocityGradient(i, 5) + velocityGradient(i, 9))
     end select

     ! Get velocity of associated object
     if (ibm_move) then
        n = objectIndex(i)
        objectVelocity = object(n)%velocity(1:nDimensions)
     else
        objectVelocity = 0.0_WP
     end if

     ! Drag production
     dragP(j) = dragP(j) + uDoublePrime(1) * pressure(i,1) * gradI(i,1) * primitiveGridNorm(i,1)
     dragV(j) = dragV(j) - uDoublePrime(1) * sum(stressTensor(i,1:nDimensions) * gradI(i,:)) * primitiveGridNorm(i,1)

  end do
  call parallel_sum(meanUU); meanUU = meanUU / meanRho / volume
  call parallel_sum(meanVV); meanVV = meanVV / meanRho / volume
  call parallel_sum(meanWW); meanWW = meanWW / meanRho / volume
  call parallel_sum(meanUV); meanUV = meanUV / meanRho / volume
  call parallel_sum(meanUW); meanUW = meanUW / meanRho / volume
  call parallel_sum(meanVW); meanVW = meanVW / meanRho / volume
  call parallel_sum(meanUprime); meanUprime = meanUprime / volume
  call parallel_sum(pU); pU = pU / volume
  call parallel_sum(uSigma); uSigma = uSigma / volume
  call parallel_sum(rhoUUU); rhoUUU = rhoUUU / volume
  call parallel_sum(pDivU); pDivU = pDivU / volume
  call parallel_sum(sigmaDivU); sigmaDivU = sigmaDivU / volume
  call parallel_sum(dragP); dragP = dragP / totalVolume
  call parallel_sum(dragV); dragV = dragV / totalVolume

!!$  ! Drag production
!!$  call ibm_integrate_forces
!!$  dragP = 0.0_WP; dragV = 0.0_WP
!!$  do n = 1, nObjects
!!$     procHasObject = .true.
!!$     do i = 1, nDimensions
!!$        if (object(n)%position(i).ge.minval(coordinates(:,i)) .and.                          &
!!$             object(n)%position(i).le.maxval(coordinates(:,i)) .and. procHasObject) then
!!$        else
!!$           procHasObject = .false.
!!$        end if
!!$     end do
!!$     if (procHasObject) then
!!$        j = 1 + nint((object(n)%position(1) - xmin) / dx)
!!$        dragP(j) = dragP(j) - (object(n)%velocity(1)-favreVelocity(j,1)) * object(n)%pForce(1)
!!$        dragV(j) = dragV(j) - (object(n)%velocity(1)-favreVelocity(j,1)) * object(n)%vForce(1)
!!$     end if
!!$  end do
!!$  call parallel_sum(dragP); dragP = dragP / totalVolume
!!$  call parallel_sum(dragV); dragV = dragV / totalVolume

  ! Root process writes
  if (iRank .eq. iRoot) then

     ! ------- WRITE THE ASCI FILE ------
     filename = "stat/stat-ibm"
     write(filename, '(2A,I8.8,A)') trim(filename), "-", timestep, ".txt"
     iunit = iopen()
     open (iunit, file=adjustl(trim(filename)), form="formatted",iostat=ierror)
     write(iunit,'(A,1ES20.12)') "t=",time
     write(iunit,'(5a20)') 'X','D','Fx','Fy','Fz'
     F = 0.0_WP
     do n = 1, nObjects
        D = (2.0_WP * real(nDimensions, WP) * object(n)%volume / pi)                         &
             ** (1.0_WP / real(nDimensions, WP))
        F(1:nDimensions) =  object(n)%pForce(1:nDimensions) +  object(n)%vForce(1:nDimensions)
        write(iunit,'(10000ES20.12)') object(n)%position(1), D, F(1), F(2), F(3)
     end do
     close(iclose(iunit)) 
     dkdt = (TKE_NEW - TKE_OLD) / timeStepSize
     filename= "stat/stat-tke"
     write(filename, '(2A,I8.8,A)') trim(filename), "-", timestep, ".txt"
     iunit = iopen()
     open (iunit, file=adjustl(trim(filename)), form="formatted",iostat=ierror)
     write(iunit,'(A,1ES21.12)') "t=",time
     write(iunit,'(22a23)') "X", "alpha", "Volume", "<U>", "<V>",  "<U'U'>", "<V'V'>",       & ! '<W>',
           "<U'V'>", "dragP", "dragV", "<rho>", "<rhoUUU>",                                  & ! "<W'W'>", "<U'W'>", "<V'W'>",
           "pU", "uSigma", "pdivU", "sigmaDivU", "<P>", "<Sigma>","<u''>","<dUdx>",          &
           "Utilde", "Vtilde", "dkdt"
     do i = 1, nStatPoints
        write(iunit,'(10000ES21.12E3)') grid1D(i), alpha(i), totalVolume(i),                 &
             meanVelocity(i,1), meanVelocity(i,2), meanUU(i), meanVV(i), meanUV(i),          & ! meanVelocity(i,3),  meanWW(i),  meanUW(i),   
              dragP(i), dragV(i), meanRho(i), rhoUUU(i), pU(i), uSigma(i), pDivU(i),         & !meanVW(i),
              sigmaDivU(i), meanP(i), meanSigma(i,1), meanUPrime(i), meandUdx(i),            &
              favreVelocity(i,1), favreVelocity(i,2), dkdt(i)
     end do
     close(iclose(iunit))
     
  end if

  ! Log
  call monitor_log("IBM STATISTICS FILE WRITTEN")
  
  return
end subroutine stat_ibm_1d_write
