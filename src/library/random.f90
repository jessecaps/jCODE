module random

  ! External modules
  use precision

  implicit none
  public
  
contains
  
  ! ------------------------------------------------------------------ !
  ! Normal distribution                                                !
  ! Adapted from the following Fortran 77 code                         !
  !      ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.                 !
  !      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE, !
  !      VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.                  !
  !  The function random_normal() returns a normally distributed       !
  !  pseudo-random number with zero mean and unit variance.            !
  !  The algorithm uses the ratio of uniforms method of A.J. Kinderman !
  !  and J.F. Monahan augmented with quadratic bounding curves.        !
  ! ------------------------------------------------------------------ !
  real(WP) function random_normal(m, sd)

    implicit none

    ! Arguments
    real(WP), intent(in), optional :: m
    real(WP), intent(in), optional :: sd

    ! Local variables
    real(WP) :: s = 0.449871_WP
    real(WP) :: t = -0.386595_WP
    real(WP) :: a = 0.19600_WP
    real(WP) :: b = 0.25472_WP
    real(WP) :: r1 = 0.27597_WP
    real(WP) :: r2 = 0.27846_WP
    real(WP) :: u,v,x,y,q
    
    ! Generate P = (u,v) uniform in rectangle enclosing acceptance region
    do    
       call random_number(u)
       call random_number(v)
       v=1.7156_WP*(v-0.5_WP)
       ! Evaluate the quadratic form
       x=u-s
       y=abs(v)-t
       q=x**2+y*(a*y-b*x)
       ! Accept P if inside inner ellipse
       if (q<r1) exit
       ! Reject P if outside outer ellipse
       if (q>r2) cycle
       ! Reject P if outside acceptance region
       if (v**2<-4.0_WP*log(u)*u**2) cycle
    end do
    
    ! Return ratio of P's coordinates as the normal deviate
    random_normal = v/u
    
    ! Modify to give correct mean and standard deviation
    if (present(sd)) random_normal = random_normal*sd
    if (present(m))  random_normal = random_normal+m
    
    return
  end function random_normal
  
  ! ---------------------------------------------------------------------------- !
  ! Log-normal distribution                                                      !
  ! If X has a lognormal distribution, then log(X) is normally distributed.      !
  ! Here the logarithm is the natural logarithm, that is to base e, sometimes    !
  ! denoted as ln.  To generate random variates from this distribution, generate !
  ! a random deviate from the normal distribution with mean and variance equal   !
  ! to the mean and variance of the logarithms of X, then take its exponential.  !
  !                                                                              !
  ! Relationship between the mean & variance of log(X) and the mean & variance   !
  ! of X, when X has a lognormal distribution.                                   !
  ! Let m = mean of log(X), and s^2 = variance of log(X)                         !
  ! Then                                                                         !
  ! mean of X     = exp(m + 0.5s^2)                                              !
  ! variance of X = (mean(X))^2.[exp(s^2) - 1]                                   !
  !                                                                              !
  ! In the reverse direction                                                     !
  ! variance of log(X) = log[1 + var(X)/(mean(X))^2]                             !
  ! mean of log(X)     = log(mean(X) - 0.5var(log(X))                            !
  ! ---------------------------------------------------------------------------- !
  real(WP) function random_lognormal(m, sd)
    implicit none

    ! Arguments
    real(WP), intent(in) :: m
    real(WP), intent(in) :: sd

    ! Local variables
    real(WP) :: x,mlog,sdlog
    
    sdlog = sqrt(log(1.0_WP+(sd/m)**2))
    mlog  = log(m)-0.5_WP*sdlog**2
    x     = random_normal(mlog,sdlog)
    random_lognormal = exp(x)
    
    return
  end function random_lognormal
    
end module random


! ------------------------- !
! Initialization of the RNG !
! Seeded based on parallel  !
! partitioning              !
! ------------------------- !
subroutine random_setup

  ! Internal modules
  use random

  ! External modules
  use parallel

  implicit none

  ! Local variables
  integer :: i, n, clock
  integer, dimension(:), allocatable :: seed

  call random_seed(size = n)
  allocate(seed(n))
  call system_clock(count=clock)
  seed = clock + 37 * (/ (i - 1, i = 1, n) /)
  seed = seed * (irank + 2)
  call random_seed(put = seed)
  deallocate(seed)
  
  return
end subroutine random_setup
