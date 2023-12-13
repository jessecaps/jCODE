module precision

  implicit none
  public

  integer, parameter  :: SP = kind(1.0)
  integer, parameter, private :: DP = kind(1.0d0)
  integer, parameter  :: WP = DP
  real(WP), private   :: sample_real_at_WP
  real(WP), parameter :: MAX_REAL_WP = HUGE(sample_real_at_WP)
  integer, private    :: sample_int
  integer, parameter  :: MAX_INTEGER = HUGE(sample_int)

end module precision
