module infnan
  implicit none

  private
  public :: isnan, isinf, isposinf, isneginf

  ! Kind numbers for single and double precision integer containers
  integer, parameter :: Single = selected_int_kind(precision(1.e0))
  integer, parameter :: Double = selected_int_kind(precision(1.d0))

  ! Single precision IEEE values
  ! Leading bit for sNegInf should be 1. But does not compile!
  integer(Single), parameter :: sNaN    = int(Z"7FC00000")
  integer(Single), parameter :: sPosInf = int(Z"7F800000")
  integer(Single), parameter :: sNegInf = int(Z"7F800000")

  ! Double precision IEEE values
  ! Leading bit for dNegInf should be 1. But does not compile!
  integer(Double), parameter :: dNaN    = int(Z"7FF8000000000000")
  integer(Double), parameter :: dPosInf = int(Z"7FF0000000000000")
  integer(Double), parameter :: dNegInf = int(Z"7FF0000000000000")

  ! Location of single and double precision sign bit (Intel)
  ! Subtract one because bit numbering starts at zero
  integer, parameter :: SPSB = bit_size(sNaN) - 1
  integer, parameter :: DPSB = bit_size(dNaN) - 1

  interface isnan
     module procedure sisnan
     module procedure disnan
  end interface isnan

  interface isinf
     module procedure sisinf
     module procedure disinf
  end interface isinf

  interface isposinf
     module procedure sisposinf
     module procedure disposinf
  end interface isposinf

  interface isneginf
     module procedure sisneginf
     module procedure disneginf
  end interface isneginf

contains    

  ! Single precision test for NaN
  elemental function sisnan(x) result(res)
    real(kind(1.e0)), intent(in) :: x
    logical :: res
    res = ieor(ibclr(transfer(x,sNaN),SPSB), sNaN) == 0
  end function sisnan

  ! Double precision test for NaN
  elemental function disnan(d) result(res)
    real(kind(1.d0)), intent(in) :: d
    logical :: res
    res = ieor(ibclr(transfer(d,dNaN),DPSB), dNaN) == 0
  end function disnan

  ! Single precision test for Inf
  elemental function sisinf(x) result(res)
    real(kind(1.e0)), intent(in) :: x
    logical :: res
    res = ieor(ibclr(transfer(x,sPosInf),SPSB), sPosInf) == 0
  end function sisinf

  ! Double precision test for Inf
  elemental function disinf(d) result(res)
    real(kind(1.d0)), intent(in) :: d
    logical :: res
    res = ieor(ibclr(transfer(d,dPosInf),DPSB), dPosInf) == 0
  end function disinf

  ! Single precision test for +Inf
  elemental function sisposinf(x) result(res)
    real(kind(1.e0)), intent(in) :: x
    logical :: res
    res = ieor(transfer(x,sPosInf), sPosInf) == 0
  end function sisposinf

  ! Double precision test for +Inf
  elemental function disposinf(d) result(res)
    real(kind(1.d0)), intent(in) :: d
    logical :: res
    res = ieor(transfer(d,dPosInf), dPosInf) == 0
  end function disposinf

  ! Single precision test for -Inf
  elemental function sisneginf(x) result(res)
    real(kind(1.e0)), intent(in) :: x
    logical :: res
    integer(Single) ix
    ix = transfer(x,sNegInf)
    res = ((ieor(ibclr(ix,SPSB), sNegInf) == 0) .and. btest(ix,SPSB))
  end function sisneginf

  ! Double precision test for -Inf
  elemental function disneginf(d) result(res)
    real(kind(1.d0)), intent(in) :: d
    logical :: res
    integer(Double) id
    id = transfer(d,dNegInf)
    res = ((ieor(ibclr(id,DPSB), dNegInf) == 0) .and. btest(id,DPSB))
  end function disneginf

end module infnan
