module string

  implicit none
  public

  integer, parameter :: str_short  = 8
  integer, parameter :: str_medium = 64
  integer, parameter :: str_long   = 8192

   ! num2str functions
  interface num2str
     module procedure int2str
     module procedure real2str
  end interface num2str

contains

  ! Changes a string to upper case
  ! ------------------------------
  pure function str2upper (str_in) result (str_out)

    implicit none

    ! Arguments
    character(*), intent(in) :: str_in

    ! Result
    character(len(str_in))      :: str_out

    ! Local variables
    integer :: ic, i
    character(26), parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    character(26), parameter :: low = 'abcdefghijklmnopqrstuvwxyz'

    ! Capitalize each letter if it is lowecase
    str_out = str_in
    do i = 1, LEN_TRIM(str_in)
        ic = INDEX(low, str_in(i:i))
        if (ic > 0) str_out(i:i) = cap(ic:ic)
    end do

  end function str2upper

  ! Changes a string to lower case
  ! ------------------------------
  pure function str2lower (str_in) result (str_out)

    implicit none

    ! Arguments
    character(*), intent(in) :: str_in

    ! Result
    character(len(str_in)) :: str_out

    ! Local variables
    integer :: ic, i
    character(26), parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    character(26), parameter :: low = 'abcdefghijklmnopqrstuvwxyz'

    ! Capitalize each letter if it is lowecase
    str_out = str_in
    do i = 1, LEN_TRIM(str_in)
        ic = INDEX(cap, str_in(i:i))
        if (ic > 0) str_out(i:i) = low(ic:ic)
    end do

  end function str2lower

  ! Convert an integer to a string
  ! ------------------------------
  character(len=128) function int2str(num)

    implicit none

    ! Arguments
    integer :: num
    
    write(int2str, *) num
    int2str = adjustl(trim(int2str))
    
  end function int2str

  ! Convert a real to a string
  ! --------------------------
  character(len=128) function real2str(num)

    ! Externaml modules
    use precision

    implicit none

    ! Arguments
    real(WP) :: num
    
    write(real2str, *) num
    real2str = adjustl(trim(real2str))
    
  end function real2str

end module string
 
