module parser

  ! External modules
  use precision

  implicit none
  
  ! Input values storage
  integer, parameter :: tag_length = 128
  integer, parameter :: line_length = 40960
  integer :: nfields

  type, private :: t_EntryType
     character(tag_length) :: tag
     character(line_length) :: value
  end type t_EntryType

  type(t_EntryType), pointer, dimension(:), public :: entries
  
  ! Define interface for reading
  interface parser_read
     module procedure parser_readlogical
     module procedure parser_readint
     module procedure parser_readintarray
     module procedure parser_readfloat
     module procedure parser_readfloatarray
     module procedure parser_readfloatarray2D
     module procedure parser_readchar
     module procedure parser_readchararray
  end interface
  
contains
  
  ! Pack the structure -----------------------------------------------------
  subroutine parser_pack(myfield)
    implicit none
    integer,intent(in) :: myfield
    type(t_EntryType), pointer, dimension(:) :: entries_new
    
    allocate(entries_new(nfields-1))
    entries_new(1:myfield-1) = entries(1:myfield-1)
    entries_new(myfield:nfields-1) = entries(myfield+1:nfields)
    deallocate(entries)
    nullify(entries)
    entries => entries_new
    nfields = nfields-1
    
    return
  end subroutine parser_pack
  
  ! Spread the structure ---------------------------------------------------
  subroutine parser_spread
    implicit none
    type(t_EntryType), pointer, dimension(:) :: entries_new
    
    if (nfields .ne. 0) then 
       allocate(entries_new(nfields+1))
       entries_new(1:nfields) = entries(1:nfields)
       deallocate(entries)
       nullify(entries)
       entries => entries_new
    else
       allocate(entries(1))
    end if
    nfields = nfields+1
    
    return
  end subroutine parser_spread
  
  ! Add a new entry in the structure ----------------------------------------
  subroutine parser_newentry(mytag,myvalue)
    implicit none
    character(*),intent(in) :: mytag
    character(*),intent(in) :: myvalue
    integer :: ifield
    logical :: isdef 
    
    call parser_fieldfortag(mytag,ifield,isdef)
    if (.not. isdef) then
       call parser_spread()  
       entries(nfields)%tag=mytag
       entries(nfields)%value=myvalue
    else
       entries(ifield)%value=myvalue
    end if
    
    return
  end subroutine parser_newentry
  
  ! Get field number from tag ----------------------------------------------
  subroutine parser_fieldfortag(mytag,myfield,isdef)
    implicit none
    character(*),intent(in) :: mytag
    integer,intent(out) :: myfield
    logical,optional,intent(out) :: isdef
    integer :: ifield
    
    isdef = .false.
    do ifield=1,nfields
       if (entries(ifield)%tag==mytag) then
          myfield=ifield
          isdef = .true.
          return
       end if
    end do
    
    return
  end subroutine parser_fieldfortag
  
  ! Check whether the field is defined -------------------------------------
  subroutine parser_is_defined(mytag,isdef)
    implicit none
    
    character(*),intent(in) :: mytag
    logical,intent(out) :: isdef
    integer :: ifield
    
    call parser_fieldfortag(mytag,ifield,isdef)
    
    return
  end subroutine parser_is_defined
  
  ! Read logicals ---------------------------------------------
  subroutine parser_readlogical(mytag,value,default)
    implicit none
    
    character(*),intent(in) :: mytag
    logical,intent(out) :: value
    logical,optional,intent(in) :: default
    integer :: ifield
    logical :: isdef
    integer :: conv
    
    call parser_fieldfortag(mytag,ifield,isdef)
    if (.not.isdef .AND. present(default)) then   
       value = default
    else if (.not.isdef .AND. .not.present(default)) then
       print*,'Parser : '// mytag //' not defined'
       stop
    else
       if (len_trim(entries(ifield)%value).eq.1) then          
          read(entries(ifield)%value,*) conv
          value = (conv.eq.1)
       else
          read(entries(ifield)%value,*) value
       end if
    end if
    
    return
  end subroutine parser_readlogical
  
  ! Read integers ---------------------------------------------
  subroutine parser_readint(mytag,value,default)
    implicit none
    
    character(*),intent(in) :: mytag
    integer,intent(out) :: value
    integer,optional,intent(in) :: default
    integer :: ifield
    logical :: isdef

    call parser_fieldfortag(mytag,ifield,isdef)
    if (.not.isdef .AND. present(default)) then   
       value = default
    else if (.not.isdef .AND. .not.present(default)) then
       print*,'Parser : '// mytag //' not defined'
       stop
    else
       read(entries(ifield)%value,*) value
    end if
    
    return
  end subroutine parser_readint
  
  ! Read floats ---------------------------------------------
  subroutine parser_readfloat(mytag,value,default)
    implicit none
    
    character(*),intent(in) :: mytag
    real(WP),intent(out) :: value
    real(WP),optional,intent(in) :: default
    integer :: ifield
    logical :: isdef
    
    call parser_fieldfortag(mytag,ifield,isdef)
    if (.not.isdef .AND. present(default)) then   
       value = default
    else if (.not.isdef .AND. .not.present(default)) then
       print*,'Parser : '// mytag //' not defined'
       stop
    else
       read(entries(ifield)%value,*) value
    end if
    
    return
  end subroutine parser_readfloat

  ! Read characters ---------------------------------------------
  subroutine parser_readchar(mytag,value,default)
    implicit none
    
    character(*),intent(in) :: mytag
    character(len=*),intent(out) :: value
    character(len=*),optional,intent(in) :: default
    integer :: ifield
    logical :: isdef
    
    call parser_fieldfortag(mytag,ifield,isdef)
    if (.not.isdef .AND. present(default)) then   
       value = default
    else if (.not.isdef .AND. .not.present(default)) then
       print*,'Parser : '// mytag //' not defined'
       stop
    else
       read(entries(ifield)%value,'(a)') value
    end if
    
    return
  end subroutine parser_readchar
  
  ! Count size of the arrays ----------------------------------------
  subroutine parser_getsize(mytag,numb)
    implicit none
    
    character(*),intent(in) :: mytag
    integer,intent(out) :: numb
    integer :: ifield
    logical :: isdef
    integer :: i
    integer, dimension(line_length) :: counter
    
    ! Read it now                                                                                    
    call parser_fieldfortag(mytag,ifield,isdef)
    if (.not. isdef) then
       print*,'Parser : '// mytag //' not defined'
       stop
    end if
    
    ! Count the number of entries                                   
    counter = 0
    do i=1,len_trim(entries(ifield)%value)
       if (entries(ifield)%value(i:i).EQ.' ') counter(i)=1
    end do
    do i=1+1,len_trim(entries(ifield)%value)
       if (counter(i).EQ.1 .AND. counter(i-1).EQ.1) counter(i-1)=0
    end do
    numb = sum(counter)+1
    
    return
  end subroutine parser_getsize
  
  ! Read integer arrays ---------------------------------------------
  subroutine parser_readintarray(mytag,value)
    implicit none
    
    character(*),intent(in) :: mytag
    integer,dimension(:),intent(out) :: value
    integer :: ifield
    logical :: isdef
    
    ! Read them
    call parser_fieldfortag(mytag,ifield,isdef)
    if (.not. isdef) then
       print*,'Parser : '// mytag //' not defined'
       stop
    end if
    read(entries(ifield)%value,*) value
    
    return
  end subroutine parser_readintarray

  ! Read float arrays ---------------------------------------------
  subroutine parser_readfloatarray(mytag,value)
    implicit none

    character(*),intent(in) :: mytag
    real(WP),dimension(:),intent(out) :: value
    integer :: ifield
    logical :: isdef
    
    ! Read them
    call parser_fieldfortag(mytag,ifield,isdef)
    if (.not. isdef) then
       print*,'Parser : '// mytag //' not defined'
       stop
    end if         
    read(entries(ifield)%value,*) value

    return
  end subroutine parser_readfloatarray

  ! Read float arrays ---------------------------------------------
  subroutine parser_readfloatarray2D(mytag,value)
    implicit none
    
    character(*),intent(in) :: mytag
    real(WP),dimension(:,:),intent(out) :: value
    integer :: ifield
    logical :: isdef
    
    ! Read them
    call parser_fieldfortag(mytag,ifield,isdef)
    if (.not. isdef) then
       print*,'Parser : '// mytag //' not defined'
       stop
    end if
    read(entries(ifield)%value,*) value
    
    return
  end subroutine parser_readfloatarray2D
  
  ! Read character arrays -----------------------------------------
  subroutine parser_readchararray(mytag,value)
    implicit none

    character(*),intent(in) :: mytag
    character(*),dimension(:),intent(out) :: value
    integer :: ifield
    logical :: isdef
    
    ! Read them
    call parser_fieldfortag(mytag,ifield,isdef)
    if (.not. isdef) then
       print*,'Parser : '// mytag //' not defined'
       stop
    end if                                                                   
    read(entries(ifield)%value,*) value

    return
  end subroutine parser_readchararray
  
end module parser

! Initialize the parser ---------------------------------------------
subroutine parser_init

  ! Internal modules
  use parser

  ! External modules
  use fileio

  implicit none
  
  nfields = 0
  if (associated(entries)) then
     deallocate(entries)
     nullify(entries)
  end if
  
  return
end subroutine parser_init

! Read & parse the input file ---------------------------------------------
subroutine parser_parsefile(input)

  ! Internal modules
  use parser

  ! External modules
  use fileio
  
  implicit none
  integer :: iunit, ierror, limiter, nlines, i, j, ntags, comment
  integer, dimension(:), allocatable :: limit,line
  character(len=line_length) :: buffer
  character(len=line_length), dimension(:), allocatable :: file
  character(len=line_length) :: value
  character(len=tag_length) :: tag
  character(len=*) :: input
  
  ! Open the file
  ierror = 0
  iunit = iopen()
  open (iunit, file = input, form='formatted', status = 'old', iostat = ierror)
  if (ierror .ne. 0) stop 'Parser : unable to open the input file.'
  
  ! Count the number of lines in the file
  ierror = 0
  nlines = 0
  do while (ierror .eq. 0)
     read(iunit,'(a)',iostat = ierror) buffer
     nlines = nlines + 1
  end do
  rewind(iunit)
  
  ! Allocate to the right size
  allocate(file(nlines+1), limit(nlines+1), line(nlines+1))
  
  ! Read everything in the buffer
  ierror = 0
  nlines = 0
  loop: do while (ierror .eq. 0)
     read(iunit,'(a)',iostat=ierror) buffer
     if (ierror .ne. 0) exit loop
     ! Remove the tabs
     do j=1,line_length
        if (ichar(buffer(j:j)).EQ.9) buffer(j:j)=' '
     end do
     ! Find comments
     comment = scan(buffer,'!#%')
     ! Remove them
     if (comment.NE.0) buffer(comment:) = ''
     ! Trim
     buffer = adjustl(buffer)
     ! Add line
     if (len_trim(buffer).NE.0) then
        nlines = nlines + 1
        file(nlines) = buffer
     end if
  end do loop
  
  ! Close de file
  close(iunit)
  ierror = iclose(iunit)
  
  ! Get the tags
  ntags = 0
  do i=1,nlines
     limiter = index(file(i),':')
     if (limiter.NE.0) then
        ntags = ntags + 1
        line(ntags) = i
        line(ntags+1) = nlines+1
        limit(ntags) = limiter
     end if
  end do
  
  ! Read everything now
  do i=1,ntags
     buffer = ''
     do j=line(i),line(i+1)-1
        if (j==line(i)) then
           buffer = trim(buffer) // trim(file(j))
        else
           buffer = trim(buffer) // ' ' // trim(file(j))
        end if
     end do
     read(buffer(1:limit(i)-1),'(a)') tag
     read(buffer(limit(i)+1:),'(a)') value
     if (len_trim(value).NE.0) then
        value = adjustl(value)
        call parser_newentry(tag,value)
     end if
  end do
  
  return
end subroutine parser_parsefile
