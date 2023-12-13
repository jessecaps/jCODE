program stl2levelset

  ! External modules
  use precision
  use string
  use parser
  use parallel
  use math
  use compgeom
  use simulation_flags
  use solver_options
  use geometry
  use grid
  use grid_levelset

  implicit none

  ! General
  integer :: ierr,iunit,i,n
  real(WP), parameter :: eps=1.0e-9_WP

  ! STL file
  character(len=str_medium) :: stlFile
  character(len=80) :: cbuff
  character(len= 2) :: padding
  integer  :: nt,nt2,nn
  real(SP) :: fltbuf
  real(WP) :: tmp
  type t_Triangle
     real(WP), dimension(3) :: norm
     real(WP), dimension(3) :: v1
     real(WP), dimension(3) :: v2
     real(WP), dimension(3) :: v3
  end type t_Triangle
  type(t_Triangle), dimension(:), allocatable :: t,t2
  real(WP) :: xmin,xmax,ymin,ymax,zmin,zmax,Lx,Ly,Lz
  real(WP), dimension(3) :: shift,scaling
  logical :: swap_xy,swap_yz,swap_zx
  logical :: swap_yx,swap_zy,swap_xz
  logical :: limit_extents,mirror_x,mirror_y,mirror_z

  ! Levelset
  integer :: count
  character(len = str_medium) :: gridFile
  real(WP) :: mydist,newdist
  real(WP), dimension(3) :: c,mynorm,myproj,newproj

  ! Input processing
  character(len=str_medium) :: input

  ! Progress monitoring
  integer :: togo,prog
  integer :: iratio_new,iratio_old

  ! Levelset file
  character(len=str_medium) :: filename

  ! Initialize parallel environment and parse the input file
  call parallel_init
  call parallel_get_inputname(input)
  call parser_init
  call parser_parsefile(input)
  disableManualDecomp = .true.

  ! Read information from input file

  if (iRank .eq. iRoot) then
     write (*,*)
     write (*,*) '======================================'
     write (*,*) '| jCODE - STL to levelset converter  |'
     write (*,*) '======================================'
     write (*,*)
  end if

  ! Read stl file
  call parser_readchar('STL file',stlFile)

  ! Open the STL file to read
  call BINARY_FILE_OPEN(iunit,trim(stlFile),"r",ierr)

  ! Read sizes
  call BINARY_FILE_READ(iunit,cbuff,80,kind(cbuff),ierr)
  call BINARY_FILE_READ(iunit,nt,1,kind(nt),ierr)
  print*,'STL name :',trim(cbuff)
  print*,'# of triangles :',nt

  ! Read triangles
  allocate(t(nt))
  do n=1,nt
     call BINARY_FILE_READ(iunit,fltbuf,1,kind(fltbuf),ierr); t(n)%norm(1)=real(fltbuf,WP)
     call BINARY_FILE_READ(iunit,fltbuf,1,kind(fltbuf),ierr); t(n)%norm(2)=real(fltbuf,WP)
     call BINARY_FILE_READ(iunit,fltbuf,1,kind(fltbuf),ierr); t(n)%norm(3)=real(fltbuf,WP)

     call BINARY_FILE_READ(iunit,fltbuf,1,kind(fltbuf),ierr); t(n)%v1(1)=real(fltbuf,WP)
     call BINARY_FILE_READ(iunit,fltbuf,1,kind(fltbuf),ierr); t(n)%v1(2)=real(fltbuf,WP)
     call BINARY_FILE_READ(iunit,fltbuf,1,kind(fltbuf),ierr); t(n)%v1(3)=real(fltbuf,WP)

     call BINARY_FILE_READ(iunit,fltbuf,1,kind(fltbuf),ierr); t(n)%v2(1)=real(fltbuf,WP)
     call BINARY_FILE_READ(iunit,fltbuf,1,kind(fltbuf),ierr); t(n)%v2(2)=real(fltbuf,WP)
     call BINARY_FILE_READ(iunit,fltbuf,1,kind(fltbuf),ierr); t(n)%v2(3)=real(fltbuf,WP)

     call BINARY_FILE_READ(iunit,fltbuf,1,kind(fltbuf),ierr); t(n)%v3(1)=real(fltbuf,WP)
     call BINARY_FILE_READ(iunit,fltbuf,1,kind(fltbuf),ierr); t(n)%v3(2)=real(fltbuf,WP)
     call BINARY_FILE_READ(iunit,fltbuf,1,kind(fltbuf),ierr); t(n)%v3(3)=real(fltbuf,WP)

     call BINARY_FILE_READ(iunit,padding,2,kind(padding),ierr)
  end do

  ! Output extents
  xmin=+huge(1.0_WP);ymin=+huge(1.0_WP);zmin=+huge(1.0_WP)
  xmax=-huge(1.0_WP);ymax=-huge(1.0_WP);zmax=-huge(1.0_WP)
  do n=1,nt
     if (t(n)%v1(1).lt.xmin) xmin=t(n)%v1(1)
     if (t(n)%v2(1).lt.xmin) xmin=t(n)%v2(1)
     if (t(n)%v3(1).lt.xmin) xmin=t(n)%v3(1)

     if (t(n)%v1(2).lt.ymin) ymin=t(n)%v1(2)
     if (t(n)%v2(2).lt.ymin) ymin=t(n)%v2(2)
     if (t(n)%v3(2).lt.ymin) ymin=t(n)%v3(2)

     if (t(n)%v1(3).lt.zmin) zmin=t(n)%v1(3)
     if (t(n)%v2(3).lt.zmin) zmin=t(n)%v2(3)
     if (t(n)%v3(3).lt.zmin) zmin=t(n)%v3(3)

     if (t(n)%v1(1).gt.xmax) xmax=t(n)%v1(1)
     if (t(n)%v2(1).gt.xmax) xmax=t(n)%v2(1)
     if (t(n)%v3(1).gt.xmax) xmax=t(n)%v3(1)

     if (t(n)%v1(2).gt.ymax) ymax=t(n)%v1(2)
     if (t(n)%v2(2).gt.ymax) ymax=t(n)%v2(2)
     if (t(n)%v3(2).gt.ymax) ymax=t(n)%v3(2)

     if (t(n)%v1(3).gt.zmax) zmax=t(n)%v1(3)
     if (t(n)%v2(3).gt.zmax) zmax=t(n)%v2(3)
     if (t(n)%v3(3).gt.zmax) zmax=t(n)%v3(3)
  end do
  print*,'Original extents - x :',xmin,'-->',xmax
  print*,'Original extents - y :',ymin,'-->',ymax
  print*,'Original extents - z :',zmin,'-->',zmax

  ! Read geometric transform info
  call parser_read('Rescale',scaling)
  call parser_read('Translate',shift)

  ! Transform the STL geometry
  print*,'Rescaling and translating...'
  do n=1,nt
     t(n)%v1=t(n)%v1*scaling+shift
     t(n)%v2=t(n)%v2*scaling+shift
     t(n)%v3=t(n)%v3*scaling+shift
     t(n)%norm=normalize(t(n)%norm*scaling)
  end do

  ! Direction swapping: x->y
  call parser_read('Swap x->y',swap_xy,.false.)
  call parser_read('Swap y->x',swap_yx,.false.)
  if (swap_xy.or.swap_yx) then
     print*,'Swapping x <=> y...'
     do n=1,nt
        tmp=t(n)%v1(1);t(n)%v1(1)=t(n)%v1(2);t(n)%v1(2)=tmp
        tmp=t(n)%v2(1);t(n)%v2(1)=t(n)%v2(2);t(n)%v2(2)=tmp
        tmp=t(n)%v3(1);t(n)%v3(1)=t(n)%v3(2);t(n)%v3(2)=tmp
        tmp=t(n)%norm(1);t(n)%norm(1)=t(n)%norm(2);t(n)%norm(2)=tmp
     end do
  end if

  ! Direction swapping: y->z
  call parser_read('Swap y->z',swap_yz,.false.)
  call parser_read('Swap z->y',swap_zy,.false.)
  if (swap_yz.or.swap_zy) then
     print*,'Swapping y <=> z...'
     do n=1,nt 
        tmp=t(n)%v1(2);t(n)%v1(2)=t(n)%v1(3);t(n)%v1(3)=tmp
        tmp=t(n)%v2(2);t(n)%v2(2)=t(n)%v2(3);t(n)%v2(3)=tmp
        tmp=t(n)%v3(2);t(n)%v3(2)=t(n)%v3(3);t(n)%v3(3)=tmp
        tmp=t(n)%norm(2);t(n)%norm(2)=t(n)%norm(3);t(n)%norm(3)=tmp
     end do
  end if

  ! Direction swapping: z->x
  call parser_read('Swap x->z',swap_xz,.false.)
  call parser_read('Swap z->x',swap_zx,.false.)
  if (swap_xz.or.swap_zx) then
     print*,'Swapping z <=> x...'
     do n=1,nt
        tmp=t(n)%v1(1);t(n)%v1(1)=t(n)%v1(3);t(n)%v1(3)=tmp
        tmp=t(n)%v2(1);t(n)%v2(1)=t(n)%v2(3);t(n)%v2(3)=tmp
        tmp=t(n)%v3(1);t(n)%v3(1)=t(n)%v3(3);t(n)%v3(3)=tmp
        tmp=t(n)%norm(1);t(n)%norm(1)=t(n)%norm(3);t(n)%norm(3)=tmp
     end do
  end if

  ! Output extents
  call parser_read('Limit extents',limit_extents,.false.)
  if (limit_extents) then
     call parser_read('xmin',xmin)
     call parser_read('ymin',ymin)
     call parser_read('zmin',zmin)
     call parser_read('xmax',xmax)
     call parser_read('ymax',ymax)
     call parser_read('zmax',zmax)
     nt2=0
     do n=1,nt
        if ( t(n)%v1(1).ge.xmin .and. t(n)%v2(1).ge.xmin .and. t(n)%v3(1).ge.xmin .and.&
             t(n)%v1(2).ge.ymin .and. t(n)%v2(2).ge.ymin .and. t(n)%v3(2).ge.ymin .and.&
             t(n)%v1(3).ge.zmin .and. t(n)%v2(3).ge.zmin .and. t(n)%v3(3).ge.zmin .and.&
             t(n)%v1(1).le.xmax .and. t(n)%v2(1).le.xmax .and. t(n)%v3(1).le.xmax .and.&
             t(n)%v1(2).le.ymax .and. t(n)%v2(2).le.ymax .and. t(n)%v3(2).le.ymax .and.&
             t(n)%v1(3).le.zmax .and. t(n)%v2(3).le.zmax .and. t(n)%v3(3).le.zmax) nt2=nt2+1
     end do
     print*,'New # of triangles :',nt2
     allocate(t2(nt2))
     nt2=0
     do n=1,nt
        if ( t(n)%v1(1).ge.xmin .and. t(n)%v2(1).ge.xmin .and. t(n)%v3(1).ge.xmin .and.&
             t(n)%v1(2).ge.ymin .and. t(n)%v2(2).ge.ymin .and. t(n)%v3(2).ge.ymin .and.&
             t(n)%v1(3).ge.zmin .and. t(n)%v2(3).ge.zmin .and. t(n)%v3(3).ge.zmin .and.&
             t(n)%v1(1).le.xmax .and. t(n)%v2(1).le.xmax .and. t(n)%v3(1).le.xmax .and.&
             t(n)%v1(2).le.ymax .and. t(n)%v2(2).le.ymax .and. t(n)%v3(2).le.ymax .and.&
             t(n)%v1(3).le.zmax .and. t(n)%v2(3).le.zmax .and. t(n)%v3(3).le.zmax) then
           nt2=nt2+1
           t2(nt2)=t(n)
        end if
     end do
     deallocate(t)
     nt=nt2
     allocate(t(nt2))
     do n=1,nt
        t(n)=t2(n)
     end do
     deallocate(t2)
  else
     xmin=+huge(1.0_WP);ymin=+huge(1.0_WP);zmin=+huge(1.0_WP)
     xmax=-huge(1.0_WP);ymax=-huge(1.0_WP);zmax=-huge(1.0_WP)
     do n=1,nt
        if (t(n)%v1(1).lt.xmin) xmin=t(n)%v1(1)
        if (t(n)%v2(1).lt.xmin) xmin=t(n)%v2(1)
        if (t(n)%v3(1).lt.xmin) xmin=t(n)%v3(1)

        if (t(n)%v1(2).lt.ymin) ymin=t(n)%v1(2)
        if (t(n)%v2(2).lt.ymin) ymin=t(n)%v2(2)
        if (t(n)%v3(2).lt.ymin) ymin=t(n)%v3(2)

        if (t(n)%v1(3).lt.zmin) zmin=t(n)%v1(3)
        if (t(n)%v2(3).lt.zmin) zmin=t(n)%v2(3)
        if (t(n)%v3(3).lt.zmin) zmin=t(n)%v3(3)

        if (t(n)%v1(1).gt.xmax) xmax=t(n)%v1(1)
        if (t(n)%v2(1).gt.xmax) xmax=t(n)%v2(1)
        if (t(n)%v3(1).gt.xmax) xmax=t(n)%v3(1)

        if (t(n)%v1(2).gt.ymax) ymax=t(n)%v1(2)
        if (t(n)%v2(2).gt.ymax) ymax=t(n)%v2(2)
        if (t(n)%v3(2).gt.ymax) ymax=t(n)%v3(2)

        if (t(n)%v1(3).gt.zmax) zmax=t(n)%v1(3)
        if (t(n)%v2(3).gt.zmax) zmax=t(n)%v2(3)
        if (t(n)%v3(3).gt.zmax) zmax=t(n)%v3(3)
     end do
  end if

  print*,'New extents - x :',xmin,'-->',xmax
  print*,'New extents - y :',ymin,'-->',ymax
  print*,'New extents - z :',zmin,'-->',zmax

  ! Mirror a direction
  nt2=nt
  call parser_read('Mirror in x',mirror_x,.false.)
  if (mirror_x) nt2=2*nt2
  call parser_read('Mirror in y',mirror_y,.false.)
  if (mirror_y) nt2=2*nt2
  call parser_read('Mirror in z',mirror_z,.false.)
  if (mirror_z) nt2=2*nt2
  if (mirror_x.or.mirror_y.or.mirror_z) then
     print*,'New # of triangles :',nt2
     allocate(t2(nt2))
     do n=1,nt
        t2(n)=t(n)
     end do
     nn=nt

     ! Mirror in x
     if (mirror_x) then
        call parser_read('Mirror length in x',Lx)
        do n=1,nn
           t2(n+nn)%norm=t2(n)%norm
           t2(n+nn)%norm(1)=-t2(n)%norm(1)
           t2(n+nn)%v1(2)=t2(n)%v1(2)
           t2(n+nn)%v1(3)=t2(n)%v1(3)
           t2(n+nn)%v2(2)=t2(n)%v2(2)
           t2(n+nn)%v2(3)=t2(n)%v2(3)
           t2(n+nn)%v3(2)=t2(n)%v3(2)
           t2(n+nn)%v3(3)=t2(n)%v3(3)
           t2(n+nn)%v1(1)=2.0_WP*Lx-t2(n)%v1(1)+2.0_WP*xmin
           t2(n+nn)%v2(1)=2.0_WP*Lx-t2(n)%v2(1)+2.0_WP*xmin
           t2(n+nn)%v3(1)=2.0_WP*Lx-t2(n)%v3(1)+2.0_WP*xmin
        end do
        nn=2*nn
     end if

     ! Mirror in y
     if (mirror_y) then
        call parser_read('Mirror length in y',Ly)
        do n=1,nn
           t2(n+nn)%norm=t2(n)%norm
           t2(n+nn)%norm(2)=-t2(n)%norm(2)
           t2(n+nn)%v1(1)=t2(n)%v1(1)
           t2(n+nn)%v1(3)=t2(n)%v1(3)
           t2(n+nn)%v2(1)=t2(n)%v2(1)
           t2(n+nn)%v2(3)=t2(n)%v2(3)
           t2(n+nn)%v3(1)=t2(n)%v3(1)
           t2(n+nn)%v3(3)=t2(n)%v3(3)
           t2(n+nn)%v1(2)=2.0_WP*Ly-t2(n)%v1(2)+2.0_WP*ymin
           t2(n+nn)%v2(2)=2.0_WP*Ly-t2(n)%v2(2)+2.0_WP*ymin
           t2(n+nn)%v3(2)=2.0_WP*Ly-t2(n)%v3(2)+2.0_WP*ymin
        end do
        nn=2*nn
     end if

     ! Mirror in z
     if (mirror_z) then
        call parser_read('Mirror length in z',Lz)
        do n=1,nn
           t2(n+nn)%norm=t2(n)%norm
           t2(n+nn)%v1(1)=t2(n)%v1(1)
           t2(n+nn)%v1(2)=t2(n)%v1(2)
           t2(n+nn)%v2(1)=t2(n)%v2(1)
           t2(n+nn)%v2(2)=t2(n)%v2(2)
           t2(n+nn)%v3(1)=t2(n)%v3(1)
           t2(n+nn)%v3(2)=t2(n)%v3(2)
           t2(n+nn)%v1(3)=2.0_WP*Lz-t2(n)%v1(3)+2.0_WP*zmin
           t2(n+nn)%v2(3)=2.0_WP*Lz-t2(n)%v2(3)+2.0_WP*zmin
           t2(n+nn)%v3(3)=2.0_WP*Lz-t2(n)%v3(3)+2.0_WP*zmin
        end do
        nn=2*nn
     end if

     ! Update triangles
     deallocate(t)
     nt=nt2
     allocate(t(nt2))
     do n=1,nt
        t(n)=t2(n)
     end do
     deallocate(t2)
  end if

  ! Close the file
  call BINARY_FILE_CLOSE(iunit,ierr)
  ! ============================================================



  ! ============================================================
  ! Set up the grid and stencil operators
  call simulation_flags_setup
  call get_dimensions(gridFile)
  call geometry_setup
  call operator_setup
  call grid_setup
  call simulation_read(IO_GRID, trim(gridFile))
  call grid_metrics_setup

  ! Allocate distance array
  allocate(levelset(nGridPoints,1)); levelset = 0.0_WP
  
  ! Prepare counter
  prog=0
  togo=nGridPoints
  iratio_old=-1
  
  ! Compute projection
  if (irank.eq.iroot) then
     print*
     print*,'Computing distance to triangles...'
  end if
  do i = 1, nGridPoints

     ! Prepare projections
     count=0
     mydist=huge(1.0_WP)

     ! Get centroid information
     c=0.0_WP
     c(1:nDimensions) = coordinates(i,1:nDimensions)

     ! Loop over triangles and check distance
     do n=1,nt

        ! Normalize triangle normal
        tmp=sqrt(dot_product(t(n)%norm,t(n)%norm))
        t(n)%norm=t(n)%norm/tmp

        call triproj(c,t(n)%v1,t(n)%v2,t(n)%v3,newproj)
        newdist=sqrt(dot_product(c-newproj,c-newproj))

        ! Check point
        if (newdist.lt.mydist-eps) then
           ! new closest point
           mydist=newdist
           myproj=newproj
           mynorm=t(n)%norm
        else if (newdist.lt.mydist+eps) then
           ! Choose better normal
           if ( sqrt(abs(dot_product(t(n)%norm,(newproj-c)/sqrt(dot_product(newproj-c,newproj-c))))) .gt. &
                sqrt(abs(dot_product(mynorm,   (myproj -c)/sqrt(dot_product(myproj -c,myproj -c))))) ) then
              mynorm=t(n)%norm
              myproj=newproj
              mydist=newdist
           end if
        end if

     end do

     ! Postprocess distance
     levelset(i,1)=mydist

     ! Get sign based on normal
     tmp=dot_product(myproj-c,mynorm)
     if (tmp.gt.0.0_WP) levelset(i,1)=-levelset(i,1)
     !levelset(i,1)=sign_switch*levelset(i,1)

     ! Add point to counter
     if (irank.eq.iroot) then
        prog=prog+1
        iratio_new=int(real(prog,WP)/real(togo,WP)*100.0_WP)
        if (iratio_new.gt.iratio_old) then
           iratio_old=iratio_new
           write(*,'(i3,x,a1)') iratio_new,'%'
        end if
     end if

  end do

  ! =========================================================
  ! Write the levelset
  ! Prefix the filename to prevent potential file locking
  call parser_read('levelset file to write', filename)
  call simulation_write(IO_LEVELSET, filename)

  ! Finalize the parallel environment
  call parallel_finalize

end program stl2levelset



! -------------------------------------------
subroutine die(errorText)

  ! External modules
  use parallel

  implicit none

  ! Arguments
  character(len = *), intent(in) :: errorText

  call parallel_kill(errorText)
  
  return
end subroutine die
