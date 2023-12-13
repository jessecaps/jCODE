program stl2ensight

  ! External modules
  use precision
  use string
  use fileio
  use parser
  use math

  implicit none

  ! General
  integer :: ierr,iunit,n
  
  ! STL file
  character(len=str_medium) :: stl_file
  character(len=80) :: cbuff
  character(len= 2) :: padding
  integer  :: nt,nt2,nn
  real(SP) :: fltbuf
  real(WP) :: tmp
  type triangle_type
     real(WP), dimension(3) :: norm
     real(WP), dimension(3) :: v1
     real(WP), dimension(3) :: v2
     real(WP), dimension(3) :: v3
  end type triangle_type
  type(triangle_type), dimension(:), allocatable :: t,t2
  real(WP) :: xmin,xmax,ymin,ymax,zmin,zmax,Lx,Ly,Lz
  real(WP), dimension(3) :: shift,scaling
  logical :: swap_xy,swap_yz,swap_zx
  logical :: swap_yx,swap_zy,swap_xz
  logical :: limit_extents,mirror_x,mirror_y,mirror_z
  
  ! Ensight file
  character(len=str_medium) :: ensight_file
  character(len=str_medium) :: ensight_norm
  character(len=str_medium) :: case_file
  character(len=80) :: cbuffer
  character(len=80) :: str
  integer  :: ibuffer
  real(SP) :: rbuffer
  
  ! Input processing
  character(len=str_medium) :: input_name
  
  ! Parse the command line
  call get_command_argument(1,input_name)
  
  ! Initialize the parser
  call parser_init
  
  ! Read the input file
  call parser_parsefile(input_name)
  
  ! Read stl file
  call parser_readchar('STL file',stl_file)
  
  ! Open the STL file to read
  call BINARY_FILE_OPEN(iunit,trim(stl_file),"r",ierr)
  
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
  ! Read ensight file
  call parser_readchar('Ensight file',ensight_file)
  ensight_norm=trim(ensight_file) // ".norm"
  
  ! Generate case file
  case_file=trim(adjustl(ensight_file)) // ".case"
  iunit = iopen()
  open(iunit,file=case_file,form="formatted",iostat=ierr,status="REPLACE")
  str='FORMAT'
  write(iunit,'(a80)') str
  str='type: ensight gold'
  write(iunit,'(a80)') str
  str='GEOMETRY'
  write(iunit,'(a80)') str
  str='model: ' // ensight_file
  write(iunit,'(a80)') str
  str='VARIABLE'
  write(iunit,'(a80)') str
  str='vector per element: normal ' // trim(ensight_norm)
  write(iunit,'(a80)') str
  close(iclose(iunit))
  
  ! Generate the geometry
  call BINARY_FILE_OPEN(iunit,trim(ensight_file),"w",ierr)
  
  ! Write header
  cbuffer = 'C Binary'
  call BINARY_FILE_WRITE(iunit,cbuffer,80,kind(cbuffer),ierr)
  cbuffer = 'Ensight Gold Geometry File'
  call BINARY_FILE_WRITE(iunit,cbuffer,80,kind(cbuffer),ierr)
  cbuffer = 'STL geometry from NGA'
  call BINARY_FILE_WRITE(iunit,cbuffer,80,kind(cbuffer),ierr)
  cbuffer = 'node id off'
  call BINARY_FILE_WRITE(iunit,cbuffer,80,kind(cbuffer),ierr)
  cbuffer = 'element id off'
  call BINARY_FILE_WRITE(iunit,cbuffer,80,kind(cbuffer),ierr)
  cbuffer = 'part'
  call BINARY_FILE_WRITE(iunit,cbuffer,80,kind(cbuffer),ierr)
  ibuffer = 1
  call BINARY_FILE_WRITE(iunit,ibuffer,1,kind(ibuffer),ierr)
  cbuffer = 'STL data from NGA'
  call BINARY_FILE_WRITE(iunit,cbuffer,80,kind(cbuffer),ierr)
  cbuffer = 'coordinates'
  call BINARY_FILE_WRITE(iunit,cbuffer,80,kind(cbuffer),ierr)
  ibuffer = 3*nt
  call BINARY_FILE_WRITE(iunit,ibuffer,1,kind(ibuffer),ierr)
  
  ! Write nodal position - x
  do n=1,nt
     rbuffer=real(t(n)%v1(1), SP)
     call BINARY_FILE_WRITE(iunit,rbuffer,1,kind(rbuffer),ierr)
     rbuffer=real(t(n)%v2(1), SP)
     call BINARY_FILE_WRITE(iunit,rbuffer,1,kind(rbuffer),ierr)
     rbuffer=real(t(n)%v3(1), SP)
     call BINARY_FILE_WRITE(iunit,rbuffer,1,kind(rbuffer),ierr)
  end do
  
  ! Write nodal position - y
  do n=1,nt
     rbuffer=real(t(n)%v1(2), SP)
     call BINARY_FILE_WRITE(iunit,rbuffer,1,kind(rbuffer),ierr)
     rbuffer=real(t(n)%v2(2), SP)
     call BINARY_FILE_WRITE(iunit,rbuffer,1,kind(rbuffer),ierr)
     rbuffer=real(t(n)%v3(2), SP)
     call BINARY_FILE_WRITE(iunit,rbuffer,1,kind(rbuffer),ierr)
  end do
  
  ! Write nodal position - z
  do n=1,nt
     rbuffer=real(t(n)%v1(3), SP)
     call BINARY_FILE_WRITE(iunit,rbuffer,1,kind(rbuffer),ierr)
     rbuffer=real(t(n)%v2(3), SP)
     call BINARY_FILE_WRITE(iunit,rbuffer,1,kind(rbuffer),ierr)
     rbuffer=real(t(n)%v3(3), SP)
     call BINARY_FILE_WRITE(iunit,rbuffer,1,kind(rbuffer),ierr)
  end do
  
  ! Write elements
  cbuffer = 'tria3'
  call BINARY_FILE_WRITE(iunit,cbuffer,80,kind(cbuffer),ierr)
  ibuffer = nt
  call BINARY_FILE_WRITE(iunit,ibuffer,1,kind(ibuffer),ierr)
  
  ! Write connectivity
  do n=1,nt
     ibuffer=3*n-2
     call BINARY_FILE_WRITE(iunit,ibuffer,1,kind(ibuffer),ierr)
     ibuffer=3*n-1
     call BINARY_FILE_WRITE(iunit,ibuffer,1,kind(ibuffer),ierr)
     ibuffer=3*n-0
     call BINARY_FILE_WRITE(iunit,ibuffer,1,kind(ibuffer),ierr)
  end do
  
  ! Close the file
  call BINARY_FILE_CLOSE(iunit,ierr)
  
  ! Generate the normal vector file
  call BINARY_FILE_OPEN(iunit,trim(ensight_norm),"w",ierr)
  
  ! Write header
  cbuffer = 'STL normals from NGA'
  call BINARY_FILE_WRITE(iunit,cbuffer,80,kind(cbuffer),ierr)
  cbuffer = 'part'
  call BINARY_FILE_WRITE(iunit,cbuffer,80,kind(cbuffer),ierr)
  ibuffer = 1
  call BINARY_FILE_WRITE(iunit,ibuffer,1,kind(ibuffer),ierr)
  cbuffer = 'tria3'
  call BINARY_FILE_WRITE(iunit,cbuffer,80,kind(cbuffer),ierr)
  
  ! Write normal vector
  do n=1,nt
     rbuffer=real(t(n)%norm(1), SP)
     call BINARY_FILE_WRITE(iunit,rbuffer,1,kind(rbuffer),ierr)
  end do
  do n=1,nt
     rbuffer=real(t(n)%norm(2), SP)
     call BINARY_FILE_WRITE(iunit,rbuffer,1,kind(rbuffer),ierr)
  end do
  do n=1,nt
     rbuffer=real(t(n)%norm(3), SP)
     call BINARY_FILE_WRITE(iunit,rbuffer,1,kind(rbuffer),ierr)
  end do
  
  ! Close the file
  call BINARY_FILE_CLOSE(iunit,ierr)
  
  ! ============================================================

end program stl2ensight


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
