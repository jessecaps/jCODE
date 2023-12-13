! ============================== !
! Computational Geometry Toolbox !
! ============================== !
module compgeom
  use precision
  use math
  implicit none
  
  ! Tolerance
  real(WP), parameter :: eps_geom=1.0e-9_WP
  
  ! Triangulation of an n-cube into n-simplices ===============
  ! Minimal triangulation of a 2-cube -------------------------
  integer, dimension(3,2) :: c2s2
  data c2s2/1,2,4,1,3,4/
  ! Minimal triangulation of a 3-cube -------------------------
  integer, dimension(4,5) :: c2s3
  data c2s3/1,2,3,5,2,3,4,8,3,5,7,8,2,3,5,8,2,5,6,8/
  ! Periodic triangulation of a 3-cube ------------------------
  integer, dimension(4,6) :: c2s3p
  data c2s3p/1,2,4,6,1,3,4,8,1,5,6,8,1,3,5,8,3,5,7,8,1,4,6,8/
  ! Minimal triangulation of a 4-cube -------------------------
  integer, dimension(5,16) :: c2s4
  data c2s4/1,2,3,5,9,5,9,13,14,15,3,9,11,12,15,3,5,7,8,15,2,&
       9,10,12,14,2,5,6,8,14,2,3,4,8,12,8,12,14,15,16,2,3,5,9,&
       15,2,5,9,14,15,2,3,9,12,15,2,3,5,8,15,2,9,12,14,15,2,5,&
       8,14,15,2,3,8,12,15,2,8,12,14,15/
  ! ===========================================================
  
contains
  
  ! TWO-DIMENSIONAL ALGORITHM =================================
  ! Calculations are performed in 3D
  ! Intersect a 2-simplex with a line
  subroutine cut_2simplex(iso,pos,val,surf,vol,com)
    implicit none
    
    real(WP) :: iso,surf,vol
    real(WP), dimension(3,3) :: pos
    real(WP), dimension(  3) :: val
    real(WP), dimension(  3) :: com
    
    integer :: i,np,nm
    logical :: swap
    real(WP), dimension(3,3) :: p
    real(WP), dimension(  3) :: v
    real(WP), dimension(  3) :: i1,i2
    real(WP), dimension(  3) :: com0
    real(WP) :: vol0
    
    ! Classify and reorder the simplex
    np=0; nm=0; swap=.false.
    do i=1,3
       if (val(i).ge.0.0_WP) then
          np=np+1; p(:,3-np+1)=pos(:,i); v(3-np+1)=val(i)
       else
          nm=nm+1; p(:,  nm  )=pos(:,i); v(  nm  )=val(i)
       end if
    end do
    if (np.gt.nm) then
       swap=.true.
       i1 =p(:,1); p(:,1)=p(:,3); p(:,3)=i1
       vol=v(  1); v(  1)=v(  3); v(  3)=vol
    end if
    
    ! Cut the simplex
    select case (min(np,nm))
    case (0)
       surf=0.0_WP
       vol =0.0_WP
       com =0.0_WP
    case (1)
       i1=vertex_interp(3,iso,p(:,1),p(:,3),v(1),v(3))
       i2=vertex_interp(3,iso,p(:,2),p(:,3),v(2),v(3))
       vol =simplex_volume(3,2,(/i1,i2,p(:,3)/))
       com =(i1+i2+p(:,3))/3.0_WP
       surf=simplex_volume(3,1,(/i1,i2/))
    end select
    
    ! Swap sign if necessary
    if (swap) then
       vol0=simplex_volume(3,2,(/p(:,1),p(:,2),p(:,3)/))
       com0=(p(:,1)+p(:,2)+p(:,3))/3.0_WP
       com =(vol0*com0-vol*com)/(vol0-vol)
       vol =vol0-vol
    end if
    
    return
  end subroutine cut_2simplex
  
  ! Intersect a 2-cube with a line
  subroutine cut_2cube(iso,pos,val,surf,vol,com)
    implicit none
    
    real(WP) :: iso,surf,vol
    real(WP), dimension(3,4) :: pos
    real(WP), dimension(  4) :: val
    real(WP), dimension(  3) :: com
    
    real(WP), dimension(3,3) :: my_pos
    real(WP), dimension(  3) :: my_val
    real(WP), dimension(  3) :: my_com
    real(WP) :: my_surf,my_vol
    integer  :: ns,i
    
    ! Initialize
    surf=0.0_WP
    vol =0.0_WP
    com =0.0_WP
    
    ! Loop over simplices
    do ns=1,2
       
       ! Form the simplex
       do i=1,3
          my_pos(:,i)=pos(:,c2s2(i,ns))
          my_val(  i)=val(  c2s2(i,ns))
       end do
       
       ! Cut the simplex
       call cut_2simplex(iso,my_pos,my_val,my_surf,my_vol,my_com)
       
       ! Increment
       surf=surf+my_surf
       vol =vol +my_vol
       com =com +my_vol*my_com
    end do
    
    ! Rescale CoM
    if (vol.gt.0.0_WP) com=com/vol
    
    return
  end subroutine cut_2cube
  ! ===========================================================
  
  
  ! THREE-DIMENSIONAL ALGORITHM ===============================
  ! Intersect a 3-simplex with a plane
  subroutine cut_3simplex(iso,pos,val,vol,com,surf,coms,nvec)
    implicit none
    
    real(WP) :: iso,surf,vol
    real(WP), dimension(3,4) :: pos
    real(WP), dimension(  4) :: val
    real(WP), dimension(  3) :: com
    real(WP), dimension(  3) :: coms
    real(WP), dimension(  3) :: nvec
    
    integer :: i,np,nm
    logical :: swap
    real(WP), dimension(3,4) :: p
    real(WP), dimension(  4) :: v
    real(WP), dimension(  3) :: i1,i2,i3,i4
    real(WP), dimension(  3) :: com0,com1,com2,com3
    real(WP) :: vol0,vol1,vol2,vol3,buf
    
    ! Classify and reorder the simplex
    np=0; nm=0; swap=.false.
    do i=1,4
       if (val(i).ge.0.0_WP) then
          np=np+1; p(:,4-np+1)=pos(:,i); v(4-np+1)=val(i)
       else
          nm=nm+1; p(:,  nm  )=pos(:,i); v(  nm  )=val(i)
       end if
    end do
    if (np.gt.nm) then
       swap=.true.
       i1 =p(:,1); p(:,1)=p(:,4); p(:,4)=i1
       vol=v(  1); v(  1)=v(  4); v(  4)=vol
    end if
    
    ! Cut the simplex
    select case (min(np,nm))
    case (0)
       vol =0.0_WP
       com =0.0_WP
       surf=0.0_WP
       coms=0.0_WP
       nvec=0.0_WP
    case (1)
       ! Create intersections from p4
       i1=vertex_interp(3,iso,p(:,1),p(:,4),v(1),v(4))
       i2=vertex_interp(3,iso,p(:,2),p(:,4),v(2),v(4))
       i3=vertex_interp(3,iso,p(:,3),p(:,4),v(3),v(4))
       ! Volume & COM of cut tet
       vol =simplex_volume(3,3,(/i1,i2,i3,p(:,4)/))
       com =0.25_WP*(i1+i2+i3+p(:,4))
       ! Area & COM of surface
       surf=simplex_volume(3,2,(/i1,i2,i3/))
       coms=(i1+i2+i3)/3.0_WP
       ! Surface normal
       nvec=get_normal(i1,i2,i3)
       if (dot_product(p(:,4)-coms,nvec).lt.0.0_WP) nvec=-nvec
    case (2)
       ! Create intersections from p3 and p4
       i1=vertex_interp(3,iso,p(:,1),p(:,3),v(1),v(3))
       i2=vertex_interp(3,iso,p(:,2),p(:,3),v(2),v(3))
       i3=vertex_interp(3,iso,p(:,1),p(:,4),v(1),v(4))
       i4=vertex_interp(3,iso,p(:,2),p(:,4),v(2),v(4))
       ! Volume & COM of cut tets
       vol1=simplex_volume(3,3,(/i1,i2,i3,p(:,3)/))
       vol2=simplex_volume(3,3,(/i2,i3,i4,p(:,4)/))
       vol3=simplex_volume(3,3,(/i2,i3,p(:,3),p(:,4)/))
       vol =vol1+vol2+vol3
       com1=0.25_WP*(i1+i2+i3+p(:,3))
       com2=0.25_WP*(i2+i3+i4+p(:,4))
       com3=0.25_WP*(i2+i3+p(:,3)+p(:,4))
       com =(vol1*com1+vol2*com2+vol3*com3)/vol
       ! Area & COM of surfaces
       vol1=simplex_volume(3,2,(/i1,i2,i3/))
       vol2=simplex_volume(3,2,(/i2,i3,i4/))
       surf=vol1+vol2
       com1=(i1+i2+i3)/3.0_WP
       com2=(i2+i3+i4)/3.0_WP
       coms=(vol1*com1+vol2*com2)/surf
       ! Surface normal
       nvec=get_normal(i1,i2,i3)
       if (dot_product(p(:,3)-com1,nvec).lt.0.0_WP) nvec=-nvec
       com1=vol1*nvec
       nvec=get_normal(i2,i3,i4)
       if (dot_product(p(:,4)-com2,nvec).lt.0.0_WP) nvec=-nvec
       nvec=vol2*nvec+com1
       buf=sqrt(dot_product(nvec,nvec))+epsilon(1.0_WP)
       nvec=nvec/buf
    end select
        
    ! Swap sign if necessary
    if (swap) then
       vol0=simplex_volume(3,3,(/p(:,1),p(:,2),p(:,3),p(:,4)/))
       com0=0.25_WP*(p(:,1)+p(:,2)+p(:,3)+p(:,4))
       com =(vol0*com0-vol*com)/(vol0-vol)
       vol =vol0-vol
       nvec=-nvec
    end if
    
    return
  end subroutine cut_3simplex
  
  ! Intersect a 3-cube with a plane
  subroutine cut_3cube(iso,pos,val,vol,com,surf,coms,nvec)
    implicit none
    
    real(WP) :: iso,surf,vol
    real(WP), dimension(3,8) :: pos
    real(WP), dimension(  8) :: val
    real(WP), dimension(  3) :: com
    real(WP), dimension(  3) :: coms
    real(WP), dimension(  3) :: nvec
    
    real(WP), dimension(3,4) :: my_pos
    real(WP), dimension(  4) :: my_val
    real(WP), dimension(  3) :: my_com
    real(WP), dimension(  3) :: my_coms
    real(WP), dimension(  3) :: my_nvec
    real(WP) :: my_surf,my_vol,buf
    integer  :: ns,i
    
    ! Initialize
    vol =0.0_WP
    com =0.0_WP
    surf=0.0_WP
    coms=0.0_WP
    nvec=0.0_WP
    
    ! Loop over simplices
    do ns=1,5
       
       ! Form the simplex
       do i=1,4
          my_pos(:,i)=pos(:,c2s3(i,ns))
          my_val(  i)=val(  c2s3(i,ns))
       end do
       
       ! Cut the simplex
       call cut_3simplex(iso,my_pos,my_val,my_vol,my_com,my_surf,my_coms,my_nvec)
       
       ! Increment
       vol =vol +my_vol
       com =com +my_vol *my_com
       surf=surf+my_surf
       coms=coms+my_surf*my_coms
       nvec=nvec+my_surf*my_nvec
    end do
    
    ! Rescale COM
    if (vol .gt.0.0_WP) com =com /vol
    if (surf.gt.0.0_WP) then
       coms=coms/surf
       buf=sqrt(dot_product(nvec,nvec))+epsilon(1.0_WP)
       nvec=nvec/buf
    end if
    
    return
  end subroutine cut_3cube
  ! ===========================================================
  
  
  ! FOUR-DIMENSIONAL ALGORITHM ================================
  ! Intersect a 4-simplex with a hyperplane
  subroutine cut_4simplex(iso,pos,val,surf,vol,com)
    implicit none
    
    real(WP) :: iso,surf,vol
    real(WP), dimension(4,5) :: pos
    real(WP), dimension(  5) :: val
    real(WP), dimension(  4) :: com
    
    integer :: i,np,nm
    logical :: swap
    real(WP), dimension(4,5) :: p
    real(WP), dimension(  5) :: v
    real(WP), dimension(  4) :: i1,i2,i3,i4,i5,i6
    real(WP), dimension(  4) :: com0,com1,com2,com3,com4
    real(WP) :: vol0,vol1,vol2,vol3,vol4
    
    ! Classify and reorder the simplex
    np=0; nm=0; swap=.false.
    do i=1,5
       if (val(i).ge.0.0_WP) then
          np=np+1; p(:,5-np+1)=pos(:,i); v(5-np+1)=val(i)
       else
          nm=nm+1; p(:,  nm  )=pos(:,i); v(  nm  )=val(i)
       end if
    end do
    if (np.gt.nm) then
       swap=.true.
       i1 =p(:,1); p(:,1)=p(:,5); p(:,5)=i1
       vol=v(  1); v(  1)=v(  5); v(  5)=vol
       i1 =p(:,2); p(:,2)=p(:,4); p(:,4)=i1
       vol=v(  2); v(  2)=v(  4); v(  4)=vol
    end if
    
    ! Cut the simplex
    select case (min(np,nm))
    case (0)
       surf=0.0_WP
       vol =0.0_WP
       com =0.0_WP
    case (1)
       i1=vertex_interp(4,iso,p(:,1),p(:,5),v(1),v(5))
       i2=vertex_interp(4,iso,p(:,2),p(:,5),v(2),v(5))
       i3=vertex_interp(4,iso,p(:,3),p(:,5),v(3),v(5))
       i4=vertex_interp(4,iso,p(:,4),p(:,5),v(4),v(5))
       vol =simplex_volume(4,4,(/i1,i2,i3,i4,p(:,5)/))
       surf=simplex_volume(4,3,(/i1,i2,i3,i4/))
       com =0.2_WP*(i1+i2+i3+i4+p(:,5))
    case (2)
       i1=vertex_interp(4,iso,p(:,1),p(:,4),v(1),v(4))
       i2=vertex_interp(4,iso,p(:,2),p(:,4),v(2),v(4))
       i3=vertex_interp(4,iso,p(:,1),p(:,5),v(1),v(5))
       i4=vertex_interp(4,iso,p(:,3),p(:,4),v(3),v(4))
       i5=vertex_interp(4,iso,p(:,2),p(:,5),v(2),v(5))
       i6=vertex_interp(4,iso,p(:,3),p(:,5),v(3),v(5))
       vol1=simplex_volume(4,4,(/i1,i2,i3,i4,p(:,5)/))
       vol2=simplex_volume(4,4,(/i2,i3,i4,i5,p(:,5)/))
       vol3=simplex_volume(4,4,(/i3,i4,i5,i6,p(:,5)/))
       vol4=simplex_volume(4,4,(/i1,i2,i4,p(:,4),p(:,5)/))
       vol =vol1+vol2+vol3+vol4
       surf=simplex_volume(4,3,(/i1,i2,i3,i4/))+&
            simplex_volume(4,3,(/i2,i3,i4,i5/))+&
            simplex_volume(4,3,(/i3,i4,i5,i6/))
       com1=0.2_WP*(i1+i2+i3+i4+p(:,5))
       com2=0.2_WP*(i2+i3+i4+i5+p(:,5))
       com3=0.2_WP*(i3+i4+i5+i6+p(:,5))
       com4=0.2_WP*(i1+i2+i4+p(:,4)+p(:,5))
       com =(vol1*com1+vol2*com2+vol3*com3+vol4*com4)/vol
    end select
    
    ! Swap sign if necessary
    if (swap) then
       vol0=simplex_volume(4,4,(/p(:,1),p(:,2),p(:,3),p(:,4),p(:,5)/))
       com0=0.2_WP*(p(:,1)+p(:,2)+p(:,3)+p(:,4)+p(:,5))
       com =(vol0*com0-vol*com)/(vol0-vol)
       vol =vol0-vol
    end if
    
    return
  end subroutine cut_4simplex
  
  ! Intersect a 4-cube with a hyperplane
  subroutine cut_4cube(iso,pos,val,surf,vol,com)
    implicit none
    
    real(WP) :: iso,surf,vol
    real(WP), dimension(4,16) :: pos
    real(WP), dimension(  16) :: val
    real(WP), dimension(   4) :: com
    
    real(WP), dimension(4,5) :: my_pos
    real(WP), dimension(  5) :: my_val
    real(WP), dimension(  4) :: my_com
    real(WP) :: my_surf,my_vol
    integer  :: ns,i
    
    ! Initialize
    surf=0.0_WP
    vol =0.0_WP
    com =0.0_WP
    
    ! Loop over simplices
    do ns=1,16
       
       ! Form the simplex
       do i=1,5
          my_pos(:,i)=pos(:,c2s4(i,ns))
          my_val(  i)=val(  c2s4(i,ns))
       end do
       
       ! Cut the simplex
       call cut_4simplex(iso,my_pos,my_val,my_surf,my_vol,my_com)
       
       ! Increment
       surf=surf+my_surf
       vol =vol +my_vol
       com =com +my_vol*my_com
       
    end do
    
    ! Rescale CoM
    if (vol.gt.0.0_WP) com=com/vol
    
    return
  end subroutine cut_4cube
  ! ===========================================================
  
  
  ! ===========================================================
  ! Surface normal calculation
  function get_normal(p1,p2,p3)
    implicit none
    real(WP), dimension(3), intent(in) :: p1,p2,p3
    real(WP), dimension(3) :: get_normal
    real(WP) :: norm
    get_normal(1)=(p3(2)-p2(2))*(p3(3)-p1(3))-(p3(3)-p2(3))*(p3(2)-p1(2))
    get_normal(2)=(p3(3)-p2(3))*(p3(1)-p1(1))-(p3(1)-p2(1))*(p3(3)-p1(3))
    get_normal(3)=(p3(1)-p2(1))*(p3(2)-p1(2))-(p3(2)-p2(2))*(p3(1)-p1(1))
    norm=sqrt(dot_product(get_normal,get_normal))+epsilon(1.0_WP)
    get_normal=get_normal/norm
  end function get_normal
  ! ===========================================================
  
  
  ! ===========================================================
  ! Linear vertex interpolation in n dimensions
  function vertex_interp(n,iso,p1,p2,v1,v2)
    implicit none
    integer,                intent(in) :: n
    real(WP), dimension(n)             :: vertex_interp
    real(WP),               intent(in) :: iso
    real(WP), dimension(n), intent(in) :: p1,p2
    real(WP),               intent(in) :: v1,v2
    real(WP)                           :: mu
    if (abs(v1-v2).lt.eps_geom) then
       vertex_interp = p1
    else
       mu = (iso-v1)/(v2-v1)
       vertex_interp = p1+mu*(p2-p1)
    end if
  end function vertex_interp
  ! ===========================================================
  
  
  ! ===========================================================
  ! Lebesgue measure of a m-simplex in n dimensions using the
  ! Cayley-Menger determinant
  function simplex_volume(n,m,S) result(v)
    implicit none
    integer, intent(in) :: n,m
    real(WP), dimension(n,m+1) :: S
    real(WP), dimension(m+2,m+2) :: D
    real(WP) :: v,fact
    integer :: i,j
    ! Determinant
    do i=1,m+1
       do j=1,m+1
          D(i+1,j+1)=dot_product(S(:,i)-S(:,j),S(:,i)-S(:,j))
       end do
    end do
    D(1,:)=1.0_WP; D(:,1)=1.0_WP; D(1,1)=0.0_WP
    v=matdet(D,m+2)
    ! Factorial
    fact=1.0_WP
    do i=2,m
       fact=fact*real(i,WP)
    end do
    ! Lebesgue measure
    v=sqrt(abs(v))/(2.0_WP**(real(m,WP)/2.0_WP)*fact)
    return
  end function simplex_volume
  ! ===========================================================
  
  
  ! ===========================================================
  ! 3D projection of a vertex onto a triangle
  subroutine triproj(myp,myt1,myt2,myt3,proj)
    implicit none
    real(WP), dimension(3), intent(in)  :: myp,myt1,myt2,myt3
    real(WP), dimension(3), intent(out) :: proj
    real(WP), dimension(3) :: v1,v2,vp
    real(WP) :: a,b,c,d,e,f
    real(WP) :: det,s,t,inv
    real(WP) :: denom,numer,tmp0,tmp1
    ! To do: check for colinearity and/or too small triangles
    
    ! Build triangle information
    v1=myt2-myt1
    v2=myt3-myt1
    vp=myt1-myp
    a=dot_product(v1,v1)
    b=dot_product(v1,v2)
    c=dot_product(v2,v2)
    d=dot_product(v1,vp)
    e=dot_product(v2,vp)
    f=dot_product(vp,vp)
    det=a*c-b*b
    s  =b*e-c*d
    t  =b*d-a*e
    
    ! Check if projection lies inside the triangle
    if (s+t.le.det) then
       if (s.lt.0.0_WP) then
          if (t.lt.0.0_WP) then
             if (d.lt.0.0_WP) then
                t=0.0_WP
                if (-d.ge.a) then
                   s=1.0_WP
                else
                   s=-d/a
                end if
             else
                s=0.0_WP
                if (e.ge.0.0_WP) then
                   t=0.0_WP
                else
                   if (-e.ge.c) then
                      t=1.0_WP
                   else
                      t=-e/c
                   end if
                end if
             end if
          else
             s=0.0_WP
             if (e.ge.0.0_WP) then
                t=0.0_WP
             else
                if (-e.ge.c) then
                   t=1.0_WP
                else
                   t=-e/c
                end if
             end if
          end if
       else
          if (t.lt.0.0_WP) then
             t=0.0_WP
             if (d.ge.0.0_WP) then
                s=0.0_WP
             else
                if (-d.ge.a) then
                   s=1.0_WP
                else
                   s=-d/a
                end if
             end if
          else
             inv=1.0_WP/det
             s=s*inv
             t=t*inv
          end if
       end if
    else
       if (s.lt.0.0_WP) then
          tmp0=b+d
          tmp1=c+e
          if (tmp1.gt.tmp0) then
             numer=tmp1-tmp0
             denom=a-2.0_WP*b+c
             if (numer.ge.denom) then
                s=1.0_WP
                t=0.0_WP
             else
                s=numer/denom
                t=1.0_WP-s
             end if
          else
             s=0.0_WP
             if (tmp1.le.0.0_WP) then
                t=1.0_WP
             else
                if (e.ge.0.0_WP) then
                   t=0.0_WP
                else
                   t=-e/c
                end if
             end if
          end if
       else
          if (t.lt.0.0_WP) then
             tmp0=b+e
             tmp1=a+d
             if (tmp1.gt.tmp0) then
                numer=tmp1-tmp0
                denom=a-2.0_WP*b+c
                if (numer.ge.denom) then
                   t=1.0_WP
                   s=0.0_WP
                else
                   t=numer/denom
                   s=1.0_WP-t
                end if
             else
                t=0.0_WP
                if (tmp1.le.0.0_WP) then
                   s=1.0_WP
                else
                   if (d.ge.0.0_WP) then
                      s=0.0_WP
                   else
                      s=-d/a
                   end if
                end if
             end if
          else
             numer=c+e-b-d
             if (numer.le.0.0_WP) then
                s=0.0_WP
                t=1.0_WP
             else
                denom=a-2.0_WP*b+c
                if (numer.ge.denom) then
                   s=1.0_WP
                   t=0.0_WP
                else
                   s=numer/denom
                   t=1.0_WP-s
                end if
             end if
          end if
       end if
    end if
    
    ! Get projection
    proj=myt1+s*v1+t*v2
    
    return
  end subroutine triproj
  ! ===========================================================
  
  
  ! ===========================================================
  !  Triangulation of cut plane on 3-cell
  subroutine triangulate_3cell(iso,pos,val,n,tri)
    implicit none
    
    ! Inputs
    real(WP), intent(in) :: iso
    real(WP), dimension(3,8), intent(in) :: pos
    real(WP), dimension(  8), intent(in) :: val
    
    ! Outputs
    integer :: n ! number of triangles
    real(WP), dimension(3,3,12) :: tri !triangle vertices

    ! Working Variables
    integer :: i,ns,p1,p2,n_intersect
    real(WP), dimension(3,4) :: my_pos
    real(WP), dimension(  4) :: my_val
    real(WP), dimension(4,3) :: my_int

    ! Initialize
    n=0; tri=0.0_WP

    ! Loop over simplices
    do ns=1,6

       ! Form the simplex
       do i=1,4
          my_pos(:,i)=pos(:,c2s3p(i,ns))
          my_val(  i)=val(  c2s3p(i,ns))
       end do

       ! Initialize counter
       n_intersect=0

       ! Loop over 6 edges of simplex and check for iso
       do i=1,6
          ! Select edge
          select case (i)
          case (1); p1=1; p2=2;
          case (2); p1=1; p2=3;
          case (3); p1=1; p2=4;
          case (4); p1=2; p2=3;
          case (5); p1=2; p2=4;
          case (6); p1=3; p2=4;
          end select
             
          ! Check for iso along edge
          if ( my_val(p1).gt.iso .and. my_val(p2).lt.iso .or. &
               my_val(p1).lt.iso .and. my_val(p2).gt.iso ) then
             n_intersect=n_intersect+1
             ! Save intersection point
             my_int(n_intersect,1)=my_pos(1,p1)+(my_pos(1,p2)-my_pos(1,p1))&
                  *(iso-my_val(p1))/(my_val(p2)-my_val(p1)) !+epsilon(1.0_WP))
             my_int(n_intersect,2)=my_pos(2,p1)+(my_pos(2,p2)-my_pos(2,p1))&
                  *(iso-my_val(p1))/(my_val(p2)-my_val(p1)) !+epsilon(1.0_WP))
             my_int(n_intersect,3)=my_pos(3,p1)+(my_pos(3,p2)-my_pos(3,p1))&
                  *(iso-my_val(p1))/(my_val(p2)-my_val(p1)) !+epsilon(1.0_WP))
          end if
       end do
          
       ! Make triangle(s)
       select case (n_intersect)
       case (3) ! 1 triangle
          n=n+1
          tri(:,1,n)=my_int(1,:)
          tri(:,2,n)=my_int(2,:)
          tri(:,3,n)=my_int(3,:)
       case (4) ! 2 triangles
          n=n+1
          tri(:,1,n)=my_int(1,:)
          tri(:,2,n)=my_int(2,:)
          tri(:,3,n)=my_int(3,:)
          n=n+1
          tri(:,1,n)=my_int(2,:)
          tri(:,2,n)=my_int(3,:)
          tri(:,3,n)=my_int(4,:)
       end select

    end do

    return
  end subroutine triangulate_3cell
  ! ===========================================================
  
  
  ! ===========================================================
  !  Calculate length of line segment on negative side of plane
  !   - returns length
  !   - line defined by end points P1, P2, where Pi=(xi,yi,zi)
  !   - plane defined by ax+by+cz+d=0 (normal on positive side)
  function intersect_line_plane(a,b,c,d,P1,P2) result(length)
    implicit none
    ! I/O
    real(WP), intent(in)  :: a,b,c,d
    real(WP), dimension(3), intent(in)  :: P1,P2
    real(WP) :: length
    ! Working variables
    real(WP) :: x1,x2,y1,y2,z1,z2
    real(WP) :: s,length_seg

    x1=P1(1); y1=P1(2); z1=P1(3);
    x2=P2(1); y2=P2(2); z2=P2(3);
    
    length_seg=sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
    
    if (a*x1+b*y1+c*z1-d.le.0.0_WP .and. a*x2+b*y2+c*z2-d.le.0.0_WP) then
       ! line segment all liquid
       length=length_seg
    else if (a*x1+b*y1+c*z1-d.gt.0.0_WP .and. a*x2+b*y2+c*z2-d.gt.0.0_WP) then
       ! line segment all gas
       length=0.0_WP
    else
       ! line segment intersects plane
       s =  -(a*x1+b*y1+c*z1-d) & 
            /(a*(x2-x1)+b*(y2-y1)+c*(z2-z1))
       length=s*length_seg
       ! Measure length of liquid side
       if (a*x1+b*y1+c*z1-d.gt.0.0_WP) length=length_seg-length
    end if

    return
  end function intersect_line_plane
  
end module compgeom
