module geodesic
  use precision
  use math

  ! Trigonometric parameters
  real(WP), parameter :: goldenRatio = 1.61803398875_WP

contains

  ! Convert spherical coordinate angles
  ! to position tangent to sphere
  ! longitude/latitude/radius => x/y/z
  function ang2pos(p,q,r) result(xyz)
    implicit none

    real(WP) :: p,q,r
    real(WP), dimension(3) :: xyz

    ! longitude p: 0     -> 2pi
    ! latitude  q: -pi/2 -> +pi/2
    ! radius    r: 0     -> infinity

    xyz(1)=r*cos(q)*p
    xyz(2)=r*q
    xyz(3)=r

    return
  end function ang2pos


  ! Convert position from spherical coordinates to Cartesian
  function sphere2cart_pos(p,q,r) result(xyz)
    implicit none

    real(WP) :: p,q,r
    real(WP), dimension(3) :: xyz

    ! longitude p: 0     -> 2pi
    ! latitude  q: -pi/2 -> +pi/2
    ! radius    r: 0     -> infinity

    xyz(1)=r*cos(q)*cos(p)
    xyz(2)=r*cos(q)*sin(p)
    xyz(3)=r*sin(q)

    return
  end function sphere2cart_pos


  ! Convert position from Cartesian coordinates to spherical
  function cart2sphere_pos(x,y,z) result(pqr)
    implicit none

    real(WP) :: x,y,z
    real(WP), dimension(3) :: pqr

    ! longitude p: 0     -> 2pi
    ! latitude  q: -pi/2 -> +pi/2
    ! radius    r: 0     -> infinity

    pqr(3)=sqrt(x**2+y**2+z**2)
    pqr(2)=asin(z/pqr(3))
    pqr(1)=arctan(y,x)

    return
  end function cart2sphere_pos


  ! Convert velocity from spherical coordinates to Cartesian
  function sphere2cart_vel(p,q,r,us,vs,ws) result(uvw)
    implicit none

    real(WP):: p,q,r,us,vs,ws
    real(WP), dimension(3) :: uvw

    ! longitude p: 0     -> 2pi
    ! latitude  q: -pi/2 -> +pi/2
    ! radius    r: 0     -> infinity

    ! zonal component:       us
    ! meridional component:  vs
    ! radial component:      ws

    uvw(1) = -us*sin(p) - vs*cos(p)*sin(q) + ws*cos(p)*cos(q)
    uvw(2) =  us*cos(p) - vs*sin(p)*sin(q) + ws*sin(p)*cos(q)
    uvw(3) =  vs*cos(q) + ws*sin(q)

    return
  end function sphere2cart_vel


  ! Convert velocity from Cartesian coordinates to spherical
  function cart2sphere_vel(x,y,z,u,v,w) result(uvw_s)
    implicit none

    real(WP) :: x,y,z,u,v,w
    real(WP), dimension(3) :: pqr,uvw_s

    ! zonal component:       uvw_s(1)
    ! meridional component:  uvw_s(2)
    ! radial component:      uvw_s(3)

    pqr=cart2sphere_pos(x,y,z)
    uvw_s(1) = -u*sin(pqr(1)) + v*cos(pqr(1))
    uvw_s(2) = -u*sin(pqr(2))*cos(pqr(1)) -v*sin(pqr(2))*sin(pqr(1)) +w*cos(pqr(2))
    uvw_s(3) =  u*cos(pqr(2))*cos(pqr(1)) +v*cos(pqr(2))*sin(pqr(1)) +w*sin(pqr(2))

    return
  end function cart2sphere_vel

  ! Compute coordinates of a geodesic grid
  ! The sphere is generated from subdivisions of the faces of an icosahedron
  function geodome(radius,freq) result(xyz)
    implicit none

    integer, intent(in) :: freq
    real(WP), intent(in) :: radius
    real(WP), dimension(2+10*freq**2,3) :: xyz

    integer :: n,f,v
    integer, dimension(12,12) :: icosA
    integer, dimension(20,3 ) :: icosF
    integer, dimension(30) :: ivec,jvec
    integer :: nVert,fm1,fm2,offset
    real(WP), allocatable, dimension(:,:) :: V1,V2,VF
    real(WP), dimension(12,3) :: icosXYZ

    if (freq.lt.1) then
       print *, 'error: Icosahedron frequency >= 1'
       stop
    end if
    call icosahedron(icosXYZ, icosA, icosF)

    ! Compute number of vertices in geodesic sphere
    fm1 = freq - 1
    fm2 = freq - 2    
    nVert = 12 + 30*fm1 + 20*fm2*fm1/2

    ! Get coordinates
    if (freq.eq. 1) then
       ! Just return the icosahedron
       xyz = icosXYZ

    else
       ! For all frequencies > 1

       ! Grab the vertices of the icosahedron
       xyz(1:12,:) = icosXYZ
       offset = 13

       ! Find the nodes that connect edges.  Note we use the upper triangular
       ! part of the adjacency matrix so that we only get each edge once
       ivec = (/ 1,3,1,3,1,3,5,2,4,2,4,7,1,2,5,7,1,2,6,8,3,4,5,7,9,3,4,6,8,10 /)
       jvec = (/ 2,4,5,5,6,6,6,7,7,8,8,8,9,9,9,9,10,10,10,10,11,11,11,11,11,12,12,12,12,12 /)

       ! Generate extra vertices on every edge of the icoshedron
       do n = 1,30
          xyz(offset:(offset+fm2) ,:) = &
               divideEdge(icosXYZ(ivec(n),:), icosXYZ(jvec(n),:), freq)
          offset = offset + fm1;
       end do

       ! Generate the extra vertices within each face of the icoshedron
       allocate(V1(freq-1,3))
       allocate(V2(freq-1,3))
       allocate(VF(freq-1,3))
       do f = 1,20 

          ! Re subdivide two of the edges of the face and get the vertices
          ! (Yes, this is wasteful but it makes code logic easier)
          V1 = divideEdge(icosXYZ(icosF(f,1),:), icosXYZ(icosF(f,2),:), freq)
          V2 = divideEdge(icosXYZ(icosF(f,1),:), icosXYZ(icosF(f,3),:), freq)            

          ! Now divide the edges that join the new vertices along the
          ! subdivided edges.
          do v = 2,fm1
             VF = divideEdge(V1(v,:), V2(v,:), v)
             xyz(offset:(offset+v-2),:) = VF
             offset = offset+v-1
          end do
       end do

    end if

    ! Adjust for radius
    xyz = xyz*radius

    ! Cleanup
    deallocate(V1,V2,VF)

    return
  end function geodome


  ! Same as geodome but returns only a subset of the total points
  ! Needed to avoid memory issues for big problems
  function geodome_subset(radius,freq,i1,i2) result(xyz)
    implicit none

    integer, intent(in) :: freq,i1,i2
    real(WP), intent(in) :: radius
    real(WP), dimension(i2-i1+1,3) :: xyz

    integer :: i,n,f,v,np
    integer, dimension(12,12) :: icosA
    integer, dimension(20,3 ) :: icosF
    integer, dimension(30) :: ivec,jvec
    integer :: nVert,fm1,fm2,offset
    real(WP), allocatable, dimension(:,:) :: V1,V2,VF,buf
    real(WP), dimension(12,3) :: icosXYZ

    if (i2.lt.i1) return

    if (freq.lt.1) then
       print *, 'error: Icosahedron frequency >= 1'
       stop
    end if
    call icosahedron(icosXYZ, icosA, icosF)

    ! Compute number of vertices in geodesic sphere
    fm1 = freq - 1
    fm2 = freq - 2    
    nVert = 12 + 30*fm1 + 20*fm2*fm1/2

    ! Compute number of particles
    np = 2+10*freq**2

    ! Get coordinates
    if (freq.eq. 1) then
       ! Just return the icosahedron
       do i = 1, np
          if (i.ge.i1 .and. i.le.i2) then
             xyz(i-i1+1,:) = icosXYZ(i,:)
          end if
       end do

    else
       ! For all frequencies > 1

       ! Grab the vertices of the icosahedron
       do i = 1, 12
          if (i.ge.i1 .and. i.le.i2) then
             xyz(i-i1+1,:) = icosXYZ(i,:)
          end if
       end do
       offset = 13

       ! Find the nodes that connect edges.  Note we use the upper triangular
       ! part of the adjacency matrix so that we only get each edge once
       ivec = (/ 1,3,1,3,1,3,5,2,4,2,4,7,1,2,5,7,1,2,6,8,3,4,5,7,9,3,4,6,8,10 /)
       jvec = (/ 2,4,5,5,6,6,6,7,7,8,8,8,9,9,9,9,10,10,10,10,11,11,11,11,11,12,12,12,12,12 /)

       ! Generate extra vertices on every edge of the icoshedron
       allocate(buf(freq-1,3))
       do n = 1,30
          buf = divideEdge(icosXYZ(ivec(n),:), icosXYZ(jvec(n),:), freq)
          do i = offset,(offset+fm2)
             if (i.ge.i1 .and. i.le.i2) then
                xyz(i-i1+1 ,:) = buf(i-offset+1,:)
             end if
          end do
          offset = offset + fm1;
       end do

       ! Generate the extra vertices within each face of the icoshedron
       allocate(V1(freq-1,3))
       allocate(V2(freq-1,3))
       allocate(VF(freq-1,3))
       do f = 1,20 

          ! Re subdivide two of the edges of the face and get the vertices
          ! (Yes, this is wasteful but it makes code logic easier)
          V1 = divideEdge(icosXYZ(icosF(f,1),:), icosXYZ(icosF(f,2),:), freq)
          V2 = divideEdge(icosXYZ(icosF(f,1),:), icosXYZ(icosF(f,3),:), freq)            

          ! Now divide the edges that join the new vertices along the
          ! subdivided edges.
          do v = 2,fm1
             VF = divideEdge(V1(v,:), V2(v,:), v)
             do i = offset,(offset+v-2)
                if (i.ge.i1 .and. i.le.i2) then
                   xyz(i-i1+1,:) = VF(i-offset+1,:)
                end if
             end do
             offset = offset+v-1
          end do
       end do

    end if

    ! Adjust for radius
    xyz = xyz*radius

    ! Cleanup
    deallocate(V1,V2,VF,buf)

    return
  end function geodome_subset


  ! Create icosahedron of radius 1
  subroutine icosahedron(xyz,A,F)
    implicit none

    ! ICOSAHEDRON   Generates vertices and graph of icosahedron
    ! Usage: [xyz, A, F] = icosahedron(radius)
    ! Argument: radius - Optional radius of icosahedron. Defaults to 1.
    ! Returns:     xyz - 12x3 matrix of vertex coordinates.
    !                A - Adjacency matrix defining connectivity of vertices.
    !                F - 20x3 matrix specifying the 3 nodes that define each face.

    real(WP), dimension(12,3) :: xyz
    real(WP), dimension(3,12) :: buf3x12,sol
    real(WP), dimension(3,3) :: buf3x3
    real(WP), dimension(3) :: X,Y,Z
    integer, dimension(12,12) :: A
    integer, dimension(20,3 ) :: F
    integer :: i,j

    ! Compute the 12 vertices
    xyz(1 ,:) = (/ 0.0_WP,  1.0_WP,  goldenRatio    /)
    xyz(2 ,:) = (/ 0.0_WP, -1.0_WP,  goldenRatio    /)
    xyz(3 ,:) = (/ 0.0_WP,  1.0_WP, -goldenRatio    /)
    xyz(4 ,:) = (/ 0.0_WP, -1.0_WP, -goldenRatio    /)
    xyz(5 ,:) = (/ 1.0_WP,  goldenRatio,     0.0_WP /)
    xyz(6 ,:) = (/-1.0_WP,  goldenRatio,     0.0_WP /)
    xyz(7 ,:) = (/ 1.0_WP, -goldenRatio,     0.0_WP /)
    xyz(8 ,:) = (/-1.0_WP, -goldenRatio,     0.0_WP /)
    xyz(9 ,:) = (/ goldenRatio,     0.0_WP,  1.0_WP /)
    xyz(10,:) = (/-goldenRatio,     0.0_WP,  1.0_WP /)
    xyz(11,:) = (/ goldenRatio,     0.0_WP, -1.0_WP /)
    xyz(12,:) = (/-goldenRatio,     0.0_WP, -1.0_WP /)

    ! Scale to radius of 1
    xyz = xyz/(sqrt(1+goldenRatio**2))

    ! Define the adjacency matrix
    A(1 ,:) = (/ 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0 /)  
    A(2 ,:) = (/ 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0 /) 
    A(3 ,:) = (/ 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1 /) 
    A(4 ,:) = (/ 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1 /) 
    A(5 ,:) = (/ 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0 /) 
    A(6 ,:) = (/ 1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1 /) 
    A(7 ,:) = (/ 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0 /) 
    A(8 ,:) = (/ 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1 /) 
    A(9 ,:) = (/ 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0 /) 
    A(10,:) = (/ 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1 /) 
    A(11,:) = (/ 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0 /) 
    A(12,:) = (/ 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0 /) 

    ! Define nodes that make up each face
    F(1 ,:) = (/ 1,  2,  9  /)
    F(2 ,:) = (/ 1,  9,  5  /)
    F(3 ,:) = (/ 1,  5,  6  /)
    F(4 ,:) = (/ 1,  6,  10 /)
    F(5 ,:) = (/ 1,  10, 2  /)
    F(6 ,:) = (/ 2,  7,  9  /)
    F(7 ,:) = (/ 9,  7,  11 /)
    F(8 ,:) = (/ 9,  11, 5  /)
    F(9 ,:) = (/ 5,  11, 3  /)
    F(10,:) = (/ 5,  3,  6  /)
    F(11,:) = (/ 6,  3,  12 /)
    F(12,:) = (/ 6,  12, 10 /)
    F(13,:) = (/ 10, 12, 8  /)
    F(14,:) = (/ 10, 8,  2  /)
    F(15,:) = (/ 2,  8,  7  /)
    F(16,:) = (/ 4,  7,  8  /)
    F(17,:) = (/ 4,  8,  12 /)
    F(18,:) = (/ 4,  12, 3  /)
    F(19,:) = (/ 4,  3,  11 /)
    F(20,:) = (/ 4,  11, 7  /)

    ! The icosahedron defined above is oriented so that an edge is at the
    ! top. The following code transforms the vertex coordinates so that a
    ! vertex is at the top to give a more conventional view.
    !
    ! Define coordinate frame where z passes through one of the vertices and
    ! transform the vertices so that one of the vertices is at the 'top'.
    Z = xyz(1,:)    ! Z passes through vertex 1
    X = xyz(2,:)    ! Choose adjacent vertex as an approximate X
    Y = cross_product(Z,X)  ! Y is perpendicular to Z and this approx X
    X = cross_product(Y,Z)  ! Final X is perpendicular to Y and Z

    ! Ensure unit vectors
    X = normalize(X)
    Y = normalize(Y)
    Z = normalize(Z)

    ! Transform points
    buf3x3(1,:) = X
    buf3x3(2,:) = Y
    buf3x3(3,:) = Z
    do i=1,3
       do j=1,12
          buf3x12(i,j) = xyz(j,i)
       end do
    end do
    sol = MATMUL(buf3x3,buf3x12)
    do i=1,12
       do j=1,3
          xyz(i,j) = sol(j,i)
       end do
    end do

    return
  end subroutine icosahedron


  ! Function to divide an edge defined between vertices V1 and V2 into nSeg
  ! segments and return the coordinates of the vertices.
  ! This function divides the *angle* between V1 and V2 rather than the distance.
  function divideEdge(V1, V2, nSeg) result(vert)
    implicit none

    real(WP), dimension(3) :: V1, V2, axis
    real(WP), dimension(nseg-1, 3) :: vert
    real(WP) :: angle,theta,Qw,Qi,Qj,Qk,t2,t3,t4,t5,t6,t7,t8,t9,t10
    real(WP), dimension(4) :: Q
    integer :: nSeg, f

    axis = cross_product(V1,V2)
    angle = atan(sqrt(axis(1)**2+axis(2)**2+axis(3)**2)/dot_product(V1,V2))
    axis = normalize(axis)

    ! Now add appropriate fractions of the edge length to the first node
    do f = 1,nSeg-1
       theta = f*angle/real(nSeg,WP)
       Q(1) = cos(0.5_WP*theta)
       Q(2:4) = sin(0.5_WP*theta)*axis

       Qw = Q(1)
       Qi = Q(2)
       Qj = Q(3)
       Qk = Q(4)

       t2 =   Qw*Qi
       t3 =   Qw*Qj
       t4 =   Qw*Qk
       t5 =  -Qi*Qi
       t6 =   Qi*Qj
       t7 =   Qi*Qk
       t8 =  -Qj*Qj
       t9 =   Qj*Qk
       t10 = -Qk*Qk

       vert(f,1) = 2.0_WP*( (t8 + t10)*V1(1) + (t6 -  t4)*V1(2) + (t3 + t7)*V1(3) ) + V1(1);
       vert(f,2) = 2.0_WP*( (t4 +  t6)*V1(1) + (t5 + t10)*V1(2) + (t9 - t2)*V1(3) ) + V1(2);
       vert(f,3) = 2.0_WP*( (t7 -  t3)*V1(1) + (t2 +  t9)*V1(2) + (t5 + t8)*V1(3) ) + V1(3);

       ! Normalize to unit length
       vert(f,:) = normalize(vert(f,:))
    end do

    return
  end function divideEdge

end module geodesic
