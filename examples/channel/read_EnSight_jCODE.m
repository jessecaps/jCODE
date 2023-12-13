% Script to read EnSight file
clear
clc

% Read in the geometry
fid = fopen('ensight-3D/geometry','r');
buf = fread(fid, 80, 'char');
buf = fread(fid, 80, 'char');
buf = fread(fid, 80, 'char');
buf = fread(fid, 80, 'char');
buf = fread(fid, 80, 'char');
buf = fread(fid, 80, 'char');
ibuf = fread(fid, 1, 'int');
buf = fread(fid, 80, 'char'); 
buf = fread(fid, 80, 'char'); 
nx = fread(fid, 1, 'int'); 
ny = fread(fid, 1, 'int'); 
nz = fread(fid, 1, 'int'); 
rbuf = fread(fid, nx*ny*nz, 'real*4');
x3d=reshape(rbuf,nx,ny,nz);
rbuf = fread(fid, nx*ny*nz, 'real*4');
y3d=reshape(rbuf,nx,ny,nz);
rbuf = fread(fid, nx*ny*nz, 'real*4');
z3d=reshape(rbuf,nx,ny,nz);
fclose(fid);

% 1D coordinates if Cartesian
x=x3d(:,1,1);
y=y3d(1,:,1);
z=z3d(1,1,:);

% Read EnSight file
fid = fopen('ensight-3D/VELOCITY/VELOCITY.000409','r');
buf = fread(fid, 80, 'char');
buf = fread(fid, 80, 'char');
ibuf = fread(fid, 1, 'int');
buf = fread(fid, 80, 'char'); 
rbuf = fread(fid, nx*ny*nz, 'real*4'); 
fclose(fid);
U=reshape(rbuf,nx,ny,nz);

% Plot
figure(1)
plot(U(1,:,1),y,'-k')
