% Script to read EnSight file
clear; clc; close all

% Analytic solution
tau = 0.2;
analytic = analytic_sod(tau);

% Directory
dir = 'ensight-3D';

% Read in the geometry
fid = fopen([dir '/geometry'],'r');
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

% Get 1D coordinates
x = x3d(:,1,1);
y = y3d(:,1,1);
z = z3d(:,1,1);

% Read EnSight file
var = 'DENSITY';
fid = fopen([dir '/' var '/' var '.000002'],'r');
buf = fread(fid, 80, 'char');
buf = fread(fid, 80, 'char');
ibuf = fread(fid, 1, 'int');
buf = fread(fid, 80, 'char'); 
rbuf = fread(fid, nx*ny*nz, 'real*4'); 
fclose(fid);
rho=reshape(rbuf,nx,ny,nz);

 var = 'BULK_VISC_SHOCK';
 fid = fopen([dir '/' var '/' var '.000002'],'r');
 buf = fread(fid, 80, 'char');
 buf = fread(fid, 80, 'char');
 ibuf = fread(fid, 1, 'int');
 buf = fread(fid, 80, 'char');
 rbuf = fread(fid, nx*ny*nz, 'real*4');
 fclose(fid);
diss=reshape(rbuf,nx,ny,nz);

% Plot density
figure()
plot(x-0.5,rho(:,1,1),'-b','linewidth',4)
hold all
plot(analytic.x,analytic.rho,'--k','linewidth',1)
hold off
%title('$\tau=0.2$','interpreter','latex','fontsize',24)
xlabel('$x$','interpreter','latex','fontsize',24)
ylabel('$\rho$','interpreter','latex','fontsize',24)
set(get(gca, 'yLabel'), 'Rotation',0);
set(gca, ...
    'box', 'on',...
    'tickdir', 'in',...
    'ticklength',[.015 .015],...
    'xminortick','on',...
    'yminortick','on',...
    'linewidth',1,...
    'fontsize',24,'FontName','Times New Roman');
xlim([-0.5 0.5])
ylim([0 1.2])
set(gcf, 'Color',[1,1,1]);
l=legend('jCODE','Analytic','location','NorthEast');

set(l,'Interpreter','latex','EdgeColor',[1 1 1],'fontsize',16,'color','none')
print('-dpng', 'density');

%Plot dissipation
figure()
diss0=max(diss(:,1,1));
plot(x-0.5,diss(:,1,1),'-k','linewidth',2)
title('$\tau=0.2$','interpreter','latex','fontsize',24)
xlabel('$x$','interpreter','latex','fontsize',24)
ylabel('Artificial viscosity','interpreter','latex','fontsize',24)
set(get(gca, 'yLabel'), 'Rotation',90);
xlim([-0.5 0.5])
%ylim([0 1.2])
set(gcf, 'Color',[1,1,1]);