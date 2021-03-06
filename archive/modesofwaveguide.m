%% Modes of waveguide
clear;clc;close all;
W      = 10.5;         % Width of the waveguide
L      = 50;           % Length of the waveguide
lambda = 10.0;         % Wavelength
k0     = 2*pi/lambda;  % Wavenumber
[x,y]  = meshgrid(linspace(0,L,L*10+1),linspace(0,W,W*10+1));
 
% Closed (reflecting) boundary conditions (along y axis)
Nc  = floor(W/(lambda/2));   % Number of waveguide modes
n   = 1;                     % Allowed range 1 <= n <= Nc
ky  = pi/W*n;
kx  = sqrt(k0^2-ky^2);
psi = sqrt(2/W)*sin(ky*y).*exp(1i*kx*x);
 
% Periodic boundary conditions (along y axis)
Np  = 2*floor(W/lambda)+1;   % Number of waveguide modes
n   = 1;                     % Allowed range -(Np-1)/2 <= n <= (Np-1)/2
ky  = 2*pi/W*n;
kx  = sqrt(k0^2-ky^2);
psi = sqrt(1/W)*exp(1i*ky*y- 5*i).*exp(1i*kx*x);
figure;pcolor(real(psi));axis equal tight;shading flat;colormap(jet);colorbar;
