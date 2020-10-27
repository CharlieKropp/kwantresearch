%% Modes of waveguide
clear;clc;close all;
W      = 10;         % Width of the waveguide
L      = 30;           % Length of the waveguide
lambda = 7.0125;        % Wavelength
k0     = 2*pi/lambda;  % Wavenumber
[x,y]  = meshgrid(linspace(0,L,L*10+1),linspace(0,W,W*10+1));

% Closed (reflecting) boundary conditions (along y axis)
Nc  = floor(W/(lambda/2));   % Number of waveguide modes
n   = 1;                     % Allowed range 1 <= n <= Nc
ky  = pi/W*n;
kx  = sqrt(k0^2-ky^2);
psi = sqrt(2/W)*sin(-ky*y).*exp(1i*kx*x);
figure(1);
subplot(3,1,1);pcolor(real(psi)   );axis equal tight;shading flat;colorbar;colormap('jet');
subplot(3,1,2);pcolor(imag(psi)   );axis equal tight;shading flat;colorbar;colormap('jet');
subplot(3,1,3);pcolor( abs(psi).^2);axis equal tight;shading flat;colorbar;colormap('jet');
 
% Periodic boundary conditions (along y axis)
Np  = 2*floor(W/lambda)+1;   % Number of waveguide modes

n1   = -1;                     % Allowed range -(Np-1)/2 <= n <= (Np-1)/2
ky1  = 2*pi/W*n1;
kx1  = sqrt(k0^2-ky1^2);
psi1 = sqrt(1/W)*exp(1i*ky1*y).*exp(1i*kx1*x);

figure(2);
subplot(3,1,1);pcolor(real(psi1)   );axis equal tight;shading flat;colorbar;colormap('jet');
subplot(3,1,2);pcolor(imag(psi1)   );axis equal tight;shading flat;colorbar;colormap('jet');
subplot(3,1,3);pcolor(real(psi1)   );axis equal tight;shading flat;colorbar;colormap('jet');
