clear;clc;close all;
load('kx.mat')

# System parameters
E = 3.45;
Ny = 20;
# Phi values
steps = 100;
phi = linspace(0,1,steps);
delta_phi = [];
# Momentum
ky = 2*pi/Ny;
kxs = [];

hold on
for n = -10:10
  for i=1:steps
    delta_phi(i) = (2*pi / Ny) * phi(i);
    kxs(i) = acos(E/2 - cos(ky * n + delta_phi(i)));
  end
  plot(phi, real(kxs), '-')
end
#plot(phis', kx, '--')
hold off