clear;clc;close all;

# System parameters
E = .5;
Ny = 140;
# Phi values
steps = 50;
phi = linspace(0,1,steps);
delta_phi = [];
# Momentum
ky = 2*pi/Ny;
kx = [];

hold on
for n = -54
  for x=1:steps
    delta_phi(x) = (2*pi/ Ny) * phi(x);
    kx(x) = acos(E/2 - cos(ky * n + delta_phi(x)));
  end
  plot(phi,kx)
end
hold off