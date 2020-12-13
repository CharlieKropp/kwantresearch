clear;clc;close all;
 
load('wf.mat');

W = 100;
n = 17;

wf0modes = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
actualmodes = [0, 1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 6, -6, 7, -7]
modemap = containers.Map(actualmodes, wf0modes)


figure;pcolor(real([psi1;psi1]));axis equal tight;shading flat;colormap(jet);colorbar;
figure;pcolor(real([psi2;psi2]));axis equal tight;shading flat;colormap(jet);colorbar;
figure;pcolor( abs([psi1;psi1]));axis equal tight;shading flat;colormap(jet);colorbar;
figure;pcolor( abs([psi2;psi2]));axis equal tight;shading flat;colormap(jet);colorbar;
%--------------------------------------------------------------------------
input('Type to close')