clear;clc;close all;

load('wf.mat');

W = 100
n = 17

figure;pcolor(real(te));axis equal tight;shading flat;colormap(jet);colorbar;
figure;pcolor(real(tp));axis equal tight;shading flat;colormap(jet);colorbar;

for x=1:17  
psi1=reshape(wfe(x,:),W,[]);
figure;pcolor(real([psi1;psi1]));axis equal tight;shading flat;colormap(jet);colorbar;
endfor