clear;clc;close all;
 
load('wfW48.mat');
%W is width. Width must be divisible by 4 and 8.
%n is desired mode to plot
W = 48;
n = 3;

wf0modes = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
actualmodes = [0, 1, -1, 2, -2, 3, -3, 4, -4, 5, -5]
modemap = containers.Map(actualmodes, wf0modes)

psi1=reshape(wf0(modemap(n),:),W,[]);psi2=reshape(wf0(modemap(n)+1,:),W,[]);

figure;pcolor(real([psi1;psi1]));axis equal tight;shading flat;colormap(jet);colorbar;
figure;pcolor(real([psi2;psi2]));axis equal tight;shading flat;colormap(jet);colorbar;
figure;pcolor( abs([psi1;psi1]));axis equal tight;shading flat;colormap(jet);colorbar;
figure;pcolor( abs([psi2;psi2]));axis equal tight;shading flat;colormap(jet);colorbar;
%--------------------------------------------------------------------------
W=size(psi1,1);

if (mod(W, 2) == 0) && (mod(W,4) == 0) && (mod(n,2) == 1)
    p11=psi1(1,1);p12=psi1(1+W/4,1);p21=psi2(1,1);p22=psi2(1+W/4,1)
elseif (mod(W, 2) == 0) && (mod(W,4) == 0) && (mod(n,2) == 0)
    p11=psi1(1,1);p12=psi1(1+W/8,1);p21=psi2(1,1);p22=psi2(1+W/8,1)
end

phi     = atan(abs( (p11-p12/1i) / (p11+p12/1i) )); 
theta   = angle(    (p11-p12/1i) / (p11+p12/1i) );
theta1  = angle(    (p11+p12/1i)/2 );
theta2  = angle(   -(p21+p22/1i)/2 );
psi1p   = ( psi1*exp(-1i*theta1)*cos(phi) - psi2*exp(-1i*theta2)*sin(phi) ); 
psi2p   = ( psi1*exp(-1i*theta1)*sin(phi) + psi2*exp(-1i*theta2)*cos(phi) )*exp(-1i*theta);

figure;pcolor(real([psi1p;psi1p]));axis equal tight;shading flat;colormap(jet);colorbar;
figure;pcolor(real([psi2p;psi2p]));axis equal tight;shading flat;colormap(jet);colorbar;
figure;pcolor( abs([psi1p;psi1p]));axis equal tight;shading flat;colormap(jet);colorbar;
figure;pcolor( abs([psi2p;psi2p]));axis equal tight;shading flat;colormap(jet);colorbar;