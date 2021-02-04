clear;clc;close all;
folder = 'home/charlie/kwantresearch/';
load('wf.mat');
W = 20;
n = rows(wf0);

for x=1:n
psi1=reshape(wf0(x,:),W,[]);
h= figure;pcolor(real([psi1;psi1]));axis equal tight;shading flat;colormap(jet);colorbar;

    
  % Then save it
saveas(h, sprintf('/home/charlie/kwantresearch/realB/fig%d.png',x))
end