clear;clc;close all;
load('wfbranch.mat');
W = 20
n = rows(wf0)

for x=1:n
psi1=reshape(wf0(x,:),W,[]);
row1 = psi1(1, :)
psi1(1, :) = []
psi1(W, :) = row1
h=figure;pcolor(real([psi1;psi1]));axis equal tight;shading flat;colormap(jet);colorbar;
saveas(h, sprintf('/home/charlie/kwantresearch/branchcut/fig%d.png',x))
endfor
