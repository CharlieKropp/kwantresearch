clear;clc;close all;
load('tau_ring.mat');
load('cvsb_ring.mat');
n = 3
phis = 100
phis_c = linspace(0,1,100)
N = 108
allcvsb = reshape(allcvsb, [phis, nrlz])
taus1 = reshape(taus, [nL,phis, 5]);
for i=1:5
figure;pcolor(log10(taus1(:,:,i)));axis equal tight;shading flat;colormap(jet);colorbar;
endfor
figure;plot(phis_c,allcvsb(:,n))
ipr = squeeze(mean(sum(taus1,1),3)).^2./squeeze(mean(sum(taus1.^2,1),3));
ipr_n = sum(taus1(1:nL,:,:),1).^2./sum(taus1(1:nL,:,:).^2,1);
#figure;semilogy(squeeze(mean(taus1(:,1,:),3)));
if n != 0
figure;plot(phis_c, ipr_n(:,:,n))
figure;plot(phis_c, allcvsb(:,n));

endif