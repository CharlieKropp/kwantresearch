clear;clc;close all;
load('tau.mat');
load('cvsb.mat');
n = 0
phis = 100
phis_c = linspace(0,1,100)
N = 108
taus1 = reshape(taus, [nL+1,phis, nrlz]);
x1 = acosh(taus1.^(-.5));
allcvsb = reshape(allcvsb, [phis,nrlz]);
###
#for i=1:10
#figure;pcolor(log10(taus1(:,:,i)));axis equal tight;shading flat;colormap(jet);colorbar;
#endfor
##

ipr = squeeze(mean(sum(taus1,1),3)).^2./squeeze(mean(sum(taus1.^2,1),3))./1.5;
ipr_n = sum(taus1(1:108,:,:),1).^2./sum(taus1(1:108,:,:).^2,1);
mean_allcvsb = mean(allcvsb,2)
#figure;plot(phis_c, ipr);
figure;plot(phis_c,mean_allcvsb(:,1))
if n != 0
figure;plot(phis_c, ipr_n(:,:,n));
figure;plot(phis_c, allcvsb(:,n));
endif
avg_rlz = squeeze(mean(taus1(1:27,:,:),3));
avg = mean(avg_rlz(:,:),2);
for i=1:100
  avg_rlz(:,i) = avg_rlz(:,i) - avg;
endfor
figure;plot(avg_rlz')
#figure;plot(squeeze(mean(taus1(1:20,:,:),3))')

#figure;plot(squeeze(mean(x1(:,:,:),3))')
#figure;semilogy(squeeze(mean(taus1(:,1,:),3)));