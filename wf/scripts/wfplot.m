f1 = squeeze(wf1(1,:));
f1_2d = reshape(f1,10,30);
col_zeros = zeros(13,1);
row_zeros = zeros(1,30);
f1_2d = [row_zeros; f1_2d; row_zeros; row_zeros;];
f1_2d = [f1_2d,col_zeros];
figure;pcolor(real(f1_2d));colormap('jet');axis equal tight;shading flat;colorbar;
