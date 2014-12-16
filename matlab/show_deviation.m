function show_deviation(X, limits, mask, cmap)

[N_f, N_r] = size(X);

imagesc(0:N_f, 0:N_r, X, limits);
axis square;
colorbar;
xlabel('N_r');
ylabel('N_f');
set(gca, 'YDir', 'normal');

return

[N_f, N_r] = size(X);

Xgray = mat2gray(X, [1-max_dev 1+max_dev]);
Xi = gray2ind(Xgray, 256);
Xrgb = ind2rgb(Xi, cmap);
Xoverlay = imoverlay(Xrgb, ~mask, [1 1 1]);
image(0:N_f, 0:N_r, Xoverlay);

axis square;
%caxis([1-max_dev 1+max_dev]);
colorbar('CLim', [1-max_dev 1+max_dev]);

xlabel('N_r');
ylabel('N_f');

set(gca, 'YDir', 'normal');


return
