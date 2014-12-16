% Plot
figure(2);
clf;

colormap(hot);
cmap = [linspace(1,0,100)' linspace(0,0,100)' linspace(0,0,100)'; linspace(0,0,100)' linspace(0,0,100)' linspace(0,1,100)'];
colormap(cmap);

% find combinations of Nf,Nr that we've tested
indices = find(mask == 1);

% one sigma
c = 1;
ci = cis(c); % expected fraction that should fall in this confidence interval

% find scale limits
X = squeeze(Pfrs_asymptotic(:,:,c) / ci);
Y = squeeze(Pfrs_bayesian(:,:,c) / ci);
max_dev = max(max([X(indices)-1 Y(indices)-1 1-X(indices) 1-Y(indices)]))

subplot(2,3,1);
show_deviation(X, [1-max_dev 1+max_dev], mask, cmap);
title(sprintf('ABAR rel fraction for %.2f CI', ci));

subplot(2,3,4);
show_deviation(Y, [1-max_dev 1+max_dev], mask, cmap);
title(sprintf('BBAR rel fraction for %.2f CI', ci));

% TWO SIGMA
c = 2;
ci = cis(c); % expected fraction that should fall in this confidence interval

% find scale limits
X = squeeze(Pfrs_asymptotic(:,:,c) / ci);
Y = squeeze(Pfrs_bayesian(:,:,c) / ci);
max_dev = max(max([X(indices)-1 Y(indices)-1 1-X(indices) 1-Y(indices)]))

subplot(2,3,2);
show_deviation(X, [1-max_dev 1+max_dev], mask, cmap);
title(sprintf('ABAR rel fraction for %.2f CI', ci));

subplot(2,3,5);
show_deviation(Y, [1-max_dev 1+max_dev], mask, cmap);
title(sprintf('BBAR rel fraction for %.3f CI', ci));

% BIAS

% find scale limits
X = squeeze(BAR_ML_bias_fr / abs(true_df));
Y = squeeze(BAR_mean_bias_fr / abs(true_df));
max_dev = max(max([X(indices)-1 Y(indices)-1 1-X(indices) 1-Y(indices)]))
max_dev = 2;

subplot(2,3,3);
show_deviation(X, [-max_dev +max_dev], mask, cmap);
title(sprintf('ABAR rel bias of ML'));

subplot(2,3,6);
show_deviation(Y, [-max_dev +max_dev], mask, cmap);
title(sprintf('BBAR rel bias of posterior mean'));

% print
filename = '../plots/bar-data.eps';
print('-depsc', filename);
unix(sprintf('epstopdf %s', filename));
