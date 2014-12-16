% Generate reference 'confidence curve validation plot' showing what confidence verification
% plot would look like with various amounts of under- or over-estimation of uncertainty.

% Clear plot
clf;

% Generate plot box.
axis([0 1 0 1]);
axis square;
hold on;

% Generate normal deviations.
%sigmas = [1, 2, 10, 100]; % standard deviations to plot
sigmas = [100, 10, 2, 1, 1/2, 1/10, 1/100]; % standard deviations to plot
npoints = 1000; % number of points in each line
desired_cis = linspace(0, 1, npoints); % desired confidence intervals
handles = [];
for sigma = sigmas
  % Generate plot for under- and over-estimates of normal posterior by 'sigma' standard deviations.

  actual_cis = zeros([1,npoints]);
  thetas_lower = norminv(0.5 - desired_cis/2);
  thetas_upper = norminv(0.5 + desired_cis/2);
  actual_cis = normcdf(thetas_upper, 0, sigma) - normcdf(thetas_lower, 0, sigma);

  handle = plot(desired_cis, actual_cis, 'k-');
  handles = [handles handle];

  if sigma == 1
    set(handle, 'LineWidth', 2);
  end
end

%legend('100 \sigma', '10 \sigma', '2 \sigma', '\sigma', '1/2 \sigma', '1/10 \sigma', '1/100 \sigma');

% Label curves.
text(0.0250, 0.9701, '100 \sigma');
text(0.1506, 0.8562, '10 \sigma');
text(0.3520, 0.6152, '2 \sigma');
text(0.4370, 0.5044, '\sigma');
text(0.5096, 0.3974, '1/2 \sigma');
text(0.6646, 0.1418, '1/10 \sigma');
text(0.7319, 0.0395, '1/100 \sigma');

% Label plot
xlabel('desired confidence level');
ylabel('actual confidence level');
set(gca, 'XTick', [0 0.5 1]);
set(gca, 'YTick', [0 0.5 1]);
axis square;
axis([0 1 0 1]);
set(gca, 'box', 'on');
%set(gca, 'FontSize', fontsize);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write plot.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filename = sprintf('../plots/confidence-level-reference.eps');
%print('-depsc', filename);
exportfig(gcf, filename, 'width', 3.25)
unix(sprintf('epstopdf %s', filename));

