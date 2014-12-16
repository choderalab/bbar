figure(1);
set(gcf,'position',[0 0 10.0 3.25]);
clf;
nstd = 3;
npoints = 1000;
xmin = min(x_0 - nstd*sigma_0, x_1 - nstd*sigma_1);
xmax = max(x_0 + nstd*sigma_0, x_1 + nstd*sigma_1);
xpoints = linspace(xmin, xmax, npoints);
% Plot potential
subplot(1,3,1);
U_0 = (K_0/2)*(xpoints-x_0).^2;
U_1 = (K_1/2)*(xpoints-x_1).^2;
plot(xpoints, U_0, 'r-', xpoints, U_1, 'b-');
xlabel('x');
ylabel('U(x) / kT');
title('potential');
oldaxis = axis;
axis([xmin xmax oldaxis(3:4)]);
% plot probability distribution
subplot(1,3,2);
P_0 = (sqrt(2*pi)*sigma_0)^(-1) * exp(-(xpoints-x_0).^2/(2*sigma_0^2));
P_1 = (sqrt(2*pi)*sigma_1)^(-1) * exp(-(xpoints-x_1).^2/(2*sigma_1^2));
plot(xpoints, P_0, 'r-', xpoints, P_1, 'b-');
xlabel('x');
ylabel('p(x)');
title('pdf');
oldaxis = axis;
axis([xmin xmax oldaxis(3:4)]);
% plot work distribution by sampling
subplot(1,3,3);
samplesize = 100000;
nbins = 100;
x_f = sigma_0 * randn([samplesize,1]) + x_0;
w_f = WF(x_f);
x_r = sigma_1 * randn([samplesize,1]) + x_1;
w_r = WR(x_r);
wmin = min([w_f; -w_r]);
wmax = max([w_f; -w_r]);
wbins = linspace(wmin, wmax, nbins);
pw_f = hist(w_f, wbins);
pw_r = hist(-w_r, wbins);
% color in minimum
woverlap = min(pw_f, pw_r);
hold on
fill([wbins(1) wbins], [0 woverlap], 0.8 * [1 1 1]);
plot(wbins, pw_f, 'r-', 'LineWidth', 2);
plot(wbins, pw_r, 'b-', 'LineWidth', 2);
oldaxis = axis;
axis([wmin wmax oldaxis(3:4)]);
%hold on
%fill(xpw_f, pw_f, 0.8*[1 1 1]);
%fill(xpw_r, pw_r, 0.6*[1 1 1]);
xlabel('work / kT');
ylabel('p(w)');
title('work distributions');
legend('overlap', 'W_f', '-W_r');
clear w_f w_r;

filename = '../plots/work-distribution.eps';
print('-depsc', filename);
unix(sprintf('epstopdf %s', filename));

