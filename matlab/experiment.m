function [w_f, w_r] = experiment(N_f, N_r)
% Conduct a given number of realizations of the forward and backward nonequilibrium pulling experiments.
%
% This function uses instantaneous switching work measurements for a pair of one-dimensional harmonic oscillators.
%
% ARGUMENTS
%  N_f - number of measurements of the forward process to collect
%  N_r - number of measurements of the reverse process to collect
%
% RETURNS
%  W_f - array of N_f measured work values of the forward process
%  W_r - array of N_r measured work values of the reverse process

% PARAMETERS
K_0 = 1.0; % spring constant of initial state
K_1 = 2.0; % spring constant of final state
x_0 = 0.0; % equilibrium spring position for initial state
x_1 = 0.5; % equilibrium spring position for final state
beta = 1.0; % inverse temperature

% Harmonic oscillator is given by potential
%   U(x) = (K / 2) (x - x_0)^2
% probability density function
%   p(x) = (1/Z) exp[-\beta U(x)]
%
% Z = \int dx \exp[-\beta U(x)]
%   = \int dx \exp[-(x - x_0)^2 / (2 sigma^2)]
%   = sqrt(2 pi) sigma
% sigma^2 = (\beta K)^{-1}

% Draw samples from stationary distributions at inverse temperature beta.
x_f = randnorm(N_f, 1) / sqrt(beta * K_0);
x_r = randnorm(N_r, 1) / sqrt(beta * K_1);

% Compute instantaneous work measurement values.
w_f = (K_1/2)*(x_f-x_1).^2 - (K_0/2)*(x_f-x_0).^2;
w_r = (K_0/2)*(x_r-x_0).^2 - (K_1/2)*(x_r-x_1).^2;

% return work measurements
return

