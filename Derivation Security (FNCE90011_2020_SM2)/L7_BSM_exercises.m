% L7_BSM_exercises.m
%
% Solutions to exercises Chapter 19
% 
%   $ Author: Thijs van der Heijden $  
%   $ Revision: 1.0.0 $  
%   $ Date: 2016/08/27 $

%% Housekeeping
clear;
fprintf('\n\n*** In %s ***\n\n',mfilename)
% format long; % use this is you want to display lots of decimal places
%#ok<*NOPTS>

%% 19.6 & 19.7
fprintf('\n** Exercises 19.6 & 19.7 **\n')
S0 = 59
K = 60
T = 44/365
sigma = 0.3
r = 0.03

[c, p] = blsprice(S0, K, r, T, sigma, 0)

%% 19.13
fprintf('\n** Exercise 19.13 **\n')

delta_c = blsdelta(S0, K, r, T, sigma, 0)
gamma_c = blsgamma(S0, K, r, T, sigma, 0)


%% 19.18
fprintf('\n** Exercise 19.18 **\n')
S0 = 11057
K = 11000
T = 45/365
sigma = 0.16
r = 0.05
delta = 0.015

[c, p] = blsprice(S0, K, r, T, sigma, delta)

p_pcp = c + K * exp(-r * T) - S0 * exp(-delta * T)