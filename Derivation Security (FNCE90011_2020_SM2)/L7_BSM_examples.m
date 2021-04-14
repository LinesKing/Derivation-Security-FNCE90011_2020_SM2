% L7_BSM_examples.m
%
% Examples of applying the Black-Scholes-Merton option pricing formula
% 
%   $ Author: Thijs van der Heijden $  
%   $ Revision: 1.0.0 $  
%   $ Date: 2016/09/02 $

%% Housekeeping
clear;
fprintf('\n\n*** In %s ***\n\n',mfilename)
% format long; % use this is you want to display lots of decimal places
%#ok<*NOPTS>

% Plot set up
width = 8;     % Width in inches
height = 4.5;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 11;      % Fontsize
lw = 1.8;      % LineWidth
msz = 8;       % MarkerSize

%% Distribution of future stock prices in a binomial tree 
% Increasing the number of time-steps in the tree while keeping the total
% horizon constant. This illustrates how the binomial tree converges to a normal
% distribution.
fprintf('\n** Distribution of future stock prices in binomial tree **\n')
sigma = 0.2
T = 1
dt = [1 1/12 1/52 1/252]
r = 0.05

figure(1);
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties
for i = 1:length(dt)
    nT = T / dt(i)
    U = exp(sigma * sqrt(dt(i)))
    D = 1 / U
    % risk-neutral probability of up move in the tree
    pie = (exp(r * dt(i)) - D) / (U - D) 
    % Potential stock prices as number of upward steps taken
    log_return = (nT:-2:-nT) * log(U);
    S = exp(log_return);
    prob = binopdf(nT:-1:0, nT, pie);
    
    subplot(2, 2, i)
    %stem(log_return, prob, 'LineWidth', lw);
    %xlabel('Cumulative log return')
    stem(S, prob, 'LineWidth', lw);
    xlabel('Stock price at maturity')
    ylabel('Probability');
    title(['Number of steps per year: ', num2str(1/dt(i)) ', \pi: ', ...
        num2str(pie)])
    axis([0 3 0 inf])
end

set(gcf, 'InvertHardcopy', 'on');
set(gcf, 'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf, 'PaperPosition', myfiguresize);

saveas(gcf, './figures/binomial_tree_growth.png')
close(gcf)


%% Call option price
fprintf('\n** Call option price in Black-Scholes-Merton model **\n')
S0 = 65
K = 70
T = 3/12
r = 0.0953
sigma = 0.35

d1 = (log(S0/K) + (r + 0.5 * sigma^2) * T) / (sigma * sqrt(T))
d2 = d1 - sigma * sqrt(T)

Nd1 = normcdf(d1)
Nd2 = normcdf(d2)

c = S0 * Nd1 - K * exp(-r*T) * Nd2

p = K * exp(-r*T) * normcdf(-d2) - S0 * normcdf(-d1)

p_pcp = c + K * exp(-r*T) - S0



%% Put price
fprintf('\n** Put price example for both Black-Scholes and binomial tree **\n')
% Black-Scholes
K = 105
S0 = 100
r = 0.05
T = 9/12
sigma = 0.15

d1 = (log(S0/K) + (r + 0.5 * sigma^2) * T) / (sigma * sqrt(T))
d2 = d1 - sigma * sqrt(T)

p = K * exp(-r * T) * normcdf(-d2) - S0 * normcdf(-d1)

% Binomial tree
% Even though we use the CRR tree here, we still need to use the 'binpriceJR' function
% as we are considering a European put option, which is not implemented in the standard
% 'binprice' function.
dt = 1/4
U = exp(sigma * sqrt(dt))
D = 1 / U
R = exp(r * dt) - 1
pie = (1 + R - D) / (U - D)

[pr, p] = binpriceJR(S0, K, r, T, dt, sigma, 2, 0, 0, 0, U, D)



%% Practice question 1
fprintf('\n** Practice question 1 **\n')
S0 = 10
K = 10.50
sigma = 0.2
T = 3/12
r = 0.02

d1 = (log(S0/K) + (r + 0.5 * sigma^2) * T) / (sigma * sqrt(T))
d2 = d1 - sigma * sqrt(T)

Nd1 = normcdf(d1)
Nd2 = normcdf(d2)

[c, p] = blsprice(S0, K, r, T, sigma, 0)

%% Practice question 2, 3 & 4
fprintf('\n** Practice questions 2, 3 & 4 **\n')
S0 = 10
K = 10.50
sigma = 0.2
T = 3/12
r = 0.02
div = 0.50
div_t = 2/12

S0t = S0 - div * exp(-r * div_t)
d1 = (log(S0t/K) + (r + 0.5 * sigma^2) * T) / (sigma * sqrt(T))
d2 = d1 - sigma * sqrt(T)

Nd1 = normcdf(d1)
Nd2 = normcdf(d2)

[c, p] = blsprice(S0t, K, r, T, sigma, 0)

p_pcp = c + K * exp(-r * T) - S0t



dt = 1/12
U = exp(sigma * sqrt(dt))
D = 1/U
pie = (exp(r * dt) - D) / (U - D)

% First part of the stock price tree, before dividend payment
S = triu(S0 * U.^[0:2; 0 0 1; 0 0 0] .* D.^[0 0 0; 0 1 1;0 0 2])

% Second part of the stock price tree, dividend payment and after.
% We use the kronecker product here to expand the tree to make room for future nodes now that the 
% tree no longer recombines.
S2 = [kron(S(:,3)-div, [1;0]) kron(S(:,3)-div, [U; D])]

% European put option in second part of tree. First determine final payoff. Then move backwards in time.
% The second last period, just after the dividend payment, there are three values again, but next period,
% because the tree does not recombine anymore, there will be six possible values. The tree is now set up
% such that we take the top two elements of the final payout of 'c' and multiply them with pie and 1-pie,
% respectively, then discount using the risk free rate to get the call value one period earlier in the 
% top node (S = 170). Similarly for the other nodes.
p2 = max(K - S2(:,2), 0);
p2 = [zeros(6, 1), p2];
p2(1:2:5, 1) = kron(eye(3), [pie 1-pie]) * p2(:,2) * exp(-r * dt)

% Move back in tree of part one using the option values in the second part
p = [zeros(3, 2), p2(1:2:5,1)];
for n = 2:-1:1
    k = 1:n;
    p(k, n) = (pie*p(k,n+1) + (1-pie)*p(k+1,n+1))*exp(-r*dt);
end
p

% American put option. Early exercise may be optimal directly after the dividend has been
% paid.
P = [zeros(3, 2), max(p2(1:2:5, 1), K - S2(1:2:5, 1))];
for n = 2:-1:1
    k = 1:n;
    discopt = (pie*P(k,n+1) + (1-pie)*P(k+1,n+1))*exp(-r*dt);
    P(k, n) = max(discopt, K - S(k, n));
end
P


%% Practice question 5
fprintf('\n** Practice question 5 **\n')
S0 = 5550
K = 5500
sigma = 0.25
T = 4/12
r = 0.02
q = 0.03 % dividend yield

F = S0 * exp((r - q) * T) % Forward price
d1 = (log(F/K) + 0.5 * sigma^2 * T) / (sigma * sqrt(T))
d2 = d1 - sigma * sqrt(T)

Nd1 = normcdf(d1)
Nd2 = normcdf(d2)

% Compute option price using alternative formula using forward price directly
c_alt = F * exp(-r*T) * Nd1 - K*exp(-r*T) * Nd2

% Compute option prices using Black-Scholes formula with dividend yield argument
[c, p] = blsprice(S0, K, r, T, sigma, q)

p_pcp = c + K * exp(-r*T) - S0*exp(-q*T)




