% PracticeQuestions.m
%
% For the practice questions of FNCE90011, 2016 S2.
%
%
%
%   $ Author: Thijs van der Heijden $  
%   $ Revision: 1.0 $  
%   $ Date: 2016/10/15 $

%% Housekeeping
clear;
fprintf('\n\n*** In %s ***\n\n',mfilename)
% format long; % use this is you want to display lots of decimal places
%#ok<*NOPTS>
mkdir('./figures')
rng(0); % initialize random number generator for reproducibility

%% Plot set up
% Reference:
% - https://dgleich.github.io/hq-matlab-figs/

% Defaults for this blog post
width = 8;     % Width in inches
height = 4.5;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 11;      % Fontsize
lw = 1.8;      % LineWidth
msz = 8;       % MarkerSize

%% Question 1 (Arbitrage / option strategies)
fprintf('\n\n* Question 1 (Arbitrage/option strategies) * \n')
S0 = 40;
Opt0 = [30	1 	10.25	0.1;
30	2	10.5	0.2;
35	1	5.5	0.25;
35	2	5.75	0.50;
40	1	1.5	1.5;
40	2	2.25	2;
45	1	0.25	5;
50	1	0.25	10]; % option prices [strike maturity (months) call price put price]

nOpt = size(Opt0,1);

% b)
SB.w = -[0 0 1 0 -2 0 1 0]; % weights of the call options in the short butterfly

ST = 20:1:50;
SB.Payoff = SB.w*max(repmat(ST,nOpt,1)-repmat(Opt0(:,1),1,length(ST)),0);


figure(1);
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties

plot(ST, SB.Payoff, 'LineWidth', lw);
xlabel('Stock price at maturity')
ylabel('Payoff of short butterfly strategy')

% Here we preserve the size of the image when we save it.
set(gcf, 'InvertHardcopy', 'on');
set(gcf, 'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf, 'PaperPosition', myfiguresize);

saveas(gcf, './figures/Q1b.png')
close(gcf);



SB.PriceC = SB.w*Opt0(:,3);
SB.PriceP = SB.w*Opt0(:,4)

% h)
BearSpread.w3040 = [0 -1 0 0 0 1 0 0];
BearSpread.w3540 = [0 0 0 -1 0 1 0 0];
% BearSpread.wC3040 = [0 -1 0 0 0 1 0 0];
% BearSpread.wC3540 = [0 0 0 -1 0 1 0 0];

BearSpread.P3040Payoff = BearSpread.w3040*max(repmat(Opt0(:,1),1,length(ST))-repmat(ST,nOpt,1),0);
BearSpread.P3540Payoff = BearSpread.w3540*max(repmat(Opt0(:,1),1,length(ST))-repmat(ST,nOpt,1),0);
BearSpread.C3040Payoff = BearSpread.w3040*max(repmat(ST,nOpt,1)-repmat(Opt0(:,1),1,length(ST)),0);
BearSpread.C3540Payoff = BearSpread.w3540*max(repmat(ST,nOpt,1)-repmat(Opt0(:,1),1,length(ST)),0);

BearSpread.P3040Cost = BearSpread.w3040*Opt0(:,4);
BearSpread.P3540Cost = BearSpread.w3540*Opt0(:,4);
BearSpread.C3040Cost = BearSpread.w3040*Opt0(:,3);
BearSpread.C3540Cost = BearSpread.w3540*Opt0(:,3)


figure(1);
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties

hh = plot(ST, [BearSpread.P3040Payoff - BearSpread.P3040Cost; ...
    BearSpread.P3540Payoff - BearSpread.P3540Cost; ...
    BearSpread.C3040Payoff - BearSpread.C3040Cost; 
    BearSpread.C3540Payoff - BearSpread.C3540Cost], 'LineWidth', lw);
hh(2).LineStyle = '--';
hh(3).LineStyle = '-.';
hh(4).LineStyle = ':';
legend('P3040', 'P3540', 'C3040', 'C3540')
xlabel('Stock price at maturity')
ylabel('Payoff of bearish vertical spread')

% Here we preserve the size of the image when we save it.
set(gcf, 'InvertHardcopy', 'on');
set(gcf, 'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf, 'PaperPosition', myfiguresize);

saveas(gcf, './figures/Q1h.png')
close(gcf);





%% Question 2 (Structured Product)
fprintf('\n\n* Question 2 (Structured Product) * \n')
SP.G0 = 1300;
SP.rf = 0.04;
SP.sig = 0.3;
SP.T = 2;
SP.KC = [1300 1500]';
SP.wC = [1 -1];
SP.KP = 1000;
SP.wP = 1;
SP.GT = 900:1600;

nC = length(SP.KC);
nP = length(SP.KP);
nGT = length(SP.GT);

SP.Payout = SP.G0 + SP.wC*max(repmat(SP.GT,nC,1)-repmat(SP.KC,1,nGT),0) + SP.wP*max(repmat(SP.KP,1,nGT)-repmat(SP.GT,nP,1),0);

figure(1);
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties

plot(SP.GT, SP.Payout, 'LineWidth', lw);
xlabel('Gold Price at Maturity ($/oz)')
ylabel('Note Payout ($)')

% Here we preserve the size of the image when we save it.
set(gcf, 'InvertHardcopy', 'on');
set(gcf, 'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf, 'PaperPosition', myfiguresize);

saveas(gcf, './figures/Q2payout.png')
close(gcf);



SP.PrBond = SP.G0*exp(-SP.rf*SP.T);
[~,SP.PrPsig25] = blsprice(SP.G0,SP.KP,SP.rf,SP.T,0.25,0);

SP.PrC = blsprice(SP.G0,SP.KC,SP.rf,SP.T,SP.sig,0);
[~,SP.PrP] = blsprice(SP.G0,SP.KP,SP.rf,SP.T,SP.sig,0);

d1 = (log(SP.G0/SP.KP) + (SP.rf + 0.5*0.25^2)*SP.T)./(0.25*sqrt(SP.T));
d2 = round(100*(d1 - 0.25*sqrt(SP.T)))/100;
d1 = round(100*d1)/100;

PutPriceRounded = SP.KP*exp(-SP.rf*SP.T)*normcdf(-d2) - SP.G0*normcdf(-d1)

SP.PricePF = SP.PrBond + SP.wC*SP.PrC + SP.wP*SP.PrP;

[SP.PrCTable,SP.PrPTable] = blsprice(SP.G0,900:100:1600,SP.rf,SP.T,SP.sig,0);
SP



%% Question 3 (Binomial model) (Hull 2011, 12.16)
fprintf('\n\n* Question 3 (Binomial model) * \n')

%% Question 4 (Black-Scholes/Greeks) (Hull 2011, 14.29 + addition)
fprintf('\n\n* Question 4 (Black-Scholes/Greeks) * \n')
Q4.S0 = 30;
Q4.K = 29;
Q4.rf = 0.05;
Q4.sig = 0.25;
Q4.T = 4/12;

[Q4.PrC, Q4.PrP] = blsprice(Q4.S0,Q4.K,Q4.rf,Q4.T,Q4.sig,0);

[Q4.deltaC, Q4.deltaP] = blsdelta(Q4.S0,Q4.K,Q4.rf,Q4.T,Q4.sig,0);

[Q4.ThetaC, Q4.ThetaP] = blstheta(Q4.S0,Q4.K,Q4.rf,Q4.T,Q4.sig,0);

d1 = (log(Q4.S0/Q4.K) + (Q4.rf + 0.5*Q4.sig^2)*Q4.T)/(Q4.sig*sqrt(Q4.T))

thetaC = -Q4.S0*normpdf(d1)*Q4.sig/(2*sqrt(Q4.T)) - Q4.rf*Q4.K*exp(-Q4.rf*Q4.T)*normcdf(d1-Q4.sig*sqrt(Q4.T))

% Put option used to gamma hedge
Q4.KP = 25;
Q4.TP = 1/12;

Q4.gammaC = blsgamma(Q4.S0, Q4.K, Q4.rf, Q4.T, Q4.sig, 0);
Q4.gammaP_hedge = blsgamma(Q4.S0, Q4.KP, Q4.rf, Q4.TP, Q4.sig, 0);
[~, Q4.deltaP_hedge] = blsdelta(Q4.S0, Q4.KP, Q4.rf, Q4.TP, Q4.sig, 0);

Q4.n_gamma_hedge = -Q4.gammaC / Q4.gammaP_hedge; % long C so we should be short P to gamma hedge
Q4.net_delta = Q4.deltaC + Q4.n_gamma_hedge * Q4.deltaP_hedge % net delta which needs to be hedged using stocks

%% Q5
% Note: since the input numbers below are rounded, some of the calculations don't yield
% the correct result. The risk-neutral probabilities will be 50% exactly in reality, but
% here they are slightly different.
fprintf('\n\n* Question 5 (Defaultable bond) * \n')


% Inputs
Q5.lambda = 0.01;  % risk-neutral probability of default per period
Q5.delta = 0.7;    % recovery rate 

Q5.B01 = 0.97;     % one-period bond price

Q5.B02 = 0.93;     % two-period bond price
Q5.B12u = 0.9620;  % two-period bond price in up state
Q5.B12d = 0.9555;  % two-period bond price in down state

Q5.B03 = 0.88;     % three-period bond price
Q5.B13u = 0.9167;  % three-period bond price in up state
Q5.B13d = 0.8977;  % three-period bond price in down state
Q5.B23 = [0.9555; % three-period bond price in up-up state
            0.9504; % three-period bond price in up-down state
            0.9440; % three-period bond price in down-up state
            0.9349]; % three-period bond price in down-down state

Q5.r2_23 = 1 ./ Q5.B23 - 1; % one-period spot interest rate at time 2



% Risk-neutral probabilities
Q5.u02 = Q5.B12u / Q5.B02;  % size of up move in first period
Q5.d02 = Q5.B12d / Q5.B02; % size of down move in first period
Q5.pi0 = (1 / Q5.B01 - Q5.d02) / (Q5.u02 - Q5.d02);    % risk-neutral probability of up move in tree in first period

Q5.B02check = (Q5.pi0 * Q5.B12u + (1 - Q5.pi0) * Q5.B12d) * Q5.B01; % current price of two-period default free bond
Q5.r0_02 = sqrt(1 / Q5.B02) - 1;       % annualized yield on two-period default free bond


% Second-period risk-neutral probability
Q5.u13 = Q5.B23(1) / Q5.B13u;
Q5.d13 = Q5.B23(2) / Q5.B13u;
Q5.pi1u = (1 / Q5.B12u - Q5.d13) / (Q5.u13 - Q5.d13);

% a) caplet
C2 = max( (Q5.r2_23 - 0.05) ./ (1 + Q5.r2_23), 0)
C1 = [Q5.B12u * (Q5.pi1u * C2(1) + (1 - Q5.pi1u) * C2(2));
        Q5.B12d * (Q5.pi1u * C2(3) + (1 - Q5.pi1u) * C2(4))]
C0 = Q5.B01 * (Q5.pi0 * C1(1) + (1 - Q5.pi0) * C1(2))


% b) floorlet
P2 = max( (0.05 - Q5.r2_23) ./ (1 + Q5.r2_23), 0)
P1 = [Q5.B12u * (Q5.pi1u * P2(1) + (1 - Q5.pi1u) * P2(2));
        Q5.B12d * (Q5.pi1u * P2(3) + (1 - Q5.pi1u) * P2(4))]
P0 = Q5.B01 * (Q5.pi0 * P1(1) + (1 - Q5.pi0) * P1(2))


% c) put-call parity
r0_23 = Q5.B02 / Q5.B03 - 1 % one-period forward rate from time 2 to time 3
PV_r0_23 = r0_23 * Q5.B03 % present value of forward rate (discount for 3 periods since it will be paid at the end)
PV_K = 0.05 * Q5.B03 % present value of the strike

P0_pcp = C0 - PV_r0_23 + PV_K


% d) bond option
KB = 0.9524

BO2 = max( (Q5.B23 - KB), 0)
BO1 = [Q5.B12u * (Q5.pi1u * BO2(1) + (1 - Q5.pi1u) * BO2(2));
        Q5.B12d * (Q5.pi1u * BO2(3) + (1 - Q5.pi1u) * BO2(4))]
BO0 = Q5.B01 * (Q5.pi0 * BO1(1) + (1 - Q5.pi0) * BO1(2))

   
    
% e) CDS valuation
E_PVpremiums = @(c) c * (1 + (1 - Q5.lambda) * Q5.B01) 
E_PVpayout = (Q5.lambda * Q5.B01 + (1 - Q5.lambda) * Q5.lambda * Q5.B02) * (1 - Q5.delta)

Q5.CDSpremium = fzero(@(c) E_PVpremiums(c) - E_PVpayout, [0 1]);



% f) One-period defaultable bond
Q5.f = Q5.B01 * (Q5.lambda * Q5.delta + (1 - Q5.lambda) * 1)



%% Question 7
fprintf('\n\n* Question 7 (Super option) * \n')
K = 210
S = [200    220     242     266.2;
     0      180     198     217.8;
     0      0       162     178.2;
     0      0        0      145.8]
 
rf = 0.01 % interest rate per month with monthly compounding

% a) Replicating portfolio
U = S(1,2) / S(1,1)
D = S(2,2) / S(1,1)

pie = (1 + rf - D) / (U - D)

C(:,4) = (S(:,4) - K).^2 
C(1:3, 3) = (pie * C(1:3, 4) + (1 - pie) * C(2:4, 4)) ./ (1 + rf)
C(1:2, 2) = (pie * C(1:2, 3) + (1 - pie) * C(2:3, 3)) ./ (1 + rf)
C(1, 1) = (pie * C(1, 2) + (1 - pie) * C(2, 2)) ./ (1 + rf)

m = triu((C(1:3, 2:4) - C(2:4, 2:4)) ./ (S(1:3, 2:4) - S(2:4, 2:4)))

b = triu((U * C(2:4, 2:4) - D * C(1:3, 2:4)) ./ ((U - D) * (1 + rf)))

% b)
% Note: the portfolio of 3 options constructed here is not completely correct. The correct
% portfolio solves the problem in Bakshi, Kapadia & Madan, as long as a continuum of 
% strike prices is available.
vS = (140:0.1:280)'; % to plot continuous payout function
vK = [210 230 250]  % strikes to use in plot to show how the payoff of the super-option
                    % can be replicated using a portfolio of standard options
Seval = [S(2, 4) 240 S(1, 4)]   % Stock values at which to make sure the portfolio of options
                                % matches the payout
Swghts(1) = 2 * (Seval(1) - 210);
Swghts(2) = 2 * (Seval(2) - 210) - Swghts(1);
Swghts(3) = 2 * (Seval(3) - 210) - Swghts(2)

bwghts(1) = (Seval(1) - 210)^2 - Swghts(1) * (Seval(1) - vK(1));
bwghts(2) = (Seval(2) - 210)^2 - Swghts(1) * (Seval(2) - vK(1)) - ...
    Swghts(2) * (Seval(2) - vK(2)) %+ bwghts(1)
bwghts(3) = (Seval(3) - 210)^2 - Swghts * (Seval(3) - vK)' %+ bwghts(2)

pf = cumsum(repmat(Swghts, length(vS), 1) .* ...
    max(repmat(vS, 1, length(vK)) - repmat(vK, length(vS), 1), 0), 2) + ...
    repmat(bwghts, length(vS), 1);

for i = size(pf,2):-1:1
   % loop over columns to retain only the part of each portfolio that is really new
   % compared to the previous portfolio with one fewer option in it.
   pf(vS < vK(i), i) = NaN;
end
% pf(pf(:, 1) < 0, 1) = NaN;

figure(1);
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties

stem(S(:,4), C(:,4), 'LineWidth', lw);
hold on;

plot(vS, [(vS - K).^2, pf], 'LineWidth', lw);
hold off;
xlabel('Stock price at option maturity ($)')
ylabel('Super option payout ($)')

% Here we preserve the size of the image when we save it.
set(gcf, 'InvertHardcopy', 'on');
set(gcf, 'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf, 'PaperPosition', myfiguresize);

saveas(gcf, './figures/Q7b.png')
close(gcf);
    
% c)
format long;
C

format short
