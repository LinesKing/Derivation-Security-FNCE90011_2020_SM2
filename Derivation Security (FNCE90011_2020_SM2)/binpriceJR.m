function [pr,opt] = binpriceJR(so,x,r,t,dt,sig,flag,q,div,exdiv, varargin) 
% THIS VERSION ADAPTED TO REFLECT THE JARROW-RUDD TREE INSTEAD OF COX-ROSS-RUBINSTEIN.
% Also extended to European options (flag == 3 for calls, flag == 2 for puts).
%
% Also allows to specify U and D directly.
%
%BINPRICE Binomial put and call pricing. 
%   [PR,OPT] = BINPRICE(SO,X,R,T,DT,SIG,FLAG,Q,DIV,EXDIV) prices an option  
%   using a binomial pricing model. SO is the underlying asset price, X is the
%   option exercise price, R is the risk-free interest rate, T is the option's
%   time until maturity in years and DT is the time increment within T.  
%   DT will be adjusted so that the length of each interval is consistent with
%   the maturity time of the option. SIG is the asset's volatility, FLAG 
%   specifies whether the option is a call (flag = 1) or a put (flag = 0),
%   Q is the dividend rate, DIV is the dividend payment at an ex-dividend date,
%   EXDIV.  EXDIV is specified in number of periods.  All inputs to this  
%   function are scalar values except DIV and EXDIV which are 1-by-n vectors.
%   For each dividend payment, there must be a corresponding ex-dividend date.
%   By default q,div, and EXDIV equal 0.  If a value is entered for the 
%   dividend rate q, DIV and EXDIV should equal 0 or not be entered. If values
%   are entered for DIV and EXDIV, set Q = 0. 
% 
%   [P,O] = binprice(52,50,.1,5/12,1/12,.4,0,0,2.06,3.5) returns the asset
%   price and option value at each node of the binary tree. 
% 
%   P = 
% 
%    52.0000   58.1367   65.0226   72.7494   79.3515   89.0642 
%          0   46.5642   52.0336   58.1706   62.9882   70.6980 
%          0         0   41.7231   46.5981   49.9992   56.1192 
%          0         0         0   37.4120   39.6887   44.5467 
%          0         0         0         0   31.5044   35.3606 
%          0         0         0         0         0   28.0688 
% 
%   O = 
% 
%     4.4404    2.1627    0.6361         0         0         0 
%          0    6.8611    3.7715    1.3018         0         0 
%          0         0   10.1591    6.3785    2.6645         0 
%          0         0         0   14.2245   10.3113    5.4533 
%          0         0         0         0   18.4956   14.6394 
%          0         0         0         0         0   21.9312 
% 
%   See also BLSPRICE. 
% 
%   Reference: Options, Futures, and Other Derivative Securities, 
%              2nd Edition, Hull, Chapter 14. 
 
%       Copyright 1995-2012 The MathWorks, Inc.
 
if nargin < 8  
  q = 0; 
end 

if nargin < 9 
  div = 0; 
  exdiv = 0; 
end 

if nargin < 7 
  error(message('finance:binprice:missingInputs')) 
end 

if flag ~= 0 && flag ~= 1 && flag ~= 2 && flag ~= 3
  error(message('finance:binprice:invalidFlag')) 
end 

if q ~= 0 && div ~= 0 
  		error(message('finance:binprice:invalidDivExDiv')) 
end 

% Calculate the number of periods, nper, and the length of each interval, dt.
% Make sure that the length of each interval is consistent with 
% the maturity time of the option.
nper = round(t/dt);             % Number of periods after time zero 
dt = t/nper;
npp = nper+1;                   % Number of periods including time zero

% Make sure there are no dividend payments after instrument maturity:
if(exdiv > nper)
    error(message('finance:binprice:invalidDivPayment'))
end

% Calculate the probability of an upward price movement 
if nargin > 10
    u = varargin{1};
    d = varargin{2};    
else 
    u = exp((r - 0.5*sig^2)*dt + sig.*sqrt(dt));
    d = exp((r - 0.5*sig^2)*dt - sig.*sqrt(dt));
end
a = exp((r-q).*dt); 
p = (a-d)./(u-d); 

% jspan = -fix(nper*.5);          % j-th node offset number 
% ispan = rem(round(t/dt),2);      % i-th node offset number 
i = 0:nper;         % i-th node numbers 
% j = (jspan:(nper+jspan))';      % j-th node numbers 
% jex = j(:,ones(size(i')));      % expand i and j to eliminate for loop 
iex = i(ones(nper+1, 1),:); 
 
pvdiv = div.*exp(-exdiv.*dt.*r);  % Find present value of all dividends 
so = so-sum(pvdiv(:));            % Find current price - div present values 
 
% Asset price at nodes, matrix is flipped so tree appears correct visually 
% pr = triu(rot90(so.*u.^jex.*d.^(iex-jex),2));   
nu = (gallery('circul', nper + 1) - 1);
nd = iex - nu;

pr = triu(so.*u.^nu.*d.^nd);


if div ~= 0                 % Present value of future dividends at nodes 
  lendiv = length(div(:)); 
  lenexdiv = length(exdiv(:)); 
  if lendiv ~= lenexdiv 
    error(message('finance:binprice:mismatchDivExDiv')) 
  end 
 
  dpvtot = zeros(npp);          % Preallocate matrix 
  for y = 1:lenexdiv 
    z = (exdiv(y):-1:0);        % Create vector from 0 to ex-div date 
    dpv = div(y)*exp(-z*dt*r);  % Discount dividends nodes 
    dpvmat = [dpv(ones(npp,1),:) zeros(npp,npp-length(dpv))]; % Expand matrix 
    dpvtot = dpvtot + dpvmat;   % Add next discounted dividend to total 
  end 
  m = find(pr~=0);              % Find nodes where option will have value 
  pr(m) = pr(m)+dpvtot(m);      % combine div pv's and prices to get new prices 
end 

[NegRow, NegCol] = find(pr<0);
if any(NegRow)
   error(message('finance:binprice:divPaymentExceedsPrice',sprintf('%g',dt*(NegCol(1)-1))));
end

opt = zeros(size(pr));    
if flag == 3                        % Option is a European call 
  opt(:,npp) = max(pr(:,npp)-x,0);  % Determine option values from prices 
  for n = nper:-1:1 
    k = 1:n; 
    % Probable option values discounted back one time step 
    discopt = (p*opt(k,n+1)+(1-p)*opt(k+1,n+1))*exp(-r*dt); 
    % Option value is max of current price - X or discopt 
    opt(:,n) = [discopt; zeros(npp-n,1)]; 
  end 
elseif flag == 2                    % Option is a European put 
  opt(:,npp) = max(x-pr(:,npp),0);  % Determine option values from prices 
  for n = nper:-1:1 
    k = 1:n; 
    % Probable option values discounted back one time step 
    discopt = (p*opt(k,n+1)+(1-p)*opt(k+1,n+1))*exp(-r*dt); 
    % Option value is max of X - current price or discopt 
    opt(:,n) = [discopt ;zeros(npp-n,1)]; 
  end 
elseif flag == 1                        % Option is an American call 
  opt(:,npp) = max(pr(:,npp)-x,0);  % Determine option values from prices 
  for n = nper:-1:1 
    k = 1:n; 
    % Probable option values discounted back one time step 
    discopt = (p*opt(k,n+1)+(1-p)*opt(k+1,n+1))*exp(-r*dt); 
    % Option value is max of current price - X or discopt 
    opt(:,n) = [max(pr(1:n,n)-x,discopt);zeros(npp-n,1)]; 
  end 
elseif flag == 0                    % Option is an American put 
  opt(:,npp) = max(x-pr(:,npp),0);  % Determine option values from prices 
  for n = nper:-1:1 
    k = 1:n; 
    % Probable option values discounted back one time step 
    discopt = (p*opt(k,n+1)+(1-p)*opt(k+1,n+1))*exp(-r*dt); 
    % Option value is max of X - current price or discopt 
    opt(:,n) = [max(x-pr(1:n,n),discopt);zeros(npp-n,1)]; 
  end 
end
