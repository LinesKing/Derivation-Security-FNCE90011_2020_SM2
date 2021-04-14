function [c, p] = blsprice(S0t, K, r, T, sigma, div)
d1 = (log(S0t/K) + (r + 0.5 * sigma^2) * T) / (sigma * sqrt(T));
d2 = d1 - sigma * sqrt(T);

Nd1 = normcdf(d1);
Nd2 = normcdf(d2);

c = S0t * exp(-div*T) * Nd1 - K * exp(-r*T) * Nd2;

p = K * exp(-r*T) * normcdf(-d2) - S0t * exp(-div*T) * normcdf(-d1);
end

