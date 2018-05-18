%% Heston model prices using calibrated parameters
function [out,H_price] = Minhestonprices(kappa,eta,theta,rho,v0,S_0,K,r,T,market_price,flag)
% Take unique values of variables
T_unique = unique(T); % take unique T values
r_unique = unique(r); % take unique r values
H_price = [];

% calculate heston prices
for i=1:length(T_unique)
    Ti = T_unique(i);
    ri = r_unique(i);
    strikes = K(T==Ti);
    flags = flag(T==Ti);
    H_price = [H_price; Priceheston(S_0,kappa,eta,theta,rho,v0,ri,Ti,strikes,flags)];
end

out = sum(abs(market_price-H_price)./(market_price))/length(T);

end 
