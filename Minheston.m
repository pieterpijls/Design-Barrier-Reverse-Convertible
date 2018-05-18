%% calibration to option market prices using Heston model
function [out] = Minheston(kappa,eta,theta,rho,v0,S0,K,r,T,market_price,flag)
% Take unique values of variables
T_unique = unique(T); % take unique T values
r_unique = unique(r); % take unique r values
H_price = [];

% calculate Heston price with function Priceheston
for i=1:length(T_unique) % loop over unique Ts
    Ti = T_unique(i);
    ri = r_unique(i);
    strikes = K(T==Ti); % Take strike corresponding with unique T
    flags = flag(T==Ti); %  Use flag to distinguish between call and put options
    H_price = [H_price; Priceheston(S0,kappa,eta,theta,rho,v0,ri,Ti,strikes,flags)]; % calculate price
end

%minimize sum of absolute errors
out = sum(abs(market_price-H_price)./(market_price))/length(T); % take average of the sum absolute relative deviations 

end 

