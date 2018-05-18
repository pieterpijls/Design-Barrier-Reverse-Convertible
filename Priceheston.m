%% Fast Frontier Transform pricing under Heston model using Car-Madan Formula
function [price,KK,C] = Priceheston(S_0,kappa,eta,theta,rho,v0,r,T,K,flag)

% variabled d and g for Heston char function
d = @(u)((rho*theta*u*1i-kappa).^2-theta.^2*(-1i*u-u.^2)).^(1/2);
g = @(u)(kappa-rho*theta*u*1i-d(u))./(kappa-rho*theta*u*1i+d(u));

q = 0.0159589; % dividend yield from put call parity

% Charasteristics function for Heston with no jumps
phi = @(u)(exp(1i*u.*(log(S_0)+(r - q).*T)) ...
      .*exp(eta*kappa*theta^(-2).*((kappa-rho*theta*u*1i-d(u)).*T-2*log((1-g(u).*exp(-d(u).*T))./(1-g(u))))) ... 
      .*exp(v0*theta^(-2)*(kappa-rho*theta*u*1i-d(u)).*(1-exp(-d(u).*T))./(1-g(u).*exp(-d(u).*T))));

  
% Inputs
N = 4096; % grid points
a = 1.5; %1.75?
eta2 = 0.25; 
lambda = 2*pi/N/eta2;
b = lambda*N/2;
k = [(-b):lambda:(b-lambda)]; %grid of strikes
KK = exp(k); 
v2 = [0:eta2:(N-1)*eta2]; % treshhold
sw =(3+(-1).^(1:1:N)); % simpsons rule
sw(1) = 1;
sw = sw/3;
rho = exp(-r.*T).*phi(v2-(a+1)*1i)./(a^2+a-v2.^2+1i*(2*a+1)*v2);
C = exp(-a*k)/pi.*real(fft(rho.*exp(1i*v2*b)*eta2.*sw));

%Calculate Put options using the Put-Call Parity
price = spline(KK,C,K); %calculate the price using the spline function and grid 
for i=1:length(K)
    if flag(i)==1 % if option
        price(i) = -exp(-q.*T)*S_0+ K(i).*exp(-r.*T)+price(i); % put call parity
    end 
end
end

