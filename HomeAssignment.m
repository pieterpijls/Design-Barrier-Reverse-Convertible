close all; % close all figures 
clear; % clear workspace
rng(037894); % use random random seed 037894 to reproduce results

%% Read in Option Data
% We start with loading the market data 
option = xlsread('MSFToptions', 'options');
stock = xlsread('MSFToptions', 'stock');
interest_rate = xlsread('MSFToptions', 'interest'); 

S0 = stock;
q = 0.0159589; % from put-call parity
T = option(:,1);
K = option(:,2);
bid = option(:,3);
ask = option(:,4);
flag = option(:,5);
market_price = option(:,6); % mean of bid-ask or mid price

% computation of the interest rate for each maturity 
% Find corresponding interest rates with given maturities
r = spline(interest_rate(:,1), interest_rate(:,2)/100, T);
T = T/365; % maturity in years

% remove bid<0.2
K = K(bid>0.2);
flag = flag(bid>0.2);
ask = ask(bid>0.2);
market_price = market_price(bid>0.2);
r = r(bid>0.2);
T = T(bid>0.2);
bid = bid(bid>0.2);

% Remove options with T close to maturity % these options prices will deviate a lot
K = K(T>0.1);
flag = flag(T>0.1);
bid = bid(T>0.1);
ask = ask(T>0.1);
market_price = market_price(T>0.1);
r = r(T>0.1);
T = T(T>0.1);

%figure()
%scatter(K,market_price,'filled')

%% HESTON CALIBRATION

%Specify Heston parameters and give initial values, lower/upper boundary
% give random values to the parameters
kappa = rand(); % speed of mean reversion (>0)
eta = rand(); % reversion mean (>0)
theta = rand(); % vol of var (>0)
rho = rand(); % correlation volatility-stock price MSFT (1>rho>-1)
v0 = rand(); % initial volatility >0
x0 = [0.2542 0.9998 0.8096 -0.7 0.0119]; % initial values to speed up algo

lb = [0 0 0 -1 0]; % lower boundary
ub = [20 1 5 0 1]; % upper boundary

fun = @(x)(Minheston(x(1),x(2),x(3),x(4),x(5),S0,K,r,T,market_price,flag));

% create empty objects
A = [];
b = [];
Aeq = [];
beq = [];

% FIND OPTIMAL PARAMETER VECTOR 
% feller constraint included
x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,@fellerconstraint); % algo stops minimizing when last step size is smaller as size tolerance

% Heston prices and minimum ARPE
[minerror,prices] = Minhestonprices(x(1),x(2),x(3),x(4),x(5),S0,K,r,T,market_price,flag);

% ARPE of market prices vs prices devided by number of prices
error = sum(abs(market_price-prices)./(market_price))/length(prices); % ARPE of 10,63%

% Figure illustrate Heston fit to market prices
fig=figure();
scatter(K,market_price, 'o')
hold on
scatter(K,prices,'r*')
title('Calibration Heston')
legend('Market price','Heston: (\sigma_0,\theta,\kappa,\eta,\rho)')
xlabel('Strike')
ylabel('Option Price')
%saveas(fig,'heston.png')

%% BARRIER REVERSE CONVERTIBLE
% Calculate the fair value of the BRC

N = 1000; % invested amount in note

T_exotic_day=365; % Maturity 1 year (do not use 252 days for interest rates)
T_exotic=1;
K_exotic=S0; % strike barrier 
H=0.85*S0; % 15% lower = 0,85*94,01 = 79,91

% computation of the interest rate
r_exotic = spline(interest_rate(:,1), interest_rate(:,2)/100, T_exotic_day);

%% 1% Zero Coupon Bond
B = N/(1+r_exotic); %  calculate value of ZCB which will grow to 1,000 in one year

%% 2% Down-and-In Put Option
m=10000; % number of paths
n=T_exotic_day; % days to maturity 

%% HESTON STOCHASTIC VOLATILITY MODEL
n = floor(252*T_exotic); % step length has to be integer => use floor
m = 100000; % number of scenarios
dt = T_exotic/n; % number time steps

% Heston parameters to simulate volatility and stock prices
kappa = 5.0587; %x(1) % speed of mean reversion (>0)
eta = 0.075; %x(2)   % reversion mean (>0)
theta = 0.8243; %x(3) % vol of var (>0)
rho = -0.5325; %x(4)   % correlation vol-stock (1>rho>-1)
v0 = 0.0502; %x(5)   % init. vol (>0)

% create empty matric for the stock prices and volatility paths
S = zeros(m,n+1);
v = zeros(m,n+1);

% simulate standard normal random numbers
eps = normrnd(0,1,m,n);
epsS = normrnd(0,1,m,n);
eps1 = eps; 

% incorperate the correlation rho in the following way
eps2 = eps*rho + sqrt(1-rho^2)*epsS;
S(:,1) = S0;
v(:,1) = v0;

% Simulate paths for Microsoft Stock using exotic sigma
for j=2:n+1
    % we avoid negative volatility by using the absolute values
        v(:,j) = abs(v(:,j-1)+ (kappa*(eta-v(:,j-1)))*dt+theta*sqrt(v(:,j-1))*sqrt(dt).*eps2(:,j-1));
    % avoid negative stock prices 
        %S(:,j) = S(:,j-1).*(1+(r_exotic-q)*dt+sqrt(v(:,j-1))*sqrt(dt).*eps1(:,j-1));
         S(:,j) = S(:,j-1).*exp(((r_exotic-q)-v(:,j-1)/2)*dt+sqrt(dt*v(:,j-1)).*eps1(:,j-1));
end

fig=figure()
for i=1:10
hold on
plot(S(i,:)) % plot of one random path MSFT stock
xlim([0 252])
end
saveas(fig,'paths.png')

% barrier option
% knocks is when the barrier of 0,85*SO is hit and pays out when S(T) is lower as K=S0
DIBP_path = max((H-min(S'))./abs(H-min(S')),0)'.*max((K_exotic-S(:,n+1)),0)*exp(-r_exotic*T_exotic);
DIBP = mean(DIBP_path);

% We sell 10,6 DIBP to create the payoff of the note if N=1,000 USD
DIBP_amount = 10.64*DIBP; % around 91.40

% Number of stocks: N/S0 = 1mio/94.01 = 10,637 stocks

%% 3% Present value of all future coupons AND Margin of our BRC

% what is left for the coupons; receive N and amount of DIBP we sell and
% invest in ZCB
N+DIBP_amount-B;

% Frequency of pay outs
% Monthly coupon of 8.5 $ for 12 months => in total 102 $
i_coupons = zeros(1,12); % create empty matrices
coupons = zeros(1,12);

% Calculate the PV of the coupons
for i = 1:12
    i_coupons(:,i)  = spline(interest_rate(:,1), interest_rate(:,2)/100, 365*i/12);
    coupons(:,i)    = 8.5/(1+i_coupons(:,i)); % monthly coupon equal to 8.5
end

coupons_total = sum(coupons); % PV is equal to 99,97 $

% Margin is equal to 
margin = (N - (B-DIBP_amount+coupons_total))/N; % we keep a 1.139% margin

% Annualized coupon rate is 9,193% => monthly coupon rate is 0,73565%
% 1000*e^rT = 1000+110,5 <=> r=9,193% with T=1,14 with 110,5 = 13*8.5

%% DELTA HEDGING OF THE BRC

S0hedge = [94,94.01]; % vector of initial values 
DIBPprices = zeros(1,2); % empty matrix to save DIBP prices

for i=1:2
%for i=1:61
   S(:,1) = S0hedge(:,i);
   v(:,1) = v0; 
   
for j=2:n+1
    % we avoid negative volatility by using the absolute values
        v(:,j) = abs(v(:,j-1)+ (kappa*(eta-v(:,j-1)))*dt+theta*sqrt(v(:,j-1))*sqrt(dt).*eps2(:,j-1));
    % avoid negative stock prices by taking the log of the SDE
         S(:,j) = S(:,j-1).*exp(((r_exotic-q)-v(:,j-1)/2)*dt+sqrt(dt*v(:,j-1)).*eps1(:,j-1));
end

DIBP_path = max((H-min(S'))./abs(H-min(S')),0)'.*max((K_exotic-S(:,n+1))*exp(-r_exotic*T_exotic),0);
DIBPhedge = mean(DIBP_path);

DIBPprices(:,i)= DIBPhedge;  
end

% We only calculate one delta of option when S_0 is 94,01
delta = (DIBPprices(:,2)-DIBPprices(:,1))/(S0hedge(2)-S0hedge(1)); %equal to -0,4237

%% DELTA FOR DIFFERENT S0

% to create figure delta vs. S0
S0hedge = 60:0.5:110; 
DIBPprices = zeros(1,101);
deltahedge = zeros(1,101);

% allows us to calculate multiple deltas for different initial values S0
for i=1:101
%for i=1:61
   S(:,1) = S0hedge(:,i);
   v(:,1) = v0; 
   
for j=2:n+1
    % we avoid negative volatility by using the absolute values
        v(:,j) = abs(v(:,j-1)+ (kappa*(eta-v(:,j-1)))*dt+theta*sqrt(v(:,j-1))*sqrt(dt).*eps2(:,j-1));
    % avoid negative stock prices by taking the log of the SDE
         S(:,j) = S(:,j-1).*exp(((r_exotic-q)-v(:,j-1)/2)*dt+sqrt(dt*v(:,j-1)).*eps1(:,j-1));
end

DIBP_path = max((H-min(S'))./abs(H-min(S')),0)'.*max((K_exotic-S(:,n+1))*exp(-r_exotic*T_exotic),0);
DIBPhedge = mean(DIBP_path);

DIBPprices(:,i)= DIBPhedge;  

   if i>1
      deltahedge(:,i)= (DIBPprices(:,i)-DIBPprices(:,i-1))/(S0hedge(:,i)-S0hedge(:,i-1));
   else
      deltahedge(:,i)=-1;
   end
end

fig=figure();
plot(S0hedge(2:101),deltahedge(2:101));
xlabel('Stock Price S_0')
ylabel('Delta \Delta')
line([94.01 94.01],[-1,0],'Color','black','LineStyle','--');
line([H H],[-1,0],'Color','red','LineStyle','-.');
%saveas(fig,'delta.png')