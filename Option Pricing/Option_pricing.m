% clc;
% clear all;
% format long

% Pricing a Barrier option using Monte Carlo simulations

S0 = 100;     % spot price of the underlying stock today
K = 105;      % strike at expiry
mu = 0.05;    % expected return
sigma = 0.2;  % volatility
r = 0.05;     % risk-free rate
T = 1.0;      % years to expiry
Sb = 110;     % barrier


% Define variable numSteps to be the number of steps for multi-step MC
% numPaths - number of sample paths used in simulations
numPaths = 10000;
[paths_1] = GRWPaths(S0, mu, sigma, T, 1, numPaths);
[paths_252] = GRWPaths(S0, mu, sigma, T, 252, numPaths);
% Implement your Black-Scholes pricing formula
[call_BS_European_Price, putBS_European_Price] = BS_european_price(S0, K, T, r, sigma);

% Implement your one-step Monte Carlo pricing procedure for European option

[callMC_European_Price_1_step, putMC_European_Price_1_step] = MC_european_price(paths_1,S0, K, T, r, mu, sigma, 1, numPaths);
numSteps = 252;
% Implement your multi-step Monte Carlo pricing procedure for European option
[callMC_European_Price_multi_step, putMC_European_Price_multi_step] = MC_european_price(paths_252,S0, K, T, r, mu, sigma, numSteps, numPaths);

% % Implement your one-step Monte Carlo pricing procedure for Barrier option
numSteps = 1;
[callMC_Barrier_Knockin_Price_1_step, putMC_Barrier_Knockin_Price_1_step] = MC_barrier_knockin_price(paths_1,S0, Sb, K, T, r, mu, sigma, numSteps, numPaths);
% 
% % Implement your multi-step Monte Carlo pricing procedure for Barrier option
numSteps = 252;
[callMC_Barrier_Knockin_Price_multi_step, putMC_Barrier_Knockin_Price_multi_step] =MC_barrier_knockin_price(paths_252,S0, Sb, K, T, r, mu, sigma, numSteps, numPaths);
% 
disp(['Black-Scholes price of an European call option is ',num2str(call_BS_European_Price)])
disp(['Black-Scholes price of an European put option is ',num2str(putBS_European_Price)])
disp(['One-step MC price of an European call option is ',num2str(callMC_European_Price_1_step)])
disp(['One-step MC price of an European put option is ',num2str(putMC_European_Price_1_step)])
disp(['Multi-step MC price of an European call option is ',num2str(callMC_European_Price_multi_step)])
disp(['Multi-step MC price of an European put option is ',num2str(putMC_European_Price_multi_step)])
disp(['One-step MC price of an Barrier call option is ',num2str(callMC_Barrier_Knockin_Price_1_step)])
disp(['One-step MC price of an Barrier put option is ',num2str(putMC_Barrier_Knockin_Price_1_step)])
disp(['Multi-step MC price of an Barrier call option is ',num2str(callMC_Barrier_Knockin_Price_multi_step)])
disp(['Multi-step MC price of an Barrier put option is ',num2str(putMC_Barrier_Knockin_Price_multi_step)])


% Plot results
figure(1);
%%%%%%%%%%% Insert your code here %%%%%%%%%%%%
figure ;
set(gcf, 'color', 'white');
plot(0:1, paths_1', 'Linewidth', 2);
title('Geometric Random Walk Paths of the Stock Price', 'FontWeight', 'bold');
xlabel('Time')
ylabel('Price')
figure(2);
%%%%%%%%%%% Insert your code here %%%%%%%%%%%%
figure ;
set(gcf, 'color', 'white');
plot(0:252, paths_252', 'Linewidth', 2);
title('Geometric Random Walk Paths of the Stock Price', 'FontWeight', 'bold');
xlabel('Time')
ylabel('Price')
find the Optimal numsteps
% Find the optimal numsteps.
% pricing_error = [];
% i = 1;
% for numSteps = 1:6:252
%     paths = GRWPaths(S0, mu, sigma, T, numSteps, numPaths);
%     [call_price1, put_price1] = BS_european_price(S0, K, T, r, sigma);
%     [call_price2, put_price2] = MC_european_price(paths,S0, K, T, r, mu, sigma, numSteps, numPaths);
%     pricing_error(i) = (abs(call_price2 -call_price1) + abs(put_price2 -put_price1))/2;
%     i=i+1;
%     display(numSteps);
% end
% fprintf('optimal number of steps =', 6*find(pricing_error == min(pricing_error)))
% figure;
% set(gcf, 'color', 'white');
% plot([1:6:252],pricing_error,'r');
% title('Pricing error for different step numbers', 'FontWeight', 'bold');
% xlabel('Number of steps')
% ylabel('Pricing error ($)')

% Pricing a European option using Black-Scholes formula and Monte Carlo simulations
function [call_price, put_price] = BS_european_price(S, K, T, r, sigma)
t = 0;
d1 = 1/(sigma*sqrt(T-t))*(log(S/K)+(r+(sigma)^(2)/2)*(T-t));
d2 = d1 - sigma*sqrt(T-t);
call_price = normcdf(d1)*S-normcdf(d2)*K*exp(-r*(T-t));
put_price = normcdf(-d2)*K*exp(-r(T-t))-normcdf(-d1)*S;
end

function paths = GRWPaths(initPrice, mu, sigma, T, numSteps, numPaths)
    % Computes numPaths random paths for a geometric random walk
    % mu is the annual drift, sigma the annual volatility
    % T is the total length of time for the path (in years)
    % dT is the time increment (in years)
       
    paths = zeros(numSteps+1, numPaths);
    dT = T/numSteps;
    
    % Vector of paths will store realizations of the asset price
    % First asset price is the initial price
    paths(1,:) = initPrice;
 
    % Generate paths
    for iPath = 1:numPaths
        for iStep = 1:numSteps
            paths(iStep+1, iPath) = paths(iStep, iPath) * exp( (mu - 0.5*sigma^2)*dT + sigma*sqrt(dT)*normrnd(0,1) );
        end
    end 
            
%     % Plot paths
%     figure;
%     set(gcf, 'color', 'white');
%     plot(0:numSteps, paths', 'Linewidth', 2);
%     title('Geometric Random Walk Paths of the Stock Price', 'FontWeight', 'bold');
%     xlabel('Time')
%     ylabel('Price')

end

% Pricing a European option using Monte Carlo simulations
function [call_price, put_price] = MC_european_price(path,S0, K, T, r, mu, sigma, stepNum, numPaths)

% Calculate the payoff for each path for a Put
PutPayoffT = max(K-(path(end,:)),0);

% Calculate the payoff for each path for a Call
CallPayoffT = max((path(end,:))-K,0);

% Discount back
put_price = mean(PutPayoffT)*exp(-r*T);
call_price = mean(CallPayoffT)*exp(-r*T);

end

function [call_price, put_price] = MC_barrier_knockin_price(path,S0, Sb, K, T, r, mu, sigma, numSteps, numPaths)


call_price = [];
put_price = [];
for i = 1: numPaths
    for j = 1: numSteps+1
        if path(j,i) >= Sb
            call_price(i) = max(path(end,i)-K,0);
            put_price(i) = max(K-path(end,i),0);
            break
        else 
            call_price(i) = 0;
            put_price(i) = 0;
        end
    end
end

call_price = mean(call_price)*exp(-r*T);
put_price = mean(put_price)*exp(-r*T);
end
            
        

