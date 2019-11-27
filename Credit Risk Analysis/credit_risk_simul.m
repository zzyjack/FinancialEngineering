clear all;
clc
format long;

Nout  = 100000; % number of out-of-sample scenarios
Nin   = 5000;   % number of in-sample scenarios
Ns    = 5;      % number of idiosyncratic scenarios for each systemic

C = 8;          % number of credit states

% Filename to save out-of-sample scenarios
filename_save_out  = 'scen_out';

% Read and parse instrument data
instr_data = dlmread('instrum_data.csv', ',');
instr_id   = instr_data(:,1);           % ID
driver     = instr_data(:,2);           % credit driver
beta       = instr_data(:,3);           % beta (sensitivity to credit driver)
recov_rate = instr_data(:,4);           % expected recovery rate
value      = instr_data(:,5);           % value
prob       = instr_data(:,6:6+C-1);     % credit-state migration probabilities (default to A)
exposure   = instr_data(:,6+C:6+2*C-1); % credit-state migration exposures (default to A)
retn       = instr_data(:,6+2*C);       % market returns

K = size(instr_data, 1); % number of  counterparties

% Read matrix of correlations for credit drivers
rho = dlmread('credit_driver_corr.csv', '\t');
sqrt_rho = (chol(rho))'; % Cholesky decomp of rho (for generating correlated Normal random numbers)

disp('======= Credit Risk Model with Credit-State Migrations =======')
disp('============== Monte Carlo Scenario Generation ===============')
disp(' ')
disp(' ')
disp([' Number of out-of-sample Monte Carlo scenarios = ' int2str(Nout)])
disp([' Number of in-sample Monte Carlo scenarios = ' int2str(Nin)])
disp([' Number of counterparties = ' int2str(K)])
disp(' ')

% Find credit-state for each counterparty
% 8 = AAA, 7 = AA, 6 = A, 5 = BBB, 4 = BB, 3 = B, 2 = CCC, 1 = default
[Ltemp, CS] = max(prob, [], 2);
clear Ltemp

% Account for default recoveries
exposure(:, 1) = (1-recov_rate) .* exposure(:, 1);

% Compute credit-state boundaries
CS_Bdry = norminv( cumsum(prob(:,1:C-1), 2) );

% -------- Insert your code here -------- %
    % generate matrix for y with 50 drivers and 100000 scenarios
    y = zeros(Nout,50);
    % generate vector matrix for counterparties corresponding y of size 100000x100
    yk = zeros(Nout,100);
    % generate matrix for y with 100 drivers and 100000 scenarios
    w = zeros(Nout,K);
    % Compute sigma for each counterparty j
    sigma = sqrt(1-beta.^2);
    % generate losses vector for 100000 scenarios
    Losses_out = zeros(Nout,K);
 %if(~exist('scenarios_out.mat','file'))
    % -------- Insert your code here -------- %
    z = randn(Nout,K);
    for s = 1:Nout
        % -------- Insert your code here -------- %
        % generate y's using correlation matrix
        y(s,:) = (sqrt_rho*randn(50,1))';
        % extract corresponding y's for each counterparty
        for i = 1: K
            yk(s,i) = y(s,driver(i));
            w(s,i) = yk(s,i)*beta(i)+sigma(i)*z(s,i);
            index = find(sort([w(s,i),CS_Bdry(i,:)])==w(s,i));
            Losses_out(s,i)=exposure(i,index);
        end
    end
    
    % Calculated out-of-sample losses (100000 x 100)
    % Losses_out
    % save('scenarios_out', 'Losses_out')
 %else
    %load('scenarios_out', 'Losses_out')
%end

% Normal approximation computed from out-of-sample scenarios
 mu_l = mean(Losses_out)';
 var_l = cov(Losses_out);
% 
% Compute portfolio weights
portf_v = sum(value);     % portfolio value
w0{1} = value / portf_v;  % asset weights (portfolio 1)
w0{2} = ones(K, 1) / K;   % asset weights (portfolio 2)
x0{1} = (portf_v ./ value) .* w0{1};  % asset units (portfolio 1)
x0{2} = (portf_v ./ value) .* w0{2};  % asset units (portfolio 2)

% Quantile levels (99%, 99.9%)
alphas = [0.99 0.999];

% Compute VaR and CVaR (non-Normal and Normal) for 100000 scenarios
for(portN = 1:2)
    for(q=1:length(alphas))
        alf = alphas(q);
        % -------- Insert your code here -------- %
        x=x0{portN};
        w=w0{portN};
        port_loss = sort(Losses_out*x);
        mean_loss = mean(port_loss);
        std_loss = std(port_loss);
        VaRout(portN,q)  = port_loss(ceil(Nout*alf));
        VaRinN(portN,q)  = mean_loss + norminv(alf,0,1)*std_loss;
        CVaRout(portN,q) = (1/(Nout*(1-alf)))*((ceil(Nout*alf)-Nout*alf)*VaRout(portN,q)+sum(port_loss(ceil(Nout*alf)+1:end)));
        CVaRinN(portN,q) = mean_loss+(normpdf(norminv(alf,0,1))/(1-alf))*std_loss;
        % -------- Insert your code here -------- %        
 end
end


% Perform 100 trials
N_trials = 100;

for(tr=1:N_trials)
    
    % Monte Carlo approximation 1

    % -------- Insert your code here -------- %
    y1=zeros(Nin/Ns,50);
    y1k = zeros(Nin/Ns,K);
    w1=zeros(Nin,K);
    Losses_inMC1=zeros(Nin,K);
    
    for s = 1:ceil(Nin/Ns) % systemic scenarios
        % -------- Insert your code here -------- %
        y1(s,:) = (sqrt_rho*randn(50,1))';
        z1 = randn(Ns,K);
        for si = 1:Ns % idiosyncratic scenarios for each systemic
            % -------- Insert your code here -------- %
            for i = 1: K
            y1k(s,i) = y1(s,driver(i));
            w1(Ns*(s-1)+si,i) = y1k(s,i)*beta(i)+sigma(i)*z1(si,i);
            index = find(sort([w1(Ns*(s-1)+si,i),CS_Bdry(i,:)])==w1(Ns*(s-1)+si,i));
            Losses_inMC1(Ns*(s-1)+si,i)=exposure(i,index);
            end
        end
    end
    
    % Calculated losses for MC1 approximation (5000 x 100)
    % Losses_inMC1
    
    % Monte Carlo approximation 2
    
    % -------- Insert your code here -------- %
    y2=zeros(Nin,50);
    y2k = zeros(Nin,K);
    w2=zeros(Nin,K);
    Losses_inMC2=zeros(Nin,K);
    z2 = randn(Nin,K);
    for s = 1:Nin % systemic scenarios (1 idiosyncratic scenario for each systemic)
        % -------- Insert your code here -------- %
        y2(s,:) = (sqrt_rho*randn(50,1))';
        % extract corresponding y's for each counterparty
        for i = 1: K
            y2k(s,i) = y2(s,driver(i));
            w2(s,i) = y2k(s,i)*beta(i)+sigma(i)*z2(s,i);
            index = find(sort([w2(s,i),CS_Bdry(i,:)])==w2(s,i));
            Losses_inMC2(s,i)=exposure(i,index);
        end
    end
        
    % Calculated losses for MC2 approximation (5000 x 100)
    % Losses_inMC2
    
    % Compute VaR and CVaR
    for(portN = 1:2)
        for(q=1:length(alphas))
            alf = alphas(q);
            % -------- Insert your code here -------- %   
            x = x0{portN};
            w = w0{portN};
            % Compute portfolio loss 
            portf_loss_inMC1 = sort(Losses_inMC1*x);
            portf_loss_inMC2 = sort(Losses_inMC2*x);
            mu_MCl = mean(Losses_inMC1)';
            var_MCl = cov(Losses_inMC1);
            mu_MC2 = mean(Losses_inMC2)';
            var_MC2 = cov(Losses_inMC2);
            % Compute portfolio mean loss mu_p_MC1 and portfolio standard deviation of losses sigma_p_MC1
            mu_p_MC1 = mean(portf_loss_inMC1);
            sigma_p_MC1 = std(portf_loss_inMC1);
            % Compute portfolio mean loss mu_p_MC2 and portfolio standard deviation of losses sigma_p_MC2
            mu_p_MC2 = mean(portf_loss_inMC2);
            sigma_p_MC2 = std(portf_loss_inMC2);
            % Compute VaR and CVaR for the current trial
            VaRinMC1{portN,q}(tr) = portf_loss_inMC1(ceil(Nin*alf));
            VaRinMC2{portN,q}(tr) = portf_loss_inMC2(ceil(Nin*alf));
            VaRinN1{portN,q}(tr) = mu_p_MC1 + norminv(alf,0,1)*sigma_p_MC1;
            VaRinN2{portN,q}(tr) = mu_p_MC2 + norminv(alf,0,1)*sigma_p_MC2;
            CVaRinMC1{portN,q}(tr) = (1/(Nin*(1-alf))) * ( (ceil(Nin*alf)-Nin*alf) * VaRinMC1{portN,q}(tr) + sum(portf_loss_inMC1(ceil(Nin*alf)+1:Nin)) );
            CVaRinMC2{portN,q}(tr) = (1/(Nin*(1-alf))) * ( (ceil(Nin*alf)-Nin*alf) * VaRinMC2{portN,q}(tr) + sum(portf_loss_inMC2(ceil(Nin*alf)+1:Nin)) );
            CVaRinN1{portN,q}(tr) = mu_p_MC1 + (normpdf(norminv(alf,0,1))/(1-alf))*sigma_p_MC1;
            CVaRinN2{portN,q}(tr) = mu_p_MC2 + (normpdf(norminv(alf,0,1))/(1-alf))*sigma_p_MC2;
            % -------- Insert your code here -------- %
        end
    end
end

% Display portfolio VaR and CVaR
for(portN = 1:2)
fprintf('\nPortfolio %d:\n\n', portN)    
 for(q=1:length(alphas))
    alf = alphas(q);
    fprintf('Out-of-sample: VaR %4.1f%% = $%6.2f, CVaR %4.1f%% = $%6.2f\n', 100*alf, VaRout(portN,q), 100*alf, CVaRout(portN,q))
    fprintf('In-sample MC1: VaR %4.1f%% = $%6.2f, CVaR %4.1f%% = $%6.2f\n', 100*alf, mean(VaRinMC1{portN,q}), 100*alf, mean(CVaRinMC1{portN,q}))
    fprintf('In-sample MC2: VaR %4.1f%% = $%6.2f, CVaR %4.1f%% = $%6.2f\n', 100*alf, mean(VaRinMC2{portN,q}), 100*alf, mean(CVaRinMC2{portN,q}))
    fprintf(' In-sample No: VaR %4.1f%% = $%6.2f, CVaR %4.1f%% = $%6.2f\n', 100*alf, VaRinN(portN,q), 100*alf, CVaRinN(portN,q))
    fprintf(' In-sample N1: VaR %4.1f%% = $%6.2f, CVaR %4.1f%% = $%6.2f\n', 100*alf, mean(VaRinN1{portN,q}), 100*alf, mean(CVaRinN1{portN,q}))
    fprintf(' In-sample N2: VaR %4.1f%% = $%6.2f, CVaR %4.1f%% = $%6.2f\n\n', 100*alf, mean(VaRinN2{portN,q}), 100*alf, mean(CVaRinN2{portN,q}))
 end
end

%Plot results
figure(1);
% -------- Insert your code here -------- %
set(gcf, 'color', 'white');
[frequencyCounts, binLocations] = hist(Losses_out*x0{2}, 150);
bar(binLocations, frequencyCounts,'DisplayName', 'Out of Sample Distribution');
hold on;
line([VaRout(2,1) VaRout(2,1)], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--','DisplayName','VaRout 99%');
line([VaRout(2,2) VaRout(2,2)], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--','DisplayName','VaRout 99.9%');
hold on;
normf = ( 1/(std(Losses_out * x0{2})*sqrt(2*pi)) ) * exp( -0.5*((binLocations-mean(Losses_out * x0{2}))/std(Losses_out * x0{2})).^2 );
normf = normf * sum(frequencyCounts)/sum(normf);
plot(binLocations, normf, 'r', 'LineWidth', 1,'DisplayName','Normal Distribution');
hold on;
line([VaRinN(2,1) VaRinN(2,1)], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '-.','DisplayName','VaRN 99%');
line([VaRinN(2,2) VaRinN(2,2)], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '-.','DisplayName','VaRN 99.9%');
hold off;
text(1*VaRout(2,1), max(frequencyCounts)/1.9, 'VaRout 99%')
text(1*VaRout(2,2), max(frequencyCounts)/1.9, 'VaRout 99.9%')
text(0.7*VaRinN(2,1), max(frequencyCounts)/1.9, 'VaRN 99%')
text(0.9*VaRinN(2,2), max(frequencyCounts)/1.9, 'VaRN 99.9%')
xlabel('Distribution of Losses')
ylabel('Frequency')
legend('show')
title('Out of Sample Distribution of Losses for Portfolio 2')

figure(2);
% -------- Insert your code here -------- %
set(gcf, 'color', 'white');
[frequencyCounts, binLocations] = hist(Losses_out*x0{1}, 150);
bar(binLocations, frequencyCounts,'DisplayName', 'Out of Sample Distribution');
hold on;
line([VaRout(1,1) VaRout(1,1)], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--','DisplayName','VaRout 99%');
line([VaRout(1,2) VaRout(1,2)], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--','DisplayName','VaRout 99.9%');
hold on;
normf = ( 1/(std(Losses_out * x0{1})*sqrt(2*pi)) ) * exp( -0.5*((binLocations-mean(Losses_out * x0{1}))/std(Losses_out * x0{1})).^2 );
normf = normf * sum(frequencyCounts)/sum(normf);
plot(binLocations, normf, 'r', 'LineWidth', 1,'DisplayName','Normal Distribution');
hold on;
line([VaRinN(1,1) VaRinN(1,1)], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '-.','DisplayName','VaRN 99%');
line([VaRinN(1,2) VaRinN(1,2)], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '-.','DisplayName','VaRN 99.9%');
hold off;
text(1*VaRout(1,1), max(frequencyCounts)/1.9, 'VaRout 99%')
text(1*VaRout(1,2), max(frequencyCounts)/1.9, 'VaRout 99.9%')
text(0.7*VaRinN(1,1), max(frequencyCounts)/1.9, 'VaRN 99%')
text(1*VaRinN(1,2), max(frequencyCounts)/1.9, 'VaRN 99.9%')
xlabel('Distribution of Losses')
ylabel('Frequency')
legend('show')
title('Out of Sample Distribution of Losses for Portfolio 1')

figure(3)
% -------- Insert your code here -------- %
set(gcf, 'color', 'white');
[frequencyCounts, binLocations] = hist(Losses_inMC1*x0{1}, 50);
bar(binLocations, frequencyCounts,'DisplayName', 'In Sample Distribution');
hold on;
line([mean(VaRinMC1{1,1}) mean(VaRinMC1{1,1})], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--','DisplayName','VaRinMC1 99%');
line([mean(VaRinMC1{1,2}) mean(VaRinMC1{1,2})], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--','DisplayName','VaRinMC1 99.9%');
hold on;
normf = ( 1/(std(Losses_inMC1*x0{1})*sqrt(2*pi)) ) * exp( -0.5*((binLocations-mean(Losses_inMC1*x0{1}))/std(Losses_inMC1*x0{1})).^2 );
normf = normf * sum(frequencyCounts)/sum(normf);
plot(binLocations, normf, 'r', 'LineWidth', 1,'DisplayName','Normal Distribution');
hold on;
line([mean(VaRinN1{1,1}) mean(VaRinN1{1,1})], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '-.','DisplayName','VaRinN1 99%');
line([mean(VaRinN1{1,2}) mean(VaRinN1{1,2})], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '-.','DisplayName','VaRinN1 99.9%');
line([VaRout(1,1) VaRout(1,1)], [0 max(frequencyCounts)/2], 'Color', 'b', 'LineWidth', 1, 'LineStyle', '-','DisplayName','VaRout 99%');
line([VaRout(1,2) VaRout(1,2)], [0 max(frequencyCounts)/2], 'Color', 'b', 'LineWidth', 1, 'LineStyle', '--','DisplayName','VaRout 99.9%');
hold off;
text(0.95*mean(VaRinMC1{1,1}), max(frequencyCounts)/1.9, 'VaRinMC1 99%')
text(0.9*mean(VaRinMC1{1,2}), max(frequencyCounts)/1.9, 'VaRinMC1 99.9%')
text(0.7*mean(VaRinN1{1,1}), max(frequencyCounts)/1.9, 'VaRinN1 99%')
text(0.9*mean(VaRinN1{1,2}), max(frequencyCounts)/1.9, 'VaRinN1 99.9%')
text(1*VaRout(1,1), max(frequencyCounts)/2.2,'VaRout 99%','color','blue')
text(1*VaRout(1,2), max(frequencyCounts)/2.2,'VaRout 99.9%','color','blue')
legend('show')
xlabel('Distribution of Losses')
ylabel('Frequency')
title('Monte Carlo Simulation 1 for Portfolio 1')

figure(4)
% -------- Insert your code here -------- %
set(gcf, 'color', 'white');
[frequencyCounts, binLocations] = hist(Losses_inMC1*x0{2}, 50);
bar(binLocations, frequencyCounts,'DisplayName', 'In Sample Distribution');
hold on;
line([mean(VaRinMC1{2,1}) mean(VaRinMC1{2,1})], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--','DisplayName','VaRinMC1 99%');
line([mean(VaRinMC1{2,2}) mean(VaRinMC1{2,2})], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--','DisplayName','VaRinMC1 99.9%');
hold on;
normf = ( 1/(std(Losses_inMC1*x0{2})*sqrt(2*pi)) ) * exp( -0.5*((binLocations-mean(Losses_inMC1*x0{2}))/std(Losses_inMC1*x0{2})).^2 );
normf = normf * sum(frequencyCounts)/sum(normf);
plot(binLocations, normf, 'r', 'LineWidth', 1,'DisplayName','Normal Distribution');
hold on;
line([mean(VaRinN1{2,1}) mean(VaRinN1{2,1})], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '-.','DisplayName','VaRinN1 99%');
line([mean(VaRinN1{2,2}) mean(VaRinN1{2,2})], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '-.','DisplayName','VaRinN1 99.9%');
line([VaRout(2,1) VaRout(2,1)], [0 max(frequencyCounts)/2], 'Color', 'b', 'LineWidth', 1, 'LineStyle', '-','DisplayName','VaRout 99%');
line([VaRout(2,2) VaRout(2,2)], [0 max(frequencyCounts)/2], 'Color', 'b', 'LineWidth', 1, 'LineStyle', '--','DisplayName','VaRout 99.9%');
hold off;
text(1*mean(VaRinMC1{2,1}), max(frequencyCounts)/1.9, 'VaRinMC1 99%')
text(0.95*mean(VaRinMC1{2,2}), max(frequencyCounts)/1.9, 'VaRinMC1 99.9%')
text(0.6*mean(VaRinN1{2,1}), max(frequencyCounts)/1.9, 'VaRinN1 99%')
text(0.9*mean(VaRinN1{2,2}), max(frequencyCounts)/1.9, 'VaRinN1 99.9%')
text(1*VaRout(2,1), max(frequencyCounts)/2.2,'VaRout 99%','color','blue')
text(1*VaRout(2,2), max(frequencyCounts)/2.2,'VaRout 99.9%','color','blue')
xlabel('Distribution of Losses')
ylabel('Frequency')
legend('show')
title('Monte Carlo Simulation 1 for Portfolio 2')

figure(5)
% -------- Insert your code here -------- %
set(gcf, 'color', 'white');
[frequencyCounts, binLocations] = hist(Losses_inMC2*x0{1}, 50);
bar(binLocations, frequencyCounts,'DisplayName', 'In Sample Distribution');
hold on;
line([mean(VaRinMC2{1,1}) mean(VaRinMC2{1,1})], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--','DisplayName','VaRinMC2 99%');
line([mean(VaRinMC2{1,2}) mean(VaRinMC2{1,2})], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--','DisplayName','VaRinMC2 99.9%');
hold on;
normf = ( 1/(std(Losses_inMC2*x0{1})*sqrt(2*pi)) ) * exp( -0.5*((binLocations-mean(Losses_inMC2*x0{1}))/std(Losses_inMC2*x0{1})).^2 );
normf = normf * sum(frequencyCounts)/sum(normf);
plot(binLocations, normf, 'r', 'LineWidth', 1,'DisplayName','Normal Distribution');
hold on;
line([mean(VaRinN2{1,1}) mean(VaRinN2{1,1})], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '-.','DisplayName','VaRinN2 99%');
line([mean(VaRinN2{1,2}) mean(VaRinN2{1,2})], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '-.','DisplayName','VaRinN2 99.9%');
line([VaRout(1,1) VaRout(1,1)], [0 max(frequencyCounts)/2], 'Color', 'b', 'LineWidth', 1, 'LineStyle', '-','DisplayName','VaRout 99%');
line([VaRout(1,2) VaRout(1,2)], [0 max(frequencyCounts)/2], 'Color', 'b', 'LineWidth', 1, 'LineStyle', '--','DisplayName','VaRout 99.9%');
hold off;
text(1*mean(VaRinMC2{1,1}), max(frequencyCounts)/1.9, 'VaRinMC2 99%')
text(0.95*mean(VaRinMC2{1,2}), max(frequencyCounts)/1.9, 'VaRinMC2 99.9%')
text(0.7*mean(VaRinN2{1,1}), max(frequencyCounts)/1.9, 'VaRinN2 99%')
text(0.95*mean(VaRinN2{1,2}), max(frequencyCounts)/1.9, 'VaRinN2 99.9%')
text(1*VaRout(1,1), max(frequencyCounts)/2.2,'VaRout 99%','color','blue')
text(1*VaRout(1,2), max(frequencyCounts)/2.2,'VaRout 99.9%','color','blue')
xlabel('Distribution of Losses')
ylabel('Frequency')
legend('show')
title('Monte Carlo Simulation 2 for Portfolio 1')

figure(6)
% -------- Insert your code here -------- %
set(gcf, 'color', 'white');
[frequencyCounts, binLocations] = hist(Losses_inMC2*x0{2}, 50);
bar(binLocations, frequencyCounts,'DisplayName', 'In Sample Distribution');
hold on;
line([mean(VaRinMC2{2,1}) mean(VaRinMC2{2,1})], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--','DisplayName','VaRinMC2 99%');
line([mean(VaRinMC2{2,2}) mean(VaRinMC2{2,2})], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--','DisplayName','VaRinMC2 99.9%');
hold on;
normf = ( 1/(std(Losses_inMC2*x0{2})*sqrt(2*pi)) ) * exp( -0.5*((binLocations-mean(Losses_inMC2*x0{2}))/std(Losses_inMC2*x0{2})).^2 );
normf = normf * sum(frequencyCounts)/sum(normf);
plot(binLocations, normf, 'r', 'LineWidth', 1,'DisplayName','Normal Distribution');
hold on;
line([mean(VaRinN2{2,1}) mean(VaRinN2{2,1})], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '-.','DisplayName','VaRinN2 99%');
line([mean(VaRinN2{2,2}) mean(VaRinN2{2,2})], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '-.','DisplayName','VaRinN2 99.9%');
line([VaRout(2,1) VaRout(2,1)], [0 max(frequencyCounts)/2], 'Color', 'b', 'LineWidth', 1, 'LineStyle', '-','DisplayName','VaRout 99%');
line([VaRout(2,2) VaRout(2,2)], [0 max(frequencyCounts)/2], 'Color', 'b', 'LineWidth', 1, 'LineStyle', '--','DisplayName','VaRout 99.9%');
hold off;
text(1*mean(VaRinMC2{2,1}), max(frequencyCounts)/1.9, 'VaRinMC2 99%')
text(0.9*mean(VaRinMC2{2,2}), max(frequencyCounts)/1.9, 'VaRinMC2 99.9%')
text(0.7*mean(VaRinN2{2,1}), max(frequencyCounts)/1.9, 'VaRinN2 99%')
text(0.9*mean(VaRinN2{2,2}), max(frequencyCounts)/1.9, 'VaRinN2 99.9%')
text(1*VaRout(2,1), max(frequencyCounts)/2.2,'VaRout 99%','color','blue')
text(1*VaRout(2,2), max(frequencyCounts)/2.2,'VaRout 99.9%','color','blue')
xlabel('Distribution of Losses')
ylabel('Frequency')
legend('show')
title('Monte Carlo Simulation 2 for Portfolio 2')