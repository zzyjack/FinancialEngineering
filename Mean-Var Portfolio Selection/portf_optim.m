clc;
clear all;
format long

% Input files
input_file_prices  = 'Daily_closing_prices.csv';

% Read daily prices
if(exist(input_file_prices,'file'))
  fprintf('\nReading daily prices datafile - %s\n', input_file_prices)
  fid = fopen(input_file_prices);
     % Read instrument tickers
     hheader  = textscan(fid, '%s', 1, 'delimiter', '\n');
     headers = textscan(char(hheader{:}), '%q', 'delimiter', ',');
     tickers = headers{1}(2:end);
     % Read time periods
     vheader = textscan(fid, '%[^,]%*[^\n]');
     dates = vheader{1}(1:end);
  fclose(fid);
  data_prices = dlmread(input_file_prices, ',', 1, 1);
else
  error('Daily prices datafile does not exist')
end

% Convert dates into array [year month day]
format_date = 'mm/dd/yyyy';
dates_array = datevec(dates, format_date);
dates_array = dates_array(:,1:3);

% Find the number of trading days in Nov-Dec 2014 and
% compute expected return and covariance matrix for period 1
day_ind_start0 = 1;
day_ind_end0 = length(find(dates_array(:,1)==2014));
cur_returns0 = data_prices(day_ind_start0+1:day_ind_end0,:) ./ data_prices(day_ind_start0:day_ind_end0-1,:) - 1;
mu = mean(cur_returns0)';
Q = cov(cur_returns0);

% Remove datapoints for year 2014
data_prices = data_prices(day_ind_end0+1:end,:);
dates_array = dates_array(day_ind_end0+1:end,:);
dates = dates(day_ind_end0+1:end,:);

% Initial positions in the portfolio
init_positions = [5000 950 2000 0 0 0 0 2000 3000 1500 0 0 0 0 0 0 1001 0 0 0]';

% Initial value of the portfolio
init_value = data_prices(1,:) * init_positions;
fprintf('\nInitial portfolio value = $ %10.2f\n\n', init_value);

% Initial portfolio weights
w_init = (data_prices(1,:) .* init_positions')' / init_value;

% Number of periods, assets, trading days
N_periods = 6*length(unique(dates_array(:,1))); % 6 periods per year
N = length(tickers);
N_days = length(dates);

% Annual risk-free rate for years 2015-2016 is 2.5%
r_rf = 0.025;

% Number of strategies
strategy_functions = {'strat_buy_and_hold' 'strat_equally_weighted' 'strat_min_variance' 'strat_max_Sharpe'};
strategy_names     = {'Buy and Hold' 'Equally Weighted Portfolio' 'Mininum Variance Portfolio' 'Maximum Sharpe Ratio Portfolio'};
N_strat = 1; % comment this in your code
% N_strat = length(strategy_functions); % uncomment this in your code
fh_array = cellfun(@str2func, strategy_functions, 'UniformOutput', false);

for (period = 1:N_periods)
   % Compute current year and month, first and last day of the period
   if(dates_array(1,1)==15)
       cur_year  = 15 + floor(period/7);
   else
       cur_year  = 2015 + floor(period/7);
   end
   cur_month = 2*rem(period-1,6) + 1;
   day_ind_start = find(dates_array(:,1)==cur_year & dates_array(:,2)==cur_month, 1, 'first');
   day_ind_end = find(dates_array(:,1)==cur_year & dates_array(:,2)==(cur_month+1), 1, 'last');
   fprintf('\nPeriod %d: start date %s, end date %s\n', period, char(dates(day_ind_start)), char(dates(day_ind_end)));

   % Prices for the current day
   cur_prices = data_prices(day_ind_start,:);

   % Execute portfolio selection strategies
   for(strategy = 1:4)

      % Get current portfolio positions
      if(period==1)
         curr_positions = init_positions;
         curr_cash = 0;
         portf_value{strategy} = zeros(N_days,1);
      else
         curr_positions = x{strategy,period-1};
         curr_cash = cash{strategy,period-1};
      end

      % Compute strategy
      [x{strategy,period} cash{strategy,period} weight{strategy,period} tran_cost{strategy,period}] = fh_array{strategy}(curr_positions, curr_cash, mu, Q, cur_prices);

      % Verify that strategy is feasible (you have enough budget to re-balance portfolio)
      % Check that cash account is >= 0
      % Check that we can buy new portfolio subject to transaction costs
      % Transaction costs and cash accumulated are computed using the following while loop, the cash is accumulated
      % until the transaction costs can be paid off then resulting a positive cash account.
      cash_deposit = 1;
          while cash{strategy,period} < 0
              %calculate the shares to buy and sell if some of the cash are
              %to be deposited in cash account to pay off the transaction
              %costs.
              x{strategy,period} = (cur_prices*curr_positions-cash_deposit)*weight{strategy,period}./cur_prices';
              x{strategy,period} = round(x{strategy,period},0);
              x_rebalance = x{strategy,period}-curr_positions;
              %calculate the transaction costs
              tran_cost{strategy,period} = 0.005*cur_prices*abs(x_rebalance);
              %calculate the cash left after transaction costs are paid
              cash{strategy,period} = curr_cash + cur_prices*(curr_positions - x{strategy,period})-tran_cost{strategy,period};
              %add increments to the cash deposit if cash account is still
              %negative
              cash_deposit = cash_deposit+0.1;
      end
           if cash{strategy,period}<0
           fprintf('insufficient budget')
           break
      end
      


      % Compute portfolio value
      portf_value{strategy}(day_ind_start:day_ind_end) = data_prices(day_ind_start:day_ind_end,:) * x{strategy,period} + cash{strategy,period};

      fprintf('   Strategy "%s", value begin = $ %10.2f, value end = $ %10.2f, cash = $ %4.2f£¬transaction cost = $ %4.2f\n', char(strategy_names{strategy}), portf_value{strategy}(day_ind_start), portf_value{strategy}(day_ind_end), cash{strategy,period},tran_cost{strategy,period});

   end
      
   % Compute expected returns and covariances for the next period
   cur_returns = data_prices(day_ind_start+1:day_ind_end,:) ./ data_prices(day_ind_start:day_ind_end-1,:) - 1;
   mu = mean(cur_returns)';
   Q = cov(cur_returns);
   
end
 figure()
 N_days = [1:504];
for i =1:4
    plot(N_days,portf_value{1,i});
    hold on
end
title('Portfolio value over time')
legend(strategy_names{1,1:4});
xlabel('days')
ylabel('value')

% manipulate the weight data structure for plot
for m = 1:4
    for i =1:12
         new_weight{1,m}(i,:)= weight{m,i};
    end
    new_weight{1,m}=[w_init';new_weight{1,m}];
end
stocks =[1:20];
period=[0:12];
for f = 1:4
    figure()
    for i = 1:20
    plot(period, new_weight{1,f}(:,i))
    hold on 
    end
    title(sprintf('Change of Weights over period for %s',strategy_names{1,f}))
    xlabel('period')
    ylabel('weight')
    legend(string(stocks))
    hold off
end
  
