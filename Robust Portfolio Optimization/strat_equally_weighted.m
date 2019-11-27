
function [x_optimal cash_optimal w_optimal tran_cost borrowed_fund] = strat_equally_weighted(x_init, cash_init, mu, Q, cur_prices)
         
% Equally weighted portfolio implies wi's is 1/20 in this case
w_optimal = ones(20,1)*1/20;
%initialize x_optimal
x_optimal = (cur_prices*x_init)*w_optimal./cur_prices';
%Calculate the transaction costs
tran_cost = 0.005*cur_prices*abs(x_optimal-x_init);
cash_optimal = cash_init + cur_prices*(x_init - x_optimal)-tran_cost;
borrowed_fund =0;
end

    

