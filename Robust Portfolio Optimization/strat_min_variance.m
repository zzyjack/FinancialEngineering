

function [x_optimal cash_optimal w_optimal tran_cost borrowed_fund] = strat_min_variance(x_init, cash_init, mu, Q, cur_prices)
% Optimization problem data
borrowed_fund =0;
addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio128\cplex\matlab\x64_win64');
lb = zeros(20,1);
ub = inf*ones(20,1);
A = ones(1,20);
b = 1;
% Compute minimum variance portfolio
cplex1 = Cplex('min_Variance');
cplex1.addCols(zeros(20,1), [], lb, ub);
cplex1.addRows(b, A, b);
cplex1.Model.Q = 2*Q;
cplex1.Param.qpmethod.Cur = 6; % concurrent algorithm
cplex1.Param.barrier.crossover.Cur = 1; % enable crossover
cplex1.DisplayFunc = []; % disable output to screen
cplex1.solve();
% Display minimum variance portfolio
w_optimal = cplex1.Solution.x;

x_optimal = (cur_prices*x_init)*w_optimal./cur_prices';%calculate the shares for the new weights
tran_cost = 0.005*cur_prices*abs(x_optimal-x_init);%calculate the transaction costs
cash_optimal = cash_init + cur_prices*(x_init - x_optimal)-tran_cost;%calculate the cash remaining


end



