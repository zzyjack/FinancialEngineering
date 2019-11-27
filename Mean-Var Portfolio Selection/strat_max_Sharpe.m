
function  [x_optimal cash_optimal w_optimal tran_cost] = strat_max_Sharpe(x_init, cash_init, mu, Q, cur_prices)

    addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio128\cplex\matlab\x64_win64');
    % total number of stocks
    n = 20; 
    % risk-free rate
    r_rf = 0.025;
    
    % modify the Q matrix to 21*21
    Q=[Q;zeros(1,n)];
    Q=[Q zeros(n+1,1)];
   
    % optimization problem data
    lb = zeros(n+1,1);
    ub = inf * ones(n+1,1);
    A = [[(mu'-r_rf/252) 0];[ones(1,n) -1]; [eye(n-1,n) -1*ones(n-1,1)]];
    lhs = [1; 0; -inf * ones(n-1,1)];
    rhs = [1; 0; zeros(n-1,1)];
    
    y = ones(n+1,1);
     
    % compute maximum sharpe ratio portfolio
    cplex1 = Cplex('max_Sharpe');
    cplex1.Model.sense = 'minimize';
    cplex1.addCols(zeros(n+1,1), [], lb, ub);
    cplex1.addRows(lhs, A, rhs);
    cplex1.Model.Q = 2 * Q;
    cplex1.Param.qpmethod.Cur = 6; % concurrent algorithm
    cplex1.Param.barrier.crossover.Cur = 1;% enable crossover
    cplex1.DisplayFunc = []; % disable output to screen
    cplex1.solve();
    
    y = cplex1.Solution.x(1:n); % re-balanced weights of each asset
    k = cplex1.Solution.x(n+1);
    w_optimal = y / k;
    
    x_optimal = (cur_prices*x_init)*w_optimal./cur_prices';%calculate the shares for the new weights
    tran_cost = 0.005*cur_prices*abs(x_optimal-x_init);%calculate the transaction costs
    cash_optimal = cash_init + cur_prices*(x_init - x_optimal)-tran_cost;%calculate the remaining cash
    




end


