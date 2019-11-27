function  [x_optimal cash_optimal w_optimal interest_payment] = strat_lever_equal_risk_contr(curr_positions, curr_cash, mu, Q, cur_prices, w_init, period)
global Q
global A_ineq
global A_eq


% total number of stocks
n =20;

% determine which risk-free rate to use

r_rf = 0.045;



curr_portf_value = cur_prices * curr_positions;
borrowed_money = curr_portf_value;
curr_portf_value = curr_portf_value + borrowed_money;

interest_payment = borrowed_money * r_rf/6;


% Equality constraints
A_eq = ones(1,n);
b_eq = 1;

% Inequality constraints
A_ineq = [];
b_ineql = [];
b_inequ = [];

% Define initial portfolio 
w0 = w_init;

options.lb = zeros(1,n);       % lower bounds on variables
options.lu = ones (1,n);       % upper bounds on variables
options.cl = [b_eq' b_ineql']; % lower bounds on constraints
options.cu = [b_eq' b_inequ']; % upper bounds on constraints

% Set the IPOPT options
options.ipopt.jac_c_constant        = 'yes';
options.ipopt.hessian_approximation = 'limited-memory';
options.ipopt.mu_strategy           = 'adaptive';
options.ipopt.tol                   = 1e-10;
options.ipopt.print_level           = 0;

% The callback functions
funcs.objective         = @computeObjERC;
funcs.constraints       = @computeConstraints;
funcs.gradient          = @computeGradERC;
funcs.jacobian          = @computeJacobian;
funcs.jacobianstructure = @computeJacobian;

% Run IPOPT
[wsol info] = ipopt(w0',funcs,options);

% Make solution a column vector
if(size(wsol,1)==1)
    w_erc = wsol';
else
    w_erc = wsol;
end


% compute variance and asset risk contributions for the ERC portfolio
std_ERC = sqrt(w_erc' * Q * w_erc);
RC_ERC = (w_erc .* (Q * w_erc)) / sqrt(w_erc' * Q * w_erc);

w_optimal = w_erc; % re-balanced weights of each asset
curr_portf_value = cur_prices * curr_positions;
x_optimal = w_optimal * curr_portf_value ./ cur_prices'; % current position
x_optimal = round(x_optimal,0); % round the x_optimal to nearest integer
borrowed_money = cur_prices * x_optimal; 
total_value = borrowed_money + curr_portf_value + curr_cash;
interest_payment = borrowed_money * (r_rf/6);
trade_volume = abs(x_optimal - curr_positions); % trade volume
trans_cost = 0.005 * cur_prices * trade_volume; % transaction cost
% the remaining goes to cash account
cash_optimal = curr_cash + borrowed_money + cur_prices * (curr_positions - total_value * w_optimal./ cur_prices') - (trans_cost + interest_payment); 

end

