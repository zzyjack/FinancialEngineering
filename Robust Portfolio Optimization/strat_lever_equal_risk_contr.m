function [x_optimal cash_optimal w_optimal tran_cost borrowed_fund] = strat_lever_equal_risk_contr(x_init, cash_init, mu, M, cur_prices)
warning('off','all')
global Q A_ineq A_eq period

Q = M;
addpath('C:\Program Files\MatlabInterface.site');
r = 0.045;
n=20;
w0 = (cur_prices'.*x_init)/(cur_prices*x_init);

% Equality constraints
A_eq = ones(1,n);
b_eq = 1;

% Inequality constraints
A_ineq = [];
b_ineql = [];
b_inequ = [];

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

  
%% Run IPOPT
[wsol info] = ipopt(w0',funcs,options);

% Make solution a column vector
if(size(wsol,1)==1)
    w_optimal = wsol';
else
    w_optimal = wsol;
end

% Compute return, variance and risk contribution for the ERC portfolio
ret_ERC = dot(mu, w_optimal);
var_ERC = w_optimal'*Q*w_optimal;
RC_ERC = (w_optimal .* ( Q*w_optimal )) / sqrt(w_optimal'*Q*w_optimal);
% fprintf('\n\nAsset risk contributions for ERC: \n')
% [RC_ERC]
    
x_optimal = (cur_prices*x_init)*w_optimal./cur_prices';%calculate the shares for the new weights
borrowed_fund = cur_prices*x_optimal; % leverage 200% on the shares calculated for this period
total_value = borrowed_fund + (cur_prices*x_init) + cash_init;
interest = borrowed_fund * (r/6); %calculate the interest
tran_cost = 0.005*cur_prices*abs(total_value*w_optimal./cur_prices'-x_init);%calculate the transaction costs
cash_optimal = cash_init + borrowed_fund + cur_prices*(x_init - total_value*w_optimal./cur_prices')-tran_cost-interest;%calculate the cash remaining
% x_optimal= x_optimal - borrowed_fund*w_optimal; %calculate the shares left after deducting the fund borrowed.

end
