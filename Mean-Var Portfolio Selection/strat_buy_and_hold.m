function  [x_optimal cash_optimal w_optimal tran_cost] = strat_buy_and_hold(x_init, cash_init, mu, Q, cur_prices)

   x_optimal = x_init;
   cash_optimal = cash_init;
  w_optimal = x_optimal.*cur_prices'./(cur_prices*x_optimal);
   tran_cost = 0;
end