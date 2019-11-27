function gval = computeGradERC (x)
global Q
n = size(Q,1) ;  
  if(size(x,1)==1)
     x = x';
  end
  % Insert your gradiant computations here
  % You can use finite differences to check the gradient
  gval = zeros(n,1); 
  % risk measure
  f = x.* (Q * x);
  for i = 1:n
      for j = 1:n
          f1 = Q(i,:) * x + Q(i,i) * x(i);
          f2 = Q(i,j) * x(i);
          delta = (f(i)-f(j)) * (f1 - f2);
          gval(i,:) = gval(i,:) + delta;
      end
      gval(i,:) = 2*2 *  gval(i,:);
  end
end
