


x0(1) = 0.2

for i = 1:5
    x1(i) = x0(i)-(7-1/x0(i))/((x0(i))^(-2))
    y(i)=(7-1/x0(i))/((x0(i))^(-2))
   
    x0(i+1) = x1(i)
    
end
