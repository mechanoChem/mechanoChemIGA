
function [N] = N_ip (i,p, x,X)

if (p == 0)
  N = (X(i) <= x & x < X(i+1));

else
  N = (x - X(i)).*divZeroSafe(N_ip(i,p-1,x,X),X(i+p) - X(i)) + ...
        (X(i+p+1) - x).*divZeroSafe(N_ip(i+1,p-1,x,X),(X(i+p+1) - X(i+1)));

end
  
end
