
function [dNdx] = dN_dx (i,p, x,X)

if (p == 0)
  dNdx = 0;

else
  dNdx = p*divZeroSafe(N_ip(i,p-1,x,X),X(i+p) - X(i)) - ...
        p*divZeroSafe(N_ip(i+1,p-1,x,X),(X(i+p+1) - X(i+1)));

end
  
end
