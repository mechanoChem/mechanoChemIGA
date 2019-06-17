
function [dkNdxk] = dkN_dxk (k,i,p,x,X)

if (k == 0)
  dkNdxk = N_ip(i,p,x,X);
elseif (p == 0)
  dkNdxk = 0;
else
  dkNdxk = p*divZeroSafe(dkN_dxk(k-1,i,p-1,x,X),X(i+p) - X(i)) - ...
        p*divZeroSafe(dkN_dxk(k-1,i+1,p-1,x,X),(X(i+p+1) - X(i+1)));
end
  
end
