
function [f] = dz_dy (x,y,p,X1,X2,C)

f = 0;
for k = 1:(length(X1) - p - 1)
  for l = 1:(length(X2) - p - 1)

    f = f + C(k,l)*N_ip(k,p,x,X1).*dkN_dxk(1,l,p,y,X2);
    
  end
end

end
