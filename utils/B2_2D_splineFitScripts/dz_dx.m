
function [f] = dz_dx (x,y,p,X1,X2,C)

f = 0;
for j = 1:(length(X1) - p - 1)
  for l = 1:(length(X2) - p - 1)

    f = f + C(j,l)*dkN_dxk(1,j,p,x,X1).*N_ip(l,p,y,X2);
    
  end
end

end
