function [y] = curveInt(x,coeff,a)
%Evaluate integral

% a - 1 or 2: 1 means domain from 0 to 1; 2 means domain from 0 to 1/2 and
% R-K polynomial is scaled accordingly

l = length(x);

n = length(coeff) - 1;

A = zeros(l,n+1);

for i = 0:n-1
    if (a==1)
        A(:,n-i) = x.*(1-x).*(2*x-1).^i; %unscaled
    elseif (a==2)
        A(:,n-i) = x.*(0.5-x).*(0.5-2*x).^i; %%scaled 
    end
end
A(:,n+1) = x; %(g1 - g0)

y =  A*coeff;

end

