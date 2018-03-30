function [y] = curveDer(x,coeff,a)
%Evaluate curve fit with given coefficients

% a - 1 or 2: 1 means domain from 0 to 1; 2 means domain from 0 to 1/2 and
% R-K polynomial is scaled accordingly

k_B = 8.6173324e-5; %Boltzmann's constant, Ev per K
T = 800; %Kelvin

l = length(x);

n = length(coeff) - 1;

A = zeros(l,n+1);

for i = 0:n-1
    if (a==1)
        A(:,n-i) = -(4*i+2)*(2*x - 1).^i + 4*i*(i-1)*x.*(1-x).*(2*x - 1).^(i-2); %unscaled
    elseif (a==2)
        A(:,n-i) = -2*(2*i+1)*(0.5 - 2*x).^i + 4*i*(i-1)*x.*(0.5-x).*(0.5-2*x).^(i-2); %%scaled 
    end
end
A(:,n+1) = zeros(l,1); %(g1 - g0)

y =  A*coeff;

end

