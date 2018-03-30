function coeff = curveFit(x,y,n,a,varargin)
%Weighted curve fit according to Redlich-Kister polynomials (n - the polynomial degree - must be at least 2 -
%quadratic)

% x - x-coordinates of data to be fitted

% y - y-coordinates of data to be fitted

% n - polynomial degree

% a - 1 or 2: 1 means domain from 0 to 1; 2 means domain from 0 to 1/2 and
% R-K polynomial is scaled accordingly

% w - (optional) vector of weights, must be equal in lengh to vectors x & y

if length(varargin) == 0
    w = ones(length(y),1);
elseif length(varargin) == 1
    w = varargin{1};
    if(length(w) ~= length(y))
        error('w should be the same length as y');
    end
end

l = length(x);

A = zeros(l,n+1);

for i = 0:n-1
    if (a==1)
        A(:,n-i) = -(2*x - 1).^(i+1) + 2*i*x.*(1-x).*(2*x-1).^(i-1); %%unscaled
    elseif (a==2)
        A(:,n-i) = (0.5 - 2*x).^(i+1) - 2*i*x.*(0.5-x).*(0.5-2*x).^(i-1); %%scaled 
    end
end
A(:,n+1) = ones(l,1); %(g1 - g0)

%Weighting
A = diag(w)*A;
y = diag(w)*y;

[q,r,p] = qr(A,0);
coeff = r\(q'*y);
coeff(p) = coeff;

end

