
function [c] = divZeroSafe (a,b)
% Assuming a can be a vector/matrix, but b is a scalar

if (b==0)
    c = zeros(size(a));
else
    c = a./b;
    
% b = b + (b == 0);
% c = a./b;
  
end
