function [MSE,MSE_val,C] = fit_dz (x,y,dzdx,dzdy,p,X1,X2)

x = x(:);
y = y(:);
dzdx = dzdx(:);
dzdy = dzdy(:);

%%% Split into training and validation data
rng(1);
shuffled = randperm(length(x));
split = round(0.75*length(x));
i_train = shuffled(1:split);
i_val = shuffled(split+1:end);

x_train = x(i_train);
y_train = y(i_train);
dzdx_train = dzdx(i_train);
dzdy_train = dzdy(i_train);

x_val = x(i_val);
y_val = y(i_val);
dzdx_val = dzdx(i_val);
dzdy_val = dzdy(i_val);

%%% Create the matrix X, each row corresponds to a data point
m = (length(X1) - p - 1);
n = (length(X2) - p - 1);
Ax_train = zeros(length(x_train),m*n);
Ay_train = Ax_train;
Ax_val = zeros(length(x_val),m*n);
Ay_val = Ax_val;
for i = 1:m
  for j = 1:n
      I = (n*(i-1)+(j-1))+1;
      %Ax(:,I) = dN_dx(i,p,x,X1).*N_ip(j,p,y,X2);
      %Ay(:,I) = N_ip(i,p,x,X1).*dN_dx(j,p,y,X2); 
      Ax_train(:,I) = dN_dx(i,p,x_train,X1).*N_ip(j,p,y_train.^2,X2); %To impose symmetry on y
      Ay_train(:,I) = 2*y_train.*N_ip(i,p,x_train,X1).*dN_dx(j,p,y_train.^2,X2); %To impose symmetry on y
      Ax_val(:,I) = dN_dx(i,p,x_val,X1).*N_ip(j,p,y_val.^2,X2); %To impose symmetry on y
      Ay_val(:,I) = 2*y_val.*N_ip(i,p,x_val,X1).*dN_dx(j,p,y_val.^2,X2); %To impose symmetry on y
  end
end

%%% Solve the least squares problem
c = (Ax_train'*Ax_train + Ay_train'*Ay_train + 1.e-4*eye(m*n))\(Ax_train'*dzdx_train + Ay_train'*dzdy_train); %A little regularization to improve the matrix condition number
C = reshape(c,[n m]);
C = C';

%%% Return the MSE as well
MSE = ((Ax_train*c - dzdx_train)'*(Ax_train*c - dzdx_train)+(Ay_train*c - dzdy_train)'*(Ay_train*c - dzdy_train))/size(x_train,1);

%%% Return the MSE as well
MSE_val = ((Ax_val*c - dzdx_val)'*(Ax_val*c - dzdx_val)+(Ay_val*c - dzdy_val)'*(Ay_val*c - dzdy_val))/size(x_val,1);

end