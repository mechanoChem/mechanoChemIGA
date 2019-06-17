% Surface using B-splines

clc
clear

%%% Load B2 data and fit
bcc = importdata('bcc_lattice_monte_carlo.txt');
eta = bcc.data(:,2:3);
cp = bcc.data(:,4:5);

p = 3;
all_MSE_train = zeros(n,1);
all_MSE_val = all_MSE_train;

n = 13; %number of interior knots
for k = 1:n
    
    % Knot vectors: 3 options
    
    %%% Use uniform
    X1 = [0 0 0 0:sqrt(2)/(k+1):sqrt(2) sqrt(2) sqrt(2) sqrt(2)];
	X2 = [0 0 0 0:.5/(k+1):.5 .5 .5 .5];
    
    %%% Use Chebyshev
%     X1 = [0 0 0 0 ChebyshevRoots(k,"Tn",[0 sqrt(2)]) sqrt(2) sqrt(2) sqrt(2) sqrt(2)];
%     X2 = [0 0 0 0 ChebyshevRoots(k,"Tn",[0 .5]) .5 .5 .5 .5];

    %%% Use genetic algorithm to optimize knot locations
%     fun = @(x) fit_dz(eta(:,1),eta(:,2),cp(:,1),cp(:,2),p,[0 0 0 0 x(1:13) sqrt(2) sqrt(2) sqrt(2) sqrt(2)],[0 0 0 0 x(14:end) .5 .5 .5 .5]);
%     cons = zeros(28,26);
%     for i=1:9
%         cons(i,i) = -1;
%         cons(i+1,i) = 1;
%         cons(i+14,i+13) = -1;
%         cons(i+15,i+13) = 1;
%     end
%     consb = zeros(20,1);
%     consb(14) = sqrt(2);
%     consb(28) = 0.5;
%     options = optimoptions('ga','MaxGenerations',100,'PlotFcn', @gaplotbestf,'Display','iter');
%     x = ga(fun,26,cons,consb,[],[],[],[],[],options);
%     X1 = [0 0 0 0 x(1:13) sqrt(2) sqrt(2) sqrt(2) sqrt(2)];
%     X2 = [0 0 0 0 x(14:end) .5 .5 .5 .5];

    [MSE,MSE_val,C] = fit_dz(eta(:,1),eta(:,2),cp(:,1),cp(:,2),p,X1,X2);
    all_MSE_train(k) = MSE;
    all_MSE_val(k) = MSE_val;
    
end

%%% Plot hyperparamater results
figure(8)
clf
semilogy(1:n,all_MSE_train,'-o')
hold on
semilogy(1:n,all_MSE_val,'-o')
legend('Training','Validation')
xlabel('Interior knots')
ylabel('MSE')

% Save bspline parameters
fid = fopen('Cpts.txt','w');
fprintf(fid,'# %d %d %d\n',p,size(C,1),size(C,2));
fclose(fid);
dlmwrite('Cpts.txt',C,'-append','delimiter',' ','precision',15);

fid = fopen('knots.txt','w');
fprintf(fid,'# %d %d %d\n',p,size(X1,2),size(X2,2));
fclose(fid);
dlmwrite('knots.txt',X1,'-append','delimiter',' ','precision',15);
dlmwrite('knots.txt',X2,'-append','delimiter',' ','precision',15);

% Plot fits
[x1,x2] = meshgrid(0:0.01:1,0:0.01:1);
e1 = (x1+x2)./sqrt(2);
e2 = (x1-x2)./sqrt(2);
Z = z(e1,e2.^2,p,X1,X2,C);
CP1 = dz_dx(e1,e2.^2,p,X1,X2,C);
CP2 = 2*e2.*dz_dy(e1,e2.^2,p,X1,X2,C);

figure(1)
clf
surf(e1,e2,Z,'EdgeColor','interp')
xlabel('c')
ylabel('\eta')
title('g(c,\eta) (B-splines)')
set(gca,'FontSize',20)
xlim([0 2/sqrt(2)])
ylim([-1/sqrt(2) 1/sqrt(2)])
view(90,0)
% view(0,0)

figure(2)
clf
surf(e1,e2,CP1,'EdgeColor','interp')
hold on
p = scatter3(eta(:,1),eta(:,2),cp(:,1),5,[0.7 0.7 0.7],'filled');
p.MarkerFaceAlpha = 0.4;
xlabel('c')
ylabel('\eta')
title('\partial{g}/\partial{c} (B-splines)')
xlim([0 2/sqrt(2)])
ylim([-1/sqrt(2) 1/sqrt(2)])
zlim([-1.,1.])
set(gca,'FontSize',20)
view(-40,20)
colormap default

figure(3)
clf
surf(e1,e2,CP2,'EdgeColor','interp')
hold on
p = scatter3(eta(:,1),eta(:,2),cp(:,2),5,[0.7 0.7 0.7],'filled');
p.MarkerFaceAlpha = 0.4;
xlabel('c')
ylabel('\eta')
title('\partial{g}/\partial{\eta} (B-splines)')
xlim([0 2/sqrt(2)])
ylim([-1/sqrt(2) 1/sqrt(2)])
zlim([-.25 .25])
set(gca,'FontSize',20)
view(-40,20)
colormap default
