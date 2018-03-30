%Greg Teichert, University of Michigan
%Matlab script to fit chemical potential data with Redlich-Kister
%polynomial expansion.

clear
addpath('RedlichKisterFit');

%Load chemical potential data (1st column: composition, 2nd column:
%chemical potential
chemPotential = load('inputData/TiO_chem_pot.txt');
%Load free energy data - this was numerically integrated from the chemical
%potential data and is only used for a comparison in the plots
free_energy = load('inputData/TiO_free_energy.txt');

%Polynomial degree
order = 15;

%Domain over [0,1] or scaled to [0,1/2]
a = 2; %1 or 2, used in the log term: k_B*T*(x.*log(x) + 1/a*(1 - a*x).*log(1 - a*x))

k_B = 8.6173324e-5; %Boltzmann's constant, Ev per K
T = 800; %Kelvin

xp = chemPotential(:,1);
yp = chemPotential(:,2);

xf = free_energy(:,1);
yf = free_energy(:,2);

%Cut off ends
yp = yp(xp>=0.001 | xp<=0.499);
xp = xp(xp>=0.001 | xp<=0.499);

%Do running average on both x and y for numerical derivative
smoothX = smooth(xp);
smoothY = smooth(yp);

l = length(smoothX);

%Find numerical derivative (for plotting)
der_x = 0.5*(smoothX(1:l-1) + smoothX(2:l));
der_y = (smoothY(2:l) - smoothY(1:l-1))./(smoothX(2:l) - smoothX(1:l-1));

%Weight data points in the two phase region (specify manually the range of
%the two phase regions)
w = ones(l,1);
for i = 1:l
    if((smoothX(i) > 0.11 && smoothX(i) < 0.15) || (smoothX(i) > 0.27 && smoothX(i) < 0.29) || (smoothX(i) > 0.45 && smoothX(i) < 0.465))
        w(i) = 10;
    end
end

fitX = xp;
fitY = yp - k_B*T*(log(a*fitX) - log(1 - a*fitX));

%Use Redlich-Kister polynomials
x = 0:0.0001:0.5;
pL = curveFit(fitX, fitY,order,a,w);
valL = curveVal(x,pL,a) + k_B*T*(log(a*x') - log(1 - a*x'));

%Differentiate Redlich-Kister
valL_der = curveDer(x,pL,a) + k_B*T./(x'.*(1 - a*x'));

%Integrate Redlich-Kister (plot w.r.t. end members at x = 0 and 0.5
valL_int = curveInt(x,pL,a);
valL_end = valL_int(end);
if(a==1)
    valL_int = valL_int' + 1/a*k_B*T*(a*x.*log(a*x) + (1 - a*x).*log(1 - a*x))...
         - 2*x*(valL_int(end) + 1/a*k_B*T*(a*0.5*log(a*0.5) + (1 - a*0.5).*log(1 - a*0.5)));
elseif(a==2)
    valL_int = valL_int' + 1/a*k_B*T*(a*x.*log(a*x) + (1 - a*x).*log(1 - a*x)) - 2*x*(valL_int(end));
end

%Write out polynomial coefficients (for chemical potential)
fileID = fopen('outputData/polyCoeffs.txt','w');
fprintf(fileID,'%.10e\n',round(pL,11,'significant'));
fclose(fileID);

%Plot polynomial fit
figure(1)
clf
subplot(1,3,1)
scatter(xp,yp,'.','r')
hold on
hold on
plot(x,valL,'linewidth',2,'color','b')
axis([0 0.5 -1.5 0.5])
legend('Data','Curve Fit')
title('Chemical potential')
xlabel('Composition(N(O)/N(Ti))')
ylabel('eV/N(O)')
grid on

%Plot derivative of chemical potential
subplot(1,3,2)
scatter(der_x,der_y,'.','r')
hold on
hold on
plot(x,valL_der,'linewidth',2,'color','b')
title('Derivative of chemical potential')
xlabel('Composition(N(O)/N(Ti))')
ylabel('eV N(Ti)/N(O)^2')
grid on
legend('Data','Curve fit')
axis([0 0.5 -10 10])

%Plot free energy
subplot(1,3,3)
scatter(xf,yf,'.','r')
hold on
%plot(x,val_int,'linewidth',2,'color','b')
hold on
plot(x,valL_int,'linewidth',2,'color','b')%,'color',[0 0.75 0])
axis([0 0.5 -0.06 0])
dim = [.5 .6 .3 .3];
legend('Data','Curve fit')
title('Free energy')
xlabel('Composition(N(O)/N(Ti))')
ylabel('eV/N(Ti)')
grid on
