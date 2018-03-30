%Greg Teichert, University of Michigan
%Octave script to fit chemical potential data with a spline

clear
clc

function [smoothx] = smoothO (x)

%Find the running average, default of 5 (2 points on either side).

l = length(x);

smoothx = x;

smoothx(2) = 1/3*(x(1) + x(2) + x(3));
smoothx(l-1) = 1/3*(x(l-2) + x(l-1) + x(l));
smoothx(3:l-2) = 1/5*(x(3:l-2) + x(2:l-3) + x(1:l-4) + x(4:l-1) + x(5:l));

endfunction

chemPotential = load('inputData/TiO_chem_pot.txt');
freeEnergy = load('inputData/TiO_free_energy.txt');

degree = 3; %Note that Octave defines the spline order to be equal to the polynomial degree

a = 2; %1 or 2, used in the log term: k_B*T*(x.*log(x) + 1/a*(1 - a*x).*log(1 - a*x))

k_B = 8.6173324e-5; %Boltzmann's constant, Ev per K
T = 800; %Kelvin

xp = chemPotential(:,1);
yp = chemPotential(:,2);

xf = freeEnergy(:,1);
yf = freeEnergy(:,2);

%Cut off noisy end data
yp = yp(xp>=0.001 & xp<=0.499);
xp = xp(xp>=0.001 & xp<=0.499);

%Smooth with running average for numerical derivative
smoothX = smoothO(xp);
smoothY = smoothO(yp);

%Find numerical derivative (central difference)
l = length(xp);
der_x = 0.5*(smoothX(1:l-1) + smoothX(2:l));
der_y = (smoothY(2:l) - smoothY(1:l-1))./(smoothX(2:l) - smoothX(1:l-1));

%Initially, use a uniform grid with many points to approximate the 2nd derivative
%Spline fit of initial grid
fitX = 0:0.005:0.5;
pp = splinefit(xp,yp-k_B*T*(log(a*xp) - log(1 - a*xp)),fitX,"order",degree);

%Find the 2nd derivative of the spline fit
x = 0:0.004:0.5;
ppd2 = ppder(pp,2);
vals_der2 = ppval(ppd2,x);

%Use somewhat sparse but uniformly distributed points...
%Plus points with a high 2nd derivative
fitX = [0:0.03:0.5 x(abs(vals_der2)>600)];

%Fit data with these knots
pp = splinefit(xp,yp-k_B*T*(log(a*xp) - log(1 - a*xp)),fitX,"order",degree);

%Evaluate and plot function, derivatives, and integral (including log term)
x1 = 0:0.0001:0.5;
vals = ppval(pp,x1) + k_B*T*(log(a*x1) - log(1 - a*x1));

%Find the derivative of the spline
ppd = ppder(pp);
vals_der = ppval(ppd,x1) + k_B*T./(x1.*(1 - a*x1));

%Find the integral (with respect to 0 and 1/2)
ppi = ppint(pp);
vals_int = ppval(ppi,x1);
vals_int = vals_int + 1/a*k_B*T*(a*x1.*log(a*x1) + (1 - a*x1).*log(1 - a*x1))...
 - 2*x1*(vals_int(end));
if(a ~= 2)
  vals_int = vals_int - 2*x1*(1/a*k_B*T*(a*0.5*log(a*0.5) + (1 - a*0.5).*log(1 - a*0.5)));
end
  
%Write spline breaks and coefficients to file
fileID = fopen('outputData/splineCoeffs.txt','w');
fprintf(fileID,'%.10e %.10e %.10e %.10e\n',pp.coefs');
fclose(fileID);

fileID = fopen('outputData/splineBreaks.txt','w');
fprintf(fileID,'%.10e\n',pp.breaks);
fclose(fileID);
 
%Plot spline fit
figure(1)
clf

%Plot the free energy
subplot(1,3,3)
scatter(xf,yf,'r','.')
hold on
plot(x1,vals_int,'linewidth',2)
title('Free energy')
xlabel('Composition(N(O)/N(Ti))')
ylabel('(eV/N(Ti))')
grid on
legend('Data','Spline')
axis([0 0.5 -0.05 0])

%Plot derivative of chemical potential
%For plotting purposes, since Octave plots out-of-range points, remove points
%with der_y > 10
subplot(1,3,2)
scatter(der_x(abs(der_y) <= 10),der_y(abs(der_y) <= 10),'r','.')
hold on
plot(x1,vals_der,'linewidth',2)
title('Derivative')
xlabel('Composition(N(O)/N(Ti))')
ylabel('(eV*N(Ti)/N(O)^2)')
grid on
legend('Data','Spline')
axis([0 0.5 -10 10])

%Plot spline fit
subplot(1,3,1)
title("Chemical potential")
xlabel("Composition(N(O)/N(Ti))")
ylabel("mu_O(eV/N(O))")
hold on
scatter(xp,yp,'r','.')
hold on
plot(x1,vals,"linewidth",3)
grid on
legend('Data','Spline')
axis([0 0.5 -2 0.5])
