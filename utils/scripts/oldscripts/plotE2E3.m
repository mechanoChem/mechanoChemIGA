close all; clear; clc;
s=0.1;
d=1;
x=-1.1*s:s/20:1.1*s;
% 
E4=3*d/(2*s^4); 
E3=-d/s^3; 
E2=-3*d/(2*s^2); 
%generate surface
intervals=100;
r=linspace(0,s,intervals);
theta=linspace(0,2*pi,intervals);
[c,r,t]=meshgrid(linspace(-s,s,intervals),linspace(-s,s,intervals),linspace(-s,s,intervals));
xx=r; yy=t;
%xx=r.*cos(t);
%yy=r.*sin(t);

%energy
hump=2;
alpha=(s+c)/(2*s);
f=hump*((d/s^4).*c.^4+(-2*d/s^2).*c.^2) + E2*alpha.*(xx.^2+yy.^2)+ E3*alpha.*yy.*(yy.^2-3*xx.^2) + E4*(xx.^2+yy.^2).^2;
contourslice((c+s)/(2*s),xx,yy,f,[0,0.5,1],[],[],35);
xlabel('c'); ylabel('e2'); zlabel('e3');

