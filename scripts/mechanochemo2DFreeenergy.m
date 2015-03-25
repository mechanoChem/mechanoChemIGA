close all; clear; figure; 
%Chemical free energy parameters
Cs=1.0; Cd=-2.0;
C4=(-16*Cd/Cs^4);
C3=(32*Cd/Cs^3);
C2=(-16*Cd/Cs^2);
%Mechanics free energy parameters
Es=0.3; Ed=-1.0;
E4=(-Ed/Es^4);
E3=0.0;
E2=(2*Ed/Es^2);

%Plot chemical part of the free energy
samplePts=20;
fac=1.1;
c=(1-fac)*Cs:Cs/samplePts:fac*Cs;
Fc=C4*c.^4+C3*c.^3+C2*c.^2;
plot (c,Fc);

%Plot mechanical part of the free energy
samplePts=20;
fac=1.1;
e=-fac*Es:Es/samplePts:fac*Es;
plusC=1;
Fu=E4.*e.^4 + plusC*E2.*e.^2;
hold on;
plot (e,Fu);

%composite manifold
scaleF=10;
samplePts=250; 
facc=1.01; %fac>1
face=1.0;
[c, e]=meshgrid((facc-1)*Cs:Cs/samplePts:facc*Cs,-face*Es:Es/samplePts:face*Es);
figure;
%convex manifold ploted at an offset
plusC=(2*0-Cs)./Cs;
offset=5.0;
F=C2*(c-0.5).^2/5 + E4.*e.^4 + plusC*E2.*e.^2 + offset;
F=F/scaleF;
h1=surf(c,e,F,'facecolor','interp','LineStyle','none');
hold on;
%non convex manifold
plusC=(2*c-Cs)./Cs;
F=C4*c.^4+C3*c.^3+C2*c.^2 + E4.*e.^4 + plusC*E2.*e.^2;
F=F/scaleF;
h2=surf(c,e,F,'facecolor','interp','LineStyle','none');
xlabel('c'); 
ylabel('e_2');
zlabel('F'); set(gca, 'ZTick', []);
set(h1,'CData',get(h1,'CData'));
set(h2,'CData',get(h2,'CData')+(offset+1)/scaleF);
axis equal; grid off;
minF1=zeros(size(c,2),3);
minF2=zeros(size(c,2),3);
alpha(.7)
for i=1:size(c,2)
    %path 1
    [val,ind]=min(F(:,i)+e(:,i)*0.01);
    minF1(i,1)=c(ind,i);
    minF1(i,2)=e(ind,i);
    minF1(i,3)=val;
    %path 2
    [val,ind]=min(F(:,i)-e(:,i)*0.01);
    minF2(i,1)=c(ind,i);
    minF2(i,2)=e(ind,i);
    minF2(i,3)=val;
end
plot3(minF1(:,1),minF1(:,2),minF1(:,3),'LineWidth',3);
plot3(minF2(:,1),minF2(:,2),minF2(:,3),'LineWidth',3);

