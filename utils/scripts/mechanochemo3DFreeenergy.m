close all; clear; clc;
Es=0.1; Cs=1.0;
Ed=0.1;   Cd=-1.0;
% 
C4=(-16*Cd/Cs^4);
C3=(32*Cd/Cs^3);
C2=(-16*Cd/Cs^2);
%
E4=3*Ed/(2*Es^4); 
E3=Ed/Es^3; 
E2=-3*Ed/(2*Es^2); 
%generate surface
intervals=100;
[c, r, t]=meshgrid(linspace(-Cs*0.0,Cs*1.0,intervals),linspace(-Es,Es,intervals),linspace(-Es,Es*1.1,intervals));
e2=r; e3=t;
%xx=rr.*cos(tt);
%yy=rr.*sin(tt);

%energy
hump=1;
alpha1=c/Cs; alpha2=(2*c-Cs)/Cs;
f=hump*(C4*c.^4+C3*c.^3+C2.*c.^2) + E2*alpha2.*(e2.^2+e3.^2)+ E3*alpha1.*e3.*(e3.^2-3*e2.^2) + E4*(e2.^2+e3.^2).^2;
contourslice(c,e2,e3,f,[0,1],[0],[],35);
xlabel('c'); ylabel('e2'); zlabel('e3');


%add minimum energy line plots
minF1=zeros(size(c,1),3);
minF2=zeros(size(c,1),3);
minF3=zeros(size(c,1),3);
%alpha(.7)
for i=1:size(c,1)
    %path 1
    tmp1=squeeze(f(:,i,:)+e3(:,i,:)*0.01);
    [val,ind]=min(tmp1(:));
    [row,col]=ind2sub(size(tmp1),ind);
    minF1(i,1)=c(row,i,col);
    minF1(i,2)=e2(row,i,col);
    minF1(i,3)=e3(row,i,col);
    %path 2
    tmp1=squeeze(f(:,i,:)-e2(:,i,:)*0.01);
    [val,ind]=min(tmp1(:));
    [row,col]=ind2sub(size(tmp1),ind);
    minF2(i,1)=c(row,i,col);
    minF2(i,2)=e2(row,i,col);
    minF2(i,3)=e3(row,i,col);
    %path 3
    tmp1=squeeze(f(:,i,:)+e2(:,i,:)*0.01);
    [val,ind]=min(tmp1(:));
    [row,col]=ind2sub(size(tmp1),ind);
    minF3(i,1)=c(row,i,col);
    minF3(i,2)=e2(row,i,col);
    minF3(i,3)=e3(row,i,col);
end
hold on;
plot3(minF1(:,1),minF1(:,2),minF1(:,3),'LineWidth',3);
plot3(minF2(:,1),minF2(:,2),minF2(:,3),'LineWidth',3);
plot3(minF3(:,1),minF3(:,2),minF3(:,3),'LineWidth',3);
