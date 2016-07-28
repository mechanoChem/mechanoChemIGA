close all; clc; clear all;
load('values2DCE2.mat')
screen_size = get(0, 'ScreenSize');
f1 = figure(1);
set(f1, 'Position', [0 0 screen_size(3) screen_size(4) ] );
%Chemical free energy parameters
Cs=1.0; Cd=-2.0;
C4=(-16*Cd/Cs^4);
C3=(32*Cd/Cs^3);
C2=(-16*Cd/Cs^2);
%Mechanics free energy parameters
Es=0.1; Ed=-0.1;
E4=(-Ed/Es^4);
E3=0.0;
E2=(2*Ed/Es^2);
%
samplePts=15; 
fac=1.1; %fac>1
[c, e]=meshgrid((1-fac)*Cs:Cs/samplePts:fac*Cs,-fac*Es:Es/samplePts:fac*Es);
plusC=(2*c-Cs)./Cs;
F=C4*c.^4+C3*c.^3+C2*c.^2 + E4.*e.^4 + plusC*E2.*e.^2;
%F=c.^2 + e.^2
%surf(c,e,F); zlabel ('F');
%contour(c,e,F);
%surf(c,e,F,F,'facecolor','interp'); 
%set(gca,'clim',[-0.001,0.001])
[px,py] = gradient(F);
v=(1-fac)*Cs:Cs/samplePts:fac*Cs;
w=-fac*Es:Es/samplePts:fac*Es;
%
writerObj = VideoWriter('manifold.avi');
writerObj.Quality=80;
writerObj.FrameRate=3;
open(writerObj);
set(0,'DefaultAxesFontSize', 14)
skip=50;
maxStep=500;
for t=0:10:maxStep
    eval(['x=T', num2str(t), ';']);
    subplot(121);
    c2=x(1:skip:end,1); e2=x(1:skip:end,2);
    plusC2=(2*c2-Cs)./Cs;
    F2=C4*c2.^4+C3*c2.^3+C2*c2.^2 + E4.*e2.^4 + plusC2*E2.*e2.^2;
    scatter(c2,e2,'k'); hold on;
    contour(v,w,F,20); hold off;
    suptitle(['Time: ' num2str(t/maxStep)]); 
    xlabel('c', 'FontSize', 20); ylabel('e_2', 'FontSize', 20);
    xlim([(1-fac)*Cs fac*Cs]); ylim([-fac*Es fac*Es]);
    subplot(122);
    surf(c,e,F,'facecolor','interp','LineStyle','none'); hold on;
    scatter3(c2,e2,F2,'fill','k'); hold off;
    %xlim([(fac-1)*Cs fac*Cs]); ylim([-fac*Es fac*Es]); zlim([(Cd+Ed) -Cd*0.5]);
    xlabel('c', 'FontSize', 20); ylabel('e_2', 'FontSize', 20); zlabel('F(c,e_2)', 'FontSize', 20);
    view([127.5 30])
    frame = getframe(gcf);   
    writeVideo(writerObj,frame);
    pause(0.01)
end
close(writerObj);