close all; clc; clear; 
screen_size = get(0, 'ScreenSize');
f1 = figure(1);
set(f1, 'Position', [0 0 screen_size(3) screen_size(4) ] );
lVal=1.0;
%Chemical free energy parameters
Cd=1.0; Cg=100*lVal^2;
C4=16*Cd; C3=-32*Cd; C2=16*Cd;
%Mechanical free energy parameters
Ed=1.0e-3; Es=0.01;
E4=1.5*Ed/Es^4;
E3=-Ed/Es^3;
E2=-1.5*Ed/Es^2;
%
%generate surface
intervals=40;
[c,e2,e3]=meshgrid(linspace(0,1,intervals),linspace(-Es,Es,intervals),linspace(-Es,Es,intervals));

%free energy slices
E2c=E2*(5*c-2.0)/3.0;
E3c=E3*c;
E4c=E4*c;
f=(C4.*c.^4 + C3.*c.^3 + C2.*c.^2) + E4c.*(e2.^2+e3.^2).^2+E3c.*e3.*(e3.^2-3*e2.^2)+E2c.*(e2.^2+e3.^2);
%f=hump*((Ed/Es^4).*c.^4+(-2*Ed/Es^2).*c.^2) + E2*alpha.*(xx.^2+yy.^2)+ E3*alpha.*yy.*(yy.^2-3*xx.^2) + E4*(xx.^2+yy.^2).^2;
%contourslice(c,e2,e3,f,[0,0.5,1],[],[],35);
xlabel('c','FontSize',24,'FontWeight','bold'); 
ylabel('e_2','FontSize',24,'FontWeight','bold'); 
zlabel('e_3','FontSize',24,'FontWeight','bold');
%axis tight;
view([-38 30])
set(gca, 'FontSize', 18)

%import data
load values3DE23.mat

%video 
writerObj = VideoWriter('manifold.avi');
writerObj.Quality=80;
writerObj.FrameRate=2;
open(writerObj);
%plot scatter points on the slices
maxStep=numIncs*100;
skipPoints=500;
for t=[100:1*200:maxStep maxStep]
    cla; 
    eval(['x=T', num2str(t), ';']);
    hold on; 
    h = contourslice(c,e2,e3,f,[0,1],[],[],25);
    set(h, 'LineWidth',2.0);
    c_2=x(1:skipPoints:end,1); e2_2=x(1:skipPoints:end,2); e3_2=x(1:skipPoints:end,3);
    c_2(c_2>1.0)=1.0;
    e2_2(e2_2>Es)=Es; e2_2(e2_2<-Es)=-Es;
    e3_2(e3_2>Es)=Es; e3_2(e3_2<-Es)=-Es;
    h = scatter3(c_2,e2_2,e3_2,'MarkerEdgeColor','k'); %,'fill','k');
    suptitle(['Time: ' num2str(double(t)/double(maxStep),2)]); 
    hold off;
    axis([0 1.0 -Es Es -Es Es])
    frame = getframe(gcf);   
    writeVideo(writerObj,frame);
    pause(.1); 
end
close(writerObj);

% samplePts=15; 
% fac=1.2; %fac>1
% [c, e]=meshgrid(-fac*ac:ac/samplePts:fac*ac,-fac*ae:ae/samplePts:fac*ae);
% minusC=(ac+3*c)./(4*ac);
% plusC=(ac+c)./(2*ac);
% F=Ac*c.^4+Bc*c.^2 + plusC.*Ae.*e.^4 + minusC*Be.*e.^2;
% %F=c.^2 + e.^2
% %surf(c,e,F); zlabel ('F');
% %contour(c,e,F);
% %surf(c,e,F,F,'facecolor','interp'); 
% %set(gca,'clim',[-0.001,0.001])
% [px,py] = gradient(F);
% v=-fac*ac:ac/samplePts:fac*ac;
% w=-fac*ae:ae/samplePts:fac*ae;
% %
% writerObj = VideoWriter('manifold.avi');
% writerObj.Quality=80;
% writerObj.FrameRate=3;
% open(writerObj);
% set(0,'DefaultAxesFontSize', 14)
% skip=5;
% maxStep=406;
% for t=1:4:maxStep
%     eval(['x=T', num2str(t), ';']);
%     subplot(121);
%     c2=x(1:skip:end,1); e2=x(1:skip:end,3);
%     minusC2=(ac+3*c2)./(4*ac);
%     plusC2=(ac+c2)./(2*ac);
%     F2=Ac*c2.^4+Bc*c2.^2 + plusC2.*Ae.*e2.^4 + minusC2*Be.*e2.^2;
%     scatter(c2,e2,'k'); hold on;
%     contour(v,w,F,20); hold off;
%     suptitle(['Time: ' num2str(t/maxStep)]); 
%     xlabel('c', 'FontSize', 20); ylabel('e_2', 'FontSize', 20);
%     xlim([-1*fac 1*fac]); ylim([-0.01*1.2 1.2*0.01]);
%     subplot(122);
%     surf(c,e,F,'facecolor','interp','LineStyle','none'); hold on;
%     scatter3(c2,e2,F2,'fill','k'); hold off;
%     xlim([-1*fac 1*fac]); ylim([-0.01*1.2 1.2*0.01]); zlim([-(de+dc) 1])
%     xlabel('c', 'FontSize', 20); ylabel('e_2', 'FontSize', 20); zlabel('F(c,e_2)', 'FontSize', 20);
%     view([127.5 30])
%     frame = getframe(gcf);   
%     writeVideo(writerObj,frame);
%     pause(0.01)
% end
%close(writerObj);