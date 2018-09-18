function LG_plot_update(fig,state,phandles)
%LG_PLOT_INITIAL Summary of this function goes here
%   Detailed explanation goes here
%% 导入数据
% Lr = 1;
xQ1  = state(1);
yQ1  = state(3);
zQ1  = state(5);

phi1    = state(7);
theta1  = state(9);
psi1    = state(11);

% isTaut = state(23);
% alpha1 = state(13);
% beta1 = state(15);

% alpha2 = state(17);
% beta2 = state(19);

% phi2    = state(21);
% theta2  = state(23);
% psi2    = 0;%state(25);

% xP  = state(17);
% yP  = state(19);
% zP  = state(21);

% xQ2 = xP + Lr*cos(beta2).*cos(alpha2);
% yQ2 = yP + Lr*cos(beta2).*sin(alpha2);
% zQ2 = zP + Lr*sin(beta2);
%% 数据预处理
R1 = euler2rotMat(phi1, theta1, psi1);
% R2 = euler2rotMat(phi2, theta2, psi2);

ux1 = R1(1,1);% each column vector of R is a direction vector
vx1 = R1(2,1);
wx1 = R1(3,1);
uy1 = R1(1,2);
vy1 = R1(2,2);
wy1 = R1(3,2);
uz1 = R1(1,3);
vz1 = R1(2,3);
wz1 = R1(3,3); 

% ux2 = R2(1,1);% each column vector of R is a direction vector
% vx2 = R2(2,1);
% wx2 = R2(3,1);
% uy2 = R2(1,2);
% vy2 = R2(2,2);
% wy2 = R2(3,2);
% uz2 = R2(1,3);
% vz2 = R2(2,3);
% wz2 = R2(3,3); 
%% 更新机体坐标系
ratio = 0.4;
set(phandles.quivXhandle, 'xdata', xQ1, 'ydata', yQ1, 'zdata', zQ1,...
    'udata', ratio*ux1, 'vdata', ratio*vx1, 'wdata', ratio*wx1);
set(phandles.quivYhandle, 'xdata', xQ1, 'ydata', yQ1, 'zdata', zQ1,...
    'udata', ratio*uy1, 'vdata', ratio*vy1, 'wdata', ratio*wy1);
set(phandles.quivZhandle, 'xdata', xQ1, 'ydata', yQ1, 'zdata', zQ1,...
    'udata', ratio*uz1, 'vdata', ratio*vz1, 'wdata', ratio*wz1);
%% 更新第一架四旋翼
    len = 0.34; % distance from the center of rotor 1 to rotor 3; 
    x1 = [xQ1,yQ1,zQ1] + 0.5*len*[ux1,vx1,wx1];
    x3 = [xQ1,yQ1,zQ1] - 0.5*len*[ux1,vx1,wx1];
    x2 = [xQ1,yQ1,zQ1] + 0.5*len*[uy1,vy1,wy1];
    x4 = [xQ1,yQ1,zQ1] - 0.5*len*[uy1,vy1,wy1];
    Q1X1 = [x1(1),x3(1)];
    Q1X2 = [x1(2),x3(2)];
    Q1X3 = [x1(3),x3(3)];
    Q1Y1 = [x2(1),x4(1)];
    Q1Y2 = [x2(2),x4(2)];
    Q1Y3 = [x2(3),x4(3)];
    
    eb1 = [ux1,vx1,wx1];
    eb2 = [uy1,vy1,wy1];
    eb3 = [uz1,vz1,wz1];
    temp = 0:(2*pi)/15 :2*pi;
    r0 = 0.04;
    x1r1 = r0*cos(temp')*eb1 + r0*sin(temp')*eb2 + ones(size(temp'))*x1;
    x1r2 = r0*cos(temp')*eb1 + r0*sin(temp')*eb2 + ones(size(temp'))*x2;
    x1r3 = r0*cos(temp')*eb1 + r0*sin(temp')*eb2 + ones(size(temp'))*x3;
    x1r4 = r0*cos(temp')*eb1 + r0*sin(temp')*eb2 + ones(size(temp'))*x4;
set(phandles.lineQXhandle,'XData',Q1X1,'YData',Q1X2,'ZData',Q1X3);
set(phandles.lineQYhandle,'XData',Q1Y1,'YData',Q1Y2,'ZData',Q1Y3);
set(phandles.circle1handle,'XData',x1r1(:,1),'YData',x1r1(:,2),'ZData',x1r1(:,3));
set(phandles.circle2handle,'XData',x1r2(:,1),'YData',x1r2(:,2),'ZData',x1r2(:,3));
set(phandles.circle3handle,'XData',x1r3(:,1),'YData',x1r3(:,2),'ZData',x1r3(:,3));
set(phandles.circle4handle,'XData',x1r4(:,1),'YData',x1r4(:,2),'ZData',x1r4(:,3));
% %% 更新第二架四旋翼
%     len = 0.34; % distance from the center of rotor 1 to rotor 3; 
%     x1 = [xQ2,yQ2,zQ2] + 0.5*len*[ux1,vx1,wx1];
%     x3 = [xQ2,yQ2,zQ2] - 0.5*len*[ux1,vx1,wx1];
%     x2 = [xQ2,yQ2,zQ2] + 0.5*len*[uy1,vy1,wy1];
%     x4 = [xQ2,yQ2,zQ2] - 0.5*len*[uy1,vy1,wy1];
%     Q2X1 = [x1(1),x3(1)];
%     Q2X2 = [x1(2),x3(2)];
%     Q2X3 = [x1(3),x3(3)];
%     Q2Y1 = [x2(1),x4(1)];
%     Q2Y2 = [x2(2),x4(2)];
%     Q2Y3 = [x2(3),x4(3)];
%     
%     eb1 = [ux2,vx2,wx2];
%     eb2 = [uy2,vy2,wy2];
%     eb3 = [uz2,vz2,wz2];
%     temp = 0:(2*pi)/15 :2*pi;
%     r0 = 0.04;
%     x2r1 = r0*cos(temp')*eb1 + r0*sin(temp')*eb2 + ones(size(temp'))*x1;
%     x2r2 = r0*cos(temp')*eb1 + r0*sin(temp')*eb2 + ones(size(temp'))*x2;
%     x2r3 = r0*cos(temp')*eb1 + r0*sin(temp')*eb2 + ones(size(temp'))*x3;
%     x2r4 = r0*cos(temp')*eb1 + r0*sin(temp')*eb2 + ones(size(temp'))*x4;
% set(phandles.lineQ2Xhandle,'XData',Q2X1,'YData',Q2X2,'ZData',Q2X3);
% set(phandles.lineQ2Yhandle,'XData',Q2Y1,'YData',Q2Y2,'ZData',Q2Y3);
% set(phandles.Q2circle1handle,'XData',x2r1(:,1),'YData',x2r1(:,2),'ZData',x2r1(:,3));
% set(phandles.Q2circle2handle,'XData',x2r2(:,1),'YData',x2r2(:,2),'ZData',x2r2(:,3));
% set(phandles.Q2circle3handle,'XData',x2r3(:,1),'YData',x2r3(:,2),'ZData',x2r3(:,3));
% set(phandles.Q2circle4handle,'XData',x2r4(:,1),'YData',x2r4(:,2),'ZData',x2r4(:,3));
%% 更新payload
% [sx,sy,sz] = sphere;
% SPx = 0.02*sx + xP;
% SPy = 0.02*sy + yP;
% SPz = 0.02*sz + zP;
% set(phandles.payloadhandle,'XData',SPx,'YData',SPy,'ZData',SPz,'FaceColor',[1 0 0],'EdgeColor','none');
%% 更新cable
% if(isTaut==false)
%     set(phandles.cablehandle,'Color','g');
% else
%     set(phandles.cablehandle,'Color','k');
% end
% set(phandles.cablehandle, 'XData',[xQ1,xP],'YData',[yQ1,yP],'ZData',[zQ1,zP]);
% set(phandles.cable1handle, 'XData',[xQ2,xP],'YData',[yQ2,yP],'ZData',[zQ2,zP]);
%% 设置坐标范围
LimitRatio = 1.0;AxisLength = 0.6;
axisLimChanged = 1;
if((xQ1 - AxisLength) < phandles.Xlim(1)), phandles.Xlim(1) = xQ1 - LimitRatio*AxisLength; axisLimChanged = true; end
if((yQ1 - AxisLength) < phandles.Ylim(1)), phandles.Ylim(1) = yQ1 - LimitRatio*AxisLength; axisLimChanged = true; end
if((zQ1 - AxisLength) < phandles.Zlim(1)), phandles.Zlim(1) = zQ1 - LimitRatio*AxisLength; axisLimChanged = true; end
if((xQ1 + AxisLength) > phandles.Xlim(2)), phandles.Xlim(2) = xQ1 + LimitRatio*AxisLength; axisLimChanged = true; end
if((yQ1 + AxisLength) > phandles.Ylim(2)), phandles.Ylim(2) = yQ1 + LimitRatio*AxisLength; axisLimChanged = true; end
if((zQ1 + AxisLength) > phandles.Zlim(2)), phandles.Zlim(2) = zQ1 + LimitRatio*AxisLength; axisLimChanged = true; end
if(axisLimChanged), set(gca, 'Xlim', phandles.Xlim, 'Ylim', phandles.Ylim, 'Zlim', phandles.Zlim); end
drawnow;
end


function R = euler2rotMat(phi, theta, psi)
    R(1,1) = cos(psi).*cos(theta);
    R(1,2) = -sin(psi).*cos(phi) + cos(psi).*sin(theta).*sin(phi);
    R(1,3) = sin(psi).*sin(phi) + cos(psi).*sin(theta).*cos(phi);
    
    R(2,1) = sin(psi).*cos(theta);
    R(2,2) = cos(psi).*cos(phi) + sin(psi).*sin(theta).*sin(phi);
    R(2,3) = -cos(psi).*sin(phi) + sin(psi).*sin(theta).*cos(phi);
    
    R(3,1) = -sin(theta);
    R(3,2) = cos(theta).*sin(phi);
    R(3,3) = cos(theta).*cos(phi);
end
