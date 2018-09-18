function [phandles] = LG_plot_initial(fig,state)
%LG_PLOT_INITIAL Summary of this function goes here
%   Detailed explanation goes here
%% 导入数据
% close all;
xQ  = state(1);
yQ  = state(3);
zQ  = state(5);

phi    = state(7);
theta  = state(9);
psi    = state(11);

%% 数据预处理
R = euler2rotMat(phi, theta, psi);

ux = R(1,1);% each column vector of R is a direction vector
vx = R(2,1);
wx = R(3,1);
uy = R(1,2);
vy = R(2,2);
wy = R(3,2);
uz = R(1,3);
vz = R(2,3);
wz = R(3,3); 
%% 初始化figure属性
position = [50 50 866 600];
Xlabel = 'X(m)';
Ylabel = 'Y(m)';
Zlabel = 'Z(m)';
set(fig, 'Position', position);
set(gca, 'drawmode', 'fast');
lighting phong;
set(gcf, 'Renderer', 'zbuffer');
hold on;
axis equal;
grid on;
az = 120;
el = 20;
view(az, el);
xlabel(Xlabel,'Interpreter','LaTex');
ylabel(Ylabel,'Interpreter','LaTex');
zlabel(Zlabel,'Interpreter','LaTex');
%% 初始化 机体坐标系句柄：
ratio = 0.4;
quivXhandle = quiver3(xQ, yQ, zQ, ratio*ux, ratio*vx, ratio*wx,  'r','MaxHeadSize', 0.999999, 'AutoScale', 'off');
quivYhandle = quiver3(xQ, yQ, zQ, ratio*uy, ratio*vy, ratio*wy,  'g','MaxHeadSize', 0.999999, 'AutoScale', 'off');
quivZhandle = quiver3(xQ, yQ, zQ, ratio*uz, ratio*vz, ratio*wz,  'b','MaxHeadSize', 0.999999, 'AutoScale', 'off');
%% 初始化 四旋翼句柄：
    len = 0.34; % distance from the center of rotor 1 to rotor 3; 
    x1 = [xQ,yQ,zQ] + 0.5*len*[ux,vx,wx];
    x3 = [xQ,yQ,zQ] - 0.5*len*[ux,vx,wx];
    x2 = [xQ,yQ,zQ] + 0.5*len*[uy,vy,wy];
    x4 = [xQ,yQ,zQ] - 0.5*len*[uy,vy,wy];
    QX1 = [x1(1),x3(1)];
    QX2 = [x1(2),x3(2)];
    QX3 = [x1(3),x3(3)];
    QY1 = [x2(1),x4(1)];
    QY2 = [x2(2),x4(2)];
    QY3 = [x2(3),x4(3)];
    
    
    eb1 = [ux,vx,wx];
    eb2 = [uy,vy,wy];
    eb3 = [uz,vz,wz];
    temp = 0:(2*pi)/15 :2*pi;
    r0 = 0.04;
    xr1 = r0*cos(temp')*eb1 + r0*sin(temp')*eb2 + ones(size(temp'))*x1;
    xr2 = r0*cos(temp')*eb1 + r0*sin(temp')*eb2 + ones(size(temp'))*x2;
    xr3 = r0*cos(temp')*eb1 + r0*sin(temp')*eb2 + ones(size(temp'))*x3;
    xr4 = r0*cos(temp')*eb1 + r0*sin(temp')*eb2 + ones(size(temp'))*x4;

    lineQXhandle = line('XData',QX1,'YData',QX2,'ZData',QX3,'LineWidth',2,'Color',[1 0 0]);
    lineQYhandle = line('XData',QY1,'YData',QY2,'ZData',QY3,'LineWidth',2,'Color',[1 0 0]);
    circle1handle = line('XData',xr1(:,1),...
        'YData',xr1(:,2),...
        'ZData',xr1(:,3),...
        'LineWidth',2,'Color',[1 0 0]);
    
    circle2handle = line('XData',xr2(:,1),...
    'YData',xr2(:,2),...
    'ZData',xr2(:,3),...
    'LineWidth',2,'Color',[1 0 0]);

    circle3handle = line('XData',xr3(:,1),...
    'YData',xr3(:,2),...
    'ZData',xr3(:,3),...
    'LineWidth',2,'Color',[1 0 0]);

    circle4handle = line('XData',xr4(:,1),...
    'YData',xr4(:,2),...
    'ZData',xr4(:,3),...
    'LineWidth',2,'Color',[1 0 0]);
%% 初始化 payload 句柄
% [sx,sy,sz] = sphere;
% SPx = 0.02*sx + xP;
% SPy = 0.02*sy + yP;
% SPz = 0.02*sz + zP;
% payloadhandle = surf(SPx,SPy,SPz,'FaceColor',[1 0 0],'EdgeColor','none');
%% 初始化 cable 句柄
% cablehandle = line('Xdata',[xQ,xP],...
%         'Ydata',[yQ,yP],...
%         'Zdata',[zQ,zP],...
%         'LineWidth',2,'Color','k');
%% 设置初始范围
LimitRatio = 1;AxisLength = 0.6;
Xlim = [xQ-AxisLength xQ+AxisLength] * LimitRatio;
Ylim = [yQ-AxisLength yQ+AxisLength] * LimitRatio;
Zlim = [zQ-AxisLength zQ+AxisLength] * LimitRatio;
set(gca, 'Xlim', Xlim, 'Ylim', Ylim, 'Zlim', Zlim);
%% 返回句柄结构体
phandles.quivXhandle = quivXhandle;
phandles.quivYhandle = quivYhandle;
phandles.quivZhandle = quivZhandle;
phandles.lineQXhandle = lineQXhandle;
phandles.lineQYhandle = lineQYhandle;
phandles.circle1handle = circle1handle;
phandles.circle2handle = circle2handle;
phandles.circle3handle = circle3handle;
phandles.circle4handle = circle4handle;
% phandles.payloadhandle = payloadhandle;
% phandles.cablehandle = cablehandle;
phandles.Xlim = Xlim;
phandles.Ylim = Ylim;
phandles.Zlim = Zlim;
end

function R = euler2rotMat(phi, theta, psi)
    R(1,1,:) = cos(psi).*cos(theta);
    R(1,2,:) = -sin(psi).*cos(phi) + cos(psi).*sin(theta).*sin(phi);
    R(1,3,:) = sin(psi).*sin(phi) + cos(psi).*sin(theta).*cos(phi);
    
    R(2,1,:) = sin(psi).*cos(theta);
    R(2,2,:) = cos(psi).*cos(phi) + sin(psi).*sin(theta).*sin(phi);
    R(2,3,:) = -cos(psi).*sin(phi) + sin(psi).*sin(theta).*cos(phi);
    
    R(3,1,:) = -sin(theta);
    R(3,2,:) = cos(theta).*sin(phi);
    R(3,3,:) = cos(theta).*cos(phi);
end