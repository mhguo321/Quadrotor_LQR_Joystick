% x,dx,ax,jx
% y,dy,ay,jy
% z,dz,az,jz
% psi,dpsi,apsi
clear;clc;
t_stop = 50;
%
Time = 20;
Fs = 1000;
Time = (0:Fs*Time-1)'/Fs;

tf = 10;
x0 = 0; v0 = 0; a0 = 0;
xf = 1; vf = 0; af = 0;    
[ c,x1,dx1,ax1 ] = LG_path_opt(x0,v0,a0,xf,vf,af,tf);

x0 = x1(end); v0 = dx1(end); a0 = ax1(end);
xf = 0; vf = 0; af = 0;    
[ c,x2,dx2,ax2 ] = LG_path_opt(x0,v0,a0,xf,vf,af,tf);
x = [x1;x2];
dx = [dx1;dx2];
ax = [ax1;ax2];

y0 = 0; v0 = 0; a0 = 0;
yf = 0; vf = 0.5; af = 0;
[ c,y1,dy1,ay1 ] = LG_path_opt(y0,v0,a0,yf,vf,af,tf);

y0 = y1(end); v0 = dy1(end); a0 = ay1(end);
xf = 0; vf = 0; af = 0;    
[ c,y2,dy2,ay2 ] = LG_path_opt(y0,v0,a0,yf,vf,af,tf);
y = [y1;y2];
dy = [dy1;dy2];
ay = [ay1;ay2];

X0 = [0;0;0;0;zeros(8,1)];
load('K.mat')
%%
% x = des_x;
% dx = LG_ddiff(x,Fs);
% ax = LG_ddiff(dx,Fs);
jx = LG_ddiff(ax,Fs);

% y = y;
% dy = LG_ddiff(y,Fs);
% ay = LG_ddiff(dy,Fs);
jy = LG_ddiff(ay,Fs);

z = 0*x;
dz = LG_ddiff(z,Fs);
az = LG_ddiff(dz,Fs);
jz = LG_ddiff(az,Fs);

psi = 0*x;
dpsi = 0*x;
ddpsi = 0*x;
%%
xw = [1,0,0]';
yw = [0,1,0]';
zw = [0,0,1]';
%%
g = 9.8;mQ = 0.55;
N = length(x);
zb = zeros(3,N);
yb = zeros(3,N);
xb = zeros(3,N);
xc = zeros(3,N);
adot = zeros(3,N);
R = zeros(3,3,N);
euler = zeros(3,N);

u1 = zeros(N,1);
hw = zeros(3,N);
phi = zeros(N,1);
theta = zeros(N,1);
psi = zeros(N,1);
p = zeros(N,1);
q = zeros(N,1);
r = zeros(N,1);
dphi = zeros(N,1);
dtheta = zeros(N,1);
dpsi = zeros(N,1);
Omega = zeros(3,N);
Wyit = zeros(3,3,N);

for i = 1:N;
    %% Euler angles:
    t = [ax(i),ay(i),az(i)+g]';
    zb(:,i) = t/norm(t);
    xc(:,i) = [cos(psi(i)),sin(psi(i)),0]';
    yb(:,i) = cross(zb(:,i),xc(:,i))/norm(cross(zb(:,i),xc(:,i)));
    xb(:,i) = cross(yb(:,i),zb(:,i));
    R(:,:,i) = [xb(:,i),yb(:,i),zb(:,i)]; 
    euler(:,i) = rotMat2euler(R(:,:,i));
    phi(i) = euler(1,i);
    theta(i) = euler(2,i);
    psi(i) = euler(3,i);
    Wyit(:,:,i) = [1,0,0;...
                    0,cos(phi(i)),-sin(phi(i));...
                    -sin(theta(i)),cos(theta(i))*sin(phi(i)),cos(theta(i))*cos(phi(i))];
    %% angular velocity
    u1(i) = mQ*norm(t);
    adot(:,i) = [jx(i),jy(i),jz(i)]';
    hw(:,i) = (mQ/u1(i)).*(adot(:,i)-dot(zb(:,i),adot(:,i))*zb(:,i));
    p(i) = -dot(hw(:,i),yb(:,i));
    q(i) = dot(hw(:,i),xb(:,i));
    r(i) = dpsi(i)*dot(zw,zb(:,i));
    Omega(:,i) = [p(i),q(i),r(i)]';
    temp = Wyit(:,:,i)*Omega(:,i);
    dphi(i) = temp(1);
    dtheta(i) = temp(2);
    dpsi(i) = temp(3);
end
%%
cmd_x = [Time,x];
cmd_dx = [Time,dx];
cmd_y = [Time,y];
cmd_dy = [Time,dy];
cmd_z = [Time,z];
cmd_dz = [Time,dz];

cmd_phi = [Time,phi];
cmd_dphi = [Time,dphi];
cmd_theta = [Time,theta];
cmd_dtheta = [Time,dtheta];
cmd_psi = [Time,psi];
cmd_dpsi = [Time,dpsi];