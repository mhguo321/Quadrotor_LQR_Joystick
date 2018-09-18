function [sys,x0,str,ts] = LG_quadcopter_Sfun(t,x,u,flag,x0,y0,z0) 
%%
 switch flag,
    case 0, 
    [sys,x0,str,ts]=mdlInitializeSizes(x0,y0,z0); 
    case 1, 
    sys=mdlDerivatives(t,x,u); 
    case 3, 
    sys=mdlOutputs(t,x,u); 
    case { 2, 4, 9 }, 
    sys = []; 
    otherwise 
    error(['Unhandled flag = ',num2str(flag)]); 
 end
function [sys,x0,str,ts]=mdlInitializeSizes(x0,y0,z0)
%% 
sizes = simsizes; 
sizes.NumContStates = 12; 
sizes.NumDiscStates = 0; 
sizes.NumOutputs = 12; 
sizes.NumInputs = 4; 
sizes.DirFeedthrough = 0; 
sizes.NumSampleTimes = 1; 
sys = simsizes(sizes); 
x0 = [x0,0,y0,0,z0,0,0,0,0,0,0,0];
% x0= [xP,dxP,yP,dyP,zP,dzP,...
%      phi1,dphi1,theta1,dtheta1,psi1,dpsi1]
%Air Force Parameters:
str = []; 
ts = [0 0]; 
function sys=mdlDerivatives(t,x,u) 
%% Model states and controllors:
xP = x(1);dxP = x(2);
yP = x(3);dyP = x(4);
zP = x(5);dzP = x(6);

phi = x(7);dphi = x(8);
theta = x(9);dtheta = x(10);
psi = x(11);dpsi = x(12);

Ft=u(1);Mx=u(2);My=u(3);Mz=u(4);
% Parameters:
% Quadcopter Parameters:
mQ=0.55;g=9.8;
LQ = 0.17;
kF = 2.98e-6;
kM = 1.14e-7;

Ixx=0.0023;
Iyy=0.0028;
Izz=0.0046;
inertia=[Ixx,0,0;0,Iyy,0;0,0,Izz];
% Payload parameters:
% mP = 0.05;
% Lr = 0.5;
% Body Frame Forces

%% TransForm Matrix

%% Generalized Force 
M = [mQ,0,0,0,0,0;0,mQ,0,0,0,0;0,0,mQ,0,0,0;0,0,0,Ixx,0,(-1).*Ixx.* ...
  sin(theta);0,0,0,0,Iyy.*cos(phi).^2+Izz.*sin(phi).^2,(Iyy+(-1).* ...
  Izz).*cos(phi).*cos(theta).*sin(phi);0,0,0,(-1).*Ixx.*sin(theta),( ...
  Iyy+(-1).*Izz).*cos(phi).*cos(theta).*sin(phi),Izz.*cos(phi).^2.* ...
  cos(theta).^2+Iyy.*cos(theta).^2.*sin(phi).^2+Ixx.*sin(theta).^2]; 

fdq1 = Ft.*(sin(phi).*sin(psi)+cos(phi).*cos(psi).*sin(theta));

fdq2 = Ft.*((-1).*cos(psi).*sin(phi)+cos(phi).*sin(psi).*sin(theta));

fdq3 = (-1).*g.*mQ+Ft.*cos(phi).*cos(theta);

fdq4 = Mx+dpsi.*dtheta.*(Ixx+(Iyy+(-1).*Izz).*cos(2.*phi)).*cos(theta)+( ...
  -1).*(Iyy+(-1).*Izz).*cos(phi).*(dtheta.^2+(-1).*dpsi.^2.*cos( ...
  theta).^2).*sin(phi);

fdq5 = cos(phi).*(My+2.*dphi.*dtheta.*(Iyy+(-1).*Izz).*sin(phi))+(1/2).*( ...
  (-2).*Mz.*sin(phi)+dpsi.*cos(theta).*((-1).*dpsi.*(Iyy+Izz).*sin( ...
  theta)+(-2).*Ixx.*(dphi+(-1).*dpsi.*sin(theta))+(Iyy+(-1).*Izz).* ...
  cos(2.*phi).*((-2).*dphi+dpsi.*sin(theta))));

fdq6 = dphi.*dpsi.*((-1).*Iyy+Izz).*cos(theta).^2.*sin(2.*phi)+(-1).*(Mx+ ...
  dtheta.^2.*((-1).*Iyy+Izz).*cos(phi).*sin(phi)).*sin(theta)+cos( ...
  theta).*(dphi.*dtheta.*Ixx+Mz.*cos(phi)+My.*sin(phi)+dpsi.* ...
  dtheta.*((-2).*Ixx+Iyy+Izz).*sin(theta)+(-1).*dtheta.*(Iyy+(-1).* ...
  Izz).*cos(2.*phi).*(dphi+dpsi.*sin(theta)));


fdq = [fdq1;fdq2;fdq3;fdq4;fdq5;fdq6];
dummy = M\fdq;

xPdot = dxP;
dxPdot = dummy(1);
yPdot = dyP;
dyPdot = dummy(2);
zdot = dzP;
dzdot = dummy(3);

phidot = dphi;
dphidot = dummy(4);
thetadot = dtheta;
dthetadot = dummy(5);
psidot = dpsi;
dpsidot = dummy(6);

sys = [xPdot;dxPdot;yPdot;dyPdot;zdot;dzdot;...
    phidot;dphidot;thetadot;dthetadot;psidot;dpsidot]; 
function sys=mdlOutputs(t,x,u)
%% Model states and controllors:
xP = x(1);dxP = x(2);
yP = x(3);dyP = x(4);
zP = x(5);dzP = x(6);

phi = x(7);dphi = x(8);
theta = x(9);dtheta = x(10);
psi = x(11);dpsi = x(12);

Ft=u(1);Mx=u(2);My=u(3);Mz=u(4);

%% TransForm Matrix

sys = [xP;dxP; yP;dyP; zP;dzP; phi;dphi; theta;dtheta; psi;dpsi];

function sys=mdlGetTimeOfNextVarHit(t,x,u)
sampleTime = 1;    %  Example, set the next hit to be one second later.
sys = t + sampleTime;

function sys=mdlTerminate(t,x,u)
sys = [];