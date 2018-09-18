t = 40; % total time
t_stop = 35; % back to (0,0,0)
Fs = 1000;
Time = (0:Fs*t-1)'/Fs;
%% z:
des_z = 0.0*sin(Time);
des_dz = LG_ddiff(des_z,Fs);
%% step test:
% posx = -1;
% posy = 1;
% posz = 1;
% x1 = LG_step(0,0,2,Fs);
% x12 = LG_step(0,posx,2,Fs);
% x2 = LG_step(posx,posx,8,Fs);
% 
% y1 = LG_step(0,0,2,Fs);
% y12 = LG_step(0,posy,2,Fs);
% y2 = LG_step(posy,posy,8,Fs);
% 
% z1 = LG_step(0,0,2,Fs);
% z12 = LG_step(0,posz,3,Fs);
% z2 = LG_step(posz,posz,7,Fs);
% 
% des_x = [x1,x12,x2];
% des_dx = LG_ddiff(des_x,Fs);
% des_y = [y1,y12,y2];
% des_dy = LG_ddiff(des_y,Fs);
% des_z = [z1,z12,z2];
% des_dz = LG_ddiff(des_z,Fs);
%% circel path
% f = 0.05;
% des_x = sqrt(2)/2*sin(2*pi*f*Time);
% des_dx = LG_ddiff(des_x,Fs);
% des_y = sqrt(2)/2*cos(2*pi*f*Time);
% des_dy = LG_ddiff(des_y,Fs);
%% "8" path
f1 = 0.05;
f2 = 0.1;
des_x = sin(2*pi*f1*Time);
des_dx = LG_ddiff(des_x,Fs);
des_y = sin(2*pi*f2*Time);
des_dy = LG_ddiff(des_y,Fs);
des_z = 0*sin(2*pi*f1*Time);
des_dz = LG_ddiff(des_z,Fs);
%% Spiral path
% f = 0.1;
% des_x = 0.05*Time.*sin(2*pi*f*Time);
% des_dx = LG_ddiff(des_x,Fs);
% des_y = 0.05*Time.*cos(2*pi*f*Time);
% des_dy = LG_ddiff(des_y,Fs);
%% Square path
% x1 = LG_step(1.5,1.5,5,Fs);
% x2 = LG_step(1.5,1.5,5,Fs);
% x3 = LG_step(-1.5,-1.5,5,Fs);
% x4 = LG_step(-1.5,-1.5,5,Fs);
% x5 = LG_step(0,0,10,Fs);
% 
% y1 = LG_step(0,0,5,Fs);
% y2 = LG_step(1.5,1.5,5,Fs);
% y3 = LG_step(1.5,1.5,5,Fs);
% y4 = LG_step(0,0,5,Fs);
% y5 = LG_step(0,0,10,Fs);
% 
% des_x = [x1,x2,x3,x4,x5];
% des_dx = LG_ddiff(des_x,Fs);
% des_y = [y1,y2,y3,y4,y5];
% des_dy = LG_ddiff(des_y,Fs);

%% Îå½ÇÐÇ£º
% r = 1.5;
% Q1 = deg2rad(72);
% Q2 = deg2rad(36);
% Q3 = deg2rad(36);
% Q4 = deg2rad(72);
% pos1 = [0,r];
% pos2 = [-sin(Q2)*r,-cos(Q2)*r];
% pos3 = [sin(Q4)*r,cos(Q4)*r];
% pos4 = [-sin(Q1)*r,cos(Q1)*r];
% pos5 = [sin(Q3)*r,-cos(Q3)*r];
% 
% x1 = LG_step(0,pos1(1),5,Fs);
% x2 = LG_step(pos1(1),pos2(1),5,Fs);
% x3 = LG_step(pos2(1),pos3(1),5,Fs);
% x4 = LG_step(pos3(1),pos4(1),5,Fs);
% x5 = LG_step(pos4(1),pos5(1),5,Fs);
% x6 = LG_step(pos5(1),pos1(1),5,Fs);
% 
% y1 = LG_step(0,pos1(2),5,Fs);
% y2 = LG_step(pos1(2),pos2(2),5,Fs);
% y3 = LG_step(pos2(2),pos3(2),5,Fs);
% y4 = LG_step(pos3(2),pos4(2),5,Fs);
% y5 = LG_step(pos4(2),pos5(2),5,Fs);
% y6 = LG_step(pos5(2),pos1(2),5,Fs);
% 
% des_x = [x1,x2,x3,x4,x5,x6];
% des_dx = LG_ddiff(des_x,Fs);
% des_y = [y1,y2,y3,y4,y5,y6];
% des_dy = LG_ddiff(des_y,Fs);
%%
% x1 = LG_step(0,pos1(1),5,Fs);
% x2 = LG_step(pos1(1),pos2(1),5,Fs);
% x3 = LG_step(pos2(1),pos3(1),5,Fs);
% x4 = LG_step(pos3(1),pos4(1),5,Fs);
% x5 = LG_step(pos4(1),pos5(1),5,Fs);
% x6 = LG_step(pos5(1),pos1(1),5,Fs);
% 
% y1 = LG_step(0,pos1(2),5,Fs);
% y2 = LG_step(pos1(2),pos2(2),5,Fs);
% y3 = LG_step(pos2(2),pos3(2),5,Fs);
% y4 = LG_step(pos3(2),pos4(2),5,Fs);
% y5 = LG_step(pos4(2),pos5(2),5,Fs);
% y6 = LG_step(pos5(2),pos1(2),5,Fs);

%%
% des_yaw1 = LG_step(0,-1/4*pi,10,Fs);
% des_yaw2 = LG_step(-1/4*pi,-1/4*pi,10,Fs);
% des_yaw3 = LG_step(-1/4*pi,0,10,Fs);
% des_yaw = [des_yaw1,des_yaw2,des_yaw3];
% des_dyaw = LG_ddiff(des_yaw,Fs);
% %% 
% %%%% with respect to Q1
% des_a1_1 = LG_step(deg2rad(90),deg2rad(45),10,Fs);
% des_a1_2 = LG_step(deg2rad(45),deg2rad(90),10,Fs);
% des_a1_3 = LG_step(deg2rad(90),deg2rad(90),10,Fs);
% des_a1_4 = LG_step(deg2rad(90),deg2rad(90),20,Fs);
% 
% des_b1_1 = LG_step(deg2rad(45),deg2rad(0),10,Fs);
% des_b1_2 = LG_step(deg2rad(0),deg2rad(-45),10,Fs);
% des_b1_3 = LG_step(deg2rad(-45),deg2rad(-45),30,Fs);
% 
% des_a1 = [des_a1_1,des_a1_2,des_a1_3,des_a1_4];
% des_da1 = LG_ddiff(des_a1,Fs);
% des_b1 = [des_b1_1,des_b1_2,des_b1_3] ;
% des_db1 = LG_ddiff(des_b1,Fs);
% %%%% with respect to Q2
% 
% % des_a2_1 = LG_step(deg2rad(90),deg2rad(5),10,Fs);
% % des_a2_2 = LG_step(deg2rad(5),deg2rad(90),40,Fs);
% % 
% % des_b2_1 = LG_step(deg2rad(-45),deg2rad(0),10,Fs);
% % des_b2_2 = LG_step(deg2rad(0),deg2rad(-45),40,Fs);
% 
% 
% 
% des_a2 =  deg2rad(180) - des_a1  ;
% des_da2 = LG_ddiff(des_a2,Fs);
% des_b2 = -des_b1 ;
% des_db2 = LG_ddiff(des_b2,Fs);
%%
cmd_x = [Time,des_x];
cmd_dx = [Time,des_dx];
cmd_y = [Time,des_y];
cmd_dy = [Time,des_dy];
cmd_z = [Time,des_z];
cmd_dz = [Time,des_dz];
%% attitude reference
des_phi = 0*Time;
des_dphi =  LG_ddiff(des_phi,Fs);
des_theta =  0*Time;
des_dtheta =  LG_ddiff(des_theta,Fs);
des_psi =  0*Time;
des_dpsi =  LG_ddiff(des_psi,Fs);

cmd_phi = [Time,(des_phi)];
cmd_dphi = [Time,(des_dphi)];
cmd_theta = [Time,(des_theta)];
cmd_dtheta = [Time,(des_theta)];
cmd_psi = [Time,(des_psi)];
cmd_dpsi = [Time,(des_dpsi)];