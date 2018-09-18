clc;close all;

%%
posQ = [xQ_out,yQ_out,zQ_out];
% Euler = [phi_z,theta_z,psi_z];
R = LG_euler2rotMat(phi_out,theta_out,psi_out);

SamplePlotFreq = 50;
isCreateAVI = false;
isFixView = false;
type1 = 'DotsOnly';
type2 = 'All';
Fig = LG_6DoFAnimation_Tra(posQ,R,posQ,SamplePlotFreq,type1,isCreateAVI,isFixView);
%%
% plot3(des_x,des_y,des_z,'r');