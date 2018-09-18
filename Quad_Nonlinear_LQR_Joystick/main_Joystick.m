%% Matlab 2015a
clear;clc;close
Sc_H2LQR_R % 得到LQR 控制器
sim('Plant_N_LQR_Joystick.slx'); % 开始simulink 仿真：35s后返回原点
%%
Sc_Plant_Animation % 仿真结束后，动画