%% Matlab 2015a
clear;clc;close
Sc_H2LQR_R % 得到LQR 控制器
Sc_reference % 生成"8" 轨迹
sim('Plant_N_LQR.slx'); % 开始simulink 仿真：35s后返回原点
%%
Sc_Plant_Animation % 仿真结束后，动画