Sc_Plant_Without_Linearization % 得到线性模型
%%%%%%%%%%%%%%%%
%% LQR Ricatti
Q1 = C'*2*C;
Q1(1,1) = 10;%xP
Q1(3,3) = 10;%yP
Q1(5,5) = 10;%zP
Q1(7,7) = 1;%phi
Q1(9,9) = 1;%theta
Q1(11,11) = 1;%psi

R1 = eye(4);
[P,L,G]=care(A,B,Q1,R1);
K_LQR = -R1\B'*P; % 基于线性模型和性能权重，得到LQR反馈增益
% K_yc = (R1\B')*((P*B/R1*B'-A')\C'*Q1);

