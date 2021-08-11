%% ####################################################################################################################
% Code for the paper:
% Utilizing Matrix Completion for Simulation and Optimization of Water Distribution Networks
% By Mashor Housh, Alaa Jamal
% University of Haifa, mhoush@univ.haifa.ac.il, alaa.gazi.50@gmial.com
%% ####################################################################################################################
% Developed under Matlab 2018b
%% ####################################################################################################################
function [f]=HQ(x,A,R,H0,demand,nnodes)
H=x(1:nnodes);
Q=x(nnodes+1:end);
f1=A*H+diag(R)*diag(Q)*abs(Q).^0.852-H0;
f2=A'*Q-demand;
f=[f1;f2];
end
