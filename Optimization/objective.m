%% ####################################################################################################################
% Code for the paper:
% Utilizing Matrix Completion for Simulation and Optimization of Water Distribution Networks
% By Mashor Housh, Alaa Jamal
% University of Haifa, mhoush@univ.haifa.ac.il, alaa.gazi.50@gmial.com
%% ####################################################################################################################
% Developed under Matlab 2018b
%% ####################################################################################################################
function [obj]=objective(x,L,scaling_f,K,alpha,nl)
r=x(1:nl);
R=L.*r;
D=((R/scaling_f*(130^alpha)/K./L).^(-1/4.87))/2.54;
C=1.1.*D.^1.5;
obj=sum(L.*C);
end

