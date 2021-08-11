%% ####################################################################################################################
% Code for the paper:
% Utilizing Matrix Completion for Simulation and Optimization of Water Distribution Networks
% By Mashor Housh, Alaa Jamal
% University of Haifa, mhoush@univ.haifa.ac.il, alaa.gazi.50@gmial.com
%% ####################################################################################################################
% Developed under Matlab 2018b
%% ####################################################################################################################
function [w]=Completion(z,R,Hfixed,demand,mat1,mat2)
Q=mat1*[demand;z];
w=mat2*(-diag(R)*diag(Q)*abs(Q).^0.852+Hfixed);
end
