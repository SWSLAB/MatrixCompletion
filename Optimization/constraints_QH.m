%% ####################################################################################################################
% Code for the paper:
% Utilizing Matrix Completion for Simulation and Optimization of Water Distribution Networks
% By Mashor Housh, Alaa Jamal
% University of Haifa, mhoush@univ.haifa.ac.il, alaa.gazi.50@gmial.com
%% ####################################################################################################################
% Developed under Matlab 2018b
%% ####################################################################################################################
function [Cineq,Ceq]=constraints_QH(x,A,HfixedS,L,scaling_f,alpha,nl,nn,NumScenarios)
r=x(1:nl);
Qv=x(nl+1:nl+nl*NumScenarios);
Hv=x(nl+nl*NumScenarios+1:nl+nl*NumScenarios+nn*NumScenarios);
Q=reshape(Qv,nl,NumScenarios);
H=reshape(Hv,nn,NumScenarios);
R=L.*r;
RS=repmat(R,1,NumScenarios);
dh=RS/scaling_f.*Q.*abs(Q).^(alpha-1);
con_eq=A*H-HfixedS+dh;
Ceq=con_eq(:);
Cineq=[];
end

