%% ####################################################################################################################
% Code for the paper:
% Utilizing Matrix Completion for Simulation and Optimization of Water Distribution Networks
% By Mashor Housh, Alaa Jamal
% University of Haifa, mhoush@univ.haifa.ac.il, alaa.gazi.50@gmial.com
%% ####################################################################################################################
% Developed under Matlab 2018b
%% ####################################################################################################################
function [f]=QChord(x,Ast,Achord,ST_id,chord_id,demand,N,R)
Qst=inv(Ast')*(demand-Achord'*x);
Q(ST_id,1)=Qst;
Q(chord_id,1)=x;
f=N'*(diag(R)*diag(Q)*abs(Q).^0.852);
end
