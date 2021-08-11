%% ####################################################################################################################
% Code for the paper:
% Utilizing Matrix Completion for Simulation and Optimization of Water Distribution Networks
% By Mashor Housh, Alaa Jamal
% University of Haifa, mhoush@univ.haifa.ac.il, alaa.gazi.50@gmial.com
%% ####################################################################################################################
% Developed under Matlab 2018b
%% ####################################################################################################################
function [Cineq,Ceq]=constraints_comp(x,iAmodt,iAmod,HfixedS,demS,HminS,L,scaling_f,alpha,nl,nn,ncomp,NumScenarios)
r=x(1:nl);
Zv=x(nl+1:nl+ncomp*NumScenarios);
Z=reshape(Zv,ncomp,NumScenarios);
Q=iAmodt*[demS;Z];
R=L.*r;
RS=repmat(R,1,NumScenarios);
dh=RS/scaling_f.*Q.*abs(Q).^(alpha-1);
tmp=iAmod*(HfixedS-dh);
H=tmp(1:nn,:);
w=tmp(nn+1:end,:);
                
Ceq=w(:);
Cineq=(HminS-H);
Cineq=Cineq(:);

end

