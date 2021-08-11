%% ####################################################################################################################
% Code for the paper:
% Utilizing Matrix Completion for Simulation and Optimization of Water Distribution Networks
% By Mashor Housh, Alaa Jamal
% University of Haifa, mhoush@univ.haifa.ac.il, alaa.gazi.50@gmial.com
%% ####################################################################################################################
% Developed under Matlab 2018b
%% ####################################################################################################################

%% Plotting Figures 3-6 in the paper
clc
clear all 
close all

load('OptimizationResults.mat')

%% Collect results of successful runs
for i=1:7
tmp=Hist_comp.obj(i,1:50);
min_f_comp(i)=min(tmp(Hist_comp.problem(i,1:50)==1 | Hist_comp.problem(i,1:50)==2));
mean_f_comp(i)=mean(tmp(Hist_comp.problem(i,1:50)==1 | Hist_comp.problem(i,1:50)==2));
max_f_comp(i)=max(tmp(Hist_comp.problem(i,1:50)==1 | Hist_comp.problem(i,1:50)==2));

tmp=Hist_comp.cputime(i,1:50);
t_comp(i)=mean(tmp(Hist_comp.problem(i,1:50)==1 | Hist_comp.problem(i,1:50)==2));
end

for i=1:7
tmp=Hist_HQ.obj(i,1:50);
min_f_HQ(i)=min(tmp(Hist_HQ.problem(i,1:50)==1 | Hist_HQ.problem(i,1:50)==2));
mean_f_HQ(i)=mean(tmp(Hist_HQ.problem(i,1:50)==1 | Hist_HQ.problem(i,1:50)==2));
max_f_HQ(i)=max(tmp(Hist_HQ.problem(i,1:50)==1 | Hist_HQ.problem(i,1:50)==2));

tmp=Hist_HQ.cputime(i,1:50);
t_HQ(i)=mean(tmp(Hist_HQ.problem(i,1:50)==1 | Hist_HQ.problem(i,1:50)==2));
end


%% Collect results of failed runs
for i=1:7
fail_comp(i)=sum(~((Hist_comp.problem(i,1:50)==1 | Hist_comp.problem(i,1:50)==2)));
tmp=Hist_comp.cputime(i,1:50);
t_fail_comp(i)=mean(tmp(~((Hist_comp.problem(i,1:50)==1 | Hist_comp.problem(i,1:50)==2))));
end
t_fail_comp(isnan(t_fail_comp))=0; %in case no failure


for i=1:7
fail_HQ(i)=sum(~((Hist_HQ.problem(i,1:50)==1 | Hist_HQ.problem(i,1:50)==2)));
tmp=Hist_HQ.cputime(i,1:50);
t_fail_HQ(i)=mean(tmp(~((Hist_HQ.problem(i,1:50)==1 | Hist_HQ.problem(i,1:50)==2))));
end
t_fail_HQ(isnan(t_fail_HQ))=0;  %in case no failure


%% Plot Figure 3 in the paper
for i=1:7
tmp1=Hist_comp.obj(i,1:50);
tmp2=tmp1(Hist_comp.problem(i,1:50)==1 | Hist_comp.problem(i,1:50)==2);
A(i,1:length(tmp2))=tmp2;
end
A(A==0)=nan;

for i=1:7
tmp1=Hist_HQ.obj(i,1:50);
tmp2=tmp1(Hist_HQ.problem(i,1:50)==1 | Hist_HQ.problem(i,1:50)==2);
B(i,1:length(tmp2))=tmp2;
end
B(B==0)=nan;

C={A' B'};

figure(3)
p1=boxplot(C{1}/(10^6),'Positions',[1,5,9,13,17,21,25]-0.5,'Widths',0.5,'Whisker',1000,'Colors','k');
set(p1,{'linew'},{2})
set(gca,'XTick',[1,5,9,13,17,21,25]-0.5)
hold on
p2=boxplot(C{2}/(10^6),'Positions',[1,5,9,13,17,21,25]+0.5,'Widths',0.5,'Whisker',1000,'Colors','r');
set(p2,{'linew'},{2})
set(gca,'XTick',[1,5,9,13,17,21,25]+0.5)
legend([p1(1),p2(2)],'r,z Formulation','r,H,Q Formulation');
xticklabels({'1' '5' '10' '15' '20' '25' '30'})
xticks([1:4:25])
ylim([5.5 10])
xlabel('Number of Demand Loads')
ylabel('Optimal Cost [$]')

%% Plot Figure 4 in the paper
figure(4)
bar([1 5:5:30],[50-fail_comp;50-fail_HQ]')
xlabel('Number of Demand Loads [-]')
ylabel('Number of Successful Runs Out of 50 [-]')
legend({'r,z Formulation' 'r,H,Q Formulation'})

%% Plot Figure 5 in the paper
figure(5)
bar([1 5:5:30],[t_comp;t_HQ]')
xlabel('Number of Demand Loads [-]')
ylabel('CPU Time [sec]')
legend({'r,z Formulation' 'r,H,Q Formulation'})

%% Plot Figure 6 in the paper
figure(6)
bar([1 5:5:30],[t_fail_comp;t_fail_HQ]')
xlabel('Number of Demand Loads [-]')
ylabel('CPU Time [sec]')
legend({'r,z Formulation' 'r,H,Q Formulation'})






