%% ####################################################################################################################
% Code for the paper:
% Utilizing Matrix Completion for Simulation and Optimization of Water Distribution Networks
% By Mashor Housh, Alaa Jamal
% University of Haifa, mhoush@univ.haifa.ac.il, alaa.gazi.50@gmial.com
%% ####################################################################################################################
% Developed under Matlab 2018b
%% ####################################################################################################################

%% Plotting Figures 1 and 2 in the paper
clc
clear all 
close all

Networks={'KL' 'Balerma' 'Fosspoly1' 'Hanoi' 'Jilin'};
FileName='SimulationResults.mat';
load(FileName)

for i=1:length(Networks)
    AvgErrs(i,:)=mean(reshape([ErrsQ.(Networks{i})],4,[])');
    MaxErrs(i,:)=max(reshape([ErrsQ.(Networks{i})],4,[])');
    MinErrs(i,:)=min(reshape([ErrsQ.(Networks{i})],4,[])');
    AvgCPU(i,:)=mean(reshape([T.(Networks{i})],4,[])');
    MaxCPU(i,:)=max(reshape([T.(Networks{i})],4,[])');
    MinCPU(i,:)=min(reshape([T.(Networks{i})],4,[])');
end


%% Figure 1 in the paper
c=0;
for i=1:5
    for j=1:4
        c=c+1;
        xerrbar(c)=i-0.25+0.7*0.25*(j-1);
    end
end
figure(1)
bar(AvgErrs)
ylabel({'MAE [m^3 hr^-^3]'})
set(gca, 'YScale', 'log')
hold on
MaxErrs=MaxErrs';
MinErrs=MinErrs';
e=errorbar(xerrbar,(MaxErrs(:)+MinErrs(:))*0.5,(MaxErrs(:)-MinErrs(:))*0.5,'k');
e.LineWidth=1.5;
set(gca, 'YScale', 'log')
e.LineStyle = 'none';
xticklabels(Networks)
legend({'H-Q' 'Qchord' 'Null-Space' 'Completion'})


%% Figure 2 in the paper
c=0;
for i=1:5
    for j=1:4
        c=c+1;
        xcpubar(c)=i-0.25+0.7*0.25*(j-1);
    end
end
figure(2)
bar(AvgCPU)
ylabel({'Average CPU Time [sec]'})
set(gca, 'YScale', 'log')
hold on
MaxCPU=MaxCPU';
MinCPU=MinCPU';
e=errorbar(xcpubar,(MaxCPU(:)+MinCPU(:))*0.5,(MaxCPU(:)-MinCPU(:))*0.5,'k');
e.LineWidth=1.5;
clear temp
set(gca, 'YScale', 'log')
e.LineStyle = 'none';
xticklabels(Networks)
legend({'H-Q' 'Qchord' 'Null-Space' 'Completing'})

