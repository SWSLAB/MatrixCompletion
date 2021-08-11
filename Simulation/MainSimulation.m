%% ####################################################################################################################
% Code for the paper:
% Utilizing Matrix Completion for Simulation and Optimization of Water Distribution Networks
% By Mashor Housh, Alaa Jamal
% University of Haifa, mhoush@univ.haifa.ac.il, alaa.gazi.50@gmial.com
%% ####################################################################################################################
% This code requires:
% Matlab Optimizatin toolbox
% Developed under Matlab 2018b
%% ####################################################################################################################

clc
close all
clear all
addpath(genpath(pwd))

% Loop for different random initial solutions
NumRuns=20;               
for i=1:NumRuns
    [It(i),T(i),ErrsQ(i),ErrsH(i),AllSolQ(i),AllSolH(i)]=RunExperiment();
end
FileName=sprintf('Results_%s.mat', datestr(now,'mm-dd-yyyy_HH-MM'));
save(FileName)

% Run experiment procedure
function [It,T,ErrsQ,ErrsH,AllSolQ,AllSolH]=RunExperiment()

% Get Data from INP
Networks={'KL' 'Balerma' 'Fosspoly1' 'Hanoi' 'Jilin'};
for n=1:length(Networks) 
    inpname=[Networks{n} '.inp'];
    d=epanet(inpname);
    NodesElev=d.getNodeElevations;
    L=d.getLinkLength;
    C=d.getLinkRoughnessCoeff;
    D=d.getLinkDiameter;
    basedemand=d.getNodeBaseDemands{1,1};
    factor=d.getOptionsPatternDemandMultiplier;
    Res_id=d.getNodeReservoirIndex;
    OD=d.getLinkNodesIndex;
    
    % Adapt Units
    switch  Networks{n}
        case 'KL'
            L=L*0.3048/1000;
            NodesElev=NodesElev*0.3048;
            D=D*2.54;
            demand=basedemand*0.2271;
        case {'Balerma', 'Fosspoly1' , 'Hanoi'}
            L=L/1000;
            NodesElev=NodesElev;
            D=D/10;
            demand=basedemand;
        case 'Jilin'
            L=L/1000;
            NodesElev=NodesElev;
            D=D/10;
            demand=basedemand*3.6;
    end
    demand=demand*factor;
    demand=double(demand');
    demand(Res_id)=[];
    
    % Build incidence matrix, resistance and fixed head vectors
    nlinks=length(L);
    nnodes=length(demand);
    At=zeros(nnodes,nlinks);
    for i=1:nlinks
        At(OD(i,:),i)=[-1 1];
    end
    ALLlink_fixed_head=[];
    Hfixed=zeros(nlinks,1);
    for i=1:length(Res_id)
        fixed_head=NodesElev(Res_id(i));
        [link_fixed_head1,~]=find(OD(:,1)==Res_id(i));
        [link_fixed_head2,~]=find(OD(:,2)==Res_id(i));
        link_fixed_head=[link_fixed_head1; link_fixed_head2];
        Hfixed(link_fixed_head)=fixed_head;
        ALLlink_fixed_head=[ALLlink_fixed_head;link_fixed_head];
    end
    TMP=At(:,ALLlink_fixed_head);
    TMP(TMP==-1)=1;
    At(:,ALLlink_fixed_head)=TMP;
    At(Res_id,:)=[];
    A=At';
    R=zeros(nlinks,1);
    for i=1:nlinks
        R(i,1)=15267331*L(i)/(C(i)^1.852)/(D(i)^4.871);
    end
    R=double(R);
    
    % Run methods
    QQ=[];
    HH=[];
    methods={'EPANET' 'HQ','Completion','NullSpace','QChord'};
    options = optimoptions('fsolve','Display','none');
    options.Display='iter';
    for m=1:length(methods)
        switch  methods{m}
            case 'EPANET'
                d.openHydraulicAnalysis;
                d.initializeHydraulicAnalysis;
                d.runHydraulicAnalysis;
                % Adapt Units for EPANET solution
                switch  Networks{n}
                    case 'KL'
                        Qsol=d.getLinkFlows'*0.2271;
                        Hsol=d.getNodePressure'*0.7031+d.getNodeElevations'*0.3048;
                    case {'Balerma', 'Fosspoly1' , 'Hanoi'}
                        Qsol=d.getLinkFlows';
                        Hsol=d.getNodePressure'+d.getNodeElevations';
                    case 'Jilin'
                        Qsol=d.getLinkFlows'*3.6;
                        Hsol=d.getNodePressure'+d.getNodeElevations';
                end
                Hsol(Res_id)=[];
                d.closeHydraulicAnalysis;
            case 'HQ'
                tic
                Q0 = 2*sum(demand)*(rand(nlinks,1)-0.5);
                H0 = 2*sum(demand)*(rand(nnodes,1)-0.5);
                x0= [H0;Q0];
                [xsol,~,~,output]=fsolve(@(x)HQ(x,A,R,Hfixed,demand,nnodes),x0,options);
                Qsol=xsol(nnodes+1:end);
                Hsol=xsol(1:nnodes);
                TodT=toc;
                TodIt=output.iterations;
            case 'Completion'
                tic
                nc=nlinks-nnodes;
                z0=2*sum(demand)*(rand(nc,1)-0.5);
                added_rows=rand(nc,nlinks);
                mat1=inv([A' ;added_rows]);
                matT=mat1';
                mat2=matT(nnodes+1:end,:);
                [zsol,~,~,output]=fsolve(@(z)Completion(z,R,Hfixed,demand,mat1,mat2),z0,options);
                Qsol=mat1*[demand;zsol];
                Hsol=matT(1:nnodes,:)*(-diag(R)*diag(Qsol)*abs(Qsol).^0.852+Hfixed);
                ComT=toc;
                ComIt=output.iterations;
            case 'NullSpace'
                tic
                nc=nlinks-nnodes;
                Qp = A'\demand;
                N=null(A');
                dQ0=2*sum(demand)*(rand(nc,1)-0.5);
                [dQsol,~,~,output]=fsolve(@(x)NullSpace(x,R,Qp,N),dQ0,options);
                Qsol=Qp+N*dQsol;
                Hsol=linsolve(A,(-diag(R)*diag(Qsol)*abs(Qsol).^0.852+Hfixed));
                UndT=toc;
                UndIt=output.iterations;
            case 'QChord'
                tic
                nc=nlinks-nnodes;
                N=null(A');
                [~,ST_id] = rref(A');
                chord_id=setdiff(1:nlinks,ST_id);
                Qchord0=2*sum(demand)*(rand(nc,1)-0.5);
                Achord=A(chord_id,:);
                Ast=A(ST_id,:);
                [Qchordsol,~,~,output]=fsolve(@(x)QChord(x,Ast,Achord,ST_id,chord_id,demand,N,R),Qchord0,options);
                Qstsol=inv(Ast')*(demand-Achord'*Qchordsol);
                Qsol(ST_id)=Qstsol;
                Qsol(chord_id)=Qchordsol;
                Hsol=linsolve(A,(-diag(R)*diag(Qsol)*abs(Qsol).^0.852+Hfixed));
                StCT=toc;
                StCIt=output.iterations;
        end
        QQ=[QQ Qsol];
        HH=[HH Hsol];
    end
    QQ(ALLlink_fixed_head,:)=abs(QQ(ALLlink_fixed_head,:));
    
    % Performance
    TodErr=max(abs(QQ(:,1)-QQ(:,2)));
    ComErr=max(abs(QQ(:,1)-QQ(:,3)));
    UndErr=max(abs(QQ(:,1)-QQ(:,4)));
    StCErr=max(abs(QQ(:,1)-QQ(:,5)));
    ErrsQ.(Networks{n})=[TodErr StCErr UndErr ComErr];
    TodErr=max(abs(HH(:,1)-HH(:,2)));
    ComErr=max(abs(HH(:,1)-HH(:,3)));
    UndErr=max(abs(HH(:,1)-HH(:,4)));
    StCErr=max(abs(HH(:,1)-HH(:,5)));
    ErrsH.(Networks{n})=([TodErr StCErr UndErr ComErr]);
    It.(Networks{n})=[TodIt,StCIt,UndIt,ComIt];
    T.(Networks{n})=[TodT,StCT,UndT,ComT];
    AllSolQ.(Networks{n})=QQ;
    AllSolH.(Networks{n})=HH;
end
end


