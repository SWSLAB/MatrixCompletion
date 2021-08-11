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
clear all
close all
%% Read problem data from Excel
A=xlsread('FujiwaraData.xlsx','IncidenceMatrix');
A=A(2:end,:);
A=A(:,2:end);
nn=size(A,1);
nl=size(A,2);
data=xlsread('FujiwaraData.xlsx','nodes');
elev=data(:,2);
pmin=data(:,3);
dem=data(:,4);
D=[12,16,20,24,30,40]';            % Avialable diameterts as in Fujiwara and Khang (1990)
C=1.1.*D.^1.5;                     % Capital cost function as in Fujiwara and Khang (1990)
data=xlsread('FujiwaraData.xlsx','pipes');
L=data(:,2);

%% Delete source node from topology and define Hfixed
A=A';
nn=nn-1;
nloop=nl-nn;
res_elev=elev(1);
Hfixed=-A(:,1)*res_elev;     % Assuming the link out from source is in the first column
A(:,1)=[];
dem(1)=[];
elev(1)=[];
pmin(1)=[];

%% Bounds of decision variables
Dmin=min(D);
Dmax=max(D);
scaling_f=1e8;
alpha=1.852;
K=1.526e4;
Rmin=(K/(130^alpha).*((2.54.*Dmax).^(-4.87)))*scaling_f;   % per unit length minimum resistance
Rmax=(K/(130^alpha).*((2.54.*Dmin).^(-4.87)))*scaling_f;   % per unit length maximum resistance

%% Generate random samples of demand
corrmat(1:15,1:15)=0.8;
corrmat(1:15,16:31)=-0.6;
corrmat(16:31,1:15)=-0.6;
corrmat(16:31,16:31)=0.8;
corrmat(eye(nn)==1)=1;
sig=0.1*dem;
covmat=corr2cov(sig,corrmat);
rand('seed',500)
Samples=mvnrnd(dem,covmat,1000)';
Samples=[dem Samples];

%% Solve Optimization
for COMPLETION=0:1      % 0 = HQ method , 1 = Completion
    Runs_vec=[1 5:5:30];
    for r=1:length(Runs_vec)     % Run for different number of demand scenarios
        NumScenarios=Runs_vec(r);
        % Solver options
        options = optimset('fmincon');
        options.Algorithm='sqp-legacy';
        options.FinDiffType='central';
        options.Display='iter';
        options.MaxFunEvals=1e6;
        options.MaxIter=5000;
        options.TolX=1e-10;
        NumStarts=50;
        
        for i=1:NumStarts       % Run from different initial guesses
            if COMPLETION
                % Initial guess
                rand('seed',sum(clock))
                R0=unifrnd(Rmin,Rmax,nl,1);
                Z0=unifrnd(-sum(dem),sum(dem),nloop,NumScenarios);
                x0=[R0;Z0(:)];
                
                % Completion matrices
                rand('seed',500)
                comp_col=rand(nl,nloop);
                Amod=[A comp_col];
                iAmodt=inv(Amod');
                iAmod=inv(Amod);
                
                % Define bounds
                Zmin=-inf*ones(nloop*NumScenarios,1);
                Zmax=inf*ones(nloop*NumScenarios,1);
                LB=[Rmin*ones(nl,1);Zmin];
                UB=[Rmax*ones(nl,1);Zmax];
                
                % Repeat data for scenarios and select scenarios from demand
                elevS=repmat(elev,1,NumScenarios);
                pminS=repmat(pmin,1,NumScenarios);
                HminS=elevS+pminS;
                HfixedS=repmat(Hfixed,1,NumScenarios);
                demS=Samples(:,1:NumScenarios);
                
                % Define objective function and constraints for fmincon
                NLPCON=@(x)constraints_comp(x,iAmodt,iAmod,HfixedS,demS,HminS,L,scaling_f,alpha,nl,nn,nloop,NumScenarios);
                obj=@(x)objective(x,L,scaling_f,K,alpha,nl);
                
                % Solve using fmincon solver
                tic
                [xval,objval,exitflag,output]=fmincon(obj,x0,[],[],[],[],LB,UB,NLPCON,options);
                cpu_time=toc;
            else
                % Initial guess
                rand('seed',sum(clock))
                R0=unifrnd(Rmin,Rmax,nl,1);
                Q0=unifrnd(-sum(dem),sum(dem),nl,NumScenarios);
                H0=unifrnd(elev(1)+pmin(1),Hfixed(1),nn,NumScenarios);
                x0=[R0;Q0(:);H0(:)];
                
                % Define linear constraints
                Aeq=[];
                beq=[];
                demS=Samples(:,1:NumScenarios);
                for s=1:NumScenarios
                    Aeq=blkdiag(Aeq,A');
                    beq=[beq;demS(:,s)];
                end
                Aeq=[zeros(nn*NumScenarios,nl) Aeq zeros(nn*NumScenarios,nn*NumScenarios)];
                
                % Repeat data for scenarios and define bounds
                HfixedS=repmat(Hfixed,1,NumScenarios);
                elevS=repmat(elev,1,NumScenarios);
                pminS=repmat(pmin,1,NumScenarios);
                HminS=elevS+pminS;
                Hmin=HminS(:);
                Hmax=inf*ones(nn*NumScenarios,1);
                Qmin=-inf*ones(nl*NumScenarios,1);
                Qmax=inf*ones(nl*NumScenarios,1);
                LB=[Rmin*ones(nl,1);Qmin;Hmin];
                UB=[Rmax*ones(nl,1);Qmax;Hmax];
                
                % Define objective function and constraints for fmincon
                NLPCON=@(x)constraints_QH(x,A,HfixedS,L,scaling_f,alpha,nl,nn,NumScenarios);
                obj=@(x)objective(x,L,scaling_f,K,alpha,nl);
                
                % Solve using fmincon solver
                tic
                [xval,objval,exitflag,output]=fmincon(obj,x0,[],[],Aeq,beq,LB,UB,NLPCON,options);
                cpu_time=toc;
            end
            
            % Collect results from solver
            if COMPLETION
                Rval=xval(1:nl);
                Zv=xval(nl+1:nl+nloop*NumScenarios);
                Zval=reshape(Zv,nloop,NumScenarios);
                Qval=iAmodt*[demS;Zval];
                RS=repmat(Rval,1,NumScenarios);
                dh=RS/scaling_f.*Qval.*abs(Qval).^(alpha-1);
                tmp=iAmod*(HfixedS-dh);
                Hval=tmp(1:nn,:);
                Dval=((Rval/scaling_f*(130^alpha)/K./L).^(-1/4.87))/2.54;
                
                Hist_comp.problem(r,i)=exitflag;
                Hist_comp.R{r}(:,i)=Rval;
                Hist_comp.D{r}(:,i)=Dval;
                Hist_comp.Q{r}{i}=Qval;
                Hist_comp.H{r}{i}=Hval;
                Hist_comp.obj(r,i)=objval;
                Hist_comp.cputime(r,i)=cpu_time;
            else
                Rval=xval(1:nl);
                Qv=xval(nl+1:nl+nl*NumScenarios);
                Hv=xval(nl+nl*NumScenarios+1:nl+nl*NumScenarios+nn*NumScenarios);
                Qval=reshape(Qv,nl,NumScenarios);
                Hval=reshape(Hv,nn,NumScenarios);
                Dval=((Rval/scaling_f*(130^alpha)/K./L).^(-1/4.87))/2.54;
                
                Hist_HQ.problem(r,i)=exitflag;
                Hist_HQ.R{r}(:,i)=Rval;
                Hist_HQ.D{r}(:,i)=Dval;
                Hist_HQ.Q{r}{i}=Qval;
                Hist_HQ.H{r}{i}=Hval;
                Hist_HQ.obj(r,i)=objval;
                Hist_HQ.cputime(r,i)=cpu_time;
            end
        end
    end
end



