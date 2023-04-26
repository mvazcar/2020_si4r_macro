% This file solves and simulates the SI4R-Macro Model, an extension for the model developed in Eichenbaum,
%Rebelo and Trabandt (2020),'The Macroeconomics of Epidemics'.
 
%Matlab 2019b used for calculations.
 


clear all; clc; close all; tic;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameters, calibration targets and other settings%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
betta=0.96^(1/52);  %Weekly household discount factor
pid=7*0.005/18;     %Weekly probability of dying
pir=7*1/18-pid;     %Weekly probability of recovering
phii=0.8;           %Productivity of infected people
sirmacro=1;         % =1 it computes SIR-Macro instead of SI4R.
optimalexogenous=1; % =1 sets SIR-Macro optimal muc as exogenous path for muc.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SI4R Extra Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pia=0.5;
pib=1-pia;
pit=0.5;

if sirmacro==1
    pia=0;
    pib=1-pia;
    pit=0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calibration targets for hours and income
n_target=28;         %Weekly hours
inc_target=58000/52; %weekly income
 
%Calibation targets for shares of pis-terms in T-function in SIR model
pis3_shr_target=2/3;                   %share of T_0 jump due general infections
pis1_shr_target=(1-pis3_shr_target)/2; %share of T_0 jump due to consumption-based infections
pis2_shr_target=(1-pis3_shr_target)/2; %share of T_0 jump due to work-based infections
RplusD_target=0.60;                    %total share of people infected and then either recovered or dead after epidemic
 
pop_ini=1;          %Initial population
i_ini=0.001;        %Initial infected
 
HH=250;             %Number of periods to solve and simulate the model
 
%containment policy

muc=zeros(HH,1);    %exogenous path for muc over time.
                    %if you want e.g. a containment policy of
                    %10 percent for weeks 1,...52, then set muc(1:52)=0.1;
                    %
                    %Make sure that muc=0 at the end of the solution and
                    %simulation horizon. Steady state assumes muc=0;
                    %
                    %With optimal policy, muc path will be chosen to
                    %maximize PV utility (switch below).
                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SI4R Extra Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if optimalexogenous==1 % Use Optimal muc of SIR Macro in SI4R as exogenous path for muc.
load last_solution_opt_policy; 
muc=muc_sol;    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
do_opt_policy=0;    %switch: if = 0, model is solved and simulated with
                    %given path for containment policy muc
                    %
                    %if =1, model is solved and simulated with optimal
                    %containment path muc. Path for muc set above will be
                    %overwritten with optimal path. COMPUTATIONS TAKE A
                    %WHILE IF YOU SET do_opt_taxes=1
                    
use_parallel=0;     %when optimal policy is computed, use_parallel=1 uses 
                    %parallel computing to maximize PV utility using fmincon.
 
%nonlinear solver and minimizer settings
opts_fsolve=optimoptions('fsolve','Algorithm', 'levenberg-marquardt', 'Display','iter'); %options for fsolve - MODIFIED TO USE MARQUARDT.
opts_fsolve_fmincon=optimoptions('fsolve','Algorithm', 'levenberg-marquardt', 'Display','iter'); %options for fsolve used opt. policy calcs. (fmincon)
% opts_fsolve_fmincon=optimoptions('fsolve','Display','off','TolFun',1e-9); %options for fsolve used opt. policy calcs. (fmincon)

if use_parallel==0
    opts_fmincon=optimoptions('fmincon','Display','iter','TolFun',1e-7,'MaxFunctionEvaluations',2000,'FiniteDifferenceStepSize',1e-2); %options for fmincon w/o parallel comp.
elseif use_parallel==1
    %opts_fmincon=optimoptions('fmincon','Display','iter','TolFun',1e-7,'MaxFunctionEvaluations',2000,'UseParallel',true,'FiniteDifferenceStepSize',1e-2); %options for fmincon with parallel comp.
    %opts_fmincon=optimoptions('fmincon','Display','iter','TolFun',1e-6,'MaxFunctionEvaluations',10000,'UseParallel',true,'FiniteDifferenceStepSize',1e-2); %options for fmincon with parallel comp.
    %opts_fmincon=optimoptions('fmincon','Display','iter','TolFun',1e-6,'MaxFunctionEvaluations',5000,'UseParallel',true,'FiniteDifferenceStepSize',1e-3);
    opts_fmincon=optimoptions('fmincon','Display','iter','TolFun',1e-7,'MaxFunctionEvaluations',10000,'UseParallel',true);
end
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Steady State Calculations%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta=1/n_target^2;     %calculate disutility of labor parameter theta
                        %so that pre-infection steady state labor is equal to
                        %desired target (using n=(1/theta)^(1/2) pre-infection 
                        %steady state equation)
A=inc_target/n_target;  %calculate parameter A such that income is equal to
                        %desired target, c=inc_target=A*n_target
 
%steady states
nr1ss=(1/theta)^(1/2);           %labor recovered (same as post-infection steady state)
cr1ss=A*nr1ss;                    %consumption recovered
ur1ss=log(cr1ss)-theta/2*nr1ss^2;  %utility recovered
Ur1ss=1/(1-betta)*ur1ss;          %PV utility recovered
Ur1ssConsUnits=Ur1ss*cr1ss;        %PV utility in cons. units (Ur1ss*Marg.Util.Cons); value of life

nib1ss=(1/theta)^(1/2);           %labor infected
cib1ss=phii*A*nib1ss;               %consumption infected
uib1ss=log(cib1ss)-theta/2*nib1ss^2;  %utility infected
Uib1ss=(1/(1-betta*(1-pir-pid)))*(uib1ss+betta*pir*Ur1ss);  %PV utility infected

nia1ss=(1/theta)^(1/2);           %labor infected
cia1ss=phii*A*nia1ss;               %consumption infected
uia1ss=log(cia1ss)-theta/2*nia1ss^2;  %utility infected
Uia1ss=(1/(1-betta*(1-pir-pid)))*(uia1ss+betta*pir*Ur1ss);  %PV utility infected

 
%Check level of present value utility
if Uib1ss-Ur1ss>0, error(['Error: parameterization implies Uib1ss>Ur1ss: ',num2str(Uib1ss-Ur1ss)]);end
 
%calibrate the pis's in T-function
go_calibrate_pis;
 
%initial guess for optimal muc and load last solution of optimal policy 
%allocations if you dont want to start maximization of PV utility from scratch.
if do_opt_policy==1      
    %muc_guess=zeros(HH,1);                  %initial guess for opt. cont. policy
    
    load last_solution_opt_policy    %uncomment if you want to use last
                                      %solution as initial guess
    muc_guess=muc_sol;
    if numel(muc_guess)~=HH
        error('Initial guess for optimal policy loaded from disk has different dimension than HH.');
    end
end
 
%initial guess of vectors of ns, ni and nr to solve nonlinear
%equilibrium equations
n_vec_guess=nr1ss*ones(7*HH,1); %guess of vectors for ns,ni,nr
 
%If optimal policy is desired, find optimal path for muc to maximize PV utility
if do_opt_policy==1 %optimal policy
    
    %minimize negative PV utility (i.e. max PV utility) to find opt. path for muc; nonlinear
    %model equations are solved inside getU.m using function get_err.m
    LB=muc_guess*0-2;UB=muc_guess*0+2;%lower and upper bounds
   
    muc_sol = fmincon(@getU,muc_guess,[],[],[],[],LB,UB,[],opts_fmincon,n_vec_guess,opts_fsolve_fmincon,A,theta,i_ini,pop_ini,pis1,pis2,pis3,pir,pid,betta,Uib1ss,Uia1ss,HH,cr1ss,nr1ss,Ur1ss,phii,pia,pib,pit);

    muc=zeros(HH,1);muc(1:numel(muc_sol))=muc_sol;
    save last_solution_opt_policy muc_sol;%save solution for possible use in subsequent maximization   
end
 
%Given either optimal path for muc or exogenous path for muc,
%solve nonlinear equilibrium model equations (i.e. adjust guesses ns,nia0,nia1,nib0,nib1,nr0,nr1)
[n_vec,fval,exitflag]=fsolve(@get_err,n_vec_guess,opts_fsolve,A,theta,i_ini,pop_ini,pis1,pis2,pis3,pir,pid,betta,Uib1ss,Uia1ss,HH,cr1ss,nr1ss,Ur1ss,muc,phii,pia,pib,pit);

if exitflag==1
    display('Equation solved. First Order Optimality is Small')
elseif exitflag==2
    display('Equation solved. The change is smaller than the tolerance specified')
elseif exitflag==3
    display('Equation solved. The change in the residuals is less than specified tolerance')
elseif exitflag==4
    display('Equation solved. The magnitude of the search direction is less than the specified tolerance')
else
    error('Fsolve could not solve the model')
end    
    
%get allocations given either exogenous or optimal path for muc at ns,nia0,nia1,nib0,nib1,nr0,nr1
%solution
[err,I,Ia0,Ia1,Ib0,Ib1,S,R,R0,R1,D,T,Pop,RnotSIRmacro,aggC,aggH,ns,nia0,nia1,nib0,nib1,nr0,nr1,cs,cia0,cia1,cib0,cib1,cr0,cr1,Us,Uia0,Uia1,Uib0,Uib1,Ur0,Ur1,U,probSD,probA0D,probR0D] = get_err(n_vec,A,theta,i_ini,pop_ini,pis1,pis2,pis3,pir,pid,betta,Uib1ss,Uia1ss,HH,cr1ss,nr1ss,Ur1ss,muc,phii,pia,pib,pit);

disp(['Max. abs. error in equilib. equations:',num2str(max(abs(err)))]);
disp(' ');
RnotSIRmacro;

%Save workspace
save results_baseline
 
%Plots
sir_plots;
si4r_plots;
  
%output some data used in paper and robustness table
aggCons_trough_percent=min((100*(aggC-cr1ss)/cr1ss))
aggCons_avg_first_year_percent=mean((100*(aggC(1:52)-cr1ss)/cr1ss))
terminal_one_minus_susceptibles_percent=100*(1-S(end))
peak_infection_percent=max(100*I)
terminal_death_share_percent=100*D(end)
terminal_number_deaths_US_millions=terminal_death_share_percent/100*330
 
toc;

