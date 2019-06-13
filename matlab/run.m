%% Main script to collect all results and perform simualtions

clear all
close all
clc 
rng('default')

%% INPUT

data = readData('../data/dati_Moody.csv');      % reading input data

T     = length(data.years);              % number of data
T_rho = 30;                              % length of tarashev zhu time series

DefRate = data.DR_SG;                    % default rates speculative grade
% DefRate = data.DR_AR;                    % default rates all rates

DefRate_mean = mean(DefRate);            % default rate 
DefRate_std  = std(DefRate);             % mean and standard deviation

d = norminv(DefRate);                         % default barriers
f_0 = @(D) integral(@(x) normcdf(D+sqrt(sum((D-d).^2)/(T-1))*x)...
           .*normpdf(x),-30,30)-DefRate_mean; % condition on d_hat     
d_hat = fzero(f_0,[-3,0]);               % corrected default barrier mean
d_std = sqrt(sum((d_hat-d).^2)/(T-1));   % standard deviation

RecRate_mean = mean(data.RR);            % recovery rate
RecRate_std  = std(data.RR);             % mean and standard deviation

rho_mean = 0.0924;                       % obligors correlation
rho_std  = 0.0386;                       % mean and standard deviation
rho_Basel = correlationFromBasel2(DefRate_mean);  % Basel correlation

N_obligors  = 50;                        % number of obligors
N_sim = 1e6;                             % number of simulations

CL = 0.99;                               % confidence level 99.0 %
% CL  = 0.999;                             % confidence level 99.9 %

idiosyncraticRisk = randn(N_sim,N_obligors); % simulation of idiosyncratic 
systematicRisk    = randn(N_sim,1);          % and systematic risk
        
%% A - COMPUTE Regulatory CapitalS IN THE NOMINAL MODEL
  %% A.i. Using mean correlation
  
disp('Required Capital nominal model')
RC_LHP = CapitalRequirementNominalLHP(RecRate_mean,DefRate_mean,...
         rho_mean,CL);                % For Large Homogeneous Portfolio
RC_HP  = CapitalRequirementNominalHP(RecRate_mean,DefRate_mean,...
         rho_mean,CL,N_obligors);     % For Homogeneous Portfolio

disp(table([RC_LHP;RC_HP], ...
    'RowNames',{'LHP','HP'},'VariableNames',{'Regulatory_Capital'})); 
  %% A.ii. Using basel correlation
  
disp('Regulatory Capital nominal model, basel correlation')
RC_LHP_B = CapitalRequirementNominalLHP(RecRate_mean,DefRate_mean,...
           rho_Basel,CL);             % For Large Homogeneous Portfolio
RC_HP_B  = CapitalRequirementNominalHP(RecRate_mean,DefRate_mean,...
           rho_Basel,CL,N_obligors);  % For Homogeneus Portfolio

disp(table([RC_LHP_B;RC_HP_B],...
    'RowNames',{'LHP','HP'},'VariableNames',{'Regulatory_Capital'}));
  %% A bonus - Check numerical accuracy of Monte Carlo 
 
disp('Check numerical accuracy')
RC_LHP_MC = CapitalRequirementAlternativeLHP(RecRate_mean,DefRate_mean,...
            rho_mean,CL,systematicRisk);% For Large Homogeneous Portfolio
RC_HP_MC  = CapitalRequirementAlternativeHP(RecRate_mean,DefRate_mean,...
            rho_mean,CL,systematicRisk,idiosyncraticRisk); % For Homogeneous Portfolio

disp(table([RC_LHP_MC;RC_HP_MC], ...
    'RowNames',{'LHP','HP'},'VariableNames',{'Regulatory_Capital'}));
%% B - FREQUENTIST INFERENCE
  %% B.1 Verify hypotesis on data
  
% Shapiro wilk normality tests
disp('Default barrier normality, Shapiro-Wilk test')
[H_d,pValue_d] = swtest(d);
if ~H_d
  fprintf('Normality assumption accepted with pvalue of: %.3f\n',pValue_d);
else
  fprintf('Normality assumption rejected with pvalue of: %.3f\n',pValue_d);
end

disp('Recovery rate normality, Shapiro-Wilk test')
[H_RR,pValue_RR] = swtest(data.RR);
if ~H_RR
  fprintf('Normality assumption accepted with pvalue of: %.3f\n',pValue_RR);
else
  fprintf('Normality assumption rejected with pvalue of: %.3f\n',pValue_RR);
end
  %% B.2 - B.3 Simulations required
   
d_sim_freq       = d_hat + randn(N_sim,1)*d_std; % simulated threshold 
DefRate_sim_freq = normcdf(d_sim_freq);          % simulated default rates

RecRate_sim_freq = RecRate_mean + randn(N_sim,1)*RecRate_std;  % simulated recovery 

[alpha,beta] = betaParameter(rho_mean,rho_std);% A & B beta parameters
rho_sim_freq = betarnd(alpha,beta,N_sim,1);    % simulated beta correlation
rho_sim_B_freq = correlationFromBasel2(DefRate_sim_freq);% simulated basel correlation

  %% B.2 Regulatory Capital & Add On using correlation from data

disp('Regulatory Capital & Add On: frequentist inference')

disp('Large Homogeneous Portfolio')
% Computing Regulatory Capitals simulating different parameters LHP CL
RCf_LHP_d     = CapitalRequirementAlternativeLHP(RecRate_mean,...
    DefRate_sim_freq,rho_mean,CL,systematicRisk);     % default rate
RCf_LHP_RR    = CapitalRequirementAlternativeLHP(RecRate_sim_freq,...
    DefRate_mean,rho_mean,CL,systematicRisk);         % recovery rate
RCf_LHP_rho   = CapitalRequirementAlternativeLHP(RecRate_mean,...
    DefRate_mean,rho_sim_freq,CL,systematicRisk);     % correlation
RCf_LHP_d_rho = CapitalRequirementAlternativeLHP(RecRate_mean,...
    DefRate_sim_freq,rho_sim_freq,CL,systematicRisk); % default rate and correlation
RCf_LHP_all   = CapitalRequirementAlternativeLHP(RecRate_sim_freq,...
    DefRate_sim_freq,rho_sim_freq,CL,systematicRisk); % all parameters

% Regulatory Capital alternative models, LHP, CL 
RCaltF_LHP = [RCf_LHP_d;RCf_LHP_RR;RCf_LHP_rho;RCf_LHP_d_rho;RCf_LHP_all];
% Add On alternative models, LHP, CL 
addOn_LHP  = RCaltF_LHP./RC_LHP - 1;

disp(table(RCaltF_LHP,addOn_LHP,...
            'RowNames',{'d','RR','rho','d_rho','All'},...
            'VariableNames',{'Regulatory_Capital','Add_On'}))
        
disp('Homogeneous Portfolio')
% Computing Regulatory Capitals simulating different parameters HP CL
RCf_HP_d     = CapitalRequirementAlternativeHP(RecRate_mean,...
    DefRate_sim_freq,rho_mean,CL,systematicRisk,idiosyncraticRisk);     % default rate
RCf_HP_RR    = CapitalRequirementAlternativeHP(RecRate_sim_freq,...
    DefRate_mean,rho_mean,CL,systematicRisk,idiosyncraticRisk);         % recovery rate
RCf_HP_rho   = CapitalRequirementAlternativeHP(RecRate_mean,...
    DefRate_mean,rho_sim_freq,CL,systematicRisk,idiosyncraticRisk);     % correlation
RCf_HP_d_rho = CapitalRequirementAlternativeHP(RecRate_mean,...
    DefRate_sim_freq,rho_sim_freq,CL,systematicRisk,idiosyncraticRisk); % default rate and correlation
RCf_HP_all   = CapitalRequirementAlternativeHP(RecRate_sim_freq,...
    DefRate_sim_freq,rho_sim_freq,CL,systematicRisk,idiosyncraticRisk); % all parameters 

% Regulatory Capital alternative models, HP, CL
RCaltF_HP = [RCf_HP_d;RCf_HP_RR ;RCf_HP_rho;RCf_HP_d_rho;RCf_HP_all];
% Add On alternative models, HP, CL
addOnF_HP = RCaltF_HP./RC_HP - 1;

disp(table(RCaltF_HP,addOnF_HP,...
            'RowNames',{'d','RR','rho','d_rho','All'},...
            'VariableNames',{'Regulatory_Capital','Add_On'}))
  %% B.3 Regulatory Capital & Add On using correlation from basel

disp('Regulatory Capital & Add On frequentist inference, basel correlation')

disp('Large Homogeneous Portfolio')
% Computing Regulatory Capitals simulating different parameters LHP CL
RCf_LHP_B_d    = CapitalRequirementAlternativeLHP(RecRate_mean,...
    DefRate_sim_freq,rho_sim_B_freq,CL,systematicRisk); % default rate
RCf_LHP_B_RR   = CapitalRequirementAlternativeLHP(RecRate_sim_freq,...
    DefRate_mean,rho_Basel,CL,systematicRisk);          % recovery rate
RCf_LHP_B_d_RR = CapitalRequirementAlternativeLHP(RecRate_sim_freq,...
    DefRate_sim_freq,rho_sim_B_freq,CL,systematicRisk); % default and recovery

% Regulatory Capital alternative models, LHP, CL, Basel
RCaltF_LHP_B = [RCf_LHP_B_d;RCf_LHP_B_RR;RCf_LHP_B_d_RR];
% Add On alternative models, LHP, CL, Basel
addOnF_LHP_B = RCaltF_LHP_B./RC_LHP_B - 1;

disp(table(RCaltF_LHP_B,addOnF_LHP_B,...
            'RowNames',{'d','RR','d_RR'},...
            'VariableNames',{'Regulatory_Capital','Add_On'}))

disp('Homogeneous Portfolio')
% Computing Regulatory Capitals simulating different parameters HP CL
RCf_HP_B_d    = CapitalRequirementAlternativeHP(RecRate_mean,...
    DefRate_sim_freq,rho_sim_B_freq,CL,systematicRisk,idiosyncraticRisk); % default rate
RCf_HP_B_RR   = CapitalRequirementAlternativeHP(RecRate_sim_freq,...
    DefRate_mean,rho_Basel,CL,systematicRisk,idiosyncraticRisk);          % recovery rate
RCf_HP_B_d_RR = CapitalRequirementAlternativeHP(RecRate_sim_freq,...
    DefRate_sim_freq,rho_sim_B_freq,CL,systematicRisk,idiosyncraticRisk); % default and recovery

% Regulatory Capital alternative models, HP, CL, Basel
RCaltF_HP_B = [RCf_HP_B_d;RCf_HP_B_RR ;RCf_HP_B_d_RR];
% Add on alternative models, HP, CL, Basel
addOnF_HP_B = RCaltF_HP_B./RC_HP_B - 1 ;

disp(table(RCaltF_HP_B,addOnF_HP_B,...
            'RowNames',{'d','RR','d_RR'},...
            'VariableNames',{'Regulatory_Capital','Add_On'}))
  %% B bonus check numerical accuracy using close formula

% Test numerical accuracy using close formula in case of uncertainty on d
disp('Large Homogeneous Portfolio, d-close formula')
RCf_LHP_d_close = CapitalRequirementAlternativeLHP_d(RecRate_mean,...
    DefRate_mean,d_hat,d_std,rho_mean,CL);
addOn_F_LHP_d_close = RCf_LHP_d_close/RC_LHP - 1;
disp(table(RCf_LHP_d_close,addOn_F_LHP_d_close,...
    'RowNames',{'d_close'},'VariableNames',{'Regulatory_Capital','Add_On'}))
%% C - BAYESIAN INFERENCE
  %% C.1 Posterior distributions fixing variance equal to empiric
  
     % i. Of threshold d - parameter mean
d_mean_post    = sum(d)/(T+var(d));        % posterion mean
d_std_post = std(d)/sqrt(T+var(d));        % posterior standard deviation
% d_mean_post = d_hat/(1+d_std^2/n);         % internal report way?
% d_std_post = d_std/sqrt(n+d_std^2);
       % distribution of X - default barriers      
d_vect      = linspace(-4,4,1000);                    % equispaced row vector
posterior_d = normpdf(d_vect,d_mean_post,d_std_post); % posterior evaluated on d_vect
xd_vect     = linspace(-4,4,1000);                    % vector on which evaluate posterior predictive distribution
f_x_d = @(x) normpdf(x,d_vect,d_std);                 % density of x given mean(unknown) and empirical std 
d_pdf = bayesianPrediction(xd_vect,d_vect,posterior_d,f_x_d); % density of X conditioned on data

d_sim_bayes  = samplingFromPdf(N_sim,xd_vect,d_pdf);  % simulated d from posterior predictive distribution
DefRate_sim_bayes = normcdf(d_sim_bayes);             % simulated default rate

     % ii. Of correlation rho
rho_vect  = linspace(0.005,0.995,1000);       % equispaced row vector 
[aa,bb]   = betaParameter(rho_vect,rho_std);  % A and B beta parameters   
posterior_rho = posteriorDistributionRho(rho_mean,rho_vect,aa,bb); % posterior distribution on rho_vect
       % distribution of X - correlation   
xrho_vect = linspace(0.005,0.995,1000);       % vector on which evaluate posterior predictive distribution
f_x_rho = @(x) betapdf(x,aa,bb);              % density of x given mean (unknown) and empirical std 
rho_pdf = bayesianPrediction(xrho_vect,rho_vect,posterior_rho,f_x_rho); % density of X conditioned on data

rho_sim_bayes = samplingFromPdf(N_sim,xrho_vect,rho_pdf); % simulated rho
  %% C.2 Regulatory Capital & Add On using correlation from data

disp('Capitar Requirement & Add On, Bayesian inference')

disp('Large Homogeneous Portfolio')
% Computing Regulatory Capitals simulating different parameters LHP CL
RCb_LHP_d   = CapitalRequirementAlternativeLHP(RecRate_mean,...
    DefRate_sim_bayes,rho_mean,CL,systematicRisk);      % default rate
RCb_LHP_rho = CapitalRequirementAlternativeLHP(RecRate_mean,...
    DefRate_mean,rho_sim_bayes,CL,systematicRisk);      % correlation
RCb_LHP_all = CapitalRequirementAlternativeLHP(RecRate_mean,...
    DefRate_sim_bayes,rho_sim_bayes,CL,systematicRisk); % default rate, correlation

% Regulatory Capital alternative models, LHP, CL 
RCaltB_LHP = [RCb_LHP_d;RCb_LHP_rho;RCb_LHP_all];
% Add on alternative models, LHP, CL
addOnB_LHP = RCaltB_LHP./RC_LHP - 1;

disp(table(RCaltB_LHP,addOnB_LHP,...
    'RowNames',{'d','rho','all'},...
    'VariableNames',{'Regulatory_Capital','Add_On'}))

disp('Homogeneous Portfolio')
% Computing Regulatory Capitals simulating different parameters HP CL
RCb_HP_d   = CapitalRequirementAlternativeHP(RecRate_mean,...
    DefRate_sim_bayes,rho_mean,CL,systematicRisk,idiosyncraticRisk);      % default rate
RCb_HP_rho = CapitalRequirementAlternativeHP(RecRate_mean,...
    DefRate_mean,rho_sim_bayes,CL,systematicRisk,idiosyncraticRisk);      % correlation
RCb_HP_all = CapitalRequirementAlternativeHP(RecRate_mean,....
    DefRate_sim_bayes,rho_sim_bayes,CL,systematicRisk,idiosyncraticRisk); % default rate, correlation

% Regulatory Capital alternative models, HP, CL
RCaltB_HP = [ RCb_HP_d ; RCb_HP_rho ; RCb_HP_all ];
% Add on alternative models, HP, CL
addOnB_HP = RCaltB_HP./RC_HP - 1;

disp(table(RCaltB_HP,addOnB_HP,...
    'RowNames',{'d','rho','all'},...
    'VariableNames',{'Regulatory_Capital','Add_On'}))
  %% C.3 Posterior distributions fixing variance equal to Cramer Rao
%-------------------------high computational cost-------------------------% 

     % Compute the Cramer Rao std for different value of d and rho
d_surf   = linspace(-4,4,80);                   % grid on which d_stdCR is  
rho_surf = linspace(0,0.999,20)';               % computed, rho and d

% takes ~4min
% CramerRao_surf  = CramerRao_d(rho_surf,d_surf,N_obligors,T);% Cramer Rao standard dev
CRsurface = struct('rho',rho_surf,'d',d_surf,'surf',CramerRao_surf);

     % i. Of threshold d
d_CramerRao_std = interp2(d_surf,rho_surf,CramerRao_surf,... % Cramer Rao stdev row
                  d_vect,rho_mean,'spline');  % on d_vect, given rho_mean
posterior_d_CramerRao = posteriorDistributionD(d_hat,d_vect,d_CramerRao_std); % posterior distribution d

       % distribution of X - default barriers     
xd_vect = linspace(-5,5,1000);                   % vector on which evaluate posterior prediction density
f_x_d_CramerRao = @(x) normpdf(x,d_vect,d_std);  % density of x given mean(unknown) and std empirical
d_pdf_CramerRao = bayesianPrediction(xd_vect,d_vect,posterior_d_CramerRao,f_x_d_CramerRao); % density of X conditioned on data

d_sim_bayes_CramerRao       = samplingFromPdf(N_sim,xd_vect,d_pdf_CramerRao); % simulated d from its density
DefRate_sim_bayes_CramerRao = normcdf(d_sim_bayes_CramerRao);                 % simulated default rate

   %  ii. Of correlation rho
rho_CramerRao_std = CramerRao_rho(10891,rho_vect,T_rho);       % Cramer Rao stdev 
[aaCR,bbCR]       = betaParameter(rho_vect,rho_CramerRao_std); % A and B beta parameters 
hposterior_rho_CR = posteriorDistributionRho(rho_mean,rho_vect,aaCR,bbCR); % posterior distribution on rho_vect
     % distribution of X - correlation   
xrho_vect = linspace(0.005,0.995,1000);      % vector on which evaluate posterior prediction density        
f_x_rhoCR = @(x) betapdf(x,aa,bb);           % density of x given mean (unknown) and empirical std 

rho_pdf_CramerRao = bayesianPrediction(xrho_vect,rho_vect,hposterior_rho_CR,f_x_rhoCR); % density of X conditioned on data

rho_sim_bayes_CramerRao = samplingFromPosterior(N_sim,xrho_vect,rho_pdf_CramerRao); % simulated correlation

    % iii. nested simulation
n_rho = 1e3;                          % number of rho         
n_d   = 1e3;                          % number of d for each rho
rho_tmp = samplingFromPdf(n_rho,xrho_vect,rho_pdf_CramerRao); 

[rho_nested, DR_nested]=nestedSimulation(d_hat,d_std,n_d,CRsurface,rho_tmp); % paired vector of simulated data
  %% C.4 Regulatory Capital & Add on
  
disp('Homogeneous Portfolio, Cramer Rao')
% Computing Regulatory Capitals simulating different parameters HP CL
RCb_HP_d_CR   = CapitalRequirementAlternativeHP(RecRate_mean,...
    DefRate_sim_bayes_CramerRao,rho_mean,CL,systematicRisk,...
    idiosyncraticRisk); % default rate
RCb_HP_rho_CR = CapitalRequirementAlternativeHP(RecRate_mean,...
    DefRate_mean,rho_sim_bayes_CramerRao,CL,systematicRisk,...
    idiosyncraticRisk); % correlation
RC_HP_allB_CR = CapitalRequirementAlternativeHP(RecRate_mean,...
    DR_nested,rho_nested,CL,systematicRisk,idiosyncraticRisk);  % all

% Regulatory Capital alternative models, HP, CL
RCaltB_HP = [RCb_HP_d_CR;RCb_HP_rho_CR;RC_HP_allB_CR];
% Add on alternative models, HP, CL
addOnB_HP = RCaltB_HP./RC_HP - 1 ;

disp(table(RCaltB_HP,addOnB_HP,...
    'RowNames',{'d','rho','all'},...
    'VariableNames',{'Regulatory_Capital','Add_On'}))
%% Sobol Inidices

Covariance = cov(norminv(data.DR_SG),data.RR);
Covariance(1,2)=0;Covariance(2,1)=0;
mu12 = mean([norminv(data.DR_SG),data.RR]);
[a,b] = betaParameter(rho_mean,rho_std);

sim12 = mvnrnd(mu12,Covariance,N_sim);
sim1 = sim12(:,1);              %barrier
sim2 = sim12(:,2);              %recovery 
sim3 = betarnd(a,b,N_sim,1);    %correlation
parameterSim = [sim1,sim2,sim3];

[S1,~] = SobolInidices(parameterSim);
