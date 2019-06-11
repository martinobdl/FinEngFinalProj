%%

clear all
close all
clc
rng('default')

%% INPUT

data = readData('../data/dati_Altman.csv');      % reading input data

DR = data.DR_SG;                 % data considered speculative grade
% DR = data.DR_AR;                 % data considered all rates

DR_mean = mean(DR);              % default rate
DR_std  = std(DR);               % mean and standard deviation

d_std = std(norminv(DR));                % standard deviation threshold
f_0   = @(D) integral(@(x) normcdf(D+d_std.*x).*normpdf(x),-30,30)-DR_mean;
d_hat = fzero(f_0,[-2,0]);               % unbiased default barrier mean

RR_mean = mean(data.RR);                 % recovery rate
RR_std  = std(data.RR);                  % mean and standard deviation

rho_mean = 0.0924;                       % obligors correlation
rho_std  = 0.0386;                       % mean and standard deviation
rho_B = correlationFromBasel2(DR_mean);  % Basel correlation

N_ob  = 60;                              % number of obligors
N_sim = 1e6;                             % number of simulations

% CL = 0.99;                               % confidence level 99.0 %
CL  = 0.999;                             % confidence level 99.9 %

idiosyncraticRisk = randn(N_sim,N_ob);   % simulation of idiosyncratic risk
systematicRisk    = randn(N_sim,1);      % simulation of systematic risk

%% A - COMPUTE Regulatory CapitalS IN THE NOMINAL MODEL
  %% A.i. Using mean correlation
disp('Required Capital nominal model')

    % For Large Homogeneous Portfolio
RC_LHP = CapitalRequirementNominalLHP(RR_mean,DR_mean,rho_mean,CL);

    % For Homogeneous Portfolio
RC_HP = CapitalRequirementNominalHP(RR_mean,DR_mean,rho_mean,CL,N_ob);

disp(table([RC_LHP;RC_HP], ...
    'RowNames',{'LHP','HP'},'VariableNames',{'Regulatory_Capital'}));
 %% A.ii. Using basel correlation
disp('Regulatory Capital nominal model, basel correlation')

    % For Large Homogeneous Portfolio
RC_LHP_B = CapitalRequirementNominalLHP(RR_mean,DR_mean,rho_B,CL);

    % For Homogeneus Portfolio
RC_HP_B = CapitalRequirementNominalHP(RR_mean,DR_mean,rho_B,CL,N_ob);

disp(table([RC_LHP_B;RC_HP_B],...
    'RowNames',{'LHP','HP'},'VariableNames',{'Regulatory_Capital'}));
%% A bonus - Check numerical accuracy of Monte Carlo
%
%  % For Large Homogeneous Portfolio
% RC_LHP_MC = CapitalRequirementAlternativeLHP(RR_mean,DR_mean,rho_mean,...
%                     CL,systematicRisk);
%
%     % For Homogeneous Portfolio
% RC_HP_MC = CapitalRequirementAlternativeHP(RR_mean,DR_mean,rho_mean,...
%                     CL,systematicRisk,idiosyncraticRisk);
%
% disp(table([RC_LHP_MC;RC_HP_MC], ...
%     'RowNames',{'LHP','HP'},'VariableNames',{'Regulatory Capital'}));
%% B - FREQUENTISTIC INFERENCE
  %% B.1 Verify hypotesis on data

d = norminv(DR);           % default barrier

% Shapiro wilk normality tests

disp('Default barrier normality, Shapiro-Wilk test')
[H_d,pValue_d] = swtest(d)
disp('Recovery rate normality, Shapiro-Wilk test')
[H_RR,pValue_RR] = swtest(data.RR)

 %% B.2 - B.3 Simulations required

d_sim_f  = d_hat + randn(N_sim,1)*d_std;     % simulated gaussian threshold
DR_sim_f = normcdf(d_sim_f);                 % simulated default rates

RR_sim_f = RR_mean + randn(N_sim,1)*RR_std;  % simulated gaussian recovery

[alpha,beta] = betaParameter(rho_mean,rho_std);% A & B beta parameters
rho_sim_f = betarnd(alpha,beta,N_sim,1);      % simulated beta correlation

rho_sim_B_f = correlationFromBasel2(DR_sim_f);% simulated basel correlation

%% B.2 Regulatory Capital & Add On using correlation from data

disp('Regulatory Capital & Add On: frequentist inference')

disp('Large Homogeneous Portfolio')
% Computing Regulatory Capitals simulating different parameters LHP CL
RCf_LHP_d     = CapitalRequirementAlternativeLHP(RR_mean,DR_sim_f,...
               rho_mean,CL,systematicRisk);  % default rate
RCf_LHP_RR    = CapitalRequirementAlternativeLHP(RR_sim_f,DR_mean,...
               rho_mean,CL,systematicRisk);  % recovery rate
RCf_LHP_rho   = CapitalRequirementAlternativeLHP(RR_mean,DR_mean,...
               rho_sim_f,CL,systematicRisk); % correlation
RCf_LHP_d_rho = CapitalRequirementAlternativeLHP(RR_mean,DR_sim_f,...
               rho_sim_f,CL,systematicRisk); % default rate and correlation
RCf_LHP_all   = CapitalRequirementAlternativeLHP(RR_sim_f,DR_sim_f,...
               rho_sim_f,CL,systematicRisk); % all parameters

% Regulatory Capital alternative models, LHP, CL
RCaltF_LHP = [RCf_LHP_d;RCf_LHP_RR;RCf_LHP_rho;RCf_LHP_d_rho;RCf_LHP_all];
% Add On alternative models, LHP, CL
addOn_LHP = RCaltF_LHP./RC_LHP - 1;

disp(table(RCaltF_LHP,addOn_LHP,...
            'RowNames',{'d','RR','rho','d_rho','All'},...
            'VariableNames',{'Regulatory_Capital','Add_On'}))

disp('Homogeneous Portfolio')
% Computing Regulatory Capitals simulating different parameters HP CL
RCf_HP_d     = CapitalRequirementAlternativeHP(RR_mean,DR_sim_f,rho_mean,...
      CL,systematicRisk,idiosyncraticRisk); % default rate
RCf_HP_RR    = CapitalRequirementAlternativeHP(RR_sim_f,DR_mean,rho_mean,...
      CL,systematicRisk,idiosyncraticRisk); % recovery rate
RCf_HP_rho   = CapitalRequirementAlternativeHP(RR_mean,DR_mean,rho_sim_f,...
      CL,systematicRisk,idiosyncraticRisk); % correlation
RCf_HP_d_rho = CapitalRequirementAlternativeHP(RR_mean,DR_sim_f,rho_sim_f,...
      CL,systematicRisk,idiosyncraticRisk); % default rate and correlation
RCf_HP_all   = CapitalRequirementAlternativeHP(RR_sim_f,DR_sim_f,rho_sim_f,...
      CL,systematicRisk,idiosyncraticRisk); % all parameters

% Regulatory Capital alternative models, HP, CL
RCaltF_HP = [RCf_HP_d;RCf_HP_RR ;RCf_HP_rho;RCf_HP_d_rho;RCf_HP_all];
% Add On alternative models, HP, CL
addOnF_HP = RCaltF_HP./RC_HP - 1;

disp(table(RCaltF_HP,addOnF_HP,...
            'RowNames',{'d','RR','rho','d_rho','All'},...
            'VariableNames',{'Regulatory_Capital','Add_On'}))

%% B.3

disp('Regulatory Capital & Add On frequentist inference, basel correlation')

disp('Large Homogeneous Portfolio')
% Computing Regulatory Capitals simulating different parameters LHP CL
RCf_LHP_B_d    = CapitalRequirementAlternativeLHP(RR_mean,DR_sim_f,...
                rho_sim_B_f,CL,systematicRisk);  % default rate
RCf_LHP_B_RR   = CapitalRequirementAlternativeLHP(RR_sim_f,DR_mean,...
                rho_B,CL,systematicRisk);        % recovery rate
RCf_LHP_B_d_RR = CapitalRequirementAlternativeLHP(RR_sim_f,DR_sim_f,...
                rho_sim_B_f,CL,systematicRisk);  % default and recovery

% Regulatory Capital alternative models, LHP, CL, Basel
RCaltF_LHP_B = [RCf_LHP_B_d;RCf_LHP_B_RR;RCf_LHP_B_d_RR];
% Add On alternative models, LHP, CL, Basel
addOnF_LHP_B = RCaltF_LHP_B./RC_LHP_B - 1;

disp(table(RCaltF_LHP_B,addOnF_LHP_B,...
            'RowNames',{'d','RR','d_RR'},...
            'VariableNames',{'Regulatory_Capital','Add_On'}))

disp('Homogeneous Portfolio')
% Computing Regulatory Capitals simulating different parameters HP CL
RCf_HP_B_d    = CapitalRequirementAlternativeHP(RR_mean,DR_sim_f,...
   rho_sim_B_f,CL,systematicRisk,idiosyncraticRisk); % default rate
RCf_HP_B_RR   = CapitalRequirementAlternativeHP(RR_sim_f,DR_mean,...
   rho_B,CL,systematicRisk,idiosyncraticRisk);       % recovery rate
RCf_HP_B_d_RR = CapitalRequirementAlternativeHP(RR_sim_f,DR_sim_f,...
   rho_sim_B_f,CL,systematicRisk,idiosyncraticRisk); % default and recovery

% Regulatory Capital alternative models, HP, CL, Basel
RCaltF_HP_B = [RCf_HP_B_d;RCf_HP_B_RR ;RCf_HP_B_d_RR];
% Add on alternative models, HP, CL, Basel
addOnF_HP_B = RCaltF_HP_B./RC_HP_B - 1 ;

disp(table(RCaltF_HP_B,addOnF_HP_B,...
            'RowNames',{'d','RR','d_RR'},...
            'VariableNames',{'Regulatory_Capital','Add_On'}))

%% B bonus check numerical accuracy using close formula
%
% % Test numerical accuracy using close formula in case of uncertainty on d
% disp('Large Homogeneous Portfolio, d-close formula')
% RCf_LHP_d_close = CapitalRequirementAlternativeLHP_d(RR_mean,DR_mean,...
%     d_hat,d_std,rho_mean,CL);
%
% disp(table(RCf_LHP_CL_d_close,...
%     'RowNames',{'d_close'},'VariableNames',{'Regulatory_Capital'}))
%% C - BAYESIAN INFERENCE
  %% C.1 Posterior distributions fixing variance equal to empiric

     % i. Of threshold d

n       = length(data.years);                % number of data
mu_d    = d_hat/(1+d_std^2/n);               % posterion mean
sigma_d = d_std/sqrt(1+d_std^2/n);           % posterior standard deviation
fprintf('The posterior distribution of d is gaussian\n')
fprintf('with mean %f and standard deviation %f\n', mu_d,sigma_d)

d_sim_b  = mu_d + randn(N_sim,1)*sigma_d;  % simulated d from posterior
DR_sim_b = normcdf(d_sim_b);               % simulated default rate

     % ii. Of correlation rho

rho_vect  = linspace(0.005,0.995,1980);       % equispaced vector
[aa,bb]   = betaParameter(rho_vect,rho_std);  % A and B beta parameters
h_rho     = posteriorDistributionRho(rho_mean,rho_vect,aa,bb); % posterior
mu_rho    = trapz(rho_vect,rho_vect.*h_rho);
sigma_rho = sqrt(trapz(rho_vect,rho_vect.^2.*h_rho)-mu_rho.^2);
fprintf('The posterior distribution of rho has\n')
fprintf('mean %f and standard deviation %f\n', mu_rho,sigma_rho)

rho_sim_b = samplingFromPosterior(N_sim,rho_vect,h_rho); % simulated rho
  %% C.2 Regulatory Capital & Add On using correlation from data

disp('Capitar Requirement & Add On, Bayesian inference')

disp('Large Homogeneous Portfolio')
% Computing Regulatory Capitals simulating different parameters LHP CL
RCb_LHP_d   = CapitalRequirementAlternativeLHP(RR_mean,DR_sim_b,...
    rho_mean,CL,systematicRisk);                % default rate
RCb_LHP_rho = CapitalRequirementAlternativeLHP(RR_mean,DR_mean,...
    rho_sim_b,CL,systematicRisk);               % correlation
RCb_LHP_all = CapitalRequirementAlternativeLHP(RR_mean,DR_sim_b,...
    rho_sim_b,CL,systematicRisk);               % default rate, correlation

% Regulatory Capital alternative models, LHP, CL
RCaltB_LHP = [RCb_LHP_d;RCb_LHP_rho;RCb_LHP_all];
% Add on alternative models, LHP, CL 99.0%
addOnB_LHP = RCaltB_LHP./RC_LHP - 1 ;

disp(table(RCaltB_LHP,addOnB_LHP,...
    'RowNames',{'d','rho','all'},...
    'VariableNames',{'Regulatory_Capital','Add_On'}))

% For Homogeneous Portfolio
disp('Homogeneous Portfolio')

% Computing Regulatory Capitals simulating different parameters HP CL
RCb_HP_d   = CapitalRequirementAlternativeHP(RR_mean,DR_sim_b,...
    rho_mean,CL,systematicRisk,idiosyncraticRisk);  % default rate
RCb_HP_rho = CapitalRequirementAlternativeHP(RR_mean,DR_mean,...
    rho_sim_b,CL,systematicRisk,idiosyncraticRisk); % correlation
RC_HP_allB = CapitalRequirementAlternativeHP(RR_mean,DR_sim_b,...
    rho_sim_b,CL,systematicRisk,idiosyncraticRisk); % default rate, correlation

% Regulatory Capital alternative models, HP, CL
RCaltB_HP = [ RCb_HP_d ; RCb_HP_rho ; RC_HP_allB ];
% Add on alternative models, HP, CL
addOnB_HP = RCaltB_HP./RC_HP - 1 ;

disp(table(RCaltB_HP,addOnB_HP,...
    'RowNames',{'d','rho','all'},...
    'VariableNames',{'Regulatory_Capital','Add_On'}))
  %% C.3 Posterior distributions fixing variance equal to Cramer Rao

     % i. Of threshold d

d_vect    = linspace(-4,0,50);                   % equispaced row vector
d_CRstd   = CramerRao_d(rho_mean,d_vect,N_ob,20);% Cramer Rao stdev
h_d_CR   = posteriorDistributionD(d_hat,d_vect,d_CRstd);
mu_dCR    = trapz(d_vect,d_vect.*h_d_CR);
sigma_dCR = sqrt(trapz(d_vect,d_vect.^2.*h_d_CR)-mu_dCR.^2);

d_sim_bCR  = samplingFromPosterior(N_sim,d_vect,h_d_CR); % simulated threshold
DR_sim_bCR = normcdf(d_sim_bCR);                  % simulated default rate

   %  ii. Of correlation rho

rho_vect  = linspace(0.005,0.995,1980);       % equispaced column vector
rho_CRstd = CramerRao_rho(10891,rho_vect,30);  % Cramer Rao stdev
[aa,bb]   = betaParameter(rho_vect,rho_CRstd); % A and B beta parameters
h_rho_CR  = posteriorDistributionRho(rho_mean,rho_vect,aa,bb); % posterior
mu_rhoCR    = trapz(rho_vect,rho_vect.*h_rho_CR);
sigma_rhoCR = sqrt(trapz(rho_vect,rho_vect.^2.*h_rho_CR)-mu_rhoCR.^2);

rho_sim_bCR = samplingFromPosterior(N_sim,rho_vect,h_rho_CR);

    % iii. nested simulation

d_surf   = linspace(-4,4,41);
rho_surf = linspace(0,0.999,31)';
CR_surf  = CramerRao_d(rho_surf,d_surf,N_ob,20);
n_rho = 1e3;
n_d = 1e3;

rho_tmp    = samplingFromPosterior(n_rho,rho_vect,h_rho_CR);
rho_nested = zeros(n_rho*n_d,1);
d_nested   = zeros(n_rho*n_d,1);
f = waitbar(0,'');
for i = 1 : n_rho
    rho_nested(1+(i-1)*n_d:i*n_d)=rho_tmp(i);
    d_CRstd = interp2(d_surf,rho_surf,CR_surf,d_vect,rho_tmp(i),'spline');
    h_d_CR_n = posteriorDistributionD(d_hat,d_vect,d_CRstd);
    d_nested(1+(i-1)*n_d:i*n_d) = samplingFromPosterior(n_d,d_vect,h_d_CR_n);
    waitbar((i)/(n_rho),f,sprintf('Please wait %2.2f%%',(i*100)/(n_rho)))
end
delete(f);
DR_nested = normcdf(d_nested);
  %% C.4 Regulatory Capital & Add on

disp('Homogeneous Portfolio, Cramer Rao')
% Computing Regulatory Capitals simulating different parameters HP CL
RCb_HP_d_CR = CapitalRequirementAlternativeHP(RR_mean,DR_sim_bCR,...
    rho_mean,CL,systematicRisk,idiosyncraticRisk);
RCb_HP_rho_CR = CapitalRequirementAlternativeHP(RR_mean,DR_mean,...
    rho_sim_bCR,CL,systematicRisk,idiosyncraticRisk);
RC_HP_allB_CR = CapitalRequirementAlternativeHP(RR_mean,DR_nested,...
    rho_nested,CL,systematicRisk,idiosyncraticRisk);

% Regulatory Capital alternative models, HP, CL
RCaltB_HP = [RCb_HP_d_CR;RCb_HP_rho_CR;RC_HP_allB_CR];
% Add on alternative models, HP, CL
addOnB_HP = RCaltB_HP./RC_HP - 1 ;

disp(table(RCaltB_HP,addOnB_HP,...
    'RowNames',{'d','rho','all'},...
    'VariableNames',{'Regulatory_Capital','Add_On'}))

%% Sobol Indices

Covariance = cov(norminv(data.DR_SG),data.RR);
mu12 = mean([norminv(data.DR_SG),data.RR]);
a = 5.1083;
b = 50.1766;

sim12 = mvnrnd(mu12,Covariance,N_sim);
sim1 = sim12(:,1);          %barrier
sim2 = sim12(:,2);          %recovery
sim3 = betarnd(a,b,N_sim,1); %correlation
parameterSim = [sim1,sim2,sim3];
%parameterSim = rand(N_sim,3);

[S1,V] = SobolInidices(parameterSim)
