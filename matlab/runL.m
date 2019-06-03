%% 

clear all
close all
clc 
rng('default')

%% INPUT

data = readData('../data/dati_Altman.csv');      % reading input data

DR_mean = mean(data.DR_SG);              % default rate speculative grade
DR_std  = std(data.DR_SG);               % mean and standard deviation

% DR_mean = mean(data.DR_AG);              % default rate all grade
% DR_std  = std(data.DR_AG);               % mean and standard deviation

<<<<<<< HEAD
d_std = std(norminv(data.DR_SG));        % standard deviation threshold
f_0 = @(D) integral(@(x) normcdf(D+d_std.*x).*normpdf(x),-30,30)-DR_mean;
d_hat = fzero(f_0,[-2,0]);               % unbiased default barrier mean
=======
d_std = std(norminv(data.DR_SG)); % ??? questo non sono ancora convinta

f_0 = @(D) integral(@(x) normpdf(x).*normcdf(D+d_std*x),-30,30)-DR_mean;
d_hat = fzero(f_0,[-5,1]);             % unbiased default barrier mean
>>>>>>> 5fa4f2d50d63115eb42902c309439b1a96b6c360

RR_mean = mean(data.RR);                 % recovery rate
RR_std  = std(data.RR);                  % mean and standard deviation

rho_mean = 0.0924;                       % obligors correlation
rho_std  = 0.0386;                       % mean and standard deviation
rho_B = correlationFromBasel2(DR_mean);  % Basel correlation

N_ob  = 60;                              % number of obligors
N_sim = 1e6;                             % number of simulations

CL1 = 0.99;                              % confidence level 99.0 %
CL2 = 0.999;                             % confidence level 99.9 %

idiosyncraticRisk = randn(N_sim,N_ob);   % simulation of idiosyncratic risk
systematicRisk    = randn(N_sim,1);      % simulation of systematic risk
        
%% A - COMPUTE CAPITAL REQUIREMENTS IN THE NOMINAL MODEL
  %% A.i. Using mean correlation
disp('Capital requirement nominal model')

    % For Large Homogeneous Portfolio
CR_LHP_CL1 = CapitalRequirementNominalLHP(RR_mean,DR_mean,rho_mean,CL1);
CR_LHP_CL2 = CapitalRequirementNominalLHP(RR_mean,DR_mean,rho_mean,CL2);

    % For Homogeneous Portfolio
CR_HP_CL1 = CapitalRequirementNominalHP(RR_mean,DR_mean,rho_mean,CL1,N_ob);
CR_HP_CL2 = CapitalRequirementNominalHP(RR_mean,DR_mean,rho_mean,CL2,N_ob);

disp(table([CR_LHP_CL1;CR_HP_CL1],[CR_LHP_CL2;CR_HP_CL2], ...
    'RowNames',{'LHP','HP'},'VariableNames',{'CL_099','CL_0999'}));

 
 %% A.ii. Using basel correlation
disp('Capital requirement nominal model, basel correlation')

    % For Large Homogeneous Portfolio
CR_LHP_CL1_B = CapitalRequirementNominalLHP(RR_mean,DR_mean,rho_B,CL1);
CR_LHP_CL2_B = CapitalRequirementNominalLHP(RR_mean,DR_mean,rho_B,CL2);

    % For Homogeneus Portfolio
CR_HP_CL1_B = CapitalRequirementNominalHP(RR_mean,DR_mean,rho_B,CL1,N_ob);
CR_HP_CL2_B = CapitalRequirementNominalHP(RR_mean,DR_mean,rho_B,CL2,N_ob);

disp(table([CR_LHP_CL1_B;CR_HP_CL1_B],[CR_LHP_CL2_B;CR_HP_CL2_B],...
    'RowNames',{'LHP','HP'},'VariableNames',{'CL_099','CL_0999'}));
%% A bonus - Check numerical accuracy of Monte Carlo 
% 
%  % For Large Homogeneous Portfolio
% CR_LHP_CL1_MC = CapitalRequirementAlternativeLHP(RR_mean,DR_mean,rho_mean,...
%                     CL1,systematicRisk);
% CR_LHP_CL2_MC = CapitalRequirementAlternativeLHP(RR_mean,DR_mean,rho_mean,...
%                     CL2,systematicRisk);
% 
%     % For Homogeneous Portfolio
% CR_HP_CL1_MC = CapitalRequirementAlternativeHP(RR_mean,DR_mean,rho_mean,...
%                     CL1,systematicRisk,idiosyncraticRisk);
% CR_HP_CL2_MC = CapitalRequirementAlternativeHP(RR_mean,DR_mean,rho_mean,...
%                     CL2,systematicRisk,idiosyncraticRisk);
% 
% disp(table([CR_LHP_CL1_MC;CR_HP_CL1_MC],[CR_LHP_CL2_MC;CR_HP_CL2_MC], ...
%     'RowNames',{'LHP','HP'},'VariableNames',{'CL_099','CL_0999'}));
%% B - FREQUENTISTIC INFERENCE
  %% B.1 Verify hypotesis on data

d = norminv(data.DR_SG);           % default barrier speculative grade
% d = norminv(data.DR_AR);           % default barrier all rates

% Shapiro wilk normality tests

disp('Default barrier normality, Shapiro-Wilk test')
[H_d,pValue_d] = swtest(d)  

<<<<<<< HEAD
disp('Recovery rate normality, Shapiro-Wilk test')
[H_RR,pValue_RR] = swtest(data.RR)

  %% B.2 - B.3 Simulations required
   
d_sim = d_hat + randn(N_sim,1)*d_std;      % simulated gaussian threshold 
=======
[d_mean,d_std] = normfit(d);     
d_sim = d_hat + randn(N_sim,1)*d_std;      % simulated default barriers
>>>>>>> 5fa4f2d50d63115eb42902c309439b1a96b6c360
DR_sim = normcdf(d_sim);                    % simulated default rates

RR_sim = RR_mean + randn(N_sim,1)*RR_std;   % simulated gaussian recovery 

[alpha,beta] = betaParameter(rho_mean,rho_std);% A & B beta parameters
rho_sim = betarnd(alpha,beta,N_sim,1);      % simulated beta correlation

rho_sim_B = correlationFromBasel2(DR_sim);  % simulated basel correlation

%% B.2 Capital Requirement & Add On using correlation from data

disp('Capital Requirement & Add On: frequentist inference')

% For Large Homogeneous Portfolio
disp('Large Homogeneous Portfolio')

% Computing capital requirements simulating different parameters LHP CL1
CR_LHP_CL1_d = CapitalRequirementAlternativeLHP(RR_mean,DR_sim,...
               rho_mean,CL1,systematicRisk); % default rate
CR_LHP_CL1_RR = CapitalRequirementAlternativeLHP(RR_sim,DR_mean,...
               rho_mean,CL1,systematicRisk); % recovery rate
CR_LHP_CL1_rho = CapitalRequirementAlternativeLHP(RR_mean,DR_mean,...
               rho_sim,CL1,systematicRisk);  % correlation
CR_LHP_CL1_d_rho = CapitalRequirementAlternativeLHP(RR_mean,DR_sim,...
               rho_sim,CL1,systematicRisk);  % default rate and correlation
CR_LHP_CL1_all = CapitalRequirementAlternativeLHP(RR_sim,DR_sim,...
               rho_sim,CL1,systematicRisk);  % all parameters

% Capital Requirement alternative models, LHP, CL 99.0 %
CRalt_LHP_CL1 = [ CR_LHP_CL1_d ; CR_LHP_CL1_RR ; CR_LHP_CL1_rho ; ...
                  CR_LHP_CL1_d_rho ; CR_LHP_CL1_all ];
% Add On alternative models, LHP, CL 99.0 %
addOn_LHP_CL1 = CRalt_LHP_CL1./CR_LHP_CL1 - 1;

% Computing capital requirements simulating different parameters LHP CL2
CR_LHP_CL2_d = CapitalRequirementAlternativeLHP(RR_mean,DR_sim,...
               rho_mean,CL2,systematicRisk); % default rate
CR_LHP_CL2_RR = CapitalRequirementAlternativeLHP(RR_sim,DR_mean,...
               rho_mean,CL2,systematicRisk); % recovery rate
CR_LHP_CL2_rho = CapitalRequirementAlternativeLHP(RR_mean,DR_mean,...
               rho_sim,CL2,systematicRisk);  % correlation 
CR_LHP_CL2_d_rho = CapitalRequirementAlternativeLHP(RR_mean,DR_sim,...
               rho_sim,CL2,systematicRisk);  % default rate and correlation
CR_LHP_CL2_all = CapitalRequirementAlternativeLHP(RR_sim,DR_sim,...
               rho_sim,CL2,systematicRisk);  % all parameters

% Capital Requirement alternative models, LHP, CL 99.9 %
CRalt_LHP_CL2 = [ CR_LHP_CL2_d ; CR_LHP_CL2_RR ; CR_LHP_CL2_rho ; ...
                  CR_LHP_CL2_d_rho ; CR_LHP_CL2_all ];
% Add On alternative models, LHP, CL 99.9 %
addOn_LHP_CL2 = CRalt_LHP_CL2./CR_LHP_CL2 - 1;

disp(table(CRalt_LHP_CL1,CRalt_LHP_CL2,addOn_LHP_CL1,addOn_LHP_CL2,...
            'RowNames',{'d','RR','rho','d_rho','All'},...
            'VariableNames',{'CR_99','CR_999','AddOn_99','AddOn_999'}))

% For Homogeneous Portfolio
disp('Homogeneous Portfolio')

% Computing capital requirements simulating different parameters HP CL1
CR_HP_CL1_d = CapitalRequirementAlternativeHP(RR_mean,DR_sim,rho_mean,...
      CL1,systematicRisk,idiosyncraticRisk); % default rate
CR_HP_CL1_RR = CapitalRequirementAlternativeHP(RR_sim,DR_mean,rho_mean,...
      CL1,systematicRisk,idiosyncraticRisk); % recovery rate
CR_HP_CL1_rho = CapitalRequirementAlternativeHP(RR_mean,DR_mean,rho_sim,...
      CL1,systematicRisk,idiosyncraticRisk); % correlation
CR_HP_CL1_d_rho = CapitalRequirementAlternativeHP(RR_mean,DR_sim,rho_sim,...
      CL1,systematicRisk,idiosyncraticRisk); % default rate and correlation
CR_HP_CL1_all = CapitalRequirementAlternativeHP(RR_sim,DR_sim,rho_sim,...
      CL1,systematicRisk,idiosyncraticRisk); % all parameters 

% Capital Requirement alternative models, HP, CL 99.0 %
CRalt_HP_CL1 = [CR_HP_CL1_d;CR_HP_CL1_RR ;CR_HP_CL1_rho;...
                CR_HP_CL1_d_rho;CR_HP_CL1_all];
% Add On alternative models, HP, CL 99.0 %
addOn_HP_CL1 = CRalt_HP_CL1./CR_HP_CL1 - 1;

% Computing capital requirements simulating different parameters HP CL2
CR_HP_CL2_d = CapitalRequirementAlternativeHP(RR_mean,DR_sim,rho_mean,...
      CL2,systematicRisk,idiosyncraticRisk); % default rate
CR_HP_CL2_RR = CapitalRequirementAlternativeHP(RR_sim,DR_mean,rho_mean,...
      CL2,systematicRisk,idiosyncraticRisk); % recovery rate
CR_HP_CL2_rho = CapitalRequirementAlternativeHP(RR_mean,DR_mean,rho_sim,...
      CL2,systematicRisk,idiosyncraticRisk); % correlation
CR_HP_CL2_d_rho = CapitalRequirementAlternativeHP(RR_mean,DR_sim,rho_sim,...
      CL2,systematicRisk,idiosyncraticRisk); % default rate and correlation
CR_HP_CL2_all = CapitalRequirementAlternativeHP(RR_sim,DR_sim,rho_sim,...
      CL2,systematicRisk,idiosyncraticRisk); % all parameters

% Capital Requirement alternative models, HP, CL 99.9 %
CRalt_HP_CL2 = [CR_HP_CL2_d;CR_HP_CL2_RR;CR_HP_CL2_rho;...
                CR_HP_CL2_d_rho;CR_HP_CL2_all];
% Add On alternative models, HP, CL 99.9 %
addOn_HP_CL2 = CRalt_HP_CL2./CR_HP_CL2 - 1;

disp(table(CRalt_HP_CL1,CRalt_HP_CL2,addOn_HP_CL1,addOn_HP_CL2,...
            'RowNames',{'d','RR','rho','d_rho','All'},...
            'VariableNames',{'CR_99','CR_999','AddOn_99','AddOn_999'}))

%% B.3 

disp('Capital Requirement & Add On frequentist inference, basel correlation')

% For Large Homogeneous Portfolio
disp('Large Homogeneous Portfolio')

% Computing capital requirements simulating different parameters LHP CL1
CR_LHP_CL1_B_d = CapitalRequirementAlternativeLHP(RR_mean,DR_sim,...
                rho_sim_B,CL1,systematicRisk);  % default rate
CR_LHP_CL1_B_RR = CapitalRequirementAlternativeLHP(RR_sim,DR_mean,...
                rho_B,CL1,systematicRisk);      % recovery rate
CR_LHP_CL1_B_d_RR = CapitalRequirementAlternativeLHP(RR_sim,DR_sim,...
                rho_sim_B,CL1,systematicRisk);  % default and recovery

% Capital Requirement alternative models, LHP, CL 99.0%, Basel
CRalt_LHP_CL1_B = [ CR_LHP_CL1_B_d ; CR_LHP_CL1_B_RR ; CR_LHP_CL1_B_d_RR ];
% Add On alternative models, LHP, confidence level 99.0%, Basel
addOn_LHP_CL1_B = CRalt_LHP_CL1_B./CR_LHP_CL1_B - 1;

% Computing capital requirements simulating different parameters LHP CL2
CR_LHP_CL2_B_d = CapitalRequirementAlternativeLHP(RR_mean,DR_sim,...
                rho_sim_B,CL2,systematicRisk);  % default rate
CR_LHP_CL2_B_RR = CapitalRequirementAlternativeLHP(RR_sim,DR_mean,...
                rho_B,CL2,systematicRisk);      % recovery rate
CR_LHP_CL2_B_d_RR = CapitalRequirementAlternativeLHP(RR_sim,DR_sim,...
                rho_sim_B,CL2,systematicRisk);  % default and recovery

% Capital Requirement alternative models, LHP, CL 99.9%, Basel
CRalt_LHP_CL2_B = [ CR_LHP_CL2_B_d ; CR_LHP_CL2_B_RR ; CR_LHP_CL2_B_d_RR ];
% Add On alternative models, LHP, confidence level 99.9%, Basel
addOn_LHP_CL2_B = CRalt_LHP_CL2_B./CR_LHP_CL2_B - 1;

disp(table(CRalt_LHP_CL1_B,CRalt_LHP_CL2_B,addOn_LHP_CL1_B,addOn_LHP_CL2_B,...
    'RowNames',{'d','RR','d_RR'},...
    'VariableNames',{'CR_99','CR_999','AddOn_99','AddOn_999'}))

% For Homogeneous Portfolio
disp('Homogeneous Portfolio')

% Computing capital requirements simulating different parameters HP CL1
CR_HP_CL1_B_d = CapitalRequirementAlternativeHP(RR_mean,DR_sim,...
   rho_sim_B,CL1,systematicRisk,idiosyncraticRisk);  % default rate
CR_HP_CL1_B_RR = CapitalRequirementAlternativeHP(RR_sim,DR_mean,...
   rho_B,CL1,systematicRisk,idiosyncraticRisk);      % recovery rate
CR_HP_CL1_B_d_RR = CapitalRequirementAlternativeHP(RR_sim,DR_sim,...  
   rho_sim_B,CL1,systematicRisk,idiosyncraticRisk);  % default and recovery

% Capital Requirement alternative models, HP, CL 99.0%, Basel
CRalt_HP_CL1_B = [CR_HP_CL1_B_d;CR_HP_CL1_B_RR ;CR_HP_CL1_B_d_RR];
% Add on alternative models, HP, CL 99.0%, Basel
addOn_HP_CL1_B = CRalt_HP_CL1_B./CR_HP_CL1_B - 1 ;

% Computing capital requirements simulating different parameters HP CL2
CR_HP_CL2_B_d = CapitalRequirementAlternativeHP(RR_mean,DR_sim,...
    rho_sim_B,CL2,systematicRisk,idiosyncraticRisk); % default rate
CR_HP_CL2_B_RR = CapitalRequirementAlternativeHP(RR_sim,DR_mean,...
    rho_B,CL2,systematicRisk,idiosyncraticRisk);     % recovery rate
CR_HP_CL2_B_d_RR = CapitalRequirementAlternativeHP(RR_sim,DR_sim,...
    rho_sim_B,CL2,systematicRisk,idiosyncraticRisk); % default and recovery

% Capital Requirement alternative models, HP, confidence level 99.9%, Basel
CRalt_HP_CL2_B = [CR_HP_CL2_B_d;CR_HP_CL2_B_RR ;CR_HP_CL2_B_d_RR];
% Add on alternative models, HP, confidence level 99.9%, Basel
addOn_HP_CL2_B = CRalt_HP_CL2_B./CR_HP_CL2_B - 1 ;

disp(table(CRalt_HP_CL1_B,CRalt_HP_CL2_B,addOn_HP_CL1_B,...
    addOn_HP_CL2_B, 'RowNames',{'d','RR','d_RR'},...
    'VariableNames',{'CR_99','CR_999','AddOn_99','AddOn_999'}))

%% B bonus check numerical accuracy using close formula
% 
% % Test numerical accuracy using close formula in case of uncertainty on d
% disp('Large Homogeneous Portfolio, d-close formula')
% CR_LHP_CL1_d_close = CapitalRequirementAlternativeLHP_D(RR_mean,DR_mean,...
%     d_hat,d_std,rho_mean,CL1);
% CR_LHP_CL2_d_close = CapitalRequirementAlternativeLHP_D(RR_mean,DR_mean,...
%     d_hat,d_std,rho_mean,CL2);
% 
% disp(table(CR_LHP_CL1_d_close,CR_LHP_CL2_d_close,...
%     'RowNames',{'d_close'},'VariableNames',{'CR_99','CR_999'}))
%% C - BAYESIAN INFERENCE
  %% C.1 Posterior distributions fixing variance equal to empiric
  
     % i. Of threshold d
     
n = length(data.years);                              % number of data 

mu_post = d_hat/(1+d_std^2/n);             % posterion mean
sigma_post = d_std/sqrt(1+d_std^2);        % posterior standard deviation
fprintf('The posterior distribution of d is gaussian\n')
fprintf('with mean %f and standard deviation %f\n', mu_post,sigma_post)

d_sim = mu_post + randn(N_sim,1)*sigma_post; % simulated d from posterior
DR_sim = normcdf(d_sim);                     % simulated default rate

     % ii. Of correlation rho
     
rho_vect = linspace(0.005,0.995,1980);       % equispaced vector 
[aa,bb] = betaParameter(rho_vect,rho_std);   % A and B beta parameters that
                                        % satisfy condition on mean and std  
h = posteriorDistributionRho(rho_mean,rho_vect,aa,bb);
mean_h = trapz(rho_vect,rho_vect.*h);
std_h = sqrt(trapz(rho_vect,rho_vect.^2.*h)-mean_h.^2);

rho_sim = samplingFromPosterior(N_sim,rho_vect,h);
  %% C.2 Capital Requirement & Add On using correlation from data
systematicRisk = randn(N_sim,1);
idiosyncraticRisk = randn(N_sim,N_ob);
disp('Capitar Requirement & Add On, Bayesian inference')

% For Large Homogeneous Portfolio
disp('Large Homogeneous Portfolio')

% Computing capital requirements simulating different parameters LHP CL1
CR_LHP_CL1_dB = CapitalRequirementAlternativeLHP(RR_mean,DR_sim,...
    rho_mean,CL1,systematicRisk);
CR_LHP_CL1_rhoB = CapitalRequirementAlternativeLHP(RR_mean,DR_mean,...
    rho_sim,CL1,systematicRisk);
CR_LHP_CL1_allB = CapitalRequirementAlternativeLHP(RR_mean,DR_sim,...
    rho_sim,CL1,systematicRisk);

% Capital Requirement alternative models, LHP, CL 99.0%
CRaltB_LHP_CL1 = [ CR_LHP_CL1_dB ; CR_LHP_CL1_rhoB ; CR_LHP_CL1_allB ];
% Add on alternative models, LHP, CL 99.0%
addOnB_LHP_CL1 = CRaltB_LHP_CL1./CR_LHP_CL1 - 1 ;

% Computing capital requirements simulating different parameters LHP CL2
CR_LHP_CL2_dB = CapitalRequirementAlternativeLHP(RR_mean,DR_sim,...
    rho_mean,CL2,systematicRisk);
CR_LHP_CL2_rhoB = CapitalRequirementAlternativeLHP(RR_mean,DR_mean,...
    rho_sim,CL2,systematicRisk);
CR_LHP_CL2_allB = CapitalRequirementAlternativeLHP(RR_mean,DR_sim,...
    rho_sim,CL2,systematicRisk);

% Capital Requirement alternative models, LHP, CL 99.9%
CRaltB_LHP_CL2 = [ CR_LHP_CL2_dB ; CR_LHP_CL2_rhoB ; CR_LHP_CL2_allB ];
% Add on alternative models, LHP, CL 99.9%
addOnB_LHP_CL2 = CRaltB_LHP_CL2./CR_LHP_CL2 - 1 ;

disp(table(CRaltB_LHP_CL1,CRaltB_LHP_CL2,addOnB_LHP_CL1,addOnB_LHP_CL2,...
    'RowNames',{'d','rho','all'},...
    'VariableNames',{'CR_99','CR_999','AddOn_99','AddOn_999'}))

% For Homogeneous Portfolio
disp('Homogeneous Portfolio')

% Computing capital requirements simulating different parameters HP CL1
CR_HP_CL1_dB = CapitalRequirementAlternativeHP(RR_mean,DR_sim,...
    rho_mean,CL1,systematicRisk,idiosyncraticRisk);
CR_HP_CL1_rhoB = CapitalRequirementAlternativeHP(RR_mean,DR_mean,...
    rho_sim,CL1,systematicRisk,idiosyncraticRisk);
CR_HP_CL1_allB = CapitalRequirementAlternativeHP(RR_mean,DR_sim,...
    rho_sim,CL1,systematicRisk,idiosyncraticRisk);

% Capital Requirement alternative models, HP, CL 99.0%
CRaltB_HP_CL1 = [ CR_HP_CL1_dB ; CR_HP_CL1_rhoB ; CR_HP_CL1_allB ];
% Add on alternative models, HP, CL 99.9%
addOnB_HP_CL1 = CRaltB_HP_CL1./CR_HP_CL1 - 1 ;

% Computing capital requirements simulating different parameters HP CL2
CR_HP_CL2_dB = CapitalRequirementAlternativeHP(RR_mean,DR_sim,...
    rho_mean,CL2,systematicRisk,idiosyncraticRisk);
CR_HP_CL2_rhoB = CapitalRequirementAlternativeHP(RR_mean,DR_mean,...
    rho_sim,CL2,systematicRisk,idiosyncraticRisk);
CR_HP_CL2_allB = CapitalRequirementAlternativeHP(RR_mean,DR_sim,...
    rho_sim,CL2,systematicRisk,idiosyncraticRisk);

% Capital Requirement alternative models, HP, CL 99.9%
CRaltB_HP_CL2 = [ CR_HP_CL2_dB ; CR_HP_CL2_rhoB ; CR_HP_CL2_allB ];
% Add on alternative models, HP, CL 99.9%
addOnB_HP_CL2 = CRaltB_HP_CL2./CR_HP_CL2 - 1 ;

disp(table(CRaltB_HP_CL1,CRaltB_HP_CL2,addOnB_HP_CL1,addOnB_HP_CL2,...
    'RowNames',{'d','rho','all'},...
    'VariableNames',{'CR_99','CR_999','AddOn_99','AddOn_999'}))
  %% C.3 Posterior distributions fixing variance equal to Cramer Rao
  
     % i. Of threshold d

% d_vect = linspace(-3.5,3.5,36); 
% rho_vect = linspace(0.005,0.995,31);
% d_CRstd = CRd(d_vect,rho_vect, N_Ob); % cramer rao surface
% mu_post = 0;

     % ii. Of correlation rho
     
rho_vect = linspace(0.005,0.995,1980);       % equispaced vector 
rho_CRstd = CRrho(10891,rho_vect,30);        % Cramer Rao stdev
[aa,bb] = betaParameter(rho_vect,rho_CRstd); % A and B beta parameters that
                                        % satisfy condition on mean and std  
h = posteriorDistributionRho(rho_mean,rho_vect,aa,bb);
mean_h = trapz(rho_vect,rho_vect.*h);
std_h = sqrt(trapz(rho_vect,rho_vect.^2.*h)-mean_h.^2);

<<<<<<< HEAD
rho_sim = samplingFromPosterior(N_sim,rho_vect,h);

  %% C.4 Capital Requirement & Add on
  
  % For Homogeneous Portfolio
disp('Homogeneous Portfolio')
=======
n=N_ob;
rho_hat = 0.0924;
CR_std = CRrho(n,rho_vect);
>>>>>>> 5fa4f2d50d63115eb42902c309439b1a96b6c360

% Computing capital requirements simulating different parameters HP CL1
CR_HP_CL1_dB = CapitalRequirementAlternativeHP(RR_mean,DR_sim,...
    rho_mean,CL1,systematicRisk,idiosyncraticRisk);
CR_HP_CL1_rhoB = CapitalRequirementAlternativeHP(RR_mean,DR_mean,...
    rho_sim,CL1,systematicRisk,idiosyncraticRisk);
CR_HP_CL1_allB = CapitalRequirementAlternativeHP(RR_mean,DR_sim,...
    rho_sim,CL1,systematicRisk,idiosyncraticRisk);

% Capital Requirement alternative models, HP, CL 99.0%
CRaltB_HP_CL1 = [ CR_HP_CL1_dB ; CR_HP_CL1_rhoB ; CR_HP_CL1_allB ];
% Add on alternative models, HP, CL 99.9%
addOnB_HP_CL1 = CRaltB_HP_CL1./CR_HP_CL1 - 1 ;

% Computing capital requirements simulating different parameters HP CL2
CR_HP_CL2_dB = CapitalRequirementAlternativeHP(RR_mean,DR_sim,...
    rho_mean,CL2,systematicRisk,idiosyncraticRisk);
CR_HP_CL2_rhoB = CapitalRequirementAlternativeHP(RR_mean,DR_mean,...
    rho_sim,CL2,systematicRisk,idiosyncraticRisk);
CR_HP_CL2_allB = CapitalRequirementAlternativeHP(RR_mean,DR_sim,...
    rho_sim,CL2,systematicRisk,idiosyncraticRisk);

% Capital Requirement alternative models, HP, CL 99.9%
CRaltB_HP_CL2 = [ CR_HP_CL2_dB ; CR_HP_CL2_rhoB ; CR_HP_CL2_allB ];
% Add on alternative models, HP, CL 99.9%
addOnB_HP_CL2 = CRaltB_HP_CL2./CR_HP_CL2 - 1 ;

disp(table(CRaltB_HP_CL1,CRaltB_HP_CL2,addOnB_HP_CL1,addOnB_HP_CL2,...
    'RowNames',{'d','rho','all'},...
    'VariableNames',{'CR_99','CR_999','AddOn_99','AddOn_999'}))
%%

<<<<<<< HEAD
rho_hat = 0.0924;
s = @(x) sqrt((2*(1-x).^2.*(1+(n-1).*x).^2)./(30*n.*(n-1)));
=======
rho_sim = samplingFromPosterior(N_sim,rho_vect,h);
s = @(x) sqrt((2*(1-x).^2.*(1+(n-1).*x).^2)./(n.*(n-1)));
>>>>>>> 5fa4f2d50d63115eb42902c309439b1a96b6c360
b = @(r) (1-r).*(r.*(1-r)-s(r).^2)./(s(r).^2);
a = @(r) r./(1-r).*(1-r).*(r.*(1-r)-s(r).^2)./(s(r).^2);

%X   = mcmc(0.1,@(x)0,@(x) log(betapdf(rho_hat,a(x),b(x))),0.2,N_sim/100);

% X = samplingFromPosterior(N_sim,rho_vect,h);
%[A,B]=betaParameter(0.0891,0.0383);
% figure()
% histogram(X,'Normalization','pdf')
% hold on
plot(rho_vect,h,'-')

CapitalRequirementAlternativeHP(...
             RR_mean,DR_mean,X,CL1,...
             randn(N_sim,1),randn(N_sim,N_ob))/CR_HP_CL1-1

%%

figure
plot(norminv(data.DR_AR),data.RR,'*')
xlabel('d')
ylabel('\pi')

<<<<<<< HEAD
[B,BINT,R] = regress(norminv(data.DR_AR),data.RR)
=======
[r,m,b] = regression(norminv(data.DR_AR)',data.RR','one');
plotregression(norminv(data.DR_AR),data.RR)

%%

d_sim   = mcmc(0.5,@(x) log(normpdf(x)),@(x) sum(log(normpdf(d,x,d_std))),0.2,N_sim/100);

[mu_post,std_post]=normfit(d_sim);

DR_Sim = normcdf(d_sim);


CapitalRequirementAlternativeLHP(RR_mean,DR_Sim,rho_mean,CL1,randn(N_sim/100,1))/CR_LHP_CL1-1

%% Sobol Indices

% INPUT PARAMETERS ARE NOT INDEPENDENT
% FIX CORRELATION (NO DATA FOR CORRELATION)

[S1,S2,S3] = SobolInidices(data.DR_SG,data.RR)
>>>>>>> 5fa4f2d50d63115eb42902c309439b1a96b6c360

