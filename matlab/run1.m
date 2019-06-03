%% 

clear all
close all
clc
rng('default')

%% INPUT

data = readData('../data/dati_Moody.csv');      % reading input data

DR_mean = mean(data.DR_SG);              % default rate speculative grade
DR_std  = std(data.DR_SG);               % mean and standard deviation

% DR_mean = mean(data.DR_AG);              % default rate all grade
% DR_std  = std(data.DR_AG);               % mean and standard deviation

d_std = std(norminv(data.DR_SG)); % ??? questo non sono ancora convinta

f_0 = @(D) integral(@(x) normpdf(x).*normcdf(D+d_std*x),-30,30)-DR_mean;
d_hat = fzero(f_0,[-5,1]);             % unbiased default barrier mean

RR_mean = mean(data.RR);                 % recovery rate
RR_std  = std(data.RR);                  % mean and standard deviation

rho_mean = 0.0924;                       % obligors correlation
rho_std  = 0.0386;                       % mean and standard deviation
rho_B_mean = correlationFromBasel2(DR_mean);  % Basel correlation

N_ob  = 60;                              % number of obligors
N_sim = 1e6;                             % number of simulations

CL1 = 0.99;                              % confidence level 99.0 %
CL2 = 0.999;                             % confidence level 99.9 %

idiosyncraticRisk = randn(N_sim,N_ob);   % simulation of idiosyncratic risk
systematicRisk    = randn(N_sim,1);      % simulation of systematic risk
        
%% A
  % i. Using mean correlation
disp('Capital requirement nominal model')
    % Large homogeneous portfolio
CR_LHP_CL1 = CapitalRequirementNominalLHP(RR_mean,DR_mean,rho_mean,CL1);
CR_LHP_CL2 = CapitalRequirementNominalLHP(RR_mean,DR_mean,rho_mean,CL2);
    % Homogeneous portfolio
CR_HP_CL1 = CapitalRequirementNominalHP(RR_mean,DR_mean,rho_mean,...
                                systematicRisk,idiosyncraticRisk,CL1);
CR_HP_CL2 = CapitalRequirementNominalHP(RR_mean,DR_mean,rho_mean,...
                                systematicRisk,idiosyncraticRisk,CL2);

disp(table([CR_LHP_CL1;CR_HP_CL1],[CR_LHP_CL2;CR_HP_CL2], ...
    'RowNames',{'LHP','HP'},'VariableNames',{'CL_099','CL_0999'}));

  % ii. Using basel correlation
disp('Capital requirement nominal model, basel correlation')
    % Large homogeneous portfolio
CR_LHP_CL1_B = CapitalRequirementNominalLHP(RR_mean,DR_mean,...
                                            rho_B_mean,CL1);
CR_LHP_CL2_B = CapitalRequirementNominalLHP(RR_mean,DR_mean,...
                                            rho_B_mean,CL2);
    % Homogeneus portfolio
CR_HP_CL1_B = CapitalRequirementNominalHP(RR_mean,DR_mean,...
            rho_B_mean,systematicRisk,idiosyncraticRisk,CL1);
CR_HP_CL2_B = CapitalRequirementNominalHP(RR_mean,DR_mean,...
            rho_B_mean,systematicRisk,idiosyncraticRisk,CL2);

disp(table([CR_LHP_CL1_B;CR_HP_CL1_B],[CR_LHP_CL2_B;CR_HP_CL2_B],...
    'RowNames',{'LHP','HP'},...
    'VariableNames',{'CL_099','CL_0999'}));

%% B.1

d = norminv(data.DR_SG);           % default barrier speculative grade
% d = norminv(data.DR_AR);      % default barrier all rates

% Shapiro wilk normality tests
disp('Default barrier normality')
[H_d,pValue_d] = swtest(d)   
disp('Recovery rate normality')
[H_RR,pValue_RR] = swtest(data.RR)

%% B.2 - B.3 simulations required

[d_mean,d_std] = normfit(d);     
d_sim = d_hat + randn(N_sim,1)*d_std;      % simulated default barriers
DR_sim = normcdf(d_sim);                    % simulated default rates

RR_sim = RR_mean + randn(N_sim,1)*RR_std;   % simulated recovery rates

[alpha,beta] = betaParameter(rho_mean,rho_std);% A & B beta parameters
rho_sim = betarnd(alpha,beta,N_sim,1);         % simulated correlation

rho_sim_B = correlationFromBasel2(DR_sim);  % simulated basel correlation

%% B.2

disp('Capital Requirement & Add On frequentist inference')

% LHP
disp('Large Homogeneous Portfolio')
% Computing capital requirements simulating different parameters LHP CL1
CR_LHP_CL1_d = CapitalRequirementAlternativeLHP(RR_mean,DR_sim,...
               rho_mean,CL1,systematicRisk);  % default barrier
CR_LHP_CL1_RR = CapitalRequirementAlternativeLHP(RR_sim,DR_mean,...
               rho_mean,CL1,systematicRisk);  % recovery rate
CR_LHP_CL1_rho = CapitalRequirementAlternativeLHP(RR_mean,DR_mean,...
               rho_sim,CL1,systematicRisk);   % correlation
CR_LHP_CL1_d_rho = CapitalRequirementAlternativeLHP(RR_mean,DR_sim,...
               rho_sim,CL1,systematicRisk);   % default barrier and corr
CR_LHP_CL1_all = CapitalRequirementAlternativeLHP(RR_sim,DR_sim,...
               rho_sim,CL1,systematicRisk);   % all parameters

% Capital Requirement alternative models, LHP, confidence level 99.0 %
CRalt_LHP_CL1 = [CR_LHP_CL1_d;CR_LHP_CL1_RR ;CR_LHP_CL1_rho;...
                 CR_LHP_CL1_d_rho;CR_LHP_CL1_all];
% Add On alternative models, LHP, confidence level 99.0 %
addOn_LHP_CL1 = CRalt_LHP_CL1./CR_LHP_CL1 - 1;

% Computing capital requirements simulating different parameters LHP CL2
CR_LHP_CL2_d = CapitalRequirementAlternativeLHP(RR_mean,DR_sim,...
               rho_mean,CL2,systematicRisk);  % default rate
CR_LHP_CL2_RR = CapitalRequirementAlternativeLHP(RR_sim,DR_mean,...
               rho_mean,CL2,systematicRisk);  % recovery rate
CR_LHP_CL2_rho = CapitalRequirementAlternativeLHP(RR_mean,DR_mean,...
               rho_sim,CL2,systematicRisk);   % correlation 
CR_LHP_CL2_d_rho = CapitalRequirementAlternativeLHP(RR_mean,DR_sim,...
               rho_sim,CL2,systematicRisk);   % default rate and corr
CR_LHP_CL2_all = CapitalRequirementAlternativeLHP(RR_sim,DR_sim,...
               rho_sim,CL2,systematicRisk);   % all parameters

% Capital Requirement alternative models, LHP, confidence level 99.9 %
CRalt_LHP_CL2 = [CR_LHP_CL2_d;CR_LHP_CL2_RR ;CR_LHP_CL2_rho;...
                 CR_LHP_CL2_d_rho;CR_LHP_CL2_all];
% Add On alternative models, LHP, confidence level 99.9 %
addOn_LHP_CL2 = CRalt_LHP_CL2./CR_LHP_CL2 - 1;

disp(table(CRalt_LHP_CL1,CRalt_LHP_CL2,addOn_LHP_CL1,addOn_LHP_CL2,...
            'RowNames',{'d','RR','rho','d_rho','All'},...
            'VariableNames',{'CR_99','CR_999','AddOn_99','AddOn_999'}))
% HP
disp('Homogeneous Portfolio')
% Computing capital requirements simulating different parameters HP CL1
CR_HP_CL1_d = CapitalRequirementAlternativeHP(RR_mean,DR_sim,rho_mean,...
        CL1,systematicRisk,idiosyncraticRisk); % default rate
CR_HP_CL1_RR = CapitalRequirementAlternativeHP(RR_sim,DR_mean,rho_mean,...
        CL1,systematicRisk,idiosyncraticRisk); % recovery rate
CR_HP_CL1_rho = CapitalRequirementAlternativeHP(RR_mean,DR_mean,rho_sim,...
        CL1,systematicRisk,idiosyncraticRisk);  % correlation
CR_HP_CL1_d_rho = CapitalRequirementAlternativeHP(RR_mean,DR_sim,rho_sim,...
        CL1,systematicRisk,idiosyncraticRisk);  % default rate and corr
CR_HP_CL1_all = CapitalRequirementAlternativeHP(RR_sim,DR_sim,rho_sim,...
        CL1,systematicRisk,idiosyncraticRisk);  % all parameters 

% Capital Requirement alternative models, HP, confidence level 99.0 %
CRalt_HP_CL1 = [CR_HP_CL1_d;CR_HP_CL1_RR ;CR_HP_CL1_rho;...
                CR_HP_CL1_d_rho;CR_HP_CL1_all];
% Add On alternative models, HP, confidence level 99.0 %
addOn_HP_CL1 = CRalt_HP_CL1./CR_HP_CL1 - 1;

% Computing capital requirements simulating different parameters HP CL2
CR_HP_CL2_d = CapitalRequirementAlternativeHP(RR_mean,DR_sim,rho_mean,...
        CL2,systematicRisk,idiosyncraticRisk);  % default rate
CR_HP_CL2_RR = CapitalRequirementAlternativeHP(RR_sim,DR_mean,rho_mean,...
        CL2,systematicRisk,idiosyncraticRisk);  % recovery rate
CR_HP_CL2_rho = CapitalRequirementAlternativeHP(RR_mean,DR_mean,rho_sim,...
        CL2,systematicRisk,idiosyncraticRisk);  % correlation
CR_HP_CL2_d_rho = CapitalRequirementAlternativeHP(RR_mean,DR_sim,rho_sim,...
        CL2,systematicRisk,idiosyncraticRisk);  % default rate and corr
CR_HP_CL2_all = CapitalRequirementAlternativeHP(RR_sim,DR_sim,rho_sim,...
        CL2,systematicRisk,idiosyncraticRisk);  % all parameters

% Capital Requirement alternative models, HP, confidence level 99.9 %
CRalt_HP_CL2 = [CR_HP_CL2_d;CR_HP_CL2_RR ;CR_HP_CL2_rho;...
                CR_HP_CL2_d_rho;CR_HP_CL2_all];
% Add On alternative models, HP, confidence level 99.9 %
addOn_HP_CL2 = CRalt_HP_CL2./CR_HP_CL2 - 1;

disp(table(CRalt_HP_CL1,CRalt_HP_CL2,addOn_HP_CL1,addOn_HP_CL2,...
            'RowNames',{'d','RR','rho','d_rho','All'},...
            'VariableNames',{'CR_99','CR_999','AddOn_99','AddOn_999'}))

%% B.3 

disp('Capital Requirement & Add On frequentist inference, basel correlation')
% LHP
disp('Large Homogeneous Portfolio')
% Computing capital requirements simulating different parameters LHP CL1
CR_LHP_CL1_B_d = CapitalRequirementAlternativeLHP(RR_mean,DR_sim,...
                rho_sim_B,CL1,systematicRisk);  % default rate
CR_LHP_CL1_B_RR = CapitalRequirementAlternativeLHP(RR_sim,DR_mean,...
                rho_B_mean,CL1,systematicRisk); % recovery rate
CR_LHP_CL1_B_d_RR = CapitalRequirementAlternativeLHP(RR_sim,DR_sim,...
                rho_sim_B,CL1,systematicRisk);  % default and recovery

% Capital Requirement alternative models, LHP, confidence level 99.0%,Basel
CRalt_LHP_CL1_B = [CR_LHP_CL1_B_d;CR_LHP_CL1_B_RR ;CR_LHP_CL1_B_d_RR];
% Add On alternative models, LHP, confidence level 99.0%, Basel
addOn_LHP_CL1_B = CRalt_LHP_CL1_B./CR_LHP_CL1_B - 1;

% Computing capital requirements simulating different parameters LHP CL2
CR_LHP_CL2_B_d = CapitalRequirementAlternativeLHP(RR_mean,DR_sim,...
                rho_sim_B,CL2,systematicRisk);  % default rate
CR_LHP_CL2_B_RR = CapitalRequirementAlternativeLHP(RR_sim,DR_mean,...
                rho_B_mean,CL2,systematicRisk); % recovery rate
CR_LHP_CL2_B_d_RR = CapitalRequirementAlternativeLHP(RR_sim,DR_sim,...
                rho_sim_B,CL2,systematicRisk);  % default and recovery

% Capital Requirement alternative models, LHP, confidence level 99.9%,Basel
CRalt_LHP_CL2_B = [CR_LHP_CL2_B_d;CR_LHP_CL2_B_RR ;CR_LHP_CL2_B_d_RR];
% Add On alternative models, LHP, confidence level 99.9%, Basel
addOn_LHP_CL2_B = CRalt_LHP_CL2_B./CR_LHP_CL2_B - 1;

disp(table(CRalt_LHP_CL1_B,CRalt_LHP_CL2_B,addOn_LHP_CL1_B,...
            addOn_LHP_CL2_B, 'RowNames',{'d','RR','d_RR'},...
            'VariableNames',{'CR_99','CR_999','AddOn_99','AddOn_999'}))

disp('Homogeneous Portfolio')
% Computing capital requirements simulating different parameters HP CL1
CR_HP_CL1_B_d = CapitalRequirementAlternativeHP(RR_mean,DR_sim,...
   rho_sim_B,CL1,systematicRisk,idiosyncraticRisk);  % default rate
CR_HP_CL1_B_RR = CapitalRequirementAlternativeHP(RR_sim,DR_mean,...
   rho_B_mean,CL1,systematicRisk,idiosyncraticRisk); % recovery rate
CR_HP_CL1_B_d_RR = CapitalRequirementAlternativeHP(RR_sim,DR_sim,...  
   rho_sim_B,CL1,systematicRisk,idiosyncraticRisk);  % default and recovery

% Capital Requirement alternative models, HP, confidence level 99.0%, Basel
CRalt_HP_CL1_B = [CR_HP_CL1_B_d;CR_HP_CL1_B_RR ;CR_HP_CL1_B_d_RR];
% Add on alternative models, HP, confidence level 99.0%, Basel
addOn_HP_CL1_B = CRalt_HP_CL1_B./CR_HP_CL1_B - 1 ;

% Computing capital requirements simulating different parameters HP CL2
CR_HP_CL2_B_d = CapitalRequirementAlternativeHP(RR_mean,DR_sim,...
    rho_sim_B,CL2,systematicRisk,idiosyncraticRisk); % default rate
CR_HP_CL2_B_RR = CapitalRequirementAlternativeHP(RR_sim,DR_mean,...
    rho_B_mean,CL2,systematicRisk,idiosyncraticRisk);% recovery rate
CR_HP_CL2_B_d_RR = CapitalRequirementAlternativeHP(RR_sim,DR_sim,...
    rho_sim_B,CL2,systematicRisk,idiosyncraticRisk); % default and recovery

% Capital Requirement alternative models, HP, confidence level 99.9%, Basel
CRalt_HP_CL2_B = [CR_HP_CL2_B_d;CR_HP_CL2_B_RR ;CR_HP_CL2_B_d_RR];
% Add on alternative models, HP, confidence level 99.9%, Basel
addOn_HP_CL2_B = CRalt_HP_CL2_B./CR_HP_CL2_B - 1 ;

disp(table(CRalt_HP_CL1_B,CRalt_HP_CL2_B,addOn_HP_CL1_B,...
    addOn_HP_CL2_B, 'RowNames',{'d','RR','d_RR'},...
    'VariableNames',{'CR_99','CR_999','AddOn_99','AddOn_999'}))

%% C.1
n = length(data.years);
  % i. d distribution
mu_post = d_hat/(1+d_std^2/n);
   %mu_post=-1.9107;% i risultati tornano se si fa sta cosa orrenda
sigma_post = d_std/sqrt(1+d_std^2);
d_sim = mu_post + sigma_post*randn(N_sim,1);

DR_sim = normcdf(d_sim);
CR = CapitalRequirementAlternativeLHP(RR_mean,DR_sim,rho_mean,CL1,systematicRisk)
CR/CR_LHP_CL1-1
  % ii. rho distribution

rho_vect = linspace(0.005,0.995,500);
[aa,bb] = betaParameter(rho_vect,rho_std);

h = posteriorDistributionRho(rho_mean,rho_vect,aa,bb);
rho_sim = samplingFromPosterior(N_sim,rho_vect,h);

%CapitalRequirementAlternativeLHP(RR_mean,defaultRateSGMean,sim_rho,confidenceLevel1,nSim)/base-1

%%

rho_vect = linspace(0.006,0.98,500);

n=N_ob;
rho_hat = 0.0924;
CR_std = CRrho(n,rho_vect);

[aa,bb] = betaParameter(rho_vect,CR_std);

h = posteriorDistributionRho(rho_hat,rho_vect,aa,bb);

rho_sim = samplingFromPosterior(N_sim,rho_vect,h);
s = @(x) sqrt((2*(1-x).^2.*(1+(n-1).*x).^2)./(n.*(n-1)));
b = @(r) (1-r).*(r.*(1-r)-s(r).^2)./(s(r).^2);
a = @(r) r./(1-r).*(1-r).*(r.*(1-r)-s(r).^2)./(s(r).^2);

%X   = mcmc(0.1,@(x)0,@(x) log(betapdf(rho_hat,a(x),b(x))),0.2,N_sim/100);

X = samplingFromPosterior(N_sim,rho_vect,h);

figure()
histogram(X,'Normalization','pdf')
hold on
plot(rho_vect,h,'-')

CapitalRequirementAlternativeHP(...
             RR_mean,DR_mean,X,CL1,...
             systematicRisk,idiosyncraticRisk)/CR_HP_CL1-1

%%

figure
plot(norminv(data.DR_AR),data.RR,'*')
hold on
plot(norminv(data.DR_SG),data.RR,'*')
xlabel('DefaultBarrier')
legend('AllGrade','SpeculativeGrade')
grid on
ylabel('Recovery')

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

[Default_S1, Recovery_S2, Correlation_S3] = SobolInidices(data.DR_SG,data.RR)

%%