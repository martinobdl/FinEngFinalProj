%% run with main.m

clear all
close all
clc

%% read Data

fileName = '../data/dati_Altman.csv';
data = readData(fileName);

%% Protfolios definitions

nObligors = 60;
confidenceLevel1 = 0.99;
confidenceLevel2 = 0.999;
rhoMean = 0.0924;
recoveryMean = mean(data.RR);
recoveryStdev = std(data.RR);
defaultRateSGMean = mean(data.DG_SG);
defaultBarrierSGMean = mean(norminv(data.DG_SG));
% defaultBarrierSGMean = -1.920994087;
defaultBarrierSGstdev = std(norminv(data.DG_SG));
nSim = 2e5;

%% A

CapitalRequirementNominalLHP(recoveryMean,defaultRateSGMean,rhoMean,confidenceLevel1)
CapitalRequirementNominalHP(recoveryMean,defaultRateSGMean,rhoMean,nSim,nObligors,confidenceLevel1)

% Basel
correlationBasel = correlationFromBasel2(defaultRateSGMean);

%% B1

% check distribution assumption

figure
ax1=subplot(2,1,1);
hist(data.RR)
ax2=subplot(2,1,2);
hist(mean(data.RR)+std(data.RR)*randn(100,1))
linkaxes([ax1,ax2],'x')
[H, pValue] = swtest(data.RR)

figure
ax1=subplot(2,1,1);
hist(norminv(data.DG_All))
ax2=subplot(2,1,2);
hist(mean(norminv(data.DG_All))+std(norminv(data.DG_All))*randn(100,1))
linkaxes([ax1,ax2],'x')
[H, pValue] = swtest(norminv(data.DG_All))

figure
ax1=subplot(2,1,1);
hist(norminv(data.DG_SG))
ax2=subplot(2,1,2);
hist(mean(norminv(data.DG_SG))+std(norminv(data.DG_SG))*randn(100,1))
linkaxes([ax1,ax2],'x')
[H, pValue] = swtest(norminv(data.DG_SG))

%% B2

base = CapitalRequirementNominalLHP(recoveryMean,defaultRateSGMean,rhoMean,confidenceLevel1);
nSim = 1e4;
a=5.1083;
b=50.1766;
rho_sim = betarnd(a,b,nSim,1);
z_sim = randn(nSim,1);

[mu_RR,sigma_RR] = normfit(data.RR);
RR_sim = randn(nSim,1)*sigma_RR+mu_RR;
d = norminv(data.DG_SG);
[mu_d,sigma_d] = normfit(d);

d_sim = randn(nSim,1)*sigma_d+mu_d;

DR_SGsim = normcdf(d_sim);

%% 

AddOn_d     = CapitalRequirementAlternativeLHP(recoveryMean,DR_SGsim,rhoMean,confidenceLevel1,nSim)/base-1
AddOn_RR    = CapitalRequirementAlternativeLHP(RR_sim,defaultRateSGMean,rhoMean,confidenceLevel1,nSim)/base-1
AddOn_rho   = CapitalRequirementAlternativeLHP(recoveryMean,defaultRateSGMean,rho_sim,confidenceLevel1,nSim)/base-1
AddOn_d_rho = CapitalRequirementAlternativeLHP(recoveryMean,DR_SGsim,rho_sim,confidenceLevel1,nSim)/base-1
AddOn_all   = CapitalRequirementAlternativeLHP(RR_sim,DR_SGsim,rho_sim,confidenceLevel1,nSim)/base-1

%%

rho_vect = linspace(0.0015,0.9985,500);
rho_hat = 0.0924;

sigma = 0.0386;

[aa,bb] = bayesianCorrelationPosterior(sigma, rho_vect);

h = posteriorDistributionRho(rho_hat,rho_vect,aa,bb);

sim_rho = samplingFromPosterior(nSim,rho_vect,h);
b = @(r) (1-r).*(r.*(1-r)-sigma.^2)./(sigma.^2);
a = @(r) r./(1-r).*(1-r).*(r.*(1-r)-sigma.^2)./(sigma.^2);
ss = MYmcmc(0.1,@(x)0,@(r)log(betapdf(rho_hat,a(r),b(r))),.2,nSim/100);

figure()
cdfplot(ss)
hold on
plot(rho_vect,cumtrapz(rho_vect,h),'-')

figure()
histogram(ss,'Normalization','pdf')
hold on
plot(rho_vect,h,'-')



figure()
cdfplot(sim_rho)
hold on
plot(rho_vect,cumtrapz(rho_vect,h),'-')

CapitalRequirementAlternativeLHP(recoveryMean,defaultRateSGMean,sim_rho,confidenceLevel1,nSim)/base-1

media   = trapz(rho_vect,rho_vect.*h);
varianza     = trapz(rho_vect,rho_vect.^2.*h)-media^2;

b = (1-media).*(media.*(1-media)-varianza)./(varianza);
a = media./(1-media).*b;

figure()
plot(rho_vect,h,'b-')
hold on
plot(rho_vect,betapdf(rho_vect, a, b),'r-')
plot(rho_vect,betapdf(rho_vect, 5.1083, 50.1766),'--')
legend('posterior','beta','freq')
%xlim([0,0.3])

%%

rho_sim = betarnd(a,b,nSim,1);
CapitalRequirementAlternativeLHP(recoveryMean,defaultRateSGMean,rho_sim,confidenceLevel1,nSim)/base-1

%% d Bayesian

n           = length(data.DG_SG);
n           = 1;
d           = norminv(data.DG_SG);
mu_post     = n*mean(d)/(n+var(d))
var_post    = var(d)/(n+var(d))

d_sim       = mu_post+sqrt(var_post)*randn(nSim,1);

DR_SGsim    = normcdf(d_sim);
AddOn_d     = CapitalRequirementAlternativeLHP(recoveryMean,DR_SGsim,rhoMean,confidenceLevel1,nSim)/base-1

%% test


