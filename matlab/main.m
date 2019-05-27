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
nSim = 1e6;
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

rho_vect = 0.005:0.005:0.9995;
[aa,bb] = bayesianCorrelationPosterior(0.0386, rho_vect);

h = posteriorDistributionRho(rho_hat,rho_vect,aa,bb);

cdf = cumtrapz(rho_vect,h);

u = rand(nSim,1);

[ ~, idx ] = min(abs(u - cdf),[],2);

sim_rho = rho_vect(idx)';
CapitalRequirementAlternativeLHP(recoveryMean,defaultRateSGMean,sim_rho,confidenceLevel1,nSim)/base-1

media   = trapz(rho_vect,rho_vect.*h);
var     = trapz(rho_vect,rho_vect.^2.*h)-media^2;

sigma_beta = @(a,rho) sqrt((1-rho)*rho^2/(a+rho));
b_a = @(a,rho) a*(1-rho)/rho;

fun = @(a,rho) sigma_beta(a,media)-sqrt(var);
a = fzero(@(a) fun(a,r),[0.00001,10000]);
b = b_a(a,media);

figure()
plot(rho_vect,h)
hold on
plot(rho_vect,betapdf(rho_vect, a, b))
legend('posterior','beta')
xlim([0,0.3])

%% test

RC = [];
NOb = 100:100:600;
nSim = 3e5;
for n=NOb
    n
    tmp = CapitalRequirementNominalHP(recoveryMean,defaultRateSGMean,rhoMean,nSim,n,confidenceLevel1);
    RC = [RC;tmp];
end

figure
plot(NOb,RC)
hold on
tmp = CapitalRequirementNominalLHP(recoveryMean,defaultRateSGMean,rhoMean,confidenceLevel1);
plot(NOb,ones(size(NOb))*tmp)

