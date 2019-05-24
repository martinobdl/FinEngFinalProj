%% run with main.m

clear all
close all
clc

%% read Data

fileName = '../data/dati_Altman.csv';
data = readData(fileName);

%% Protfolios definitions

nObligors = 50;
confidenceLevel1 = 0.99;
confidenceLevel2 = 0.999;
rhoMean = 0.0924;
recoveryMean = mean(data.RR);
defaultRateSGMean = mean(data.DG_SG);
nSim = 2e5;

%% A

correlationBasel = correlationFromBasel2(defaultRateSGMean);

CapitalRequirementNominalLHP(recoveryMean,defaultRateSGMean,rhoMean,confidenceLevel1)
CapitalRequirementNominalLHP(recoveryMean,defaultRateSGMean,rhoMean,confidenceLevel2)
CapitalRequirementNominalLHP(recoveryMean,defaultRateSGMean,correlationBasel,confidenceLevel1)
CapitalRequirementNominalLHP(recoveryMean,defaultRateSGMean,correlationBasel,confidenceLevel2)

%CapitalRequirementNominalHP(recoveryMean,defaultRateSGMean,rhoMean,nSim,nObligors,confidenceLevel1)
%CapitalRequirementNominalHP(recoveryMean,defaultRateSGMean,rhoMean,nSim,nObligors,confidenceLevel2)
%CapitalRequirementNominalHP(recoveryMean,defaultRateSGMean,correlationBasel,nSim,nObligors,confidenceLevel1)
%CapitalRequirementNominalHP(recoveryMean,defaultRateSGMean,correlationBasel,nSim,nObligors,confidenceLevel2)

%% B

%% test

NOb = 10:10:100;
RC = [];

for n=NOb
    tmp = CapitalRequirementNominalHP(recoveryMean,defaultRateSGMean,rhoMean,nSim,n,confidenceLevel1);
    RC = [RC;tmp];
end

figure
plot(NOb,RC)
hold on
tmp = CapitalRequirementNominalLHP(recoveryMean,defaultRateSGMean,rhoMean,confidenceLevel1);
plot(NOb,ones(size(NOb))*tmp)

