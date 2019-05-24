function [CapitalRequrement] = CapitalRequirementNominalHP(recoveryRate,defaultRate,correlation,nSim,nObligors,confidenceLevel)
% data = readData()
% Read the data file into a struct with
% year, specultaive grade default rate, total defalut rate
% and recovery rate.
%
% @inputs: None.
%
% @outputs: data: struct with year, DG_SG, DG_All, RR
%

rng(1)

exposureAtDefault   = 1/nObligors;
lossGivenDefault    = 1-recoveryRate;

systematicRisk      = randn(nSim,1);
idiosyncraticRisk   = randn(nSim,nObligors);

systematicRisk      = [systematicRisk;-systematicRisk];
idiosyncraticRisk   = [idiosyncraticRisk; -idiosyncraticRisk];

firmValues = sqrt(correlation)*systematicRisk + sqrt(1-correlation)*idiosyncraticRisk;
defaultBarrier      = norminv(defaultRate);

numberOfDefaults    = sum(firmValues < defaultBarrier,2);

loss = exposureAtDefault * lossGivenDefault * numberOfDefaults;

VaR                 = prctile(loss,confidenceLevel*100);

expectedLoss        = mean(loss);

CapitalRequrement   = VaR - expectedLoss;

end