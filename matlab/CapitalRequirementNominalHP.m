function [CapitalRequrement] = CapitalRequirementNominalHP(recoveryRate,...
                    defaultRate,correlation,systematicRisk,idiosyncraticRisk,confidenceLevel)
% data = readData()
% Read the data file into a struct with
% year, specultaive grade default rate, total defalut rate
% and recovery rate.
%
% @inputs: None.
%
% @outputs: data: struct with year, DG_SG, DG_All, RR
%

nObligors = size(idiosyncraticRisk,2);

exposureAtDefault   = 1/nObligors;
lossGivenDefault    = 1-recoveryRate;

systematicRiskVR      = [systematicRisk;-systematicRisk];
idiosyncraticRiskVR   = [idiosyncraticRisk; -idiosyncraticRisk];

firmValues = sqrt(correlation)*systematicRiskVR + sqrt(1-correlation)*idiosyncraticRiskVR;
defaultBarrier      = norminv(defaultRate);

numberOfDefaults    = sum(firmValues < defaultBarrier,2);

loss = exposureAtDefault * lossGivenDefault * numberOfDefaults;

valueAtRisk         = prctile(loss,confidenceLevel*100);

expectedLoss        = mean(loss);

CapitalRequrement   = valueAtRisk - expectedLoss;

end