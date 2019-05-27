function CapitalRequrement = CapitalRequirementNominalHP(recoveryRate,defaultRate,correlation,nSim,nObligors,confidenceLevel)
% CapitalRequrement = CapitalRequirementNominalHP(recoveryRate,defaultRate,correlation,nSim,nObligors,confidenceLevel)
% Computes the Capital Requirement in the case of a Homogeneous
% Portfolio of nObligors given the common recovery, dafalt probability and correlation
% at the required confidence level, with a MC tecnique for nSim simualtions
%
% @inputs: recoveryRate,defaultRate,correlation,nSim,nObligors,confidenceLevel
%
% @outputs: CapitalRequrement
%

rng(1)

exposureAtDefault   = 1/nObligors;
lossGivenDefault    = 1-recoveryRate;

systematicRisk      = randn(nSim,1);
idiosyncraticRisk   = randn(nSim,nObligors);

% AV
systematicRisk      = [systematicRisk;-systematicRisk];
idiosyncraticRisk   = [idiosyncraticRisk; -idiosyncraticRisk];

firmValue = sqrt(correlation)*systematicRisk + sqrt(1-correlation)*idiosyncraticRisk;
defaultBarrier      = norminv(defaultRate);

numberOfDefaults    = sum(firmValue < defaultBarrier,2);

loss = exposureAtDefault * lossGivenDefault * numberOfDefaults;

VaR                 = prctile(loss,confidenceLevel*100);

expectedLoss        = (1-recoveryRate)*defaultRate;

CapitalRequrement   = VaR - expectedLoss;

end