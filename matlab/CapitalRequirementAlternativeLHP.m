function CapitalRequirement = CapitalRequirementAlternativeLHP(recoveryRate,defaultRate,correlation,confidenceLevel,nSim)
% CapitalRequrement = CapitalRequirementNominalLHP(recoveryRate,defaultRate,correlation,confidenceLevel)
% Computes the Capital Requirement in the case of a Large Homogeneous
% Portfolio given the common recovery, dafalt probability and correlation
% at the required confidence level.
%
% @inputs: recoveryRate,defaultRate,correlation,confidenceLevel
%
% @outputs: CapitalRequrement
%

%% simulation

commonRiskFactor = randn(nSim,1);
defaultBarrier = norminv(defaultRate);

%% Loss

loss = (1-recoveryRate).*normcdf((defaultBarrier - sqrt(correlation).*commonRiskFactor)./sqrt(1-correlation));

%% VaR and CR

VaR = prctile(loss,confidenceLevel*100);
meanLoss = mean(loss);
CapitalRequirement = VaR - meanLoss;

end

