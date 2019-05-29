function CapitalRequirement = CapitalRequirementAlternativeLHP(...
recoveryRate,defaultRate,correlation,confidenceLevel,systematicRisk)
% CapitalRequrement = CapitalRequirementNominalLHP(recoveryRate,...
%defaultRate,correlation,confidenceLevel)
% Computes the Capital Requirement in the case of a Large Homogeneous
% Portfolio given the common recovery, dafualt probability and correlation
% at the required confidence level.
%
% @inputs: recoveryRate,defaultRate,correlation,confidenceLevel
%
% @outputs: CapitalRequrement
%

%% simulation

defaultBarrier = norminv(defaultRate);

%% Loss

loss = (1-recoveryRate).*normcdf((defaultBarrier - sqrt(correlation).*systematicRisk)./sqrt(1-correlation));

%% VaR and CR

VaR = prctile(loss,confidenceLevel*100);
meanLoss = mean(loss);
CapitalRequirement = VaR - meanLoss;

end

