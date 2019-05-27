function [CapitalRequrement] = CapitalRequirementNominalLHP(recoveryRate,defaultRate,correlation,confidenceLevel)
% CapitalRequrement = CapitalRequirementNominalLHP(recoveryRate,defaultRate,correlation,confidenceLevel)
% Computes the Capital Requirement in the case of a Large Homogeneous
% Portfolio given the common recovery, dafalt probability and correlation
% at the required confidence level.
%
% @inputs: recoveryRate,defaultRate,correlation,confidenceLevel
%
% @outputs: CapitalRequrement
%

VaR = normcdf((norminv(defaultRate)-sqrt(correlation)*norminv(1-confidenceLevel))/...
    sqrt(1-correlation));
lossGivenDefault = (1-recoveryRate);
expectedLoss = lossGivenDefault * defaultRate;

CapitalRequrement = lossGivenDefault*VaR - expectedLoss;

end

