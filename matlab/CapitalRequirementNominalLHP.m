function [CapitalRequirement] = CapitalRequirementNominalLHP(recoveryRate,...
    defaultRate,correlation,confidenceLevel)
% Computes capital requirement in the nominal model large homogeneous
% porftolio case
%
% @inputs:   recovery rate, 
%            probability of default, 
%            correlation
%            confidence level
% @outputs: CapitalRequirement
%

lossGivenDefault = (1-recoveryRate);
valueAtRisk = lossGivenDefault.*normcdf((norminv(defaultRate)-sqrt(correlation).*norminv(1-confidenceLevel))./sqrt(1-correlation));
expectedLoss = lossGivenDefault.*defaultRate;

CapitalRequirement = valueAtRisk - expectedLoss;

end

