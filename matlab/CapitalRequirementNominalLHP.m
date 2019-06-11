function CapitalRequirement = CapitalRequirementNominalLHP(recoveryRate,...
    defaultRate,correlation,confidenceLevel)
% CapitalRequirement = CAPITALREQUIREMENTNOMINALLHP(recoveryRate,
% defaultRate,correlation,confidenceLevel)
%
% Computes the capital requirement for a Large Homogeneous Portfolio, given 
% market parameters at a given confidence level, nominal model.
%
% @inputs:        - recoveryRate:      scalar
%                 - defaultRate:       scalar
%                 - correlation:       scalar 
%                 - confidenceLevel:   scalar
% @outputs:       - CapitalRequirement: scalar

lossGivenDefault = (1-recoveryRate);
defaultBarrier   = norminv(defaultRate);

valueAtRisk      = lossGivenDefault.*normcdf((defaultBarrier-...
                        sqrt(correlation).*norminv(1-confidenceLevel))./...
                        sqrt(1-correlation));
expectedLoss     = lossGivenDefault.*defaultRate;

CapitalRequirement = valueAtRisk - expectedLoss;

end

