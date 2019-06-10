function CapitalRequirement = CapitalRequirementAlternativeLHP_d(...
             recoveryRate,defaultRate,defaultBarrierMean,...
             defaultBarrierStd,correlation,confidenceLevel)
% CapitalRequrement = CAPITALREQUIREMENTALTERNATIVELHP_D(recoveryRate,
% defaultRAte,defaultBarrierMean,defaultBarrierStd,correlation,
% confidenceLevel)

% Computes the Capital Requirement in the case of a Large Homogeneous
% Portfolio given the recovery, default barrier parameters, at the required
% confidence level, using close formula in case of uncertainty of default
% barrier, assumed normally distributed.
%
% @inputs:        - recoveryRate:      scalar
%                 - defaultRate:       scalar
%                 - defaultBarriermean:scalar
%                 - defaultBarrierStd: scalar
%                 - correlation:       scalar 
%                 - confidenceLevel:   scalar
% @outputs:       - CapitalRequirement: scalar

lossGivenDefault = (1-recoveryRate);

valueAtRisk      = lossGivenDefault.*normcdf((defaultBarrierMean - ...
    sqrt(correlation+defaultBarrierStd^2).*norminv(1-confidenceLevel))./...
    sqrt(1-correlation));
expectedLoss     = lossGivenDefault.*defaultRate;

CapitalRequirement = valueAtRisk - expectedLoss;

end

