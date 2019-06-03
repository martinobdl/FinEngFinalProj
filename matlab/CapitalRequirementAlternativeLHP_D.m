function CapitalRequirement = CapitalRequirementAlternativeLHP_D(...
    recoveryRate,defaultRate,defaultBarrierMean,defaultBarrierStd,correlation,...
    confidenceLevel)
% CapitalRequrement = CapitalRequirementNominalLHP(recoveryRate,...
% defaultBarrierMean,defaultBarrierStd,correlation,confidenceLevel)
% Computes the Capital Requirement in the case of a Large Homogeneous
% Portfolio given the recovery, default barrier parameters, at the required
% confidence level, using close formula in case of uncertainty of default
% barrier, assumed normal distributed
%
% @inputs: recoveryRate,defaultBarrierMean,defaultBarrierStd,correlation,
%           confidenceLevel
%
% @outputs: CapitalRequrement
%

lossGivenDefault = (1-recoveryRate);
valueAtRisk = lossGivenDefault.*normcdf((defaultBarrierMean - ...
    sqrt(correlation+defaultBarrierStd^2).*norminv(1-confidenceLevel))./...
    sqrt(1-correlation));
expectedLoss = lossGivenDefault.*defaultRate;

CapitalRequirement = valueAtRisk - expectedLoss;

end

