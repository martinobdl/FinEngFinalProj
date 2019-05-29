function CapitalRequirement = CapitalRequirementAlternativeHP(...
             recoveryRate,defaultRate,correlation,confidenceLevel,...
             systematicRisk,idiosyncraticRisk)
% CapitalRequrement = CapitalRequirementAlternativeHP(recoveryRate,
% defaultRate,correlation,confidenceLevel,systematicRisk,idiosyncraticRisk)
% Computes the Capital Requirement in the case of a Large Homogeneous
% Portfolio given the recovery, default probability and correlation
% at the required confidence level. Inputs can be simulated vectors.
%
% @inputs: recoveryRate,defaultRate,correlation,confidenceLevel,
%          systematicRisk, idiosyncraticRisk.
%
% @outputs: CapitalRequrement.
%

nObligors          = size(idiosyncraticRisk,2);
exposureAtDefault  = 1/nObligors;
lossGivenDefault   = 1-recoveryRate;

defaultBarrier     = norminv(defaultRate);
firmValues         = sqrt(correlation).*systematicRisk + ...
                     sqrt(1-correlation).*idiosyncraticRisk;
numberOfDefaults   = sum(firmValues < defaultBarrier,2);
loss               = exposureAtDefault.*lossGivenDefault.*numberOfDefaults;

valueAtRisk        = prctile(loss,confidenceLevel*100);
expectedLoss       = mean(loss);
CapitalRequirement = valueAtRisk - expectedLoss;

end