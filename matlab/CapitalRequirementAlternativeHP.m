function CapitalRequirement = CapitalRequirementAlternativeHP(...
             recoveryRate,defaultRate,correlation,confidenceLevel,...
             systematicRisk,idiosyncraticRisk)
% CapitalRequrement = CAPITALREQUIREMENTALTERNATIVEHP(recoveryRate,
% defaultRate,correlation,confidenceLevel,systematicRisk,idiosyncraticRisk)
%
% Computes the Capital Requirement in the case of a Homogeneous Portfolio
% given the recovery, default probability and correlation at the required
% confidence level. Inputs can be vectors.
%
% @inputs:        - recoveryRate:      scalar or N_sim x 1 vector
%                 - defaultRate:       scalar or N_sim x 1 vector
%                 - correlation:       scalar or N_sim x 1 vector 
%                 - confidenceLevel:   scalar
%                 - systematicRisk:    N_sim x 1 vector
%                 - idiosyncraticRisk: N_sim x N_ob matrix
%
% @outputs:       - CapitalRequirement: scalar


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