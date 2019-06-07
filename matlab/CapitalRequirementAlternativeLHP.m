function CapitalRequirement = CapitalRequirementAlternativeLHP(...
             recoveryRate,defaultRate,correlation,confidenceLevel,...
             systematicRisk)
% CapitalRequrement = CAPITALREQUIREMENTALTERNATIVELHP(recoveryRate,
% defaultRate,correlation,confidenceLevel,systematicRisk)
%
% Computes the Capital Requirement in the case of a Large Homogeneous 
% Portfolio given the recovery, default probability and correlation at the
% required confidence level. Inputs can be vectors.
%
% @inputs:        - recoveryRate:      scalar or N_sim x 1 vector
%                 - defaultRate:       scalar or N_sim x 1 vector
%                 - correlation:       scalar or N_sim x 1 vector 
%                 - confidenceLevel:   scalar
%                 - systematicRisk:    N_sim x 1 vector
%
% @outputs:       - CapitalRequirement: scalar

defaultBarrier = norminv(defaultRate);

loss = (1-recoveryRate).*normcdf((defaultBarrier -...
            sqrt(correlation).*systematicRisk)./sqrt(1-correlation));

VaR = prctile(loss,confidenceLevel*100);
expectedLoss = mean(loss);

CapitalRequirement = VaR - expectedLoss;

end

