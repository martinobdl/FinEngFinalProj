function [CapitalRequirement] = CapitalRequirementNominalHP(recoveryRate,...
                    defaultRate,correlation,confidenceLevel,nObligors)
% Computes the capital requirement for a homogeneous portfolio at a given
% confidence level
%
% @inputs: recoveryRate
%
% @outputs: CapitalRequirement
%
 
exposureAtDefault   = 1/nObligors;
lossGivenDefault    = 1-recoveryRate;
defaultBarrier = norminv(defaultRate);
binomialcoef = @(N,n) factorial(N)./(factorial(n).*factorial(N-n));

p = @(y) normcdf((defaultBarrier-sqrt(correlation).*y)./(sqrt(1-correlation)));
cond_law = @(y,m) (p(y).^m).*(1-p(y)).^(nObligors-m);
prob_m = @(m) integral(@(y) normpdf(y).*cond_law(y,m).*binomialcoef(nObligors,m),-30,30,'ArrayValued',true);
M=(0:nObligors)';
P=prob_m(M);    

P_cum = cumsum(P);
idx = find(P_cum>=confidenceLevel, 1 )-1;
meanloss = exposureAtDefault*P'*M*lossGivenDefault;
CapitalRequirement = idx*exposureAtDefault*lossGivenDefault - meanloss;


end