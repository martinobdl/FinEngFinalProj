function CapitalRequirement = CapitalRequirementNominalHP(recoveryRate,...
             defaultRate,correlation,confidenceLevel,nObligors)
% CapitalRequirement = CapitalRequirementNominalHP(recoveryRate,
% defaultRate,correlation,confidenceLevel,nObligors)
%
% Computes the capital requirement for a Homogeneous Portfolio, given 
% market parameters at a given confidence level, nominal model.
%
% @inputs:        - recoveryRate:      scalar
%                 - defaultRate:       scalar
%                 - correlation:       scalar 
%                 - confidenceLevel:   scalar
%                 - nObligors:         scalar
%
% @outputs:       - CapitalRequirement: scalar
 
exposureAtDefault   = 1/nObligors;
lossGivenDefault    = 1-recoveryRate;
defaultBarrier      = norminv(defaultRate);

% binomial coefficient
binomialcoef = @(N,n) factorial(N)./(factorial(n).*factorial(N-n));

% probability of one default given market risk
p        = @(y) normcdf((defaultBarrier-sqrt(correlation).*y)./...
                                (sqrt(1-correlation)));
% probability of m defaults given market risk
cond_law = @(y,m) (p(y).^m).*(1-p(y)).^(nObligors-m).*...
                                binomialcoef(nObligors,m);
% probability of m defaults
prob_m   = @(m) integral(@(y) normpdf(y).*cond_law(y,m),-30,30,...
                                'ArrayValued',true); 
                            
M        = (0:nObligors)';  % vector of possible number of defaults
P        = prob_m(M);       % probability of m defaults for m = 0:nObligors

P_cum = cumsum(P);          % cumulative density function of P
idx = min(find(P_cum>=confidenceLevel))-1; % quantile index

valueAtRisk  = idx*exposureAtDefault*lossGivenDefault;
expectedLoss = exposureAtDefault*P'*M*lossGivenDefault;

CapitalRequirement = valueAtRisk - expectedLoss;

end