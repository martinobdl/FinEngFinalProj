function expectedLL = ExpectedSecondDerivativeLL(correlation, defaultBarrier, nObligors)


epsilon = 1e-6;

defaultBarrier_pluseps = defaultBarrier + epsilon;
defaultBarrier_mineps = defaultBarrier - epsilon;
binomialcoef = @(N,n) factorial(N)./(factorial(n).*factorial(N-n));

p = @(y) normcdf((defaultBarrier-sqrt(correlation)*y)/...
                        (sqrt(1-correlation)));
p_plus = @(y) normcdf((defaultBarrier_pluseps-sqrt(correlation)*y)/...
                        (sqrt(1-correlation)));
p_min = @(y) normcdf((defaultBarrier_mineps-sqrt(correlation)*y)/...
                        (sqrt(1-correlation)));

cond_law = @(y,m) (p(y).^m).*(1-p(y)).^(nObligors-m).*binomialcoef(nObligors,m);
cond_law_plus = @(y,m) (p_plus(y).^m).*(1-p_plus(y)).^(nObligors-m).*binomialcoef(nObligors,m);
cond_law_min = @(y,m) (p_min(y).^m).*(1-p_min(y)).^(nObligors-m).*binomialcoef(nObligors,m);

m=(0:nObligors)';

P_m = integral(@(y) normpdf(y).*cond_law(y,m),-30,30,'ArrayValued',true);
P_m_plus = integral(@(y) normpdf(y).*cond_law_plus(y,m),-30,30, 'ArrayValued',true);
P_m_min = integral(@(y) normpdf(y).*cond_law_min(y,m),-30,30,'ArrayValued',true);

LL = log(P_m);
LL_pluseps = log(P_m_plus);
LL_mineps = log(P_m_min);

first_plus = (LL_pluseps-LL)./epsilon;
first_min = (LL-LL_mineps)./epsilon;

second = (first_plus - first_min)./epsilon;

expectedLL = second'*P_m;

end
