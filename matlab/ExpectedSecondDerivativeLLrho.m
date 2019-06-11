function expectedLL = ExpectedSecondDerivativeLLrho(correlation, defaultBarrier, nObligors)

epsilon = 1e-6;

correlation_pluseps = correlation + epsilon;
correlation_mineps = correlation - epsilon;

Sigma = (correlation).*ones(nObligors);
Sigmadiag = ones(nObligors,1).*(1-correlation);
Diagg = diag(Sigmadiag)+Sigma;

Sigmaplus = (correlation_pluseps).*ones(nObligors);
Sigmadiagplus = ones(nObligors,1).*(1-correlation_pluseps);
Diaggplus = diag(Sigmadiag)+Sigma;

Sigmamin = (correlation_mineps).*ones(nObligors);
Sigmadiag_mineps = ones(nObligors,1).*(1-correlation_mineps);
Diaggmin = diag(Sigmadiag)+Sigma;
densita = @(x)(2*pi)^(-nObligors/2).*sqrt(det(Diagg).*exp(-1./2.(x)'*inv(Diagg)*x);
p = @(x) log((2*pi)^(-nObligors/2).*sqrt(det(Diagg).*exp(-1./2.(x)'*inv(Diagg)*x));
p_plus = @(x) log((2*pi)^(-nObligors/2).*sqrt(det(Diaggplus).*exp(-1./2.(x)'*inv(Diaggplus)*x));
p_min = @(x) log((2*pi)^(-nObligors/2).*sqrt(det(Diaggmin).*exp(-1./2.(x)'*inv(Diaggmin)*x));

first_plus = @(x)(p_plus(x)-p(x))./epsilon;
first_min = @(x)(p(x)-p_min(x))./epsilon;

second = @(X)(first_plus(X) - first_min(X))./epsilon;

expectedLL = integral(@(X) second(X)*,-15,

end