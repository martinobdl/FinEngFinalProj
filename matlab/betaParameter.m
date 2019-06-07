function [alpha,beta] = betaParameter(mu,sigma)
% [alpha,beta] = BETAPARAMETER(mu,sigma) 
%
% Computes alpha e beta parameters of a beta distribution given mean and
% standard deviation.
% 
% @inputs:        - mu:     mean of beta distribution
%                 - sigma:  standard deviation of beta distribution
% @outputs:       - alpha:  alpha patameter of beta distribution
%                 - beta:   beta parameter of beta distribution

% Check admissibility of inputs
if min(mu-1/2*(1-sqrt(1-4*sigma.^2)))<0 ...
    || max(mu-1/2*(1+sqrt(1-4*sigma.^2)))>0
    error("mu and sigma are not acceptable")
end

% Computations
beta  = (1-mu).*(mu.*(1-mu)-sigma.^2)./(sigma.^2);
alpha = mu./(1-mu).*beta;

end
