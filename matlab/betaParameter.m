function [ alpha , beta ] = betaParameter(mu,sigma)
% Computes alpha e beta parameters of a beta distribution given mean and
% standard deviation.
% INPUTS
% - mu:           mean of the Beta distribution
% - sigma:        standard deviation of the Beta distibution
% OUTPUTS
% - alpha:        alpha parameter of Beta distribution
% - beta:         beta parameter of Beta distribution

beta = (1-mu).*(mu.*(1-mu)-sigma.^2)./(sigma.^2);
alpha = mu./(1-mu).*beta;

end
