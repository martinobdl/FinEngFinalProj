function h_posterior = posteriorDistributionD(d,d_vect,d_std)
% h_posterior = POSTERIORDISTRIBUTIONRHO(d_hat,d_vect,d_std)
%
% Computes posterior distribution when the prior is standard normal and the
% likelihood is a normal given the observed d_hat, with mean d_vect and
% standard deviation d_std.
% 
% @inputs:    - d:  observations on which the posterior is conditioned                    
%             - d_vect: vector of points on which evaluate posterior 
%             - d_std:  vector of standard deviations on d_vect points
% @outputs:   - h_posterior: vector of posterior distribution on d_vect
%
    
tmp         = normpdf(d_vect).*prod(normpdf(d,d_vect,d_std),1); 
h_posterior = tmp./(trapz(d_vect,tmp));

end

