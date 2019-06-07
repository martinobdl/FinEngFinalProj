function h_posterior = posteriorDistributionD(d,d_vect,d_CR_std)
% h = POSTERIORDISTRIBUTIONRHO(d_hat,d_vect,d_CR_std)
%
% Computes posterior distribution when the prior is standard normal and the
% likelihood is a normal given the observed d_hat
% @inputs:    - d_hat:  observation on which the posterior is conditioned                    
%             - d_vect: vector of points on which evaluate posterior 
%
% @outputs:   - h_posterior: vector of posterior distribution on d_vect
    
tmp         = normpdf(d_vect).*normpdf(d,d_vect,d_CR_std); 
h_posterior = tmp./(trapz(d_vect,tmp));

end

