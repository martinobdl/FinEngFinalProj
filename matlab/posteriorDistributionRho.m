function h_posterior = posteriorDistributionRho(rho_hat,rho_vect,alpha,beta)
% h = POSTERIORDISTRIBUTIONRHO(rho_hat,rho_vect,alpha,beta)
%
% Computes posterior distribution when the prior is uniform and the
% likelihood is a beta given the observed rho_hat
% @inputs:    - rho_hat:  observation on which the posterior is conditioned                    
%             - rho_vect: vector of points on which evaluate posterior 
%             - alpha:    vector of alpha parameter assiciated to rho_vect       
%             - beta:     vector of beta parameter associated to rho_vect
%
% @outputs:   - h_posterior: vector of posterior distribution on rho_vect
    
tmp         = betapdf(rho_hat,alpha,beta); 
h_posterior = tmp./(trapz(rho_vect,tmp));

end

