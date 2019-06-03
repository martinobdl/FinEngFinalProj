function h = posteriorDistributionRho(rho_hat,rho_vect,aa,bb)
% Computes posterior distribution when the prior is uniform and the
% likelihood is a beta given the observed rho_hat
% @inputs:  rho_hat, rho_vect, aa, bb
%
% @outputs: h
%
    tmp = betapdf(rho_hat,aa,bb);
    h = tmp./(trapz(rho_vect,tmp));
end

