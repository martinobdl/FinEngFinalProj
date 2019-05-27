function h = posteriorDistributionRho(rho_hat,rho_vect,aa,bb)
    tmp = betapdf(rho_hat,aa,bb);
    h = tmp./(trapz(rho_vect,tmp));
end

