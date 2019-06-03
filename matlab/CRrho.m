function rho_CRstd = CRrho(Nob,rho_vect,T)
%Computes Cramer Rao bound for the correlation 
%
%@input: Nob, rho_vect
%@output: CR_std
    rho_CRstd = sqrt((2*(1-rho_vect).^2.*(1+(Nob-1).*rho_vect).^2)./(T*Nob*(Nob-1)));
end

