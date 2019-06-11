function rho_CRstd = CramerRao_rho(Nob,rho,T)
% rho_CRstd = CRAMERRAO_RHO(Nob,rho_vect,T)
% Computes Cramer Rao bound, standard deviation, for the correlation.  
%
% @input:           - Nob:      scalar, number of obligors
%                   - rho:      scalar or vector on which compute the limit
% @output:          - rho_CRstd:scalar or vector

rho_CRstd = sqrt((2*(1-rho).^2.*(1+(Nob-1).*rho).^2)./...
                        (T*Nob*(Nob-1)));
end

