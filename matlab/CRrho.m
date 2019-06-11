function rho_CRstd = CramerRao_rho(Nob,rho_vect,T)
% rho_CRstd = CRAMERRAO_RHO(Nob,rho_vect,T)
% Computes Cramer Rao bound for the correlation 
%
% @input:           - Nob: scalar, numbero of obligors
%                   - rho_vect: vector
% @output:          - CR_std: vector

rho_CRstd = sqrt((2*(1-rho_vect).^2.*(1+(Nob-1).*rho_vect).^2)./...
                        (T*Nob*(Nob-1)));
end

