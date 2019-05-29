function CR_std = CRrho(Nob,rho_vect)
    CR_std = sqrt((2*(1-rho_vect).^2.*(1+(Nob-1).*rho_vect).^2)./(Nob*(Nob-1)));
end

