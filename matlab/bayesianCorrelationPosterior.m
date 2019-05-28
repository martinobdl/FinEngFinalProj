function [aa,bb] = bayesianCorrelationPosterior(sigma, rho_vect)
% data = readData()
% Read the data file into a struct with
% year, specultaive grade default rate, total defalut rate
% and recovery rate.
%
% @inputs: None.
%
% @outputs: data: struct with year, DG_SG, DG_All, RR
%

if min(rho_vect-1/2*(1-sqrt(1-4*sigma.^2)))<0 || max(rho_vect-1/2*(1+sqrt(1-4*sigma.^2)))>0
    error("rho_vect is too small or too big")
end

bb = (1-rho_vect).*(rho_vect.*(1-rho_vect)-sigma.^2)./(sigma.^2);
aa = rho_vect./(1-rho_vect).*bb;


end

