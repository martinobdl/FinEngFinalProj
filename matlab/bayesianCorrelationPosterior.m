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

sigma_beta = @(a,rho) sqrt((1-rho)*rho^2/(a+rho));
b_a = @(a,rho) a*(1-rho)/rho;

aa = zeros(size(rho_vect));
bb = zeros(size(rho_vect));

fun = @(a,rho) sigma_beta(a,rho)-sigma;

for i=1:length(aa)
    r = rho_vect(i);
    eps = 1e-6;
    a = fzero(@(a) fun(a,r),[eps,10000]);
    b = b_a(a,r);
    aa(i) = a;
    bb(i) = b;
end

end

