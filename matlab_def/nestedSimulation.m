function [rho_nested,DR_nested]=nestedSimulation(d_hat,d_std,nSim_d,CRsurface,rho)
% [rho_nested,DR_nested]=NESTEDSIMULATION(nSim_rho,nSim_d,CRsurface,rho)
%
% Simulates paired vectors of correlations and default rates, given the
% Cramer Rao variance for the dafault barrier, and a vector of simulated 
% rho 
%
% @inputs:                 - d_hat:  value of d observed
%                          - nSim_d: number of simulation of d for each rho
%                          - CRsurface: struct containing d and rho fields,
%                            that determine a grid on which the field surf
%                            is defined
%                          - rho:    vector of given simulated correlation
% @outputs:                - rho_nested: vector of simulated rho
%                          - d_nested:   vector of simulated d

nSim_rho   = length(rho);                 % number of given rho
rho_nested = zeros(nSim_rho*nSim_d,1);    % preallocate rho_nested
d_nested   = zeros(nSim_rho*nSim_d,1);    % preallocate d_nested

d_vect = linspace(-4,4,1000);             % vector on which evaluate posterior distribution
xd_vect = linspace(-5,5,1000); 
f_x_dCR = @(x) normpdf(x,d_vect,d_std);           % density of x given mean(unknown) and std
  
for i = 1 : nSim_rho                      % cycle over the desired rho
    rho_nested(1+(i-1)*nSim_d:i*nSim_d)=rho(i);
    d_CRstd = interp2(CRsurface.d,CRsurface.rho,CRsurface.surf,...
                      d_vect,rho(i),'spline');             % CramerRao 
    h_d_CR = posteriorDistributionD(d_hat,d_vect,d_CRstd); % posterior distribution
    xd_pdfCR = bayesianPrediction(xd_vect,d_vect,h_d_CR,f_x_dCR); % density of X conditioned on data
    d_nested(1+(i-1)*nSim_d:i*nSim_d) = samplingFromPosterior(nSim_d,xd_vect,xd_pdfCR);   % simulated d from its density
end 

DR_nested = normcdf(d_nested);            % default rates

end