function X = samplingFromPosterior(nSim,x,y)
% X = SAMPLINGFROMPOSTERIOR(nSim,x,y)
%
% Simulates a vector of variables extracted from an empirical density using
% the cumulative density function inversion sampling.
%
% @inputs:               - nSim: number of variable to be simulated
%                        - x: vector of points on which density in given
%                        - y: density from which extract simulations
% @outputs:              - X: sampled vector

cdf = cumtrapz(x,y);             % cumulative density function 
[cdf_unique, idx] = unique(cdf); % solve the uniqueness problem for interp1

u = rand(nSim,1);                % uniform r.v.

X = interp1(cdf_unique,x(idx),u);%sampled vector
    
end

