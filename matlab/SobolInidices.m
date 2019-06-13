function [S1,V] = SobolInidices(parameterSimulations)
% [S1,V] = SOBOLINDICES(parameterSimulations)
%
% Computes the Sobol's Indices for the CapitalRequirementNominalLHP
%
% @inputs:            - parameterSimulations: simulation of the parameter of size: Nsim x Nparam
% @outputs:           - S1: Sobol 1st level indices
%                     - V: variances of the S1 indices

    M = parameterSimulations;
    nSim = size(M,1);

    confidenceLevel = 0.99;

    A = M(1:nSim/2,:);
    B = M(1+nSim/2:end,:);

    S = [0,0,0];
    V = [0,0,0];

    Y_A = CapitalRequirementNominalLHP(A(:,2),normcdf(A(:,1)),A(:,3),confidenceLevel);
    Y_B = CapitalRequirementNominalLHP(B(:,2),normcdf(B(:,1)),B(:,3),confidenceLevel);

    for i=1:3
        BA = B;
        BA(:,i) = A(:,i);
        Y = Y_A.*...
        (CapitalRequirementNominalLHP(BA(:,2),normcdf(BA(:,1)),BA(:,3),confidenceLevel)+...
            -Y_B);
        S(i) = mean(Y);
        V(i) = var(Y);

    end

    S1 = S/var([Y_A;Y_B]);

end
