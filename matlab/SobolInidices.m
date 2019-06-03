function [S1,S2,S3] = SobolInidices(recovery,defaultRate)
    
    nSim = 1e6;
    [mu1,s1] = normfit(recovery);
    [mu2,s2] = normfit(norminv(defaultRate));
    a = 5.1083;
    b = 50.1766;
    confidenceLevel = 0.99;
    
    sim1 = mu1 + s1*randn(nSim,1);
    sim2 = mu2 + s2*randn(nSim,1);
    sim3 = betarnd(a,b,nSim,1);
    
    sim1 = rand(nSim,1);
    sim2 = rand(nSim,1);
    sim3 = rand(nSim,1);
    
    l = 100;
    e = 1e-3;
    
    discretization_V1 = linspace(e,1-e,l);
    discretization_V2 = linspace(e,1-e,l);
    discretization_V3 = linspace(e,1-e,l);
    
    e1 = zeros(length(discretization_V1),1);
    e2 = zeros(length(discretization_V1),1);
    e3 = zeros(length(discretization_V1),1);
    
    for i=1:length(discretization_V1)
        v1 = discretization_V1(i);
        v2 = discretization_V2(i);
        v3 = discretization_V3(i);
        x1 = CapitalRequirementNominalLHP(v1,...
    normcdf(sim2),sim3,confidenceLevel);
        x2 = CapitalRequirementNominalLHP(sim1,...
    normcdf(v2),sim3,confidenceLevel);
        x3 = CapitalRequirementNominalLHP(sim1,...
    normcdf(sim2),v3,confidenceLevel);
        e1(i) = var(x1);
        e2(i) = var(x2);
        e3(i) = var(x3);
    end
    
    S = mean([e1,e2,e3])';
    S1 = S(1)/sum(S);
    S2 = S(2)/sum(S);
    S3 = S(3)/sum(S);
    
end

