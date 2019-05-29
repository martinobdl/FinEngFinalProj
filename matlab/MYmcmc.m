function X=MYmcmc(x0,loglikelihood,logprior,stepsize,nSim)

X=nan(nSim,1);

logPrior=logprior(x0);
logL=loglikelihood(x0);
X(1)=x0; 

for ii=2:nSim
    x1=x0+randn()*stepsize;
    proposed_logprior=logprior(x1);
    if log(rand)<proposed_logprior-logPrior
        proposed_logL=loglikelihood(x1);
        if log(rand)<proposed_logL-logL
            x0=x1;
            logL=proposed_logL;
            logPrior=proposed_logprior;
        end
    end
    X(ii)=x0;
end
