function X=mcmc(x0,loglikelihood,logprior,stepsize,nSim)
% X = MCMC(x0,loglikelihood,logprior,stepsize,nSimf)
%
% Samples nSim samples from a Bayesian Posterior distribution described by loglikelihood and logprior
%
%
% @inputs:               - x0: starting point of the Markov Chain
%                        - loglikelihood: handle function of the loglikelihood
%                        - logprior: handle function of the logprior
%                        - stepsize: variance of the symmetric distribution used for generating the porposals
%                        - nSim: number of simulations
% @outputs:              - X: Samples drawn
%

X=nan(nSim,1);

logPrior=logprior(x0);
logL=loglikelihood(x0);
X(1)=x0;
f = waitbar(0,'');
sim = randn(nSim-1,1);

for ii=2:nSim
    waitbar(ii/nSim,f,sprintf('Please wait: %2.2f%%',ii*100/nSim))
    x1=x0+sim(ii-1)*stepsize;
    proposed_logprior=logprior(x1);
    proposed_logL=loglikelihood(x1);
    if log(rand)<proposed_logprior-logPrior+proposed_logL-logL
        x0=x1;
        logL=proposed_logL;
        logPrior=proposed_logprior;
    end
    X(ii)=x0;
end

F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F)

end
