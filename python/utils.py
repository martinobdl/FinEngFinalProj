import numpy as np
import scipy
import scipy.interpolate

def betaParameter(mu,sigma):

    beta  = (1-mu)*(mu*(1-mu)-sigma**2)/(sigma**2)
    alpha = mu/(1-mu)*beta

    return [alpha,beta]

def correlationFromBasel2(defaultRate):
    rhoB = 0.12*(1-np.exp(-50*defaultRate))/(1-np.exp(-50))+0.24*(1-(1-np.exp(-50*defaultRate))/(1-np.exp(-50)))

    return rhoB

def samplingFromPosterior(N_sim,x,y):
    cdf = scipy.integrate.cumtrapz(y,x)
    cdf_unique, idx = np.unique(cdf, return_index=True)
    u = np.random.uniform(0,1,N_sim)

    f = scipy.interpolate.interp1d(cdf_unique, x[idx], kind='slinear', fill_value="extrapolate")

    return f(u)


def mcmc(x0,loglikelihood,logprior,stepsize,nSim):

    X=np.zeros(int(nSim))

    logPrior=logprior(x0)
    logL=loglikelihood(x0)
    X[0]=x0
    sim = np.random.normal(0,1,nSim-1)

    for ii in range(1,nSim):
        x1=x0+sim[ii-1]*stepsize
        proposed_logprior=logprior(x1)
        proposed_logL=loglikelihood(x1)
        if np.log(np.random.uniform())<proposed_logprior-logPrior+proposed_logL-logL:
            x0=x1
            logL=proposed_logL
            logPrior=proposed_logprior
        X[ii]=x0

    return X

def posteriorDistributionRho(rho_hat,rho_vect,alpha,beta):
    tmp         = scipy.stats.beta.pdf(rho_hat,alpha,beta)
    h_posterior = tmp/(np.trapz(tmp,rho_vect))

    return h_posterior

def posteriorDistributionD(d,d_vect,d_CR_std):
    tmp         = scipy.stats.norm.pdf(d_vect)*scipy.stats.norm.pdf(d,d_vect,d_CR_std)
    h_posterior = tmp/(np.trapz(tmp,d_vect))

    return h_posterior

def CramerRao_d(correlation, defaultBarrier, nObligors, T):

    epsilon = 1e-6
    X = np.arange(nObligors+1)

    binomialcoef = scipy.special.factorial(nObligors)/(scipy.special.factorial(X)*scipy.special.factorial(nObligors-X))

    p = lambda s,y: scipy.stats.norm.cdf((defaultBarrier+s*epsilon-np.sqrt(correlation)*y)/
                (np.sqrt(1-correlation)))
    integrand = lambda y,s,m: scipy.stats.norm.pdf(y)*(p(s,y)**m)*(1-p(s,y))**(nObligors-m)*binomialcoef[m]

    P = lambda s: np.array([ np.array([scipy.integrate.quad(lambda y,s,m: integrand(y,s,m)[k],-10,10,args=(s,m,))[0] for m in range(nObligors)]) for k in range(len(defaultBarrier))])

    P_central = P(0)    # without change on defaultBarrier
    P_plus = P(1)       # with indrement on defaultBarrier of epsilon
    P_minus = P(-1)     # with decrement on defaultBarrier of -epsilon

    # Logarithm of the above matrices -> LogLikelihood
    LL_central = np.log(P_central)
    LL_plus = np.log(P_plus)
    LL_minus = np.log(P_minus)

    secondDerivativeLL = (LL_plus-2*LL_central+LL_minus)/epsilon**2
    expected = np.sum(secondDerivativeLL*P_central,axis=1)

    d_CRstd = np.sqrt(-expected**-1/T)

    return d_CRstd

def CramerRao_rho(Nob,rho_vect,T):

    return np.sqrt((2*(1-rho_vect)**2*(1+(Nob-1)*rho_vect)**2)/\
        (T*Nob*(Nob-1)))
