"""
Collection of main functions for the Final Project of the Course in Fin Eng in capital requirement
model risk.
"""

import numpy as np
import scipy
import scipy.stats
import scipy.integrate
import scipy.special

def CapitalRequirementNominalLHP(recoveryRate,defaultRate,correlation,confidenceLevel):
    """
    CapitalRequirement = CAPITALREQUIREMENTNOMINALLHP(recoveryRate,
    defaultRate,correlation,confidenceLevel)

    Computes the capital requirement for a Large Homogeneous Portfolio, given
    market parameters at a given confidence level, nominal model.

    @inputs:        - recoveryRate:      scalar
                    - defaultRate:       scalar
                    - correlation:       scalar
                    - confidenceLevel:   scalar
    @outputs:       - CapitalRequirement: scalar
    """
    lossGivenDefault = (1 - recoveryRate)
    valueAtRisk = lossGivenDefault*scipy.stats.norm.cdf((scipy.stats.norm.ppf(defaultRate) - np.sqrt(correlation)* scipy.stats.norm.ppf(1 - confidenceLevel))
                                            / np.sqrt(1 - correlation))
    expectedLoss = lossGivenDefault*defaultRate

    CapitalRequirement = valueAtRisk - expectedLoss
    return CapitalRequirement

def CapitalRequirementNominalHP(recoveryRate,defaultRate,correlation,confidenceLevel,numberOfObligor):
    """
    CapitalRequirement = CapitalRequirementNominalHP(recoveryRate,
    defaultRate,correlation,confidenceLevel,nObligors)

    Computes the capital requirement for a Homogeneous Portfolio, given
    market parameters at a given confidence level, nominal model.
    @inputs:        - recoveryRate:      scalar
                    - defaultRate:       scalar
                    - correlation:       scalar
                    - confidenceLevel:   scalar
                    - nObligors:         scalar
    @outputs:       - CapitalRequirement: scalar
    """
    lossGivenDefault = (1 - recoveryRate)
    exposureAtDefault = 1/numberOfObligor
    prob_y = lambda y: scipy.stats.norm.cdf((scipy.stats.norm.ppf(defaultRate)-np.sqrt(correlation)*y)/np.sqrt(1-correlation))
    prob_m = lambda y,m: prob_y(y)**m * (1-prob_y(y))**(numberOfObligor-m)*scipy.special.binom(numberOfObligor,m)*scipy.stats.norm.pdf(y)
    prob = np.zeros([1,numberOfObligor+1])
    for m in range(numberOfObligor):
        prob[0,m] = scipy.integrate.quad(prob_m,-30,30,args=(m))[0]
    n = np.linspace(0,numberOfObligor,numberOfObligor+1)
    prob_cumul = np.cumsum(prob)
    valueAtRisk = min(np.argwhere(prob_cumul>confidenceLevel))*lossGivenDefault*exposureAtDefault
    expectedLoss = np.vdot(prob, n)*lossGivenDefault*exposureAtDefault
    CapitalRequirement = float(valueAtRisk - expectedLoss)
    return(CapitalRequirement)

def CapitalRequirementAlternativeLHP(recoveryRate,defaultRate,correlation,confidenceLevel,systematicRisk):
    """
    CapitalRequrement = CAPITALREQUIREMENTALTERNATIVELHP(recoveryRate,
    defaultRate,correlation,confidenceLevel,systematicRisk)

    Computes the Capital Requirement in the case of a Large Homogeneous
    Portfolio given the recovery, default probability and correlation at the
    required confidence level. Inputs can be vectors.

    @inputs:        - recoveryRate:      scalar or N_sim x 1 numpy array
                    - defaultRate:       scalar or N_sim x 1 numpy array
                    - correlation:       scalar or N_sim x 1 numpy array
                    - confidenceLevel:   scalar
                    - systematicRisk:    N_sim x 1 vector

    @outputs:       - CapitalRequirement: scalar
    """
    defaultBarrier = scipy.stats.norm.ppf(defaultRate)
    loss = (1-recoveryRate)*scipy.stats.norm.cdf((defaultBarrier -
        np.sqrt(correlation)*systematicRisk)/np.sqrt(1-correlation))
    VaR = np.percentile(loss,confidenceLevel*100)
    expectedLoss = np.mean(loss)
    CapitalRequirement = VaR - expectedLoss
    return CapitalRequirement

def CapitalRequirementAlternativeHP(
            recoveryRate,defaultRate,correlation,confidenceLevel,
            systematicRisk,idiosyncraticRisk):
    """
    CapitalRequirement = CapitalRequirementNominalHP(recoveryRate,
    defaultRate,correlation,confidenceLevel,nObligors)

    Computes the capital requirement for a Homogeneous Portfolio, given
    market parameters at a given confidence level, nominal model.

    @inputs:        - recoveryRate:      scalar
                    - defaultRate:       scalar
                    - correlation:       scalar
                    - confidenceLevel:   scalar
                    - nObligors:         scalar

    @outputs:       - CapitalRequirement: scalar
    """
    nObligors          = np.shape(idiosyncraticRisk)[1]
    exposureAtDefault  = 1./nObligors
    lossGivenDefault   = 1-recoveryRate

    defaultBarrier     = scipy.stats.norm.ppf(defaultRate)
    firmValues         = (np.sqrt(correlation)*systematicRisk +\
                            np.sqrt(1-correlation)*idiosyncraticRisk.T).T
    numberOfDefaults   = np.sum(firmValues.T < defaultBarrier,axis=0).T

    loss               = exposureAtDefault*lossGivenDefault*numberOfDefaults

    valueAtRisk        = np.percentile(loss,confidenceLevel*100)
    expectedLoss       = np.mean(loss)

    CapitalRequirement = valueAtRisk - expectedLoss

    return CapitalRequirement

def beyesianPrediction():


