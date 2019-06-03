import csv
from numpy import *
from scipy.stats import norm, shapiro
from scipy.integrate import *
from scipy.special import binom

def CapitalRequirementNominalLHP(recoveryRate,defaultRate,correlation,confidenceLevel):
    lossGivenDefault = (1 - recoveryRate)
    valueAtRisk = lossGivenDefault*norm.cdf((norm.ppf(defaultRate) - sqrt(correlation)* norm.ppf(1 - confidenceLevel))
                                            / sqrt(1 - correlation))
    expectedLoss = lossGivenDefault*defaultRate

    CapitalRequirement = valueAtRisk - expectedLoss
    return (CapitalRequirement)

def CapitalRequirementNominalHP(recoveryRate,defaultRate,correlation,confidenceLevel,numberOfObligor):
    lossGivenDefault = (1 - recoveryRate)
    exposureAtDefault = 1/numberOfObligor
    prob_y = lambda y: norm.cdf((norm.ppf(defaultRate)-sqrt(correlation)*y)/sqrt(1-correlation))
    prob_m = lambda y,m: prob_y(y)**m * (1-prob_y(y))**(numberOfObligor-m)*binom(numberOfObligor,m)*norm.pdf(y)
    prob =zeros([1,numberOfObligor+1])
    for m in range(numberOfObligor):
        prob[0,m] = quad(prob_m,-30,30,args=(m))[0]
    n = linspace(0,numberOfObligor,numberOfObligor+1)
    prob_cumul = cumsum(prob)
    valueAtRisk = min(argwhere(prob_cumul>confidenceLevel))*lossGivenDefault*exposureAtDefault
    expectedLoss = vdot(prob, n)*lossGivenDefault*exposureAtDefault
    CapitalRequirement = float(valueAtRisk - expectedLoss)
    return(CapitalRequirement)

def correlationFromBasel2(defaultRate):
    rhoB = 0.12*(1-exp(-50*defaultRate))/(1-exp(-50))+0.24*(1-(1-exp(-50*defaultRate))/(1-exp(-50)))
    return rhoB


with open('dati_Moody.csv') as csvfile:
    readCSV = csv.reader(csvfile, delimiter=',')
    years = []
    DR_SG = []
    DR_AR = []
    RR = []
    for row in readCSV:
        year = int(row[0])
        SG = float(row[1])
        AR = float(row[2])
        rr = float(row[3])
        years.append(year)
        DR_SG.append(SG)
        DR_AR.append(AR)
        RR.append(rr)

        

DR_mean = mean(DR_SG)
DR_std = std(DR_SG)

# DR_mean = mean(DR_AR)
# DR_std = pstdev(DR_AG)

d_mean = norm.ppf(DR_mean)
# trovare lo 0

RR_mean = mean(RR)
RR_std = std(RR)

rho_mean = 0.0924
rho_std = 0.0388
rho_B_mean = correlationFromBasel2(DR_mean)

N_ob = 50
N_sim = 5*10**5

CL1 = 0.99
CL2 = 0.999

#idiosyncraticRisk = random.normal(0,1,[N_sim ,N_ob])
#systematicRisk = random.normal(0,1,[N_sim,1])

# A.i
    # Large Homogeneous Portfolio
RC_LHP_CL1 = CapitalRequirementNominalLHP(RR_mean,DR_mean,rho_mean,CL1)
RC_LHP_CL2 = CapitalRequirementNominalLHP(RR_mean,DR_mean,rho_mean,CL2)
    # Homogeneous PortFolio
RC_HP_CL1 = CapitalRequirementNominalHP(RR_mean,DR_mean,rho_mean,CL1,N_ob)
RC_HP_CL2 = CapitalRequirementNominalHP(RR_mean,DR_mean,rho_mean,CL2,N_ob)
    # Results
print('RC LHP 99% = ', RC_LHP_CL1)
print('RC LHP 99.9% = ',RC_LHP_CL2)
print('RC HP 99% = ',RC_HP_CL1)
print('RC HP 99.9% = ',RC_HP_CL2)

# A.ii
    # Large Homogeneous Portfolio
RC_LHP_CL1_B = CapitalRequirementNominalLHP(RR_mean,DR_mean,rho_B_mean,CL1)
RC_LHP_CL2_B = CapitalRequirementNominalLHP(RR_mean,DR_mean,rho_B_mean,CL2)
    # Homogeneous PortFolio
RC_HP_CL1_B = CapitalRequirementNominalHP(RR_mean,DR_mean,rho_B_mean,CL1,N_ob)
RC_HP_CL2_B = CapitalRequirementNominalHP(RR_mean,DR_mean,rho_B_mean,CL2,N_ob)
    # Results
print('RC LHP 99% basel = ',RC_LHP_CL1_B)
print('RC LHP 99.9% basel = ',RC_LHP_CL2_B)
print('RC HP 99% basel = ',RC_HP_CL1_B)
print('RC HP 99.9% basel = ',RC_HP_CL2_B)

# B.1

d = norm.ppf(DR_SG)
    # Shapiro-Wilk test
[w_d,pValue_d]=shapiro(d)
[w_RR,pValue_RR]=shapiro(RR)
    # Results
print("pValue_d =",pValue_d)
print("pValue_RR =",pValue_RR)
