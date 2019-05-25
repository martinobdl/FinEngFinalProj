esatto = 0.051307884733499


a=5.1083;
b=50.1766;
rho_sim = betarnd(a,b,1e6,1);
z_sim = randn(1e6,1);

[mu_RR,sigma_RR] = normfit(data.RR);
RR_sim = randn(1e6,1)*sigma_RR+mu_RR;
d = norminv(data.DG_SG);
[mu_d,sigma_d] = normfit(d);

d_sim = randn(1e6,1)*sigma_d+mu_d;

DR_SGsim = normcdf(d_sim,mu_d,sigma_d);

% d simulato

loss = (1-recoveryMean).*normcdf((d_sim - sqrt(rhoMean).*z_sim)./sqrt(1-rhoMean));

VaR = prctile(loss,confidenceLevel1*100);

lossmean = mean(loss);

aCR_d = VaR-lossmean

AddOn_d = aCR_d/esatto-1

% RR simulato

loss = (1-RR_sim).*normcdf((norminv(defaultRateSGMean) - sqrt(rhoMean).*z_sim)./sqrt(1-rhoMean));

VaR = prctile(loss,confidenceLevel1*100);

lossmean = mean(loss);

aCR_RR = VaR-lossmean

AddOn_RR = aCR_RR/esatto -1
% rho simulato

loss = (1-recoveryMean).*normcdf((norminv(defaultRateSGMean) - sqrt(rho_sim).*z_sim)./sqrt(1-rho_sim));

VaR = prctile(loss,confidenceLevel1*100);

lossmean = mean(loss);

aCR_rho = VaR-lossmean

AddOn_rho = aCR_rho/esatto -1

% caso globale

loss = (1-RR_sim).*normcdf((d_sim - sqrt(rho_sim).*z_sim)./sqrt(1-rho_sim));

VaR = prctile(loss,confidenceLevel1*100);

lossmean = mean(loss);

aCR_all = VaR- lossmean

AddOn_all = aCR_all/esatto -1
