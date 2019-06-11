

rho_sim = samplingFromPosterior(1000,rho_vect,H_rho);

d_vect = -4:0.005:4;
CR_tmp = interp2(CRsurfno20.d,CRsurfno20.rho,CRsurfno20.surface,d_vect,rho_sim,'spline');
rho_post=zeros(1e6,1);
Cr = size(1e3,1);
d_sim=zeros(1e6,1);
for i = 1:1000
   rho_post(1+(i-1)*1e3:1e3*i)=rho_sim(i);
   tmp = normpdf(d_vect,0,1).*prod(normpdf(d,d_vect,CR_tmp(i,:)));
   d_H = tmp./trapz(d_vect,tmp);
   d_sim(1+(i-1)*1e3:1e3*i)=samplingFromPosterior(1e3,d_vect,d_H);
   %d_sim=samplingFromPosterior(1e3,d_vect,d_H);
   %Cr(i)=CapitalRequirementAlternativeHP(RR_mean,normcdf(d_sim),rho_sim(i),CL1,...
    %        randn(1e3,1),randn(1e3,N_ob));
end


rho_mean = 0.0924;
d_vect = -4:0.005:4;
CR_tmp = interp2(CRsurfno20.d,CRsurfno20.rho,CRsurfno20.surface,d_vect,rho_mean,'spline');

   tmp = normpdf(d_vect,0,1).*normpdf(d_hat,d_vect,CR_tmp/sqrt(20));
   d_H = tmp./trapz(d_vect,tmp);
   d_sim=samplingFromPosterior(1e6,d_vect,d_H);
   %d_sim=samplingFromPosterior(1e3,d_vect,d_H);
   %Cr(i)=CapitalRequirementAlternativeHP(RR_mean,normcdf(d_sim),rho_sim(i),CL1,...
    %        randn(1e3,1),randn(1e3,N_ob));




 

 